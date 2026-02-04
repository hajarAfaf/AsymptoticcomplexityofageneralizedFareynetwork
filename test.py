import networkx as nx
import pandas as pd
import math
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

def format_huge_number(log_val):
    """
    Utilitaire : Convertit un logarithme naturel (ln) en format scientifique A x 10^B.
    Permet d'afficher des nombres trop grands pour être calculés directement.
    """
    if log_val == 0: return "1"
    log10_val = log_val / math.log(10)
    exposant_b = int(log10_val)
    reste = log10_val - exposant_b
    coefficient_a = 10**reste
    return f"{coefficient_a:.3f} x 10^{exposant_b}"

class FareyMethodAnalyzer:
    def __init__(self, logger_func=None):
        # Utilisation de MultiGraph pour supporter les arêtes parallèles temporaires
        self.G = nx.MultiGraph()
        self.log_tau = 0.0  # On travaille en Log pour éviter l'overflow
        self.initial_node_count = 0
        self.logger = logger_func
    def log(self, message):
        """Envoie le message à l'interface si elle existe, sinon print"""
        if self.logger:
            self.logger(message)
        else:
            print(message)
    def load_snap_data(self, filepath):
        """
        Charge un fichier de données SNAP (.edges).
        Initialise tous les poids (conductances) à 1.0 comme spécifié dans la méthodologie.
        Charge un fichier. Supporte :
        - 2 colonnes : u v (Poids par défaut = 1.0)
        - 3 colonnes : u v w (Poids = w)
        """
        self.log(f"Chargement du fichier : {os.path.basename(filepath)} ...")
        try:
            try:
                # Tentative de lecture avec poids
                df = pd.read_csv(filepath, sep='\s+', header=None, engine='python')
                num_cols = len(df.columns)
            except:
                num_cols = 0

            self.G.clear()
            edges = []

            if num_cols >= 3:
                self.log(" Format détecté : Graphe PONDÉRÉ (3 colonnes)")
                # On suppose que les colonnes sont 0=u, 1=v, 2=weight
                for _, row in df.iterrows():
                    u, v = int(row[0]), int(row[1])
                    w = float(row[2])
                    if u != v and w > 0: # On ignore les poids nuls ou négatifs pour la conductance
                        edges.append((u, v, w))
            else:
                self.log("Format détecté : Graphe NON PONDÉRÉ (poids=1.0)")
                # Lecture standard u, v
                df = pd.read_csv(filepath, sep='\s+', header=None, names=['u', 'v'])
                for _, row in df.iterrows():
                    u, v = row['u'], row['v']
                    if u != v: 
                        edges.append((u, v, 1.0))
            
            self.G.add_weighted_edges_from(edges)
            
            # Initialisation
            self.initial_node_count = self.G.number_of_nodes()
            self.log_tau = 0.0 
            self.log(f" Graphe chargé : {self.G.number_of_nodes()} nœuds, {self.G.number_of_edges()} arêtes.")
            return True
            
        except Exception as e:
            self.log(f" Erreur : {e}")
            return False
    def get_first_edge_weight(self, u, v):
        """Récupère le poids sans planter si l'index n'est pas 0"""
        edge_data = self.G.get_edge_data(u, v)
        if edge_data:
            first_key = list(edge_data.keys())[0]
            return edge_data[first_key]['weight']
        return 0
    
    # --- ALGORITHME 1 : Parallel Edge ---
    def algo_parallel(self):
        """
        Algorithme 1 : Transformation des Arêtes Parallèles.
        Fusionne plusieurs arêtes reliant les deux mêmes nœuds (u, v) en une seule arête.
        Le nouveau poids est la SOMME des poids individuels (Loi des conductances en parallèle).
        """
        has_changed = False
        # On itère sur les nœuds pour trouver ceux qui ont des voisins connectés plusieurs fois
        for u in list(self.G.nodes()):
            for v in set(self.G.neighbors(u)):
                if u >= v: continue
                if self.G.number_of_edges(u, v) > 1:
                    edges = self.G.get_edge_data(u, v)
                    total_weight = sum(d['weight'] for d in edges.values())
                    self.G.remove_edges_from([(u, v, k) for k in edges])
                    self.G.add_edge(u, v, weight=total_weight)
                    has_changed = True
        return has_changed

    # --- ALGORITHME 2 : Serial Edge ---
    def algo_serial(self):
        """
        Algorithme 2 : Transformation des Arêtes en Série.
        Supprime un nœud de degré 2 et fusionne ses deux voisins.
        Si u --(a)-- node --(b)-- v, alors u --(c)-- v.
        Le nouveau poids c est calculé par la formule : (a * b) / (a + b).
        Met à jour l'invariant Tau en ajoutant ln(a + b).
        """
        has_changed = False
        nodes = list(self.G.nodes())
        for node in nodes:
            if not self.G.has_node(node): continue
            if self.G.degree(node) == 2:
                neighbors = list(self.G.neighbors(node))
                # Vérification qu'il s'agit bien de 2 voisins distincts (pas de boucle)
                if len(set(neighbors)) == 2:
                    u, v = neighbors[0], neighbors[1]
                    try:
                        wa = self.G[u][node][0]['weight']
                        wb = self.G[node][v][0]['weight']
                        denom = wa + wb
                        wc = (wa * wb) / denom
                        
                        self.log_tau += math.log(denom) # Mise à jour Tau
                        
                        self.G.remove_node(node)
                        self.G.add_edge(u, v, weight=wc)
                        has_changed = True
                    except KeyError: continue
        return has_changed

    # --- ALGORITHME 3 : Wye-Delta (Star-Mesh pour n=3) ---
    def algo_wye_delta(self):
        """
        Algorithme 3 : Transformation Wye-Delta (Étoile vers Triangle).
        Transforme un nœud de degré 3 en connectant ses 3 voisins entre eux (Triangle).
        Les nouveaux poids sont calculés selon : x = (a*b)/S, y = (a*c)/S, z = (b*c)/S.
        Où S = a + b + c.
        Met à jour Tau en ajoutant ln(S).
        """
        for node in list(self.G.nodes()):
            if not self.G.has_node(node): continue
            if self.G.degree(node) == 3:
                neighbors = list(self.G.neighbors(node))
                if len(set(neighbors)) == 3:
                    try:
                        n1, n2, n3 = neighbors
                        a = self.G[node][n1][0]['weight']
                        b = self.G[node][n2][0]['weight']
                        c = self.G[node][n3][0]['weight']
                        s = a + b + c
                        
                        self.log_tau += math.log(s) # Mise à jour Tau
                        
                        self.G.remove_node(node)
                        self.G.add_edge(n1, n2, weight=(a*b)/s)
                        self.G.add_edge(n1, n3, weight=(a*c)/s)
                        self.G.add_edge(n2, n3, weight=(b*c)/s)
                        return True # On retourne True pour re-nettoyer (Parallel)
                    except KeyError: continue
        return False

    # --- ALGORITHME 5 : Star-Mesh (LA SOLUTION POUR GRAPHES DENSES) ---
    def algo_star_mesh(self):
        """
        Algorithme 5 : Transformation Star-Mesh (Étoile vers Maillage).
        Généralisation de Wye-Delta pour les nœuds de degré k > 3.
        Supprime le nœud central et crée une clique (maillage complet) entre tous ses voisins.
        Chaque nouvelle arête (u, v) a pour poids : (w_u * w_v) / Somme_des_poids.
        Indispensable pour réduire les graphes denses (réseaux sociaux).
        """
        # Stratégie : On prend le nœud avec le PLUS PETIT degré (>3) pour limiter l'explosion des arêtes
        candidates = [n for n in self.G.nodes() if self.G.degree(n) > 3]
        if not candidates:
            return False
            
        # On trie pour traiter le cas le moins coûteux d'abord
        candidates.sort(key=lambda n: self.G.degree(n))
        node = candidates[0] # Le meilleur candidat
        
        neighbors = [n for n in set(self.G.neighbors(node)) if n != node]
        if len(neighbors) < 2: return False
        
        try:
            # 1. Calcul de la somme S (Algorithm 5, ligne 9)
            # Récupération des poids de toutes les arêtes connectées au nœud central
            weights = {}
            total_s = 0.0
            for neighbor in neighbors:
                # On prend la première arête dispo (normalement parallel edge a déjà nettoyé)
                w = self.get_first_edge_weight(node, neighbor)
                weights[neighbor] = w
                total_s += w
            
            # 2. Mise à jour de Tau (Algorithm 5, ligne 16)
            self.log_tau += math.log(total_s)
            
            # 3. Création du maillage (clique) entre les voisins (Lignes 10-14)
            # On connecte chaque voisin à tous les autres
            new_edges = []
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    u, v = neighbors[i], neighbors[j]
                    wu = weights[u]
                    wv = weights[v]
                    
                    # Formule : x = (a * b) / S
                    new_weight = (wu * wv) / total_s
                    new_edges.append((u, v, new_weight))
            
            # 4. Suppression du nœud central
            self.G.remove_node(node)
            
            # 5. Ajout des nouvelles arêtes
            self.G.add_weighted_edges_from(new_edges)
            return True
        except Exception as e:
            self.log(f"Erreur Star-Mesh: {e}")
            return False

    # --- ALGORITHME 6 : MAIN LOOP ---
    def run_analysis(self):
        """
        Exécute la boucle de réduction itérative selon l'ordre de priorité :
        1. Parallèle (Nettoyage)
        2. Série (Simplification rapide)
        3. Wye-Delta (Structure locale)
        4. Star-Mesh (Déblocage des nœuds denses)
        """
        self.log("--- Début de l'analyse par Transformations ---")
        iteration = 0
        while True:
            iteration += 1
            if iteration % 10 == 0:
                self.log(f"Iter {iteration} | Nœuds restants : {self.G.number_of_nodes()} | Arêtes : {self.G.number_of_edges()}")

            # Condition d'arrêt
            if self.G.number_of_nodes() <= 2:
                break

            # 1. Parallel (Toujours en premier pour nettoyer)
            if self.algo_parallel(): continue
            
            # 2. Serial (Le plus efficace)
            if self.algo_serial(): continue
            
            # 3. Wye-Delta (Degré 3)
            if self.algo_wye_delta(): continue
            
            # 4. Star-Mesh (Degré > 3) - C'est lui qui va sauver ton analyse sur SNAP
            # Si on arrive ici, c'est qu'il n'y a ni parallèle, ni série, ni degré 3
            if self.algo_star_mesh(): continue
            
            # Si rien ne marche
            self.log("Aucune transformation possible.")
            break
        
        # Résultat final
        final_tau = self.log_tau
        # Si on finit avec 2 nœuds et X arêtes parallèles (ou 1 arête)
        if self.G.number_of_edges() > 0:
             # On fait une passe finale de nettoyage parallèle pour être sûr
            self.algo_parallel()
            if self.G.number_of_edges() > 0:
                last_edge = list(self.G.edges(data=True))[0]
                final_tau += math.log(last_edge[2]['weight'])

        return final_tau

    def calculate_entropy(self):
        """ Calcule l'entropie asymptotique : Rho = ln(Tau) / N """
        if self.initial_node_count == 0: return 0
        return self.log_tau / self.initial_node_count

# --- 3. L'INTERFACE GRAPHIQUE (GUI) ---
class FareyApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Network Complexity Analyzer")
        self.root.geometry("950x700")
        self.root.configure(bg="#F8F9FA")
        
        self.style_setup()
        
        self.filepath = tk.StringVar()
        self.filepath.set("Aucun fichier sélectionné")
        
        self.create_widgets()

    def style_setup(self):
        style = ttk.Style()
        style.theme_use('clam')
        
        style.configure("Card.TFrame", background="#FFFFFF", relief="flat", borderwidth=0)
        style.configure("Main.TFrame", background="#F8F9FA")
        style.configure("Header.TLabel", background="#F8F9FA", foreground="#4A5568", font=("Segoe UI", 22, "bold"))
        style.configure("SubHeader.TLabel", background="#FFFFFF", foreground="#4A5568", font=("Segoe UI", 12, "bold"))
        style.configure("Path.TLabel", background="#FFFFFF", foreground="#000000", font=("Consolas", 10))
        
        # Boutons Baby Blue
        style.configure("Action.TButton", font=("Segoe UI", 10, "bold"), background="#BA68C8", foreground="white", padding=10, borderwidth=0)
        style.map("Action.TButton", background=[("active", "#F0A7DC")])
        style.configure("Browse.TButton", font=("Segoe UI", 9), background="#E0E0E0", padding=6)
    
    def create_widgets(self):
        main = ttk.Frame(self.root, style="Main.TFrame")
        main.pack(fill="both", expand=True, padx=25, pady=25)
        ttk.Label(main, text="Analyseur de Complexité de Réseaux", style="Header.TLabel").pack(anchor="w", pady=(0, 25))

        # Selection
        card = ttk.Frame(main, style="Card.TFrame")
        card.pack(fill="x", pady=(0, 20))
        inner = tk.Frame(card, bg="#FFFFFF", padx=20, pady=15)
        inner.pack(fill="x")
        
        ttk.Label(inner, text="1. Importation", style="SubHeader.TLabel").pack(anchor="w")
        row = ttk.Frame(inner, style="Card.TFrame")
        row.pack(fill="x", pady=10)
        ttk.Button(row, text="Parcourir...", style="Browse.TButton", command=self.browse_file).pack(side="left")
        tk.Label(row, textvariable=self.filepath, fg="#000000", bg="#FFFFFF").pack(side="left", padx=10)
        
        self.btn_run = ttk.Button(inner, text="▶ LANCER L'ANALYSE", style="Action.TButton", command=self.start_thread, state="disabled")
        self.btn_run.pack(anchor="e")

        # Logs
        log_c = ttk.Frame(main, style="Card.TFrame")
        log_c.pack(fill="both", expand=True, pady=(0, 20))
        ttk.Label(log_c, text="Journal", style="SubHeader.TLabel", background="white").pack(anchor="w", padx=20, pady=10)
        self.log_area = scrolledtext.ScrolledText(log_c, height=8, font=("Consolas", 10), bg="#EBE2E2", fg="#000000", relief="flat")
        self.log_area.pack(fill="both", expand=True, padx=20, pady=(0, 20))

        # Resultats
        res = ttk.Frame(main, style="Card.TFrame")
        res.pack(fill="x")
        res_in = tk.Frame(res, bg="white", padx=20, pady=15)
        res_in.pack(fill="x")
        ttk.Label(res_in, text="3. Résultats & Conclusion", style="SubHeader.TLabel").pack(anchor="w", pady=(0, 15))
        grid = tk.Frame(res_in, bg="white")
        grid.pack(fill="x")
        self.create_card(grid, 0, "Log(Tau)", "Arbres couvrants", "lbl_tau", "lbl_trees")
        self.create_card(grid, 1, "Entropie (Rho)", "Référence Farey: 0.9457", "lbl_rho", None, color="#89CFF0")
        
        self.lbl_concl = tk.Label(res_in, text="En attente...", font=("Segoe UI", 12, "bold"), 
                                  bg="#F8F9FA", fg="#666666", pady=15, width=60)
        self.lbl_concl.pack(pady=20)

    def create_card(self, parent, col, title, subtitle, attr_val, attr_sub, color="#444"):
        f = tk.Frame(parent, bg="#F8F9FA", padx=15, pady=10)
        f.grid(row=0, column=col, sticky="ew", padx=5)
        parent.grid_columnconfigure(col, weight=1)
        tk.Label(f, text=title.upper(), font=("Segoe UI", 9, "bold"), fg="#888", bg="#F8F9FA").pack(anchor="w")
        l = tk.Label(f, text="-", font=("Segoe UI", 16, "bold"), fg=color, bg="#F8F9FA")
        l.pack(anchor="w")
        setattr(self, attr_val, l)
        if attr_sub:
            l_sub = tk.Label(f, text=subtitle, font=("Segoe UI", 10), fg="#666", bg="#F8F9FA")
            l_sub.pack(anchor="w")
            setattr(self, attr_sub, l_sub)
        else:
            tk.Label(f, text=subtitle, font=("Segoe UI", 9, "italic"), fg="#888", bg="#F8F9FA").pack(anchor="w")

    def log_message(self, msg):
        """Fonction appelée par l'analyzer pour écrire dans la GUI"""
        self.root.after(0, lambda: self._write(msg))

    def _write(self, msg):
        self.log_area.config(state='normal')
        self.log_area.insert(tk.END, ">> " + msg + "\n")
        self.log_area.see(tk.END)
        self.log_area.config(state='disabled')

    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("Edges Files", "*.edges"), ("Text Files", "*.txt"), ("All Files", "*.*")])
        if filename:
            self.filepath.set(filename)
            self.btn_run.config(state="normal")
            self.log_message(f"Fichier sélectionné : {os.path.basename(filename)}")
            self.lbl_tau.config(text="-")
            self.lbl_trees.config(text="-")
            self.lbl_rho.config(text="-")
            self.lbl_concl.config(text="Prêt à lancer", fg="#6c757d", bg="#f8f9fa")
            self.reset_display()
    def reset_display(self):
        """Réinitialise l'affichage des résultats avant un nouveau calcul"""
        self.lbl_tau.config(text="-")
        self.lbl_rho.config(text="-")
        if hasattr(self, 'lbl_trees'): 
            self.lbl_trees.config(text="Arbres couvrants")
        self.lbl_concl.config(text="Prêt", bg="#F8F9FA", fg="#666")
    def start_thread(self):
        """Lance le calcul dans un thread séparé pour ne pas geler la fenêtre"""
        if not self.filepath.get(): return
        
        self.btn_run.config(state="disabled")
        self.log_area.config(state='normal')
        self.log_area.delete(1.0, tk.END) # Clear logs
        self.log_area.config(state='disabled')
        
        # Reset labels
        self.lbl_tau.config(text="Log(Tau) : Calcul...")
        self.lbl_rho.config(text="Entropie : Calcul...")
        self.lbl_concl.config(text="Conclusion : ...", fg="#DC99DC", bg="#EAF6FF")

        # Threading
        thread = threading.Thread(target=self.run_process)
        thread.daemon = True
        thread.start()

    def run_process(self):
        analyzer = FareyMethodAnalyzer(logger_func=self.log_message)
        success = analyzer.load_snap_data(self.filepath.get())
        
        if success:
            log_tau = analyzer.run_analysis()
            rho = analyzer.calculate_entropy()
            
            # Mise à jour de l'interface (depuis le thread)
            self.root.after(0, lambda: self.show_results(log_tau, rho))
        else:
            self.root.after(0, lambda: messagebox.showerror("Erreur", "Impossible de lire le fichier."))
            self.root.after(0, lambda: self.btn_run.config(state="normal"))

    def show_results(self, log_tau, rho):
        huge_num = format_huge_number(log_tau)
        self.lbl_tau.config(text=f"Log(Tau) : {log_tau:.4f}")
        if hasattr(self, 'lbl_trees'):
                self.lbl_trees.config(text=f"Arbres: {huge_num}")
        self.lbl_rho.config(text=f"Entropie (Rho) : {rho:.4f}")
        
        if rho > 0.9457:
                msg = "CONCLUSION : PLUS robuste que Farey"
                bg_col = "#C1E1C1"
                fg_col = "#2E5D3B"
        else:
                msg = "CONCLUSION : MOINS robuste que Farey"
                bg_col = "#EDB6C3" 
                fg_col = "#7F2A41" 
            
        self.lbl_concl.config(text=msg, bg=bg_col, fg=fg_col)
        self.lbl_concl.pack_forget() 
        self.lbl_concl.pack(pady=20, fill="x", padx=20)
            
        # Force la mise à jour visuelle de la fenêtre
        self.root.update_idletasks()
        self.btn_run.config(state="normal", text="▶ LANCER L'ANALYSE")
        messagebox.showinfo("Succès", "Analyse terminée !")

if __name__ == "__main__":
    root = tk.Tk()
    app = FareyApp(root)
    root.mainloop()
