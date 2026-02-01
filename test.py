import networkx as nx
import pandas as pd
import math
import os

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
    def __init__(self):
        # Utilisation de MultiGraph pour supporter les arêtes parallèles temporaires
        self.G = nx.MultiGraph()
        self.log_tau = 0.0  # On travaille en Log pour éviter l'overflow
        self.initial_node_count = 0
        self.history = []

    def load_snap_data(self, filepath):
        """
        Charge un fichier de données SNAP (.edges).
        Initialise tous les poids (conductances) à 1.0 comme spécifié dans la méthodologie.
        """
        print(f"Chargement : {filepath} ...")
        try:
            df = pd.read_csv(filepath, sep='\s+', header=None, names=['u', 'v'])
            self.G.clear()
            # Poids initial = 1.0 partout
            edges = [(row['u'], row['v'], 1.0) for _, row in df.iterrows()]
            self.G.add_weighted_edges_from(edges)
            
            self.initial_node_count = self.G.number_of_nodes()
            self.log_tau = 0.0 
            print(f"✅ Graphe chargé : {self.G.number_of_nodes()} nœuds, {self.G.number_of_edges()} arêtes.")
        except Exception as e:
            print(f"❌ Erreur : {e}")

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
        if has_changed: self.history.append("Réduction Parallèle")
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
        if has_changed: self.history.append("Réduction Série")
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
                        
                        self.history.append(f"Wye-Delta (Node {node})")
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
        
        neighbors = list(self.G.neighbors(node))
        # On s'assure que les voisins sont uniques
        neighbors = list(set(neighbors))
        degree = len(neighbors)
        
        try:
            # 1. Calcul de la somme S (Algorithm 5, ligne 9)
            # Récupération des poids de toutes les arêtes connectées au nœud central
            weights = {}
            total_s = 0.0
            for neighbor in neighbors:
                # On prend la première arête dispo (normalement parallel edge a déjà nettoyé)
                w = self.G[node][neighbor][0]['weight']
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
            
            self.history.append(f"Star-Mesh (Node {node}, Degré {degree})")
            return True
            
        except Exception as e:
            print(f"Erreur Star-Mesh sur nœud {node}: {e}")
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
        print("--- Début de l'analyse par Transformations ---")
        iteration = 0
        while True:
            iteration += 1
            if iteration % 10 == 0:
                print(f"Iter {iteration} | Nœuds restants : {self.G.number_of_nodes()} | Arêtes : {self.G.number_of_edges()}")

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
            print("Aucune transformation possible.")
            break
            
        print("--- Fin de la réduction ---")
        
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

# --- MAIN ---
if __name__ == "__main__":
    analyzer = FareyMethodAnalyzer()
    
    # FICHIER A ANALYSER
    mon_fichier = "24117694.edges"  
    
    # Création fichier test si absent
    if not os.path.exists(mon_fichier):
        with open(mon_fichier, "w") as f: f.write("1 2\n2 3\n3 1\n1 4")

    analyzer.load_snap_data(mon_fichier)
    
    # Calcul
    log_tau = analyzer.run_analysis()
    rho = analyzer.calculate_entropy()
    nombre_arbres_lisible = format_huge_number(log_tau)
    print("\n" + "="*40)
    print("      RÉSULTATS FINAUX")
    print("="*40)
    print(f"Log(Tau) : {log_tau:.4f}")
    print(f"Nombre d'arbres est.  : {nombre_arbres_lisible}")
    print(f"Entropie (Rho) : {rho:.4f}")
    print(f"Farey (Ref)    : 0.9457")
    
    if rho > 0.9457:
        print("CONCLUSION : Ton réseau est PLUS robuste que Farey.")
    else:
        print("CONCLUSION : Ton réseau est MOINS robuste que Farey.")
