# Network Complexity Analyzer 

Une application scientifique en Python pour l'analyse structurelle de réseaux complexes (sociaux, biologiques, communication).

Cet outil calcule la **complexité (Tau)** et l'**entropie structurelle (Rho)** en utilisant des méthodes de réduction de graphes inspirées de la théorie des circuits résistifs, permettant de contourner les limitations de complexité algorithmique classiques.


##  Contexte Théorique & Graph Theory

### 1. Le Problème : Compter les Arbres Couvrants
La mesure fondamentale utilisée ici est le nombre d'**Arbres Couvrants (Spanning Trees, $\tau$)**.
En théorie des graphes, un arbre couvrant est un sous-graphe qui connecte tous les sommets sans former de cycles. Le nombre total de ces arbres est un indicateur direct de la **robustesse** et de la **complexité** du réseau (plus il y a d'arbres, plus le réseau est redondant et résistant aux pannes).

### 2. Graphes Pondérés vs Non-Pondérés
L'algorithme a été conçu pour être **universel** :

* **Graphes Non-Pondérés (Topologiques)** :
    * *Cas d'usage :* Réseaux sociaux (SNAP), interactions biologiques.
    * *Traitement :* Chaque arête est considérée comme ayant un poids (conductance) unitaire $w = 1.0$.
    * *Résultat :* $\tau$ représente le nombre exact de configurations topologiques possibles.

* **Graphes Pondérés (Weighted Graphs)** :
    * *Cas d'usage :* Réseaux routiers (distance), Internet (bande passante), réseaux neuronaux (force synaptique).
    * *Traitement :* L'algorithme lit le poids $w$ fourni dans le fichier.
    * *Résultat :* $\tau$ devient une "Complexité Pondérée", reflétant non seulement la connectivité mais aussi la "facilité" de passage à travers le réseau.

---

##  De la Physique à l'Algorithme : L'Analogie Électrique

Pourquoi utilisons-nous des termes comme "Conductance", "Série" ou "Parallèle" pour des graphes abstraits ?

### L'Isomorphisme Laplacien
Il existe une identité mathématique stricte entre :
1.  Le calcul des arbres couvrants d'un graphe (Théorème Arbre-Matrice de Kirchhoff).
2.  Le calcul de la conductance équivalente d'un réseau de résistances électriques.

Classiquement, pour calculer $\tau$, il faut calculer le déterminant de la **Matrice Laplacienne** du graphe.
* **Problème :** Le calcul d'un déterminant est en $O(N^3)$. Pour un graphe de 100 000 nœuds (Twitter/Facebook), c'est impossible (la matrice ne tient même pas en mémoire RAM).

### Notre Approche : La Réduction Itérative
Au lieu d'attaquer la matrice entière, nous utilisons les **transformations locales**. Si on remplace les arêtes par des résistances électriques, on peut simplifier le graphe petit à petit sans changer sa "conductance totale" (qui est proportionnelle à $\tau$).

L'outil applique dynamiquement ces lois de conservation :

1.  **Loi des Nœuds (Transformation Étoile-Maillage / Star-Mesh)** : Un nœud central peut être supprimé si l'on reconnecte tous ses voisins entre eux en redistribuant les poids. C'est la généralisation de la transformation *Y-Δ (Wye-Delta)*.
2.  **Loi des Conductances en Série** : $w_{eq} = (w_1 \cdot w_2) / (w_1 + w_2)$.
3.  **Loi des Conductances en Parallèle** : $w_{eq} = w_1 + w_2$.

Cette méthode permet de réduire des graphes massifs que l'approche matricielle classique ne pourrait jamais traiter.

---

##  Architecture Technique

### Stack Technologique
* **Python 3.11** : Langage principal.
* **NetworkX** : Structure de données (MultiGraph) pour manipuler les nœuds et les arêtes.
* **Pandas** : Parsing haute performance des fichiers `.edges` (millions de lignes).
* **Tkinter** : Interface graphique (GUI) native et réactive (Threadée).

### Algorithme de Réduction
Le moteur d'analyse fonctionne par priorités pour optimiser la vitesse de convergence :
1.  **Nettoyage (Parallel)** : Fusion immédiate des arêtes redondantes.
2.  **Simplification Rapide (Series)** : Élimination des nœuds de passage (degré 2).
3.  **Restructuration (Star-Mesh)** : Attaque des nœuds denses (Hubs) en commençant par les degrés les plus faibles pour limiter l'explosion combinatoire des arêtes.

---

## Installation et Utilisation

### Prérequis
```bash
pip install networkx pandas
