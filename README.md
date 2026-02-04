# Network Complexity Analyzer 

A scientific Python application for the structural analysis of complex networks (social, biological, communication systems).

This project is the practical software implementation of the research work presented in the paper:
> **"Asymptotic complexity of a generalized Farey network"**
> *By: Mokhlissi Raihana, Lotfi Dounia, Joyati Debnath, Soukayna Qarboua, El marraki Mohamed.*

This tool calculates **Complexity (Tau)** and **Structural Entropy (Rho)** using graph reduction methods inspired by resistive circuit theory, effectively bypassing classical algorithmic complexity limitations to analyze massive datasets.

---

##  Theoretical Background & Graph Theory

### 1. The Core Metric: Spanning Trees ($\tau$)
The fundamental measure used in this analysis is the number of **Spanning Trees** ($\tau$).
In graph theory, a spanning tree is a subgraph that connects all vertices together without any cycles. The total number of these trees is a direct indicator of the network's **robustness**, **redundancy**, and **complexity**.

### 2. Scientific Basis (The Reference Paper)
According to the research by *Mokhlissi et al.*, the asymptotic entropy of generalized Farey networks converges to a specific constant. Our software allows researchers to verify if real-world networks (Twitter, Facebook, Power Grids) approach this **Farey Entropy** or diverge from it.

* **Unweighted Graphs (Topological):** Treated with unitary weight ($w=1$). $\tau$ represents the exact number of topological configurations.
* **Weighted Graphs:** The algorithm supports edge weights (conductance), allowing the analysis of flow capability, bandwidth, or connection strength.

---

##  From Physics to Algorithms: The 5 Transformations

Why do we use terms like "Conductance", "Series", or "Delta-Wye" for abstract graphs?

### The Laplacian Isomorphism
There is a strict mathematical identity between calculating spanning trees (Kirchhoff's Matrix Tree Theorem) and calculating the equivalent conductance of an electrical resistor network.
* **The Classical Problem:** Computing $\tau$ via the Laplacian Matrix determinant is $O(N^3)$.
* **Our Solution:** Iterative Local Reduction using 5 specific transformations.

The tool applies these laws dynamically to reduce the graph while preserving the total complexity:

1.  **Parallel Edge Transformation:**
    * *Action:* Merges multiple edges between the same two nodes.
    * *Formula:* $w_{eq} = w_1 + w_2$ (Summation of conductances).

2.  **Serial Edge Transformation:**
    * *Action:* Eliminates a node of degree 2 (a simple path).
    * *Formula:* $w_{eq} = \frac{w_1 \cdot w_2}{w_1 + w_2}$ (Harmonic mean).

3.  **Wye-Delta Transformation ($Y \to \Delta$):**
    * *Action:* Eliminates a central node of degree 3 (Star) and creates a triangle (Delta) among its neighbors.
    * *Purpose:* Reduces the node count.

4.  **Delta-Wye Transformation ($\Delta \to Y$):**
    * *Action:* Transforms a triangle loop (Delta) into a central node structure (Star).
    * *Purpose:* Used to break complex loops or create Series/Parallel opportunities that didn't exist before.

5.  **Star-Mesh Transformation:**
    * *Action:* The generalization of Wye-Delta for nodes of degree $k > 3$.
    * *Method:* The central node is removed, and a weighted clique (mesh) is created connecting all its neighbors.

---

##  Technical Architecture

### Tech Stack
* **Python 3.11**: Core logic.
* **NetworkX**: Data structures (`MultiGraph`) for efficient node/edge manipulation.
* **Pandas**: High-performance parsing of `.edges` files.
* **Tkinter**: Native GUI with **Multithreading** to prevent application freezing during intensive calculations.

### Reduction Strategy
The engine uses a priority queue system to optimize convergence speed:
1.  **Cleanup:** Parallel edges are merged immediately.
2.  **Simplification:** Serial edges are processed to reduce path lengths.
3.  **Topology Change:** $Y-\Delta$ and $\Delta-Y$ are applied to reshape local structures.
4.  **Heavy Reduction:** Star-Mesh is used as a last resort for dense hubs (high-degree nodes).

---

## Installation & Usage

### Prerequisites
```bash
pip install networkx pandas
