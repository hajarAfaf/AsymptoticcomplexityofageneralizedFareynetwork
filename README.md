# Network Complexity Analyzer

A scientific Python application for the structural analysis of complex networks (social, biological, communication systems).

This project is the practical software implementation of the research work presented in the paper:
> **"Asymptotic complexity of a generalized Farey network"**
> *By: Mokhlissi Raihana, Lotfi Dounia, Joyati Debnath, Soukayna Qarboua, El marraki Mohamed.*

This tool calculates **Complexity (Tau)** and **Structural Entropy (Rho)** using graph reduction methods inspired by resistive circuit theory, effectively bypassing classical algorithmic complexity limitations to analyze massive datasets (e.g., SNAP).

---

##  Theoretical Background & Graph Theory

### 1. The Core Metric: Spanning Trees ($\tau$)
The fundamental measure used in this analysis is the number of **Spanning Trees** ($\tau$).
In graph theory, a spanning tree is a subgraph that connects all vertices together without any cycles. The total number of these trees is a direct indicator of the network's **robustness**, **redundancy**, and **complexity**.

### 2. Scientific Basis (The Reference Paper)
According to the research by *Mokhlissi et al.*, the asymptotic entropy of generalized Farey networks converges to a specific constant. Our software allows researchers to verify if real-world networks (Twitter, Facebook, Power Grids) approach this **Farey Entropy** or diverge from it.

* **Unweighted Graphs (Topological):** Treated with unitary weight ($w=1$). Here, $\tau$ represents the exact number of topological configurations.
* **Weighted Graphs:** The algorithm supports edge weights (conductance), allowing the analysis of flow capability, bandwidth, or connection strength.

---

##  From Physics to Algorithms: The Electrical Analogy

Why do we use terms like "Conductance", "Series", or "Parallel" for abstract graphs?

### The Laplacian Isomorphism
There is a strict mathematical identity between calculating spanning trees (Kirchhoff's Matrix Tree Theorem) and calculating the equivalent conductance of an electrical resistor network.

* **The Classical Problem:** Computing $\tau$ usually requires finding the determinant of the **Laplacian Matrix**. This operation is $O(N^3)$. For a graph with $N > 10,000$ nodes, this is computationally impossible.
* **Our Solution:** Iterative Local Reduction.

The tool dynamically applies reduction transformations described in the methodology, preserving the "Total Conductance" (Complexity) at every step:

1.  **Star-Mesh Transformation (Node Elimination):** A generalization of the *Y-Δ (Wye-Delta)* transform. A central node is removed, and its neighbors form a weighted clique (mesh).
2.  **Series Edge Law:** $w_{eq} = \frac{w_1 \cdot w_2}{w_1 + w_2}$ (Harmonic mean).
3.  **Parallel Edge Law:** $w_{eq} = w_1 + w_2$ (Summation).

---

##  Technical Architecture

### Tech Stack
* **Python 3.x**: Core logic.
* **NetworkX**: Data structures (`MultiGraph`) for efficient node/edge manipulation.
* **Pandas**: High-performance parsing of `.edges` files (millions of rows).
* **Tkinter**: Native GUI with **Multithreading** to prevent application freezing during intensive calculations.

### Reduction Algorithm Priority
The engine uses a priority queue system to optimize convergence speed:
1.  **Cleanup (Parallel):** Immediate fusion of redundant edges.
2.  **Simplification (Series):** Fast elimination of path nodes (degree 2).
3.  **Restructuring (Star-Mesh):** Handling dense nodes (Hubs), prioritizing lower-degree nodes first to minimize edge explosion.

---

##  Installation & Usage

### Prerequisites
```bash
pip install networkx pandas
