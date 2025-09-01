# 🧬 tRNA-Sec Study Project 🔬

This repository provides a reproducible framework to explore the biology, structural features, and computational analysis of **tRNA-Sec** (the tRNA that incorporates the non-canonical amino acid selenocysteine) and related small RNAs.

---

## 🎯 Objectives

The main goals of this project are:

-   **Understand the biology**: Learn the essential features of selenocysteine incorporation and secondary/tertiary structure, including how non-canonical base pairs and nucleotide modifications influence conformation.
-   **Assemble and curate a dataset**: Fetch sequences from RNAcentral and annotate them with modification information from MODOMICS, consensus secondary structures from Rfam, gene context from GtRNAdb, and relevant literature.
-   **Perform bioinformatic analyses**: Align and build phylogenetic trees 🌳 of the collected sequences; analyze predicted secondary structures using tools such as ViennaRNA, NUPACK, and RNAstructure; and investigate the influence of non-canonical base pairs.
-   **Investigate three-dimensional structures**: Search the Protein Data Bank for existing tRNA-Sec structures or model them using coarse-grained or atomistic methods; describe characteristic structural motifs and discuss how modifications affect them.
-   **Apply machine learning (optional)**: Explore basic clustering and classification methods using simple descriptors derived from sequences, motifs, modifications, and predicted structures.
-   **Ensure reproducibility and communication**: Provide scripts and notebooks so that analyses can be reproduced end-to-end. Build an informative presentation 📊 summarizing the workflow and findings.

---

## 📁 Repository Structure

```bash
project/
│   README.md               # this file
│   .gitignore              # patterns for files/folders to exclude from version control
│
├── data/
│   ├── raw/                # 📥 downloaded datasets (e.g. RNACentral tRNA-Sec sequences modification table)
│   └── processed/          # 🧼 curated datasets with annotations and derived features
│
├── scripts/                # 💻 command-line scripts for data retrieval, cleaning and analysis
│
├── notebooks/              # 📝 Jupyter notebooks for exploratory analyses, alignments and visualisation
│
├── results/
│   ├── figures/            # 🖼️ figures generated during analyses
│   └── tables/             # 📈 summary tables and intermediate outputs
│
├── docs/
│   └── presentation_outline.md  # outline for the final slide deck
│
└── reports/                # 📄 longer reports or exported notebooks (e.g. PDF or markdown)
```

---

---

## 🛠️ Installation and Prerequisites

This project requires both **Python libraries** and some **external command-line tools**.  
We recommend Python **3.9+** and a virtual environment.

1.  **Clone this repository**:
    ```bash
    git clone <repository_url>
    cd project
    ```
2.  **Create a virtual environment and activate it**:
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```
3.  **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
---
##⚙️ External Tools (must be installed separately)

Some analyses require third-party bioinformatics software. Make sure these are installed and available in your $PATH:

1. Infernal (cmalign, esl-reformat) → structural alignments
2. FastTree → fast phylogenetic tree construction
3. ViennaRNA (RNAfold, RNAsubopt, RNAplot) → RNA secondary structure prediction
---

## 🚀 Usage

Jupyter notebooks in the `notebooks/` folder should be designed to run from the command line and accept relevant arguments. For instance:

-   `fetch_and_explore.ipynb`: Query RNAcentral for tRNA-Sec sequences and save them to `data/raw/`.
-   `align_sequences_prediction.ipynb`: Perform multiple sequence alignment with `Infernal` and phylogenetic tree construction. Run secondary structure predictions using `ViennaRNA`.
-   `trnasec_structural_exploration.ipynb`: Analyze secondary and tertiary structures of tRNA-Sec across domains of life, identify non-canonical motifs, and (optionally) run MD simulations.
- `rna_ml_exploration.ipynb`: Extract sequence/structural features, perform clustering of RNAs, and train baseline ML classifiers for exploratory prediction tasks.

- 🚀 Usage

Main notebooks and their purposes:

-   `fetch_and_explore.ipynb`:  ✅ → Query RNAcentral for tRNA-Sec sequences, store in data/raw/.
-   `align_sequences_prediction.ipynb`: ⏳ → Perform alignments with Infernal, build phylogenetic trees, run secondary structure predictions.
-   `trnasec_structural_exploration.ipynb`: ✅ → Analyze secondary and tertiary structures of tRNA-Sec across domains of life.
-   `rna_ml_exploration.ipynb`: 🚧 → Extract features, perform clustering, and train baseline ML classifiers. (Currently functional for basic classification; more complex features planned.)
---
📌 Current Status

-  ✅ Most notebooks are functional and reproducible.
-  🛠️ Some steps are in progress (especially advanced structural analysis and ML feature engineering).
-  📈 All notebooks can already be run end-to-end, though results will improve as new features are added.

---

## 🙏 Acknowledgements

This repository draws on open-source software and publicly available databases.  

For references to tools and libraries used in this project, see [docs/bibliography.md](docs/bibliography.md).

For a biological overview of tRNA-Sec, see [docs/trnasec_overview.md](docs/trnasec_overview.md).


For questions, feedback, or collaboration opportunities, please contact:  

📧 **Rodrigo Machado Birollo** — rodrigomachadobirollo@gmail.com


I hope you enjoy reading this!




