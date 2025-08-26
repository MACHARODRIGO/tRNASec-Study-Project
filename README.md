# 🧬 tRNA-Sec Study Project 🔬

This repository provides a reproducible framework to explore the biology, structural features, and computational analysis of **tRNA-Sec** (the tRNA that incorporates the non-canonical amino acid selenocysteine) and related small RNAs. It's organized around specific hiring assessment goals and is intended as a starting point for further work.

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

This project is designed to be executed in a **Python 3.9+** environment. To install the required packages:

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
    Recommended dependencies include: `numpy`, `pandas` for data handling; `biopython` for sequence parsing and alignments; `scikit-learn` for basic machine learning tasks; `matplotlib`, `seaborn` for plotting. Command-line tools like **ViennaRNA**, **NUPACK**, and **RNAstructure** are also required.

---

## 🚀 Usage

Scripts in the `scripts/` folder should be designed to run from the command line and accept relevant arguments. For instance:

-   `fetch_rnacentral.py`: Query RNAcentral for tRNA-Sec sequences and save them to `data/raw/`.
-   `annotate_modifications.py`: Integrate modification data from MODOMICS and other sources.
-   `align_sequences.py`: Perform multiple sequence alignment and phylogenetic tree construction.
-   `predict_structures.py`: Run secondary structure predictions using ViennaRNA, NUPACK, or RNAstructure.

Jupyter notebooks in the `notebooks/` directory can be used for exploratory analyses, visualization, and documentation of intermediate findings.

---

## ✅ Reproducibility

For reproducible results, record the versions of all external tools and databases used. Wherever possible, run scripts with deterministic random seeds and provide the input data and output in the repository. Document the workflow clearly in notebooks and commit changes regularly.

---

## ⚖️ License

This project is released under the **MIT License**. See the `LICENSE` file for details.

---

## 🙏 Acknowledgements

This repository draws on open-source software and publicly available databases. Refer to the `reports/bibliography.md` for citations of the literature and resources used.






