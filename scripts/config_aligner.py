"""
config.py - Local paths and executables for aligner functions.
⚠️ Adapt these paths to your system.
"""

# ===============================
# Input files
# ===============================
input_file = "../data/raw/200_trna_sec.csv"

# ===============================
# Processed data outputs
# ===============================
fasta_default = "../data/processed/tRNA_sequences.fasta"
sto_default = "../data/processed/tRNA_sequences.sto"
aligned_default = "../data/processed/tRNA_alignment.fasta"

# WSL paths (only used when running in WSL)
fasta_wsl = "/mnt/c/Users/ro-ma/Documents/GitHub/tRNASec-Study-Project/data/processed/tRNA_sequences.fasta"
sto_wsl = "/mnt/c/Users/ro-ma/Documents/GitHub/tRNASec-Study-Project/data/processed/tRNA_sequences.sto"
aligned_wsl = "/mnt/c/Users/ro-ma/Documents/GitHub/tRNASec-Study-Project/data/processed/tRNA_alignment.fasta"

# ===============================
# Analysis results (relative paths inside repo)
# ===============================
tree_default        = "../results/tables/tRNA_phylogenetic_tree.nhx"
tree_report_default = "../results/tables/tree_report.md"
tree_image_default  = "../results/figures/tree_radial.png"
logo_default        = "../results/figures/tRNA_logo.png"
rna_plots_dir       = "../results/figures/rna_plots_by_species"

# ===============================
# External executables (real paths)
# ===============================
cmalign_bin      = "/usr/bin/cmalign"
esl_reformat_bin = "/home/ro-ma/infernal_models/infernal/easel/miniapps/esl-reformat"
cm_model         = "/home/ro-ma/infernal_models/infernal/RF01852.cm"
fasttree_exe     = r"C:\Program Files\FastTree\FastTree.exe"



