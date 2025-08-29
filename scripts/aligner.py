"""
aligner.py - Utility functions for tRNA-Sec alignment and phylogenetic analysis

This module provides:
1. Convert RNAcentral CSV datasets to FASTA.
2. Run alignments using Infernal (cmalign + esl-reformat).
3. Build phylogenetic trees with FastTree.
4. Visualize phylogenetic trees (radial layout with ete3).
5. Interpret tree structure into Markdown reports.
6. Generate sequence logos from aligned FASTA files.
7. Predict and plot RNA secondary structures with ViennaRNA.

⚠️ Dependencies: Biopython, pandas, ete3, logomaker, matplotlib, ViennaRNA, cairosvg.
Also requires external tools (cmalign, esl-reformat, FastTree).
"""

import os
from collections import Counter
from Bio import Phylo, SeqIO


# =====================================================
# 1. Convert CSV to FASTA
# =====================================================
def convert_csv_to_fasta(csv_path, fasta_path, mapping_file="../results/tables/id_mapping.csv"):
    """
    Convert RNAcentral CSV to FASTA format.
    - Extracts a taxonomic name from description when possible.
    - Cleans species names for FastTree compatibility.
    - Stores a mapping CSV for later reference.
    """
    import pandas as pd
    import re

    try:
        df = pd.read_csv(csv_path)
        mapping = []

        with open(fasta_path, "w") as fasta_file:
            for _, row in df.iterrows():
                urs_id = str(row['URS_ID']).strip()
                description = str(row['description']).strip()

                # Extract binomial name or "Candidatus" if possible
                species_name = None
                match_binomial = re.search(r'([A-Z][a-z]+ [a-z]+)', description)
                if match_binomial:
                    species_name = match_binomial.group(1)

                match_candidatus = re.search(r'([A-Z][a-z]+ [A-Z][a-z]+(?: [a-z]+)?)', description)
                if match_candidatus:
                    species_name = match_candidatus.group(1)

                if not species_name:
                    species_name = description

                safe_species = re.sub(r"[^A-Za-z0-9_.-]", "_", species_name)
                header = f">{urs_id}|{safe_species}|{description}"
                sequence = str(row['Sequence']).strip()
                fasta_file.write(f"{header}\n{sequence}\n")

                mapping.append([urs_id, safe_species, description])

        os.makedirs(os.path.dirname(mapping_file), exist_ok=True)
        import pandas as pd
        pd.DataFrame(mapping, columns=["URS_ID", "safe_species", "description"]).to_csv(mapping_file, index=False)

        print(f"✅ FASTA file created: {fasta_path}")
        print(f"✅ Mapping saved to: {mapping_file}")

    except Exception as e:
        print(f"❌ Error while converting CSV to FASTA: {e}")


# =====================================================
# 2. Alignment pipeline (Infernal)
# =====================================================
def run_cmalign(fasta_wsl, sto_windows, cmalign_bin, cm_model):
    """Run cmalign (Infernal) through WSL."""
    import subprocess
    print("⏳ Running cmalign...")
    try:
        with open(sto_windows, "w") as sto_out:
            subprocess.run(["wsl", cmalign_bin, cm_model, fasta_wsl],
                           stdout=sto_out, check=True)
        print(f"✅ Stockholm file saved: {sto_windows}")
    except subprocess.CalledProcessError as e:
        print("❌ Error running cmalign:", e)


def run_esl_reformat(sto_wsl, aligned_windows, esl_reformat_bin):
    """Run esl-reformat to convert Stockholm to aligned FASTA."""
    import subprocess
    print("⏳ Running esl-reformat...")
    try:
        with open(aligned_windows, "w") as fasta_out:
            subprocess.run(["wsl", esl_reformat_bin, "afa", sto_wsl],
                           stdout=fasta_out, check=True)
        print(f"✅ Aligned FASTA file saved: {aligned_windows}")
    except subprocess.CalledProcessError as e:
        print("❌ Error running esl-reformat:", e)


def preview_alignment(aligned_windows, n=20):
    """Preview the first lines of the alignment."""
    if os.path.exists(aligned_windows):
        print("\n📖 Alignment preview:")
        with open(aligned_windows, "r") as f:
            for line in f.readlines()[:n]:
                print(line.strip())
    else:
        print("⚠️ Aligned file not found.")


def run_alignment_pipeline(fasta_wsl, sto_windows, sto_wsl, aligned_windows,
                           cmalign_bin, cm_model, esl_reformat_bin, preview_lines=20):
    """Full alignment pipeline: cmalign + esl-reformat + preview."""
    run_cmalign(fasta_wsl, sto_windows, cmalign_bin, cm_model)
    run_esl_reformat(sto_wsl, aligned_windows, esl_reformat_bin)
    preview_alignment(aligned_windows, n=preview_lines)


# =====================================================
# 3. Build Phylogenetic Tree with FastTree
# =====================================================
def build_phylogenetic_tree(input_alignment, output_tree, fasttree_exe):
    """Build phylogenetic tree with FastTree (nucleotide mode)."""
    import subprocess
    if not os.path.exists(input_alignment):
        print(f"❌ Alignment file not found: {input_alignment}")
        return

    print("⏳ Building phylogenetic tree with FastTree...")
    try:
        with open(output_tree, "w") as out_tree:
            subprocess.run([fasttree_exe, "-nt", input_alignment],
                           stdout=out_tree, check=True)
        print(f"✅ Tree saved: {output_tree}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running FastTree: {e}")


# =====================================================
# 4. Radial Tree Visualization
# =====================================================
def visualize_tree_radial(tree_file, out_image="../results/figures/tree_radial.png"):
    """Render radial phylogenetic tree using ete3."""
    from ete3 import Tree, TreeStyle
    from IPython.display import Image, display

    if not os.path.exists(tree_file):
        print(f"❌ Tree file not found: {tree_file}")
        return

    print("📊 Rendering radial phylogenetic tree...")
    os.makedirs(os.path.dirname(out_image), exist_ok=True)

    t = Tree(tree_file)
    ts = TreeStyle()
    ts.mode = "c"  # circular
    ts.show_leaf_name = True
    ts.show_branch_support = True
    ts.scale = 50

    t.render(out_image, tree_style=ts, w=2000, units="px")
    print(f"✅ Radial tree image saved: {out_image}")
    display(Image(filename=out_image))


# =====================================================
# 5. Tree Interpretation (Markdown report)
# =====================================================
def find_maximal_clades(clades):
    """Filter clades to keep only the maximals (remove subsets)."""
    maximal = []
    for clade in clades:
        contained = False
        for other in clades:
            if clade is not other and set(clade.get_terminals()).issubset(set(other.get_terminals())):
                contained = True
                break
        if not contained:
            maximal.append(clade)
    return maximal


def interpret_tree_markdown(tree_file, min_clade_size=4, support_cutoff=0.7, report_file="../results/tables/tree_report.md"):
    """Interpret phylogenetic tree and generate a Markdown report."""
    if not os.path.exists(tree_file):
        print(f"❌ Tree file not found: {tree_file}")
        return

    tree = Phylo.read(tree_file, "newick")
    leaves = tree.get_terminals()
    n_total = len(leaves)

    report = []
    report.append("# 📖 Automatic interpretation of the phylogenetic tree\n")
    report.append(f"- Total number of sequences: **{n_total}**\n")

    # --- Global species count
    species = [str(leaf).split("|")[1].strip() for leaf in leaves if "|" in str(leaf)]
    sp_counts = Counter(species)
    report.append("## 🧾 Most frequent species\n")
    for sp, c in sp_counts.most_common(10):
        report.append(f"- {sp}: {c} sequences")

    os.makedirs(os.path.dirname(report_file), exist_ok=True)
    with open(report_file, "w", encoding="utf-8") as f:
        f.write("\n".join(report))

    print(f"✅ Report saved: {report_file}")


# =====================================================
# 6. Generate Sequence Logo
# =====================================================
def generate_sequence_logo(alignment_fasta, out_logo="../results/figures/tRNA_logo.png", gap_threshold=0.9):
    """
    Generate a sequence logo from a multiple alignment.
    - alignment_fasta: aligned FASTA file
    - out_logo: PNG output path
    - gap_threshold: discard columns with >90% gaps
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    import logomaker

    if not os.path.exists(alignment_fasta):
        print(f"❌ Alignment file not found: {alignment_fasta}")
        return

    seqs = [str(record.seq) for record in SeqIO.parse(alignment_fasta, "fasta")]
    df = pd.DataFrame([list(seq) for seq in seqs])

    # Filter gap-rich columns
    gap_fraction = (df == "-").sum() / len(df)
    df = df.loc[:, gap_fraction < gap_threshold]
    df.replace("-", pd.NA, inplace=True)

    # Frequency matrix
    freq_matrix = pd.DataFrame(
        {i: df[i].value_counts(normalize=True) for i in df.columns}
    ).fillna(0).T

    os.makedirs(os.path.dirname(out_logo), exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(10, freq_matrix.shape[0] * 0.2), 6))
    logo = logomaker.Logo(freq_matrix, ax=ax, shade_below=.5, fade_below=.5, color_scheme="classic")
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90)
    logo.ax.set_ylabel('Frequency')

    plt.savefig(out_logo, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ Sequence logo saved: {out_logo}")


# =====================================================
# 7. RNA Secondary Structures by Species (ViennaRNA)
# =====================================================
def plot_rna_by_species(fasta_file, out_dir="../results/figures/rna_plots_by_species", max_per_fig=12):
    """
    Generate RNA secondary structure plots by species using ViennaRNA.
    Each subplot shows URS ID + minimum free energy (MFE).
    """
    import re
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import ViennaRNA as RNA

    os.makedirs(out_dir, exist_ok=True)

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    species_groups = {}
    for record in sequences:
        parts = record.id.split('|')
        description = " ".join(parts[1:]) if len(parts) > 1 else record.id
        match_binomial = re.search(r"([A-Z][a-z]+ [a-z]+)", description)
        species_name = match_binomial.group(1) if match_binomial else (parts[1] if len(parts) > 1 else record.id)

        safe_species = species_name.replace(" ", "_")
        species_groups.setdefault(safe_species, {"name": species_name, "records": []})
        species_groups[safe_species]["records"].append(record)

    print(f"📥 {len(sequences)} sequences loaded, grouped into {len(species_groups)} species.")

    for species_id, info in species_groups.items():
        species_name = info["name"]
        group = info["records"]
        print(f"\nGenerating figures for species: {species_name} ({len(group)} sequences)")

        for i in range(0, len(group), max_per_fig):
            block = group[i:i+max_per_fig]
            fig, axes = plt.subplots(3, 4, figsize=(18, 14))
            axes = axes.flatten()

            for j, record in enumerate(block):
                seq = str(record.seq).replace("T", "U")
                if not seq:
                    continue

                ss, mfe = RNA.fold(seq)
                safe_id = record.id.replace('|', '_').replace(" ", "_")
                svg_file = os.path.join(out_dir, f"{safe_id}.svg")
                RNA.svg_rna_plot(seq, ss, svg_file)

                png_file = svg_file.replace(".svg", ".png")
                try:
                    import cairosvg
                    cairosvg.svg2png(url=svg_file, write_to=png_file)
                except ImportError:
                    png_file = None

                urs_id = record.id.split('|')[0]
                ax = axes[j]
                if png_file and os.path.exists(png_file):
                    img = mpimg.imread(png_file)
                    ax.imshow(img)

                ax.set_title(f"{urs_id}\nMFE={mfe:.2f}", fontsize=9)
                ax.axis("off")

            for k in range(len(block), len(axes)):
                axes[k].axis("off")

            fig.suptitle(f"RNA structures — {species_name}", fontsize=16)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            output_fig_path = os.path.join(out_dir, f"rna_structures_{species_id}_{i//max_per_fig+1}.png")
            fig.savefig(output_fig_path, dpi=300)
            plt.close(fig)

            print(f"✅ Figure saved: {output_fig_path}")
