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
from IPython.display import display, Markdown

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
def run_cmalign(cmalign_bin, cm_model, fasta_wsl, sto_wsl, timeout=600):
    """Run cmalign correctly for unaligned sequences"""
    import subprocess, time
    
    print("⏳ Running cmalign for unaligned sequences...")
    print("⚠️ This may take several minutes depending on sequence count")
    
    start_time = time.time()
    try:
        # CORRECT ORDER: cmalign -o output.cm cm_model.cm input.fasta
        result = subprocess.run(
            ["wsl", cmalign_bin, "-o", sto_wsl, cm_model, fasta_wsl],
            capture_output=True,
            text=True,
            timeout=timeout
        )
        
        elapsed = time.time() - start_time
        print(f"⏱️ Execution time: {elapsed:.2f} seconds")
        print(f"Exit code: {result.returncode}")
        
        if result.stderr:
            print(f"stderr (first 300 chars):\n{result.stderr[:300]}...\n")
        
        if result.returncode == 0:
            print(f"✅ Success! Stockholm file generated: {sto_wsl}")
        else:
            print("❌ cmalign failed")
    
    except subprocess.TimeoutExpired:
        print(f"❌ cmalign timeout - {timeout/60} minutes was not enough")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")


def run_esl_reformat(sto_wsl, aligned_wsl, esl_reformat_bin):
    """Alternative method using stdout redirection"""
    import subprocess
    print("🔄 Trying alternative esl-reformat method...")
    
    try:
        # Método alternativo: usar stdout en lugar de -o
        result = subprocess.run(
            ["wsl", esl_reformat_bin, "afa", sto_wsl],
            capture_output=True,
            text=True,
            check=True,
            timeout=120
        )
        
        if result.returncode == 0:
            # Guardar manualmente el output
            with open(aligned_wsl.replace('/mnt/c/', 'C:/').replace('/', '\\'), 'w') as f:
                f.write(result.stdout)
            print(f"✅ Aligned FASTA generated (alternative method): {aligned_wsl}")
        else:
            print(f"❌ Alternative method failed: {result.stderr[:200]}")
            
    except Exception as e:
        print(f"❌ Alternative method error: {e}")


def preview_alignment(aligned_windows, n=20):
    """Muestra las primeras líneas del alineamiento"""
    import os
    if os.path.exists(aligned_windows):
        print("\n📖 Vista previa del alineamiento:")
        with open(aligned_windows, "r") as f:
            for line in f.readlines()[:n]:
                print(line.strip())
    else:
        print("⚠️ No existe el archivo de alineamiento.")


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
# 5. Tree Interpretation (Markdown report for Notebook)
# =====================================================
import os
from collections import Counter
from Bio import Phylo
from IPython.display import display, Markdown

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
    """Interpret phylogenetic tree and generate a Markdown report.
       Shows it nicely in Jupyter Notebook and saves it to disk.
    """
    if not os.path.exists(tree_file):
        print(f"❌ Tree file not found: {tree_file}")
        return

    tree = Phylo.read(tree_file, "newick")
    leaves = tree.get_terminals()
    n_total = len(leaves)

    report = []
    report.append("# 📖 Automatic interpretation of the phylogenetic tree\n")
    report.append(f"- **Total sequences in tree:** {n_total}\n")

    # --- Global species count
    species = [str(leaf).split("|")[1].strip() for leaf in leaves if "|" in str(leaf)]
    sp_counts = Counter(species)

    report.append("## 🧾 Most frequent species\n")
    report.append("| Rank | Species | # Sequences |")
    report.append("|------|---------|-------------|")
    for i, (sp, c) in enumerate(sp_counts.most_common(10), start=1):
        report.append(f"| {i} | {sp} | {c} |")

    # --- Largest clades
    clades = find_maximal_clades(tree.get_nonterminals())
    big_clades = [c for c in clades if len(c.get_terminals()) >= min_clade_size]

    if big_clades:
        report.append("\n## 🌳 Major clades detected\n")
        report.append(f"Clades with at least **{min_clade_size} sequences** and support > {support_cutoff}:\n")
        for i, cl in enumerate(big_clades, start=1):
            sp_list = [str(leaf).split('|')[1] if '|' in str(leaf) else str(leaf) for leaf in cl.get_terminals()]
            sp_summary = Counter(sp_list).most_common(3)
            support = getattr(cl, "confidence", None)
            report.append(f"- **Clade {i}:** {len(cl.get_terminals())} seqs | Support: {support if support else 'n/a'}")
            report.append(f"  - Main species: {', '.join([f'{s} ({c})' for s, c in sp_summary])}")

    # Save to file
    os.makedirs(os.path.dirname(report_file), exist_ok=True)
    with open(report_file, "w", encoding="utf-8") as f:
        f.write("\n".join(report))

    # 👇 Show pretty output inside Jupyter
    display(Markdown("\n".join(report)))

    print(f"✅ Report saved: {report_file}")
    
# =====================================================
# 6. Generate Sequence Logo
# =====================================================
def generate_sequence_logo(
    alignment_fasta,
    out_logo="../results/figures/tRNA_logo.png",
    gap_threshold=0.9,
    scale_factor=0.35  # 👈 controls how wide the figure will be
):
    """
    Generate a sequence logo from a multiple alignment.
    
    Parameters
    ----------
    alignment_fasta : str
        Path to the aligned FASTA file.
    out_logo : str
        Output path for the PNG logo.
    gap_threshold : float
        Discard alignment columns with more than this fraction of gaps (default = 0.9).
    scale_factor : float
        Controls figure width scaling (default = 0.35, increase if letters look cramped).
    """
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import logomaker
    from Bio import SeqIO
    from IPython.display import Image, display

    # --- Check input
    if not os.path.exists(alignment_fasta):
        print(f"❌ Alignment file not found: {alignment_fasta}")
        return

    # --- Load sequences
    seqs = [str(record.seq) for record in SeqIO.parse(alignment_fasta, "fasta")]
    df = pd.DataFrame([list(seq) for seq in seqs])

    # --- Filter gap-rich columns
    gap_fraction = (df == "-").sum() / len(df)
    df = df.loc[:, gap_fraction < gap_threshold]
    df.replace("-", pd.NA, inplace=True)

    # --- Build frequency matrix
    freq_matrix = pd.DataFrame(
        {i: df[i].value_counts(normalize=True) for i in df.columns}
    ).fillna(0).T

    # --- Ensure output directory exists
    os.makedirs(os.path.dirname(out_logo), exist_ok=True)

    # --- Plot logo
    fig_width = max(12, freq_matrix.shape[0] * scale_factor)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    logo = logomaker.Logo(
        freq_matrix,
        ax=ax,
        color_scheme="classic",
        width=1.0,     # 👈 width of each column (increase for more separation)
        vpad=0.1,
        show_spines=False
    )

    # --- Clean style
    logo.ax.set_ylabel('Frequency')
    logo.ax.set_xticks(range(0, freq_matrix.shape[0], 10))   # tick marks every 10 positions
    logo.ax.set_xticklabels(range(0, freq_matrix.shape[0], 10))
    logo.ax.grid(False)

    # --- Save and close
    plt.savefig(out_logo, dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

    # --- Show inline in Jupyter
    display(Image(filename=out_logo))
    print(f"✅ Sequence logo saved and displayed: {out_logo}")


# =====================================================
# 7. RNA Secondary Structures by Species (ViennaRNA)
# =====================================================
def plot_rna_by_species(
    fasta_file,
    out_dir="../results/figures/rna_plots_by_species",
    max_per_fig=12,
    total_figures=None   # 👈 global hard limit across ALL species
):
    """
    Generate RNA secondary structure plots by species using ViennaRNA.

    Parameters
    ----------
    fasta_file : str
        Path to input FASTA file.
    out_dir : str
        Directory to save figures.
    max_per_fig : int
        Maximum sequences per figure (default = 12).
    total_figures : int or None
        If set, stops after generating this many figures in total across all species.
        If None (default), generate all.
    """
    import os
    import re
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from Bio import SeqIO
    import ViennaRNA as RNA
    from IPython.display import display, Image

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

    global_count = 0

    for species_id, info in species_groups.items():
        species_name = info["name"]
        group = info["records"]
        print(f"\nGenerating figures for species: {species_name} ({len(group)} sequences)")

        for i in range(0, len(group), max_per_fig):
            # --- HARD STOP condition
            if total_figures is not None and global_count >= total_figures:
                print(f"⏹️ Global limit of {total_figures} figures reached. Stopping.")
                return

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

            # Hide empty subplots
            for k in range(len(block), len(axes)):
                axes[k].axis("off")

            fig.suptitle(f"RNA structures — {species_name}", fontsize=16)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            output_fig_path = os.path.join(out_dir, f"rna_structures_{species_id}_{i//max_per_fig+1}.png")
            fig.savefig(output_fig_path, dpi=300)
            plt.close(fig)

            print(f"✅ Figure saved: {output_fig_path}")

            # Show preview only for the first figure
            if global_count == 0:
                display(Image(filename=output_fig_path))

            # increment global counter
            global_count += 1