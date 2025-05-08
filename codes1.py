import streamlit as st
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from io import StringIO, BytesIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Page config
st.set_page_config(page_title="Phylogenetic Tree App", page_icon="🧬", layout="wide")

# Custom CSS with background images
st.markdown("""
    <style>
    .tab1 {
        background-image: url('https://images.unsplash.com/photo-1623487882315-bb65cb0b4b54?ixlib=rb-4.0.3&auto=format&fit=crop&w=1600&q=60');
        background-size: cover;
        padding: 20px;
        border-radius: 10px;
    }
    .tab2 {
        background: linear-gradient(135deg, #ece9e6, #ffffff);
        padding: 20px;
        border-radius: 10px;
    }
    .stMarkdown a {
        color: #1a73e8;
        text-decoration: none;
        font-weight: bold;
    }
    </style>
""", unsafe_allow_html=True)

# Tabs
tabs = st.tabs(["Phylogenetic Tree Builder", "About"])

# PHYLOGENETIC TREE BUILDER tab
with tabs[0]:
    st.markdown('<div class="tab1">', unsafe_allow_html=True)
    st.title("🧬 Phylogenetic Tree Builder & Explorer")

    fasta_input = st.text_area("Paste your DNA/Protein FASTA sequences here", height=300)

    if fasta_input:
        # Parse FASTA
        seq_records = []
        lines = fasta_input.strip().split('\n')
        current_id = ""
        current_seq = ""
        for line in lines:
            if line.startswith('>'):
                if current_id:
                    seq_records.append(SeqRecord(Seq(current_seq), id=current_id))
                current_id = line[1:].strip()
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_id:
            seq_records.append(SeqRecord(Seq(current_seq), id=current_id))

        if len(seq_records) < 2:
            st.error("Please paste at least two FASTA sequences.")
        else:
            min_len = min(len(rec.seq) for rec in seq_records)
            for rec in seq_records:
                rec.seq = rec.seq[:min_len]

            names = [rec.id for rec in seq_records]
            matrix = []
            for rec1 in seq_records:
                row = []
                for rec2 in seq_records:
                    matches = sum(a == b for a, b in zip(rec1.seq, rec2.seq))
                    identity = matches / min_len
                    distance = 1 - identity
                    row.append(distance)
                matrix.append(row)
            np_matrix = np.array(matrix)

            lower_triangle = []
            for i in range(len(names)):
                lower_triangle.append(list(np_matrix[i][:i+1]))

            dm = _DistanceMatrix(names=names, matrix=lower_triangle)

            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)

            newick_str = StringIO()
            Phylo.write(tree, newick_str, "newick")
            newick_str.seek(0)

            tree_from_file = Phylo.read(newick_str, "newick")

            # ASCII Tree View
            st.markdown("### Phylogenetic Tree (ASCII View)")
            ascii_str = StringIO()
            Phylo.draw_ascii(tree_from_file, file=ascii_str)
            ascii_output = ascii_str.getvalue()
            st.text(ascii_output)

            # Download button for ASCII tree
            ascii_bytes = BytesIO(ascii_output.encode())
            st.download_button("Download ASCII Tree (.txt)", ascii_bytes, file_name="phylo_ascii_tree.txt")

            # Matplotlib Tree View
            st.markdown("### Phylogenetic Tree (Matplotlib View)")
            fig, ax = plt.subplots(figsize=(10, 6))
            Phylo.draw(tree_from_file, axes=ax, do_show=False)

            # Save tree figure to BytesIO buffer
            tree_buffer = BytesIO()
            fig.savefig(tree_buffer, format='png', bbox_inches='tight')
            tree_buffer.seek(0)

            st.pyplot(fig)
            plt.close(fig)

            # Download button for matplotlib tree
            st.download_button("Download Tree Image (.png)", tree_buffer, file_name="phylo_tree.png", mime="image/png")

            # Heatmap View
            st.markdown("### Distance Matrix Heatmap")
            dm_df = pd.DataFrame(np_matrix, index=names, columns=names)
            fig2, ax = plt.subplots()
            sns.heatmap(dm_df, annot=True, cmap="coolwarm", ax=ax)

            # Save heatmap to BytesIO buffer
            heatmap_buffer = BytesIO()
            fig2.savefig(heatmap_buffer, format='png', bbox_inches='tight')
            heatmap_buffer.seek(0)

            st.pyplot(fig2)
            plt.close(fig2)

            # Download button for heatmap
            st.download_button("Download Heatmap (.png)", heatmap_buffer, file_name="distance_heatmap.png", mime="image/png")

    else:
        st.info("Paste your FASTA sequences above to begin.")
    st.markdown('</div>', unsafe_allow_html=True)

# ABOUT tab
with tabs[1]:
    st.markdown('<div class="tab2">', unsafe_allow_html=True)
    st.title("🏠 HOME PAGE")

    st.header("About")
    st.write("""
        Welcome to the **Phylogenetic Tree Builder & Explorer**!  
        This web application allows users to:
        - Paste DNA or protein sequences (FASTA format)
        - Calculate distance matrices
        - Build phylogenetic trees using Neighbor-Joining
        - Explore interactive visualizations and ASCII trees

        This tool is designed to help students, researchers, and bioinformatics enthusiasts easily perform phylogenetic analysis.
    """)

    st.header("Creator")
    col1, col2 = st.columns([1, 3])

    with col1:
        st.image("my_photo.jpg", width=150)  # Replace with your image file or URL

    with col2:
        st.subheader("Mrudula Tushar Talegaonkar")
        st.write("""
            **MSc Bioinformatics Student (DES Pune University)**  
            Passionate about genomics, computational biology, and developing bioinformatics tools to simplify complex analyses.

            🔗 [LinkedIn](https://www.linkedin.com/in/mrudula-talegaonkar-22b8482b4?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app)  
            💻 [GitHub](https://github.com/Mrudula-02)
        """)

    st.header("Acknowledgement")
    st.write("""
        I would like to sincerely thank **Dr. Kushagra Kashyap Sir** and **Dr. Poonam Deshpande Ma’am** for their valuable support and guidance throughout this project.
        Their encouragement, suggestions, and help were key to the successful completion of this phylogenetic tree app.
        I am truly grateful for their time and advice, which greatly improved my learning and this project. Thank you!
    """)
    st.markdown('</div>', unsafe_allow_html=True)
