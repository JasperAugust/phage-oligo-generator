import streamlit as st

st.set_page_config(layout="wide", page_title="Phage Tools")

st.title("Phage Analysis Tools")

st.markdown(
    """

This application provides two main tools for phage analysis:

### 1. Oligonucleotide Generator
- Generate sense and antisense oligonucleotides
- Optimize sequences for E. coli K12 codon usage
- Support for degenerate codons (NNK) for library construction

### 2. Phage Sequencing Analysis
- Upload and analyze .ab1 Sanger sequencing files
- Identify peptides and evaluate sequencing quality
- Generate quality metrics and visualizations

Choose a tool from the sidebar to get started.

#### Made by **Jasper August Tootsi** with ❤️
"""
)
