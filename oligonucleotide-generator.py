import streamlit as st

st.set_page_config(layout="wide")
import pandas as pd
from Bio.Seq import Seq

ecoli_codon_usage = {
    "A": "GCG",
    "R": "CGC",
    "N": "AAC",
    "D": "GAT",
    "C": "TGC",
    "Q": "CAG",
    "E": "GAA",
    "G": "GGC",
    "H": "CAT",
    "I": "ATT",
    "L": "CTG",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCG",
    "S": "TCT",
    "T": "ACC",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTG",
    "X": "NNK",  # any amino acid (used in degenerate codon scenarios)
}


def optimize_sequence(peptide_seq):
    """
    Optimizes a peptide sequence into a DNA sequence using e. coli k12 codon usage.
    Uses NNK for X (degenerate amino acids).
    """
    codons = []
    for aa in peptide_seq:
        if aa in ecoli_codon_usage:
            codons.append(ecoli_codon_usage[aa])
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    return codons


def reverse_complement_dna(seq):
    """
    Returns the reverse complement of a DNA sequence with degenerate bases.
    For NNK codons, returns MNN instead of complementing individual bases.
    """
    # Split into codons
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]

    reverse_complement_codons = []
    for codon in codons:
        if codon == "NNK":
            reverse_complement_codons.append("MNN")
        else:
            # Regular DNA complement for non-degenerate codons
            complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
            rev_comp = "".join(complement_map[base] for base in codon[::-1])
            reverse_complement_codons.append(rev_comp)

    return "".join(reverse_complement_codons[::-1])  # Reverse the order of codons


def generate_oligonucleotides(peptide_sequences):
    results = []

    for name, peptide_seq in peptide_sequences:
        codons = optimize_sequence(peptide_seq)
        peptide_dna = "".join(codons)

        # sense and antisense sequences
        sense_seq = "AATTCT" + peptide_dna + "TA"
        antisense_seq = "AGCTT" + reverse_complement_dna(peptide_dna) + "AG"

        results.append(
            {
                "Name": name,
                "Peptide Sequence": peptide_seq,
                "Sense Oligonucleotide (5'-3')": f"5'Phos-{sense_seq}",
                "Antisense Oligonucleotide (5'-3')": f"5'Phos-{antisense_seq}",
            }
        )

    return results


example_input = "RPAR\tGRPARPAR\niRGD\tCRGDKGPDC\ncx7c\tCXXXXXXXC"

st.title("Custom oligo generator")

st.subheader(
    "Paste your data below (tab or comma-separated, Name and Peptide Sequence):"
)
raw_data = st.text_area("Example input:", example_input, height=200)

run_button = st.button("Generate Oligonucleotides")

if raw_data:
    peptide_sequences = []
    try:
        for line in raw_data.strip().split("\n"):
            # Try tab first, then comma if tab fails
            parts = line.split("\t") if "\t" in line else line.split(",")
            if len(parts) == 2:
                peptide_sequences.append((parts[0].strip(), parts[1].strip()))
            else:
                st.error(
                    "Each line must contain exactly two columns separated by tab or comma."
                )
    except Exception as e:
        st.error(f"Error parsing input: {e}")

    if peptide_sequences:
        results = generate_oligonucleotides(peptide_sequences)
        df = pd.DataFrame(results)
        st.subheader("Resulting oligonucleotides:")
        st.dataframe(df)

        # download button for results
        st.download_button(
            "Download as CSV",
            data=df.to_csv(index=False),
            file_name="optimized_oligonucleotides.csv",
            mime="text/csv",
        )

st.markdown(
    """
    # **More info**

    This tool generates sense and antisense oligonucleotides for constructing phage display libraries or cloning individual peptide sequences. The process involves encoding your peptide sequences into dna oligonucleotides that are optimized for *e. coli k12* codon usage.

    ## **How it works**
    Phage libraries are constructed by inserting random or specific peptide sequences into a vector (e.g., t7select vectors from novagen). The peptide sequences are back-translated into DNA sequences, which must meet specific requirements for compatibility with the vector:

    1. **Random libraries**: 
    - if your peptide includes random amino acids (e.g., `X`), the tool uses **nnk codons**. These codons:
        - encode all 20 amino acids.
        - avoid generating stop codons like taa and tga.
    - this ensures coverage of all possible peptide variations in your library.

    2. **Individual phages**:
    - for specific peptide sequences (e.g., `CRGDKGPDC`), the tool optimizes the dna sequence for **e. coli k12 codon usage** to maximize expression efficiency.

    3. **Output design**:
    - the generated oligonucleotides are synthesized with:
        - **complementary overhangs** to be compatible with commercial *ecoRI* and *HindIII-cleaved* t7 vectors.
        - **5'-phosphorylation** to allow for ligation into the vector.
    - both **sense (5'-3')** and **antisense (5'-3')** strands are outputted, ready for ordering from a synthesis company.

    ## **How to use**
    1. paste your peptide data in the following format:
    name1 peptide_sequence1 name2 peptide_sequence2

    example input:
    ```css
    RPAR GRPARPAR 
    iRGD CRGDKGPDC 
    cx7c CXXXXXXXC
    ```
    PS! You can also use comma separated values.

    2. The tool will:
    - back-translate the peptides into optimized dna oligonucleotides.
    - generate **sense** and **antisense** oligos for each peptide.
    - display the results in a table and provide an option to download as a csv file.

    3. copy the generated oligonucleotides and submit them to a dna synthesis company for your phage library construction or peptide cloning experiments.

    ## **About this tool**
    this tool:
    - performs *e. coli k12 codon optimization* to ensure efficient peptide expression.
    - uses **nnk codons** for random amino acids in phage display libraries.
    - simplifies the process of designing oligonucleotides compatible with commercial t7 vectors.

    #### Made by **Jasper August Tootsi** with ❤️
    """
)
