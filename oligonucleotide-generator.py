import streamlit as st
import pandas as pd
from Bio.Seq import Seq

# codon optimization table for e. coli k12
ecoli_codon_usage = {
    "A": "GCT",
    "R": "CGT",
    "N": "AAT",
    "D": "GAT",
    "C": "TGC",
    "Q": "CAA",
    "E": "GAA",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "L": "CTT",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCT",
    "S": "TCT",
    "T": "ACT",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTT",
    "X": "NNK",  # variability placeholder for degenerate codons
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


def generate_oligonucleotides(peptide_sequences):
    """
    Generates sense and antisense oligonucleotides for given peptide sequences,
    optimized for e. coli codons.

    Parameters:
        peptide_sequences (list): List of tuples with (name, peptide sequence).

    Returns:
        list: List of dictionaries with oligonucleotide sequences.
    """
    results = []

    for name, peptide_seq in peptide_sequences:
        peptide_seq_with_ns = "NS" + peptide_seq

        codons = optimize_sequence(peptide_seq_with_ns)

        sense_seq = "AATTCT" + "".join(codons) + "TAA"  # add stop codon

        antisense_seq = (
            "AGCTTAGCA"
            + "".join(
                ["NNM" if aa == "X" else ecoli_codon_usage[aa] for aa in peptide_seq]
            )
            + "GCAAG"
        )

        # store results
        results.append(
            {
                "Name": name,
                "Peptide Sequence": peptide_seq,
                "Sense Oligonucleotide (5'-3')": f"5'Phos-{sense_seq}",
                "Antisense Oligonucleotide (5'-3')": f"5'Phos-{antisense_seq}",
            }
        )

    return results


st.title("Custom oligo generator")

st.subheader("Paste your data below (tab-separated, Name and Peptide Sequence):")
example_input = "RPAR\tGRPARPAR\niRGD\tCRGDKGPDC\ncx7c\tCXXXXXXXC"
raw_data = st.text_area("Example input:", example_input, height=200)

run_button = st.button("Generate Oligonucleotides")


if raw_data:
    peptide_sequences = []
    try:
        for line in raw_data.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                peptide_sequences.append((parts[0], parts[1]))
            else:
                st.error("Each line must contain exactly two tab-separated columns.")
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

    This tool generates sense and antisense oligonucleotides for constructing phage display libraries or cloning individual peptide sequences. The process involves encoding your peptide sequences into dna oligonucleotides that are optimized for *e. coli k12* codon usage. This ensures efficient expression of your peptides in bacterial systems.

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
