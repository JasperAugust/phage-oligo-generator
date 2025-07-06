import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


EXPECTED_CODONS = {
    "M": "ATG",
    "L": "TTG",
    "G": "GGT",
    "D": "GAT",
    "P": "CCA",
    "N": "AAT",
    "S": "TCT",
}


def find_peptide(
    sequence, motif="MLGDPNS", start_codon="ATG", next_gene_indicator="M"
):
    """
    Searches for the peptide motif in all reading frames. Handles cases where the sequence is
    missing a STOP codon or includes parts of the next gene.
    Returns the peptide sequence or an appropriate error message.
    """
    try:
        # Iterate through all three reading frames 
        for frame in range(3):
            subseq = sequence[frame:]  # Adjust for the reading frame
            translated = Seq(subseq).translate()

            # Search for the motif in the translated sequence
            motif_index = str(translated).find(motif)
            if motif_index != -1:
                # Extract the peptide sequence starting from the motif
                peptide_start = motif_index + len(motif)
                peptide = translated[peptide_start:]

                # Check for STOP codon
                stop_index = peptide.find("*")
                if stop_index != -1:
                    return (
                        str(peptide[:stop_index]),
                        frame,
                        None,
                    )  # Valid peptide up to the STOP codon
                else:
                    # Handle cases where STOP codon is missing
                    if next_gene_indicator in peptide:
                        peptide = peptide.split(next_gene_indicator)[0]
                    return (
                        str(peptide),
                        frame,
                        "STOP codon missing, sequence may contain the start of the next gene.",
                    )

        return None, None, "Expected sequence motif not found in any reading frame."
    except Exception as e:
        return None, None, str(e)


# Function to extract quality metrics
def extract_quality_metrics(record):
    """
    Extracts quality metrics from an ABI sequencing record.
    Returns a dictionary of metrics.
    """
    try:
        quality_scores = record.letter_annotations["phred_quality"]

        # Calculate metrics
        avg_quality = np.mean(quality_scores)
        std_quality = np.std(quality_scores)
        q20_count = sum(q >= 20 for q in quality_scores)
        q30_count = sum(q >= 30 for q in quality_scores)
        total_bases = len(quality_scores)
        percent_q20 = (q20_count / total_bases) * 100
        percent_q30 = (q30_count / total_bases) * 100
        ambiguous_bases = record.seq.count("N")

        return {
            "Avg Quality Score": avg_quality,
            "Quality Std Dev": std_quality,
            "Percent Q20": percent_q20,
            "Percent Q30": percent_q30,
            "Ambiguous Bases (N)": ambiguous_bases,
        }
    except Exception as e:
        return {"Error": str(e)}


def is_valid_ab1(file):
    """
    Validates if the uploaded file is a valid .ab1 file by checking the "ABIF" marker.
    """
    try:
        file.seek(0)  # Reset file pointer
        marker = file.read(4)  # Read the first 4 bytes
        file.seek(0)  # Reset file pointer again for SeqIO
        return marker == b"ABIF"
    except Exception:
        return False


def plot_quality_scores(quality_scores, start, end):
    """
    Plots quality scores for a specific region of the sequence.
    """
    region_scores = quality_scores[start:end]
    plt.figure(figsize=(10, 4))
    plt.plot(region_scores, label="Quality Score")
    plt.axhline(20, color="orange", linestyle="--", label="Q20 Threshold")
    plt.axhline(30, color="green", linestyle="--", label="Q30 Threshold")
    plt.title("Quality Scores for the region")
    plt.xlabel("Base Position (Relative to Region)")
    plt.ylabel("Quality Score")
    plt.legend()
    plt.grid()
    st.pyplot(plt)


st.title("T7 Phage Peptide Sequence Analyzer with Quality Metrics")
st.write(
    "Upload your .ab1 Sanger sequencing files to identify peptides and evaluate sequencing quality."
)

uploaded_files = st.file_uploader(
    "Upload .ab1 files",
    type=["ab1"],
    accept_multiple_files=True,
    key="unique_file_uploader",
)
# Add clear button
if st.button("Clear Files and Data"):
    st.session_state.clear()
    uploaded_files = None
    st.rerun()

if uploaded_files:
    results = []
    for uploaded_file in uploaded_files:
        # Validate file format
        if not is_valid_ab1(uploaded_file):
            st.warning(
                f"File '{uploaded_file.name}' is not a valid .ab1 file or is corrupted."
            )
            continue

        try:
            record = SeqIO.read(uploaded_file, "abi")
            nucleotide_seq = str(record.seq)

            # Find peptide sequence
            peptide, frame, error = find_peptide(nucleotide_seq)

            # Extract quality metrics
            quality_metrics = extract_quality_metrics(record)

            # Append results
            results.append(
                {
                    "File Name": uploaded_file.name,
                    "Peptide Sequence": peptide if peptide else "Error",
                    "Peptide Error": error if error else "None",
                    "Avg Quality Score": quality_metrics.get(
                        "Avg Quality Score", "Error"
                    ),
                    "Percent Q20": quality_metrics.get("Percent Q20", "Error"),
                    "Percent Q30": quality_metrics.get("Percent Q30", "Error"),
                    "Ambiguous Bases (N)": quality_metrics.get(
                        "Ambiguous Bases (N)", "Error"
                    ),
                    "Quality Scores": record.letter_annotations.get(
                        "phred_quality", []
                    ),
                    "Nucleotide Sequence": nucleotide_seq,
                    "Frame": frame,
                }
            )
        except Exception as e:
            st.error(f"An error occurred while processing '{uploaded_file.name}': {e}")
            continue

    # Display results table
    if results:
        df = pd.DataFrame(
            [
                {
                    "File Name": r["File Name"],
                    "Peptide Sequence": r["Peptide Sequence"],
                    "Peptide Error": r["Peptide Error"],
                    "Avg Quality Score": r["Avg Quality Score"],
                    "Percent Q20": r["Percent Q20"],
                    "Percent Q30": r["Percent Q30"],
                    "Ambiguous Bases (N)": r["Ambiguous Bases (N)"],
                }
                for r in results
            ]
        )
        st.write("### Analysis Results")
        st.write(df)

        # Allow results to be downloaded as CSV
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Download Analysis Results",
            data=csv,
            file_name="analysis_results.csv",
            mime="text/csv",
        )

        # Generate plots under expandable sections
        st.write("### Detailed Plots")
        for result in results:
            if result["Peptide Sequence"] != "Error":
                quality_scores = result["Quality Scores"]
                nucleotide_seq = result["Nucleotide Sequence"]
                frame = result["Frame"]
                motif_start = (
                    nucleotide_seq.find("ATG") + frame
                )  # Start of MLGDPNS in nucleotides
                region_start = max(0, motif_start)
                region_end = min(len(quality_scores), region_start + 100)

                with st.expander(f"Plots for {result['File Name']}"):
                    st.write(f"#### Quality Scores for {result['File Name']}")
                    plot_quality_scores(quality_scores, region_start, region_end)
