import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import copy

def data_load(file):
    #dna_data = file.read().decode('utf-8').replace('\n', '')
    dna_data = pd.read_csv(file, sep='\t', engine='python')
    dna_data['sequence'] = dna_data['sequence'].str.replace('\n', '', regex=False)
    dna_data=dna_data["sequence"].to_string()
    return dna_data

def initialize_dict(sequence):
    binary_sequence = string_to_binary(sequence)
    dictionary = {}
    for char in binary_sequence:
        if char not in dictionary:
            dictionary[char] = len(dictionary)
    return dictionary

def string_to_binary(string_sequence):
    return ''.join(format(ord(x), 'b') for x in string_sequence)

def lzw_compression(binary_value, dictionary):
    index = 0
    compression = []
    dictc = copy.deepcopy(dictionary)
    while index < len(binary_value):
        i = index + 1
        strng = binary_value[index]
        while i < len(binary_value) and (strng + binary_value[i]) in dictc:
            strng += binary_value[i]
            i += 1
        compression.append(dictc.get(strng, len(dictc)))
        if i < len(binary_value):
            dictc[strng + binary_value[i]] = len(dictc)
        index = i
    return len(compression), compression

def compare_dna(seq1, seq2):
    min_length = min(len(seq1), len(seq2))
    match_count = sum(1 for i in range(min_length) if seq1[i] == seq2[i])
    mismatch_count = min_length - match_count
    similarity = (match_count / min_length) * 100
    return match_count, mismatch_count, similarity

def calculate_lzw(sequence, window_size, step):
    dictionary = initialize_dict(sequence)
    start, end = 0, window_size
    lzw_values = []
    while end <= len(sequence):
        window = sequence[start:end]
        binary_window = string_to_binary(window)
        lzw_length, _ = lzw_compression(binary_window, dictionary)
        lzw_values.append(lzw_length / window_size)
        start += step
        end += step
    return lzw_values

def compare_lzw(seq1, seq2, window_size, step):
    lzw1 = calculate_lzw(seq1, window_size, step)
    lzw2 = calculate_lzw(seq2, window_size, step)
    min_length = min(len(lzw1), len(lzw2))
    lzw1 = lzw1[:min_length]
    lzw2 = lzw2[:min_length]
    return lzw1, lzw2

def visualize_lzw(lzw1, lzw2):
    plt.figure(figsize=(10, 5))
    plt.plot(lzw1, label="Sequence 1 LZW Complexity", color="dodgerblue",alpha=0.7, linestyle='--')
    plt.plot(lzw2, label="Sequence 2 LZW Complexity", color="crimson",alpha=0.7, linestyle='-')
    plt.xlabel("Window Position")
    plt.ylabel("Normalized LZW Complexity")
    plt.title("LZW Complexity Comparison")
    plt.grid(alpha=0.3)
    plt.legend()
    st.pyplot(plt)

def main():
    st.title("DNA Sequence LZW Comparison Tool")

    uploaded_file1 = st.file_uploader("Upload the first DNA Sequence file (Text format):", type=["txt"], key="file1")
    uploaded_file2 = st.file_uploader("Upload the second DNA Sequence file (Text format):", type=["txt"], key="file2")

    window_size = st.number_input("Enter the window size:", min_value=1, step=1, value=50)
    step = st.number_input("Enter the step size:", min_value=1, step=1, value=10)

    if uploaded_file1 is not None and uploaded_file2 is not None:
        dna_data1 = data_load(uploaded_file1)
        dna_data2 = data_load(uploaded_file2)
        st.success("Both DNA sequences loaded successfully!")

        if st.button("Compare DNA Sequences Using LZW"):
            match_count, mismatch_count, similarity = compare_dna(dna_data1, dna_data2)
            lzw1, lzw2 = compare_lzw(dna_data1, dna_data2, window_size, step)

            st.write("**LZW Complexity Visualization:**")
            visualize_lzw(lzw1, lzw2)
            st.write(f"- Number of Matches: {match_count}")
            st.write(f"- Number of Mismatches: {mismatch_count}")
            st.write(f"- Similarity: {similarity:.2f}%")

if __name__ == "__main__":
    main()
