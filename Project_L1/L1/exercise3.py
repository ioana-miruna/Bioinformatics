from collections import Counter

def read_fasta(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("File is empty!")

    sequence_lines = []
    description = None

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            description = line
            continue
        sequence_lines.append(line)

    sequence = "".join(sequence_lines)
    return description, sequence

def analyze_sequence(sequence):
    filtered = [base for base in sequence if base in "AGCT"]
    total = len(filtered)

    if total == 0:
        print("No valid nucleotide bases found.")
        return

    counts = Counter(filtered)
    print("Sequence length:", total)
    for base in "AGCT":
        freq = (counts[base] / total) * 100
        print(f"The relative frequency of {base} is {freq:.2f}%")


if __name__ == "__main__":
    fasta_file = "influenza.fna"

    description, sequence = read_fasta(fasta_file)
    analyze_sequence(sequence)