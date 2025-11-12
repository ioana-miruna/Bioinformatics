import matplotlib.pyplot as plt
from collections import Counter
import os

def read_fasta(filename):
    with open(filename) as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def find_repetitions(sequence, min_len=6, max_len=10):
    repeats = Counter()
    for size in range(min_len, max_len + 1):
        for i in range(len(sequence) - size + 1):
            fragment = sequence[i:i+size]
            repeats[fragment] += 1
    repeated = {k: v for k, v in repeats.items() if v > 1}
    return repeated

def plot_top_repetitions(repetitions, title, top_n=20, save=False):
    if not repetitions:
        print(f"No repetitions found in {title}.")
        return

    sorted_repeats = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:top_n]
    seqs, freqs = zip(*sorted_repeats)

    plt.figure(figsize=(12, 6))
    plt.bar(seqs, freqs, width=0.5)
    plt.title(f"Top {top_n} DNA Repetitions in {title} (6–10 bp)")
    plt.xlabel("DNA Fragment")
    plt.ylabel("Frequency")
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()

    if save:
        out_name = f"{title.replace(' ', '_')}_top{top_n}_repetitions.png"
        plt.savefig(out_name, dpi=300)
        print(f"Saved plot: {out_name}")

    plt.show()

if __name__ == "__main__":
    for i in range(1, 11):
        filename = f"influenza_{i}.fasta"
        if not os.path.exists(filename):
            print(f"⚠️ File not found: {filename}")
            continue

        print(f"\nProcessing {filename} ...")
        dna_seq = read_fasta(filename)
        repetitions = find_repetitions(dna_seq)
        print(f"  → Found {len(repetitions)} unique repeated fragments.")

        plot_top_repetitions(repetitions, title=f"Influenza {i}", top_n=20, save=False)
