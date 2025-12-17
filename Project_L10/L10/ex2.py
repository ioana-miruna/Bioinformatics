import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

training_sequences = [
    "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
    "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
]
bases = ['A', 'C', 'G', 'T']
seq_len = 9
null_model = 0.25

count_matrix = pd.DataFrame(0, index=bases, columns=range(1, seq_len + 1))
for seq in training_sequences:
    for i, base in enumerate(seq):
        count_matrix.loc[base, i + 1] += 1

pwm = np.log((count_matrix + 0.1) / (len(training_sequences) + 0.4) / null_model)

directory = "/Users/ioanabadica/Documents/UNI/YEAR 4/BioInformatics/L10"
genome_files = [f"influenza_{i}.fasta" for i in range(1, 11)]

for file_name in genome_files:
    file_path = os.path.join(directory, file_name)

    if not os.path.exists(file_path):
        print(f"Skipping: {file_name} (File not found)")
        continue

    with open(file_path, 'r') as f:
        genome_seq = "".join(line.strip() for line in f if not line.startswith(">")).upper()

    scores = []
    positions = []

    for i in range(len(genome_seq) - seq_len + 1):
        window = genome_seq[i : i + seq_len]
        current_score = 0
        valid_window = True
        for j, base in enumerate(window):
            if base in pwm.index:
                current_score += pwm.loc[base, j + 1]
            else:
                valid_window = False
                break

        if valid_window:
            scores.append(current_score)
            positions.append(i)

    plt.figure(figsize=(14, 5))
    plt.plot(positions, scores, color='#2c3e50', linewidth=0.7, label='Log-Likelihood Score')
    threshold = 3.0
    high_signals = [(p, s) for p, s in zip(positions, scores) if s > threshold]

    if high_signals:
        peak_pos, peak_val = zip(*high_signals)
        plt.scatter(peak_pos, peak_val, color='#e74c3c', s=15, label=f"High Signal (>{threshold})", zorder=3)

    plt.title(f"Splice Site Motif Analysis: {file_name}", fontsize=14)
    plt.xlabel("Genome Position (Nucleotides)", fontsize=12)
    plt.ylabel("PWM Score", fontsize=12)
    plt.axhline(y=threshold, color='red', linestyle='--', alpha=0.5, label='Threshold')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')
    plt.show()

    if scores:
        max_idx = np.argmax(scores)
        print(f"DONE: {file_name} | Max Score: {scores[max_idx]:.2f} at position {positions[max_idx]}")