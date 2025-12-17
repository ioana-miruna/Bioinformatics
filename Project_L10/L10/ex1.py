import numpy as np
import pandas as pd

sequences = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
bases = ['A', 'C', 'G', 'T']
seq_len = 9
n_seq = len(sequences)

count_matrix = pd.DataFrame(0, index=bases, columns=range(1, seq_len + 1))
for seq in sequences:
    for i, base in enumerate(seq):
        count_matrix.loc[base, i + 1] += 1

freq_matrix = count_matrix / n_seq

# formula: ln( P(N) / Null_Model ) where Null_Model = 0.25
epsilon = 0.0001
log_likelihood_matrix = freq_matrix.apply(lambda x: np.log((x + epsilon) / 0.25))

results = []
for i in range(len(S) - seq_len + 1):
    window = S[i : i + seq_len]
    score = 0
    for j, char in enumerate(window):
        score += log_likelihood_matrix.loc[char, j + 1]
    results.append({'Position': i, 'Window': window, 'Score': round(score, 4)})

df_results = pd.DataFrame(results).sort_values('Score', ascending=False)

print("--- 1. COUNT MATRIX ---")
print(count_matrix)

print("\n--- 2. RELATIVE FREQUENCIES MATRIX ---")
print(freq_matrix.round(3))

print("\n--- 3. LOG-LIKELIHOODS MATRIX ---")
# Note: Values near -7.8 indicate bases that never appeared (mathematical -infinity)
print(log_likelihood_matrix.round(3))

print("\n--- 4. SLIDING WINDOW SCAN RESULTS FOR S (Top 5) ---")
print(df_results.head(5).to_string(index=False))

# analysis output
best_hit = df_results.iloc[0]
print(f"\nAnalysis: The highest signal is at position {best_hit['Position']} ('{best_hit['Window']}') "
      f"with a score of {best_hit['Score']}. This strongly indicates an exon-intron border.")