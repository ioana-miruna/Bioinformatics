import random

def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

TRANSPOSONS = [
    "ATGCGTAC",      # TE1 – 8 bp
    "CGTATGCAAG",    # TE2 – 10 bp
    "TTCGAGCTA",     # TE3 – 9 bp
    "GGTACCAT"       # TE4 – 8 bp
]

def insert_transposons(seq, transposons):
    seq_list = list(seq)
    for te in transposons:
        pos = random.randint(0, len(seq) - len(te))
        seq_list[pos:pos+len(te)] = te
    return "".join(seq_list)

base_length = random.randint(200, 400)
dna_sequence = random_dna(base_length)

num_tes = random.randint(3, 4)
chosen_tes = random.sample(TRANSPOSONS, num_tes)
dna_sequence = insert_transposons(dna_sequence, chosen_tes)

def find_all_occurrences(dna, pattern):
    positions = []
    start = 0
    while True:
        idx = dna.find(pattern, start)
        if idx == -1:
            break
        positions.append((idx, idx + len(pattern)))
        start = idx + 1
    return positions

results = {}
for te in TRANSPOSONS:
    matches = find_all_occurrences(dna_sequence, te)
    if matches:
        results[te] = matches

print("\n Artificial DNA Sequence (length =", len(dna_sequence), ")")
print(dna_sequence)

print("\n Detected Transposable Elements")
for te, match_list in results.items():
    print(f"\nTransposon: {te}")
    for (start, end) in match_list:
        print(f"  → Occurrence at positions [{start}, {end})")
