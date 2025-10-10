S = "ATTGTCCCAATCTGTTG"

def calc_percentages(S, n):
    counts = {}
    total = len(S) - n + 1
    
    for i in range(total):
        fragment = S[i:i+n]
        if fragment in counts:
            counts[fragment] += 1
        else:
            counts[fragment] = 1
    
    percentages = {k: round(v / total * 100, 2) for k, v in counts.items()}
    return percentages

din_p = calc_percentages(S, 2)
print("dinucleotide percentages:")
for k, v in din_p.items():
    print(f"{k}: {v}%")

tri_p = calc_percentages(S, 3)
print("\n trinucleotide percentages:")
for k, v in tri_p.items():
    print(f"{k}: {v}%")
