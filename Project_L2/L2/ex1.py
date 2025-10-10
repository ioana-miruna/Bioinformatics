S = "ATTGTCCCAATCTGTTG"

nucs = ['A', 'C', 'G', 'T']

dins = [a + b for a in nucs for b in nucs]
tris = [a + b + c for a in nucs for b in nucs for c in nucs]

def percentage(seq, n):
    total = len(S) - n + 1
    return round(seq.count(S[:0]) / total * 100, 2)

def calc_percentages(S, combos, n):
    total = len(S) - n + 1
    result = {}
    for combo in combos:
        count = 0
        for i in range(total):
            if S[i:i+n] == combo:
                count += 1
        result[combo] = round(count / total * 100, 2)
    return result

din_p = calc_percentages(S, dins, 2)
tri_p = calc_percentages(S, tris, 3)

print("dinucleotide percentages")
for k, v in din_p.items():
    print(f"{k}: {v}%")

print("\n trinucleotide percentages")
for k, v in tri_p.items():
    print(f"{k}: {v}%")