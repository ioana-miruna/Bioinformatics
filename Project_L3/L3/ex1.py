import math

def tm(dna, na_conc=0.05):
    length = len(dna)
    A = dna.count('A')
    T = dna.count('T')
    G = dna.count('G')
    C = dna.count('C')

    tm1 = 4 * (G + C) + 2 * (A + T)
 
    gc_percent = ((G + C) / length) * 100
    tm2 = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)

    return tm1, tm2

dna = "CGACGGACTGCTGCCAACACCCAGG"

tm1, tm2 = tm(dna)

print(f"dna length: {len(dna)} bases")
print(f"GC procent: {((dna.count('G') + dna.count('C')) / len(dna)) * 100:.2f}%")
print(f"tm (using formula 1) = {tm1:.2f} °C")
print(f"tm (using formula 2 -> NA+=0.05) = {tm2:.2f} °C")
