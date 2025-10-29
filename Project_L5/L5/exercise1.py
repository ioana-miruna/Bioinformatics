import random
from Bio import Entrez, SeqIO

Entrez.email = "ioana.badica22@gmail.com"

handle = Entrez.efetch(db="nucleotide", id="NM_001200001", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq)[:2000]
with open("original_dna.txt", "w") as f:
    f.write(sequence)

samples = []
for _ in range(2000):
    start = random.randint(0, len(sequence) - 151)
    length = random.randint(100, 150)
    fragment = sequence[start:start + length]
    samples.append(fragment)

reconstructed = ''.join(samples)

with open("reconstructed_dna.txt", "w") as f:
    f.write(reconstructed)

print("Files created:")
print(" - original_dna.txt (original sequence)")
print(" - reconstructed_dna.txt (naive reconstruction)")
print(" - answers.txt (explanation)")

