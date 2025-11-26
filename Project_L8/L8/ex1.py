import re
enzymes = {
    "EcoRI": {
        "site": "GAATTC",
        "cut_pos": 1
    },
    "BamHI": {
        "site": "GGATCC",
        "cut_pos": 1
    },
    "HindIII": {
        "site": "AAGCTT",
        "cut_pos": 1
    },
    "TaqI": {
        "site": "TCGA",
        "cut_pos": 1
    },
    "HaeIII": {
        "site": "GGCC",
        "cut_pos": 2
    }
}

dna_sequence = """
ATGCGGATCCTTAACTGAATTCCTAGGTACTTCCGGAATTCATTCGAGGCCATTTAGCGTA
GCTAGCTTACGATCGATTCGAAGCTTAGGCGATCGGCCATCGATAGCTTCGATTCGGATCC
GATCGATGCTTAGCGGCCATTAAGGATCCGATTCGAACGTTAGCTTAGCTTACTTAAAGCT
AAGCTTCGTACGAATTCGATCGATCGGCCAGCTTAGGCCGAATTCGATCGTTAGCGGCCGA
""".replace("\n","").upper()

def digest_dna(enzyme_name, dna, site, cut_offset):
    cut_positions = []
    pattern = re.compile(site)

    for match in pattern.finditer(dna):
        cut_pos = match.start() + cut_offset
        cut_positions.append(cut_pos)

    # Generate fragments
    fragments = []
    if not cut_positions:
        return [], [len(dna)]

    prev = 0
    for pos in cut_positions:
        fragments.append(pos - prev)
        prev = pos
    fragments.append(len(dna) - prev)

    return cut_positions, fragments

def simulate_gel(results):
    print("\n=== Simulated Electrophoresis Gel ===")
    print("(shorter fragments migrate further ↓)\n")

    # Determine range of fragment sizes
    all_fragments = [f for r in results.values() for f in r["fragments"]]
    max_len = max(all_fragments)
    scale = max_len / 40  # compress for gel output

    for enzyme, data in results.items():
        print(f"{enzyme}:")
        for size in sorted(data["fragments"], reverse=True):
            band = int(size / scale)
            print(" " * (45 - band) + "#" * 5)
        print()


results = {}

for name, info in enzymes.items():
    cuts, fragments = digest_dna(name, dna_sequence, info["site"], info["cut_pos"])
    results[name] = {
        "cuts": cuts,
        "fragments": fragments
    }

for enzyme, data in results.items():
    print("\n========================================")
    print(f"Enzyme: {enzyme}")
    print("Recognition site:", enzymes[enzyme]["site"])
    print("Number of cuts:", len(data["cuts"]))
    print("Cut positions:", data["cuts"])
    print("Fragment lengths:", data["fragments"])

simulate_gel(results)