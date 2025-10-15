import math
import matplotlib.pyplot as plt

def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq.upper()


def calculate_tm(seq, na_conc=0.05):
    G = seq.count('G')
    C = seq.count('C')
    length = len(seq)
    if length == 0:
        return 0.0
    gc_percent = ((G + C) / length) * 100
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm


def sliding_window_tm(seq, window_size=8, step=1):
    results = []
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        tm = calculate_tm(window)
        results.append((i + 1, tm))
    return results


def main():
    fasta_file = input("Enter the path to the FASTA file: ").strip()
    sequence = read_fasta(fasta_file)
    print(f"\nLoaded sequence of length {len(sequence)} bases.")

    window_size = 8
    tm_values = sliding_window_tm(sequence, window_size)

    print(f"\nSliding Window (size={window_size}) Tm values:\n")
    print("Position\tTm (°C)")
    for pos, tm in tm_values:
        print(f"{pos}\t\t{tm:.2f}")

    positions = [pos for pos, _ in tm_values]
    temps = [tm for _, tm in tm_values]

    plt.figure(figsize=(10, 5))
    plt.plot(positions, temps, color="blue", linewidth=1.5)
    plt.title(f"Melting Temperature (Tm) Profile\nSliding Window = {window_size} bp")
    plt.xlabel("Position along sequence (bp)")
    plt.ylabel("Tm (°C)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()