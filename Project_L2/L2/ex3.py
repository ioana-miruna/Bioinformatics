import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt

def sliding_window_frequencies(sequence, window_size=30):
    nucleotides = sorted(set(sequence))
    total_windows = len(sequence) - window_size + 1
    
    freq_vectors = {nuc: [] for nuc in nucleotides}
    
    for i in range(total_windows):
        window = sequence[i:i+window_size]
        for nuc in nucleotides:
            freq_vectors[nuc].append(window.count(nuc)/window_size)
    
    return freq_vectors

def select_file():
    filepath = filedialog.askopenfilename(
        title="Upload FASTA file", 
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )
    if not filepath:
        return
    
    try:
        record = next(SeqIO.parse(filepath, "fasta"))
        sequence = str(record.seq).upper()
        
        freq_vectors = sliding_window_frequencies(sequence, window_size=30)
        
        plt.figure(figsize=(12, 6))
        for nuc, freqs in freq_vectors.items():
            plt.plot(freqs, label=nuc)
        
        plt.title("nucleotide frequencies")
        plt.xlabel("window position")
        plt.ylabel("relative frequency")
        plt.legend()
        plt.grid(True)
        plt.show()
    
    except Exception as e:
        messagebox.showerror("error", f"failure:\n{e}")

root = tk.Tk()
root.title("Analyser - FASTA files")

root.geometry("500x200")  
root.minsize(500, 200)   
root.maxsize(800, 400)   
root.resizable(True, True)

title_label = tk.Label(root, text="FASTA sliding window analyzer")
title_label.pack(pady=20)

btn_select_file = tk.Button(root, text="select FASTA file", command=select_file)
btn_select_file.pack(padx=20, pady=20)

root.mainloop()
