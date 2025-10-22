genetic_code = {
    'UUU':'Phe', 'UUC':'Phe', 'UUA':'Leu', 'UUG':'Leu',
    'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu',
    'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met',
    'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val',

    'UCU':'Ser', 'UCC':'Ser', 'UCA':'Ser', 'UCG':'Ser',
    'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
    'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
    'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',

    'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
    'CAU':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
    'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',

    'UGU':'Cys', 'UGC':'Cys', 'UGA':'Stop', 'UGG':'Trp',
    'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
    'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
    'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
}

def rna_to_protein(rna_sequence):
    start = rna_sequence.find('AUG')
    if start == -1:
        return "No start sequence found."

    protein = []
    for i in range(start, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = genetic_code.get(codon, '???')
        if amino_acid == 'Stop':
            break
        protein.append(amino_acid)

    return ','.join(protein)

rna_input = "AUGGUUUUCGAACUGAGUGA"
protein_sequence = rna_to_protein(rna_input)
print("RNA sequence:", rna_input)
print("Protein sequence:", protein_sequence)
