from collections import defaultdict

MIN_REPEAT_LENGTH = 4
MAX_REPEAT_LENGTH = 6
MIN_SPACER_DISTANCE = 5
MAX_SPACER_DISTANCE = 50
MAX_GENOME_SLICE = 10000

def parse_fasta(file_content):
    lines = file_content.split('\n')
    header = ""
    sequence = ""
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:].strip()
        elif line and not line.startswith('['):
            sequence += "".join(c for c in line.upper() if c in 'ACGT')

    return header, sequence

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def find_inverted_repeats(dna_sequence, min_len, max_len, min_spacer, max_spacer):
    results = defaultdict(list)
    n = len(dna_sequence)

    for length in range(min_len, max_len + 1):
        for i in range(n - length):
            s1 = dna_sequence[i:i + length]
            rc_s1 = reverse_complement(s1)
            start_j = i + length + min_spacer
            max_j = min(n - length, i + length + max_spacer)

            for j in range(start_j, max_j + 1):
                s2 = dna_sequence[j:j + length]
                if s2 == rc_s1:
                    spacer_distance = j - (i + length)
                    results[s1].append({
                        "length": length,
                        "s1_start": i,
                        "s1_end": i + length,
                        "s2_start": j,
                        "s2_end": j + length,
                        "spacer_length": spacer_distance
                    })
    return results

file_contents = {
    "GCA_000006945.2_ASM694v2_genomic.fna": ">AE006468.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome\nAGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATATCGCCAGCAGCACATGAACAA\nGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGGCGGGCAGATACTTTAACCAATATAGGAATACAAGACAGAC\nAAATAAAAATGACAGAGTACACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT\n...\n",
    "GCA_000195955.2_ASM19595v2_genomic.fna": ">AL123456.3 Mycobacterium tuberculosis H37Rv complete genome\nTTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGA\nCGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCAT\nTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAGCAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCC\n...\n",
    "GCA_001457635.1_NCTC7465_genomic.fna": ">LN831051.1 Streptococcus pneumoniae genome assembly NCTC7465, chromosome : 1\nAGGTCATGAAGTCAACTTTGATGTTTTGGTATCTCCAAAAGCAGCATTGAAACGATGAGTACCTCAAACAAAAAATTAGC\nAGAAGTTGAAAGCTGCTAAAGAACAGGTGGTCTTATAAACTTTGCTTCACCAGAAAGTAAAAAAGAAGCTTTCTTGAAAG\nCAATTGAAGCGCGAACAAGTGTTGAAAGACCATGAAATTAGCACCCAAGATCAAGTCAATGACCGACTTAATAAATTGAC\n...\n"
}

all_results = {}
for filename, content in file_contents.items():
    header, full_sequence = parse_fasta(content)
    sequence_slice = full_sequence[:MAX_GENOME_SLICE]

    ir_results = find_inverted_repeats(
        sequence_slice,
        MIN_REPEAT_LENGTH,
        MAX_REPEAT_LENGTH,
        MIN_SPACER_DISTANCE,
        MAX_SPACER_DISTANCE
    )

    all_results[filename] = {
        "header": header,
        "slice_length": len(sequence_slice),
        "total_matches": sum(len(matches) for matches in ir_results.values()),
        "unique_repeats": len(ir_results),
        "repeats": ir_results
    }

output_markdown = "## Inverted Repeat (Transposon Motif) Search Results\n\n"
output_markdown += f"*(Analysis performed on the first **{MAX_GENOME_SLICE} bp** of each genome slice using repeat lengths of {MIN_REPEAT_LENGTH}-{MAX_REPEAT_LENGTH} bp and a spacer distance of {MIN_SPACER_DISTANCE}-{MAX_SPACER_DISTANCE} bp.)*\n\n"

for filename, data in all_results.items():
    output_markdown += f"***\n"
    output_markdown += f"### {data['header']} ({filename})\n"
    output_markdown += f"**Total sequence analyzed:** {data['slice_length']} bp\n"
    output_markdown += f"**Total IR Pairs Found:** {data['total_matches']}\n"
    output_markdown += f"**Unique Repeat Motifs Found:** {data['unique_repeats']}\n\n"

    top_repeats = list(data['repeats'].items())[:5]

    if top_repeats:
        output_markdown += "| Repeat (S1) | Length | Reverse Complement (S2) | \n"
        output_markdown += "|:------------|:-------|:------------------------|\n"

        for s1, matches in top_repeats:
            example = matches[0]
            s2 = reverse_complement(s1)
            output_markdown += f"| `{s1}` | {example['length']} | `{s2}` | [{example['s1_start']}, {example['s1_end']}) | [{example['s2_start']}, {example['s2_end']}) | {example['spacer_length']} |\n"
    else:
        output_markdown += "No inverted repeats found under the specified constraints.\n"

print(output_markdown)