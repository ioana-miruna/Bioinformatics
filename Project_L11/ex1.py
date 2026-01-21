import math

S1_CpG_plus = "ATCGATTCGATATCATACACGTAT"
S2_CpG_minus = "CTCGACTAGTATGAAGTCCACGCTTG"
S_target = "CGCGACGCGTCG"

bases = ['A', 'C', 'G', 'T']

def get_transition_probabilities(sequence, bases):
    counts = {b1: {b2: 1 for b2 in bases} for b1 in bases}

    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i+1]
        counts[current_base][next_base] += 1

    probs = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}
    for b1 in bases:
        row_total = sum(counts[b1].values())
        for b2 in bases:
            probs[b1][b2] = counts[b1][b2] / row_total

    return probs

def print_matrix(matrix, title):
    print(f"\n--- {title} ---")
    header = "\t" + "\t".join(bases)
    print(header)
    for b1 in bases:
        row_str = f"{b1}\t" + "\t".join([f"{matrix[b1][b2]:.3f}" for b2 in bases])
        print(row_str)

prob_plus = get_transition_probabilities(S1_CpG_plus, bases)
prob_minus = get_transition_probabilities(S2_CpG_minus, bases)

print_matrix(prob_plus, "Step 1: CpG (+) Probabilities (from S1)")
print_matrix(prob_minus, "Step 2: Non-CpG (-) Probabilities (from S2)")

llr_matrix = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}

for b1 in bases:
    for b2 in bases:
        p_plus = prob_plus[b1][b2]
        p_minus = prob_minus[b1][b2]
        # formula: log2( P+(st) / P-(st) )
        llr_matrix[b1][b2] = math.log2(p_plus / p_minus)

print_matrix(llr_matrix, "Step 3: Log-Likelihood Matrix")

score = 0
print(f"\n--- Analysis of Sequence S: {S_target} ---")
print("Transitions scores:")

for i in range(len(S_target) - 1):
    curr = S_target[i]
    nxt = S_target[i+1]
    val = llr_matrix[curr][nxt]
    score += val
    print(f"  {curr}->{nxt}: {val:.3f}")

print("-" * 30)
print(f"Total Log-Likelihood Score: {score:.4f}")

if score > 0:
    print("Result: The sequence belongs to a CpG Island (Model +)")
else:
    print("Result: The sequence does NOT belong to a CpG Island (Model -)")