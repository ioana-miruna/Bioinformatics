def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
        traceback[i][0] = "UP"

    for j in range(1, m + 1):
        score[0][j] = j * gap
        traceback[0][j] = "LEFT"

    traceback[0][0] = "DONE"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = score[i - 1][j - 1] + match
            else:
                diag = score[i - 1][j - 1] + mismatch

            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = max(diag, up, left)
            score[i][j] = best

            if best == diag:
                traceback[i][j] = "DIAG"
            elif best == up:
                traceback[i][j] = "UP"
            else:
                traceback[i][j] = "LEFT"

    align1 = []
    align2 = []

    i, j = n, m
    while traceback[i][j] != "DONE":
        if traceback[i][j] == "DIAG":
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == "UP":
            align1.append(seq1[i - 1])
            align2.append("-")
            i -= 1
        elif traceback[i][j] == "LEFT":
            align1.append("-")
            align2.append(seq2[j - 1])
            j -= 1

    align1.reverse()
    align2.reverse()

    return score, "".join(align1), "".join(align2)


def alignment_stats(aln1, aln2):
    matches = sum(1 for a, b in zip(aln1, aln2) if a == b)
    length = len(aln1)
    similarity = (matches / length) * 100
    return matches, length, similarity


S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

score_matrix, aligned_s1, aligned_s2 = needleman_wunsch(
    S1, S2, match=1, mismatch=-1, gap=-1
)

matches, length, similarity = alignment_stats(aligned_s1, aligned_s2)

# Print results
print("Aligned Sequences:\n")
print("S1:", aligned_s1)
print("    ", "".join("|" if a == b else " " for a, b in zip(aligned_s1, aligned_s2)))
print("S2:", aligned_s2)

print("\nStatistics:")
print("Matches   =", matches)
print("Length    =", length)
print(f"Similarity= {similarity:.2f} %")

print("\nFinal alignment score:", score_matrix[len(S1)][len(S2)])
