seq = 'ATTTCGCCGATA'
alphabet = {}

for letter in seq:
    if (letter in alphabet):
        alphabet[letter] += 1
    else:
        alphabet[letter] = 1

print(alphabet)