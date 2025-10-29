seq = 'ATTTCGCCGATA'
alphabet = {}

for letter in seq:
    if (letter in alphabet):
        alphabet[letter] += 1
    else:
        alphabet[letter] = 1

for item in alphabet:
    freq = (alphabet[item] / len(seq)) * 100
    print('The relative frequency of ' + item + ' is ' + str(freq) + '%')