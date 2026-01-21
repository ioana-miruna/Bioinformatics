import math
import string
import matplotlib.pyplot as plt

#Mihai Eminescu - "La Steaua"
eminescu_train = """
La steaua care-a răsărit E-o cale-atât de lungă,
Că mii de ani i-au trebuit Luminii să ne-ajungă.
Poate de mult s-a stins în drum În depărtări albastre,
Iar raza ei abia acum Luci vederii noastre.
Icoana stelei ce-a murit Încet pe cer se suie:
Era pe când nu s-a zărit, Azi o vedem, și nu e.
Tot astfel când al nostru dor Pieri în noapte-adâncă,
Lumina stinsului amor Ne urmărește încă.
"""

#Nichita Stănescu - "Emoție de toamnă"
stanescu_train = """
A venit toamna, acoperă-mi inima cu ceva,
cu umbra unui copac sau mai bine cu umbra ta.
Mă tem că n-am să te mai văd, uneori,
că or să-mi crească aripi ascuțite până la nori,
că ai să te ascunzi într-un ochi străin,
și el o să se-nchidă cu-o frunză de pelin.
Și-atunci mă apropii de pietre și tac,
iau cuvintele și le-nec în mare.
Șuier luna și o răsar și o prefac
într-o dragoste mare.
"""

# Test Text (Combined): Eminescu ("Somnoroase păsărele") + Stănescu ("Leoaică tânără, iubirea")
test_text_combined = """
Somnoroase păsărele Pe la cuiburi se adună,
Se ascund în rămurele - Noapte bună!
Doar izvoarele suspină, Pe când codrul negru tace;
Dorm și florile-n grădină - Dormi în pace!
Trece lebăda pe ape Între trestii să se culce -
Fie-ți îngerii aproape, Somnul dulce!
Peste-a nopții feerie Se ridică mândra lună,
Totu-i vis și armonie - Noapte bună!

Leoaica tânără, iubirea mi-ai sărit în față.
Mă pândise-n încordare mai demult.
Colții albi mi i-a înfipt în față,
m-a mușcat leoaica, azi, de față.
Și deodata-n jurul meu, natura
se făcu un cerc, de-a-dura,
când mai larg, când mai aproape,
ca o strângere de ape.
Și privirea-n sus țișni, curcubeu tăiat în două,
și auzul o-ntâlni tocmai lângă ciocârlii.
Mi-am dus mâna la sprânceană,
la tâmplă și la bărbie, dar mâna nu le mai știe.
"""

def clean_and_tokenize(text):
    translator = str.maketrans(string.punctuation, ' ' * len(string.punctuation))
    text = text.translate(translator).lower()
    return text.split()

tokens_eminescu = clean_and_tokenize(eminescu_train)
tokens_stanescu = clean_and_tokenize(stanescu_train)
tokens_test = clean_and_tokenize(test_text_combined)

vocabulary = set(tokens_eminescu).union(set(tokens_stanescu))

def train_markov_model(tokens, vocab):
    counts = {w1: {w2: 1 for w2 in vocab} for w1 in vocab}

    for i in range(len(tokens) - 1):
        w1 = tokens[i]
        w2 = tokens[i+1]
        if w1 in vocab and w2 in vocab:
            counts[w1][w2] += 1

    probs = {w1: {} for w1 in vocab}
    for w1 in vocab:
        total_count = sum(counts[w1].values())
        for w2 in vocab:
            probs[w1][w2] = counts[w1][w2] / total_count

    return probs

prob_eminescu = train_markov_model(tokens_eminescu, vocabulary)
prob_stanescu = train_markov_model(tokens_stanescu, vocabulary)

llr_matrix = {w1: {} for w1 in vocabulary}

for w1 in vocabulary:
    for w2 in vocabulary:
        p1 = prob_eminescu[w1][w2]
        p2 = prob_stanescu[w1][w2]
        llr_matrix[w1][w2] = math.log(p1 / p2)

window_size = 15
scores = []
indices = []

for i in range(len(tokens_test) - window_size):
    window = tokens_test[i : i + window_size]
    window_score = 0

    for j in range(len(window) - 1):
        w1 = window[j]
        w2 = window[j+1]

        if w1 in llr_matrix and w2 in llr_matrix[w1]:
            window_score += llr_matrix[w1][w2]
        else:
            pass

    scores.append(window_score)
    indices.append(i)

plt.figure(figsize=(12, 6))

plt.plot(indices, scores, color='purple', linewidth=1.5, label='Sliding Window Score')
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.fill_between(indices, scores, 0, where=[s > 0 for s in scores],
                 facecolor='blue', alpha=0.3, interpolate=True, label='Classified: Eminescu')
plt.fill_between(indices, scores, 0, where=[s < 0 for s in scores],
                 facecolor='orange', alpha=0.3, interpolate=True, label='Classified: Stănescu')
plt.title("Stylometric Analysis: Eminescu vs. Stănescu (Log-Likelihood Scanner)")
plt.xlabel("Window Position (Word Index)")
plt.ylabel("Log-Likelihood Score\n(>0 Eminescu, <0 Stănescu)")
plt.legend()
plt.grid(True, alpha=0.3)

plt.text(10, max(scores)*0.8, "Actual: Eminescu\n(Somnoroase păsărele)",
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

plt.text(len(tokens_test)-40, min(scores)*0.8, "Actual: Stănescu\n(Leoaica tânără, iubirea)",
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('poetry_analysis.png')
print("Chart generated: poetry_analysis.png")