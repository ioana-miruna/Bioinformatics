import numpy as np

def predict_steps(matrix, initial_vector, steps=5):
    A = np.array(matrix)
    v = np.array(initial_vector)

    predictions = [v.copy()]
    for _ in range(steps):
        v = A @ v
        predictions.append(v.copy())

    return predictions

if __name__ == "__main__":
    A = [[0.5, 0.1],
         [0.2, 0.7]]
    v0 = [1, 0]

    predictions = predict_steps(A, v0, steps=5)

    for step, vec in enumerate(predictions):
        print(f"Step {step}: {vec}")
