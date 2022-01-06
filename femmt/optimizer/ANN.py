import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
import pandas as pd


def shuffle_in_unison(a, b):
    assert len(a) == len(b)
    shuffled_a = np.empty(a.shape, dtype=a.dtype)
    shuffled_b = np.empty(b.shape, dtype=b.dtype)
    permutation = np.random.permutation(len(a))
    for old_index, new_index in enumerate(permutation):
        shuffled_a[new_index] = a[old_index]
        shuffled_b[new_index] = b[old_index]
    return shuffled_a, shuffled_b


# Goal function
goal = np.array([0, 9.5e-06, 2.2e-05, -1.31e-05])  # [Pv, L11, L22, M]

# Load training and test data
parameters = np.load("parameters_tmp.npy")  # inputs to FEM
results = np.load("results_tmp.npy")  # raw results from FEM

# -- Normalize
# - Inputs (parameter vector)
#   parameters[:, -6:] = parameters[:, -6:] * 1000
#   parameters = parameters / np.linalg.norm(parameters)
max_par = np.abs(parameters).max(axis=0)
norm_parameters = parameters / max_par
# - Results
#   results = results[:, 0]
#   norm_results = (results-goal) / np.linalg.norm(results-goal)
norm_results = (results-goal)
max_res = np.abs(norm_results).max(axis=0)
norm_results = norm_results / max_res

# --
# Initial Shuffle
norm_parameters, norm_results = shuffle_in_unison(norm_parameters, norm_results)

# -- Split the Dataset
# - Training
train_parameters = norm_parameters[:-10, :]
train_results = norm_results[:-10]
# - Testing
test_parameters = norm_parameters[-10:, :]
test_results = norm_results[-10:, :] * max_res + goal
print(max_par, max_res)

# Save the split set
np.save("test_parameters_1", test_parameters)
np.save("test_results_1", test_results)
np.save("train_parameters_1", train_parameters)
np.save("train_results_1", train_results)

"""
print(len(train_parameters[0]))
print(len(test_parameters[0]))
print(len(train_results))
print(len(test_results))
print(train_parameters)
print(test_parameters)
print(train_results)
print(test_results)
"""

# Setup of the Multi Layer Perceptron (the Neuronal Network)
model = tf.keras.models.Sequential()
model.add(tf.keras.layers.Flatten())
model.add(tf.keras.layers.Dense(len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(10*len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(100*len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(100*len(train_parameters[0]), activation=tf.nn.relu))
#model.add(tf.keras.layers.Dense(1000*len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(100*len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(10*len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(len(train_parameters[0]), activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(4))
# model.add(tf.keras.layers.Dense(1, activation=tf.nn.relu))
# We have 4 possible outputs: Pv, L11, L22, M
# model.add(tf.keras.layers.Dense(len(results[0]), activation=tf.nn.softmax))  # classification type

model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_squared_error'])

for i in range(0, 100):
    train_parameters, train_results = shuffle_in_unison(train_parameters, train_results)
    model.fit(train_parameters, train_results, epochs=5)  # Train for some epochs

# Save the model
model.save('loss_optimizer.model')

# Save the predictions of the test features
result_estimations = model.predict(test_parameters)
result_estimations = np.asarray(result_estimations * max_res + goal)  # Unnormalize it
# result_estimations = result_estimations.flatten()


print(f"\n"
      f"{result_estimations}\n"
      f"{test_results}\n")

# Calc the Error
error = (result_estimations-test_results)/test_results


# Print Date
print(f"result_estimations: {result_estimations}\n"
      f"test_results: {[test_results]}\n"
      f"error: {error}")

# Save Data
np.save('result_estimations_1', result_estimations)
np.save('error_1', error)

# np.save('test_features', test_features)
# np.save('test_labels', test_labels)


# Visualize Results
df = pd.DataFrame(error, columns = ['Pv', 'L11', 'L22', 'M'])
df.plot()
plt.show()


