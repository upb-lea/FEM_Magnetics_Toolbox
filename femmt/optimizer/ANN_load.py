import numpy as np
import tensorflow as tf
import itertools

# create dense parameter space
N1 = [5, 6, 7]
N2 = [9, 10, 11]
window_h = [0.025, 0.030, 0.035]
window_w = [0.08, 0.01, 0.012]
core_w = [0.018, 0.02, 0.022]
air_gap_h = [0.0019, 0.0020, 0.0021]
conductor_radius1 = [0.0014, 0.0015, 0.0016, 0.0017, 0.0018]
conductor_radius2 = [0.0014, 0.0015, 0.0016, 0.0017, 0.0018]
# x_ = [N1, N2, window_h, window_w, core_w, air_gap_h, conductor_radius1, conductor_radius2]
x_dense = list(itertools.product(N1, N2, window_h, window_w, core_w, air_gap_h, conductor_radius1, conductor_radius2))


# Load Test Data
train_parameters = np.load("train_parameters.npy")
train_results = np.load("train_results.npy")
print(np.load("results_tmp.npy"))
print(np.load("parameters_tmp.npy"))

# Load the model
model = tf.keras.models.load_model('number_reader.model')


# Load test Data
test_parameters = np.load("test_parameters.npy")
test_results = np.load("test_results.npy")

# Save the predictions of the test features
result_estimations = model.predict(test_parameters)
result_estimations = np.asarray(result_estimations)
result_estimations = result_estimations.flatten()

# Calc the Error
error = (result_estimations-test_results)/test_results

# Print Date
print(f"result_estimations: {result_estimations}\n"
      f"test_results: {[test_results]}\n"
      f"error: {error}")

# Try on dense parameter set
dense_estimations = model.predict(x_dense)
best_one = np.argmin(dense_estimations)

print(x_dense[best_one], dense_estimations[best_one])