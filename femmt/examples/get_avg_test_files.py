import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# the losses are extracted and then plotted
# the plot is not precise until now
# the window size and threshold can be controlled to detect the steady state region, so it is a static code
# just a try to detect the steady state region and to get avg.
# # it is approperiate for 200 step time


def get_avg_file(file_path, threshold, loss_label):
    # Read and parse the file
    times = []
    losses = []
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        line_values = line.split()
        if len(line_values) == 2:
            times.append(float(line_values[0]))
            losses.append(float(line_values[1]))

    # Convert to a pandas DataFrame for easier manipulation
    df = pd.DataFrame({
        'Time': times,
        loss_label: losses
    })

    # Compute rolling average with a window of 27
    df['rolling_avg'] = df[loss_label].rolling(window=27).mean()

    # Calculate the differences between the rolling average at each point and the next point
    differences = df['rolling_avg'].diff().abs()

    # Find the first point where the difference is smaller than the threshold
    constant_index = np.where(differences < threshold)[0]

    if len(constant_index) > 0:
        steady_start_idx = constant_index[0]
        steady_end_idx = len(df) - 1
        print(f"Steady state starts at index {steady_start_idx} and ends at {steady_end_idx} for {file_path}")
    else:
        print(f"The rolling average for {file_path} does not reach a constant value within the given threshold.")
        return

    # Calculate the average from the start of the steady state to the end
    avg_steady = df[loss_label].iloc[steady_start_idx:].mean()
    print("avg_steady for", file_path, ":", avg_steady)

    # Calculate the running average until steady_start_idx
    running_avg_transient = np.cumsum(df[loss_label].iloc[:steady_start_idx]) / np.arange(1, steady_start_idx + 1)

    # Combine the running averages
    running_avg = np.concatenate([running_avg_transient, [avg_steady] * (len(df) - steady_start_idx)])

    # Plot the data and the running average
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time'], df[loss_label], label=loss_label)
    plt.plot(df['Time'], running_avg, label='Running average')
    plt.axvline(x=df['Time'].iloc[steady_start_idx], color='r', linestyle='--', label='Estimated start of steady state')
    plt.axvline(x=df['Time'].iloc[steady_end_idx], color='g', linestyle='--', label='Estimated end of steady state')
    plt.xlabel('Time (s)')
    plt.ylabel(loss_label)
    plt.legend()
    plt.grid(True)
    plt.title(f"Analysis for {file_path}")
    plt.show()


# Paths for the files
j2F_path = "C:\\Users\\uthmn\\PycharmProjects\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\inductor\\results\\values\\j2F_1.dat"
turn_losses_path = "C:\\Users\\uthmn\\PycharmProjects\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\inductor\\results\\values\\winding_1\\Losses_turn_1.dat"

# Analyze the files
get_avg_file(j2F_path, 0.0001, 'j2F_losses')
get_avg_file(turn_losses_path, 0.000001, 'turn_losses')
