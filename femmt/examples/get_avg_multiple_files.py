import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# the losses are extracted and then plotted
# the plot is not precise until now
# the window size and threshold can be controlled to detect the steady state region, so it is a static code
# just a try to detect the steady state region and to get avg.
# it is approperiate for 200 step time

# Process function
def get_avg(file_path, column_name, threshold):
    # Read the file
    times = []
    losses = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        line_values = line.split()
        if len(line_values) == 2:
            times.append(float(line_values[0]))
            losses.append(float(line_values[1]))

    # Convert the lists to a DataFrame
    df = pd.DataFrame({
        'Time': times,
        column_name: losses
    })

    # Calculate the rolling average
    window_size = 250
    df['rolling_avg'] = df[column_name].rolling(window_size, min_periods=1).mean()

    differences = df['rolling_avg'].diff().abs()

    # Find the first point where the difference is smaller than the threshold
    constant_index = np.where(differences < threshold)[0]

    if len(constant_index) > 0:
        steady_start = df.loc[constant_index[0], 'Time']
        steady_end = df['Time'].iloc[-1]  # Assuming steady state lasts till the end of the data
        print(f"Steady state starts at time {steady_start} and ends at {steady_end}")
    else:
        print("The rolling average does not reach a constant value within the given threshold.")

    half_period = (steady_end - steady_start) / 2

    # Get the data for steady state
    steady_df = df[(df['Time'] >= steady_start) & (df['Time'] <= steady_end)]
    full_period_avg = steady_df[column_name].mean()
    half_period_avg = steady_df[steady_df['Time'] <= steady_start + half_period][column_name].mean()

    print(f'Average in steady state for full period ({column_name}):', full_period_avg)
    print(f'Average in steady state for half period ({column_name}):', half_period_avg)

    plt.figure(figsize=(10, 6))
    plt.plot(df['Time'], df[column_name], label='Raw data')
    plt.plot(df['Time'], df['rolling_avg'], label='Rolling average')
    plt.plot(steady_df['Time'], steady_df[column_name], label='Steady state')
    plt.axvline(x=steady_start, color='r', linestyle='--', label='Start of steady state')
    plt.axvline(x=steady_end, color='g', linestyle='--', label='End of steady state')
    plt.axhline(y=full_period_avg, color='b', linestyle='--', label=f'Average for the full interval ({column_name})')
    plt.axhline(y=half_period_avg, color='m', linestyle='--', label=f'Average for the half interval ({column_name})')
    plt.xlabel('Time (s)')
    plt.ylabel(column_name)
    plt.title(f'Time vs {column_name}')
    plt.legend()
    plt.show()

# List of files with full paths, corresponding column names, and thresholds
file_data = [
    {
        "filepath": "C:\\Users\\uthmn\\PycharmProjects\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\inductor\\results\\values\\j2F_1.dat",
        "column_name": "j2F_losses_1",
        "threshold": 0.0001
    },
    {
        "filepath": "C:\\Users\\uthmn\\PycharmProjects\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\inductor\\results\\values\\winding_1\\Losses_turn_1.dat",
        "column_name": "Losses_turn_1",
        "threshold": 0.000001
    }
]

# Loop through the file_data and process each file
for data in file_data:
    get_avg(data["filepath"], data["column_name"], data["threshold"])
