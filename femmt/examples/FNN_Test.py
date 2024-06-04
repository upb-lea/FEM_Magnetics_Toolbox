import time
import torch
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
import seaborn as sns
import FNN_train


def create_acutual_and_predicted_df(actual_values, predicted_values):
    # Dummy Spaltennamen für die ersten paar Spalten
    first_columns = FNN_train.target_arr  # Passe dies entsprechend an

    # Anzahl der restlichen Spalten (nach den ersten paar Spalten)
    num_additional_columns = predicted_values.shape[1] - len(first_columns)

    # Generiere die restlichen Spaltennamen als "Turn 1", "Turn 2", etc.
    additional_columns = [f'Turn {i + 1}' for i in range(num_additional_columns)]

    # Kombiniere die ersten Spaltennamen mit den restlichen Spaltennamen
    all_columns = first_columns + additional_columns

    all_columns_actual = [f'{all_columns[i]}_actual' for i in range(len(all_columns))]
    all_columns_predicted = [f'{all_columns[i]}_predicted' for i in range(len(all_columns))]

    # Erstelle DataFrames für predicted und actual values
    predicted_df = pd.DataFrame(predicted_values_descaled, columns=all_columns_predicted)
    actual_df = pd.DataFrame(actual_values, columns=all_columns_actual)

    return actual_df, predicted_df

def actual_vs_predictet_value(actual_values_descaled, predicted_values_descaled):
    # Erstelle ein Diagramm mit Subplots
    nrows = 3
    ncols = 3

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1, 1))

    for i in range(len(FNN_train.target_arr)):
        # Scatterplot für jeden Ausgabewert
        axs[i // nrows, i % ncols].scatter(actual_values_descaled[:, i], predicted_values_descaled[:, i])

        # Diagonale Linie für den idealen Fall hinzufügen
        axs[i // nrows, i % ncols].plot([min(actual_values_descaled[:, i].min(), predicted_values_descaled[:, i].min()),
                                         max(actual_values_descaled[:, i].max(),
                                             predicted_values_descaled[:, i].max())],
                                        [min(actual_values_descaled[:, i].min(), predicted_values_descaled[:, i].min()),
                                         max(actual_values_descaled[:, i].max(),
                                             predicted_values_descaled[:, i].max())],
                                        linestyle='--', color='gray', label='Ideal Line')

        axs[i // nrows, i % ncols].set_xlabel(f'Actual Values Output {FNN_train.target_arr[i]}')
        axs[i // nrows, i % ncols].set_ylabel(f'Predicted Values Output {FNN_train.target_arr[i]}')
        axs[i // nrows, i % ncols].legend()

    plt.tight_layout(pad=10.0)
    plt.show()

def actual_vs_predicted_turn_loss(actual_values_descaled, predicted_values_descaled):
    for j in range(15):
        i = j + len(FNN_train.target_arr)
        plt.scatter(actual_values_descaled[:, i], predicted_values_descaled[:, i])

        # Diagonale Linie für den idealen Fall hinzufügen
        plt.plot([min(actual_values_descaled[:, i].min(), predicted_values_descaled[:, i].min()),
                  max(actual_values_descaled[:, i].max(),
                      predicted_values_descaled[:, i].max())],
                 [min(actual_values_descaled[:, i].min(), predicted_values_descaled[:, i].min()),
                  max(actual_values_descaled[:, i].max(),
                      predicted_values_descaled[:, i].max())],
                 linestyle='--', color='gray', label='Ideal Line')
        plt.xlabel(f'Actual Values Turn {j}')
        plt.ylabel(f'Predicted Values Turn {j}')
        plt.legend()
        plt.show()



def scatterplot_rel_err():
    for j, column1 in enumerate(FNN_train.features):
        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        axs = axs.flatten()
        for i, column2 in enumerate(percentage_errors_df_columns):
            colors = ['red' if val <0.6 else 'blue' for val in (FNN_train.X_test['Current'] * FNN_train.X_test['Turns'])]
            axs[i].scatter(FNN_train.X_test[column1], percentage_errors_df[column2], c = colors)
            axs[i].set_title(f'{column1} vs {column2}')
            axs[i].set_xlabel(f'{column1}')
            axs[i].set_ylabel(column2)

        plt.tight_layout()
        plt.show()

model = FNN_train.model

# Modellbewertung auf Testdaten
model.eval()
with torch.no_grad():
    # Normalisierung der Input-Daten für die Vorhersage
    predicted = model(FNN_train.X_test_tensor)


# Offset entfernen falls nötig
# predicted = predicted - 1e-6


# Skallierte Zielvariablen
predicted_values = predicted.numpy()
actual_values = FNN_train.y_test_scaled

# Deskalierung der Zielvariablen
predicted_values_inverse = FNN_train.scaler_target_minmax.inverse_transform(predicted_values)
predicted_values_descaled_log = np.expm1(predicted_values_inverse[:, :len(FNN_train.log_scaled_columns)])

actual_values_inverse = FNN_train.scaler_target_minmax.inverse_transform(actual_values)
actual_values_descaled_log = np.expm1(actual_values_inverse[:, :len(FNN_train.log_scaled_columns)])

# Zusammenführen der deskalierten Daten
predicted_values_descaled = np.concatenate((predicted_values_descaled_log, predicted_values_inverse [:, len(FNN_train.log_scaled_columns):]), axis=1)
actual_values_descaled = np.concatenate((actual_values_descaled_log, actual_values_inverse [:, len(FNN_train.log_scaled_columns):]), axis=1)



# TODO def für perdcentage df

# Berechne den durchschnittlichen absoluten Fehler für jeden Input
absolute_errors = abs(predicted_values_descaled - actual_values_descaled)
average_absolute_errors = absolute_errors.mean(axis=0)

# Berechne den durchschnittlichen prozentualen Fehler für jeden Input
percentage_errors = (absolute_errors[:,:len(FNN_train.target_arr)] / actual_values_descaled[:,:len(FNN_train.target_arr)]) * 100

mean_percentage_errors = np.mean(np.abs(percentage_errors), axis=0)
median_percentage_errors = np.median(np.abs(percentage_errors), axis=0)


# Erstelle DataFrame für prozentualen Fehler
percentage_errors_df_columns = [f'{i}_%-Error' for i in FNN_train.target_arr]
percentage_errors_df = pd.DataFrame(percentage_errors, columns= percentage_errors_df_columns)


# Drucke die durchschnittlichen Fehler für jeden Input
for i, variable in enumerate(FNN_train.target_arr):
    print(f'Average Absolute Error ({variable}): {average_absolute_errors[i]:.10f}')
    print(f'Average relative Percentage Error ({variable}): {mean_percentage_errors[i]:.4f}%')
    print(f'Average median Error ({variable}): {median_percentage_errors[i]:.4f}%\n')


# Erstelle DataFrame
actual_df, predicted_df = create_acutual_and_predicted_df(actual_values_descaled, predicted_values_descaled)
percentage_errors_df_columns = [f'{i}_%-Error' for i in FNN_train.target_arr]
percentage_errors_df = pd.DataFrame(percentage_errors, columns= percentage_errors_df_columns)


#Visualisierungen/Plots
actual_vs_predictet_value(actual_values_descaled, predicted_values_descaled)
scatterplot_rel_err()
actual_vs_predicted_turn_loss(actual_values_descaled, predicted_values_descaled)


