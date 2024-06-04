import time
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
import seaborn as sns


def load_or_create_df():
    df = pd.read_pickle('saved_dataframe_turns.pkl')
    return df


class MultiOutputRegressionModel(nn.Module):
    def __init__(self, input_size, output_size, dropout_rate=0.5):
        super(MultiOutputRegressionModel, self).__init__()
        self.layer1 = nn.Linear(input_size, 512)

        self.activation1 = nn.LeakyReLU(negative_slope=0.3)
        self.layer2 = nn.Linear(512, 512)
        self.activation2 = nn.LeakyReLU(negative_slope=0.3)
        self.layer3 = nn.Linear(512, 512)
        self.activation3 = nn.LeakyReLU(negative_slope=0.3)
        self.layer4 = nn.Linear(512, 512)
        self.activation4 = nn.LeakyReLU(negative_slope=0.3)
        self.layer6 = nn.Linear(512, 512)
        self.activation6 = nn.LeakyReLU(negative_slope=0.3)
        self.layer5 = nn.Linear(512, output_size)

    def forward(self, x):
        x = self.layer1(x)
        x = self.activation1(x)
        x = self.layer2(x)
        x = self.activation2(x)
        x = self.layer3(x)
        x = self.activation3(x)
        x = self.layer4(x)
        x = self.activation4(x)
        x = self.layer6(x)
        x = self.activation6(x)

        x = self.layer5(x)
        return x


class CustomLoss(nn.Module):
    def __init__(self, epsilon = 1e-6):
        super(CustomLoss, self).__init__()
        self.epsilon = epsilon

    def forward(self, predicted, target):
        target = target + self.epsilon
        absolute_error = torch.abs(predicted - target)
        percentage_error = absolute_error / target
        mape = torch.mean(percentage_error) * 100
        return mape


def load_or_create_model(x_train, y_train, x_val, y_val, input_size, output_size):
    try:
        model = MultiOutputRegressionModel(input_size, output_size)
        optimizer = optim.Adam(model.parameters(), lr=0.00001)

        checkpoint = torch.load('basic_model_turn.pth')
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        epoch = checkpoint['epoch']
        loss = checkpoint['loss']

        return model, optimizer, epoch, loss

    except:
        # Trainings-, Validierungs- und Testdaten laden
        train_dataset = torch.utils.data.TensorDataset(x_train, y_train)
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=64, shuffle=True)

        val_dataset = torch.utils.data.TensorDataset(x_val, y_val)
        val_loader = torch.utils.data.DataLoader(dataset=val_dataset, batch_size=64, shuffle=False)

        model = MultiOutputRegressionModel(input_size, output_size)

        # Verlustfunktion und Optimierer definieren
        criterion = nn.L1Loss()
        #criterion = nn.CustomLoss()

        optimizer = optim.Adam(model.parameters(), lr=0.00001)

        print("Start")

        best_val_loss = float('inf')
        patience = 50  # Anzahl der Epochen ohne Verbesserung, bevor das Training gestoppt wird
        num_epochs = 20000
        early_stop_counter = 0
        train_loss_values = []
        val_loss_values = []

        for epoch in range(num_epochs):
            # Trainingsschleife
            model.train()
            train_losses = []
            for inputs, targets in train_loader:
                outputs = model(inputs)
                loss_i = criterion(outputs, targets)
                optimizer.zero_grad()
                loss_i.backward()
                optimizer.step()
                train_losses.append(loss_i.item())

            avg_train_loss = sum(train_losses) / len(train_losses)
            train_loss_values.append(avg_train_loss)

            # Validierungsschleife
            model.eval()
            val_losses = []
            with torch.no_grad():
                for val_inputs, val_targets in val_loader:
                    val_outputs = model(val_inputs)
                    val_loss_i = criterion(val_outputs, val_targets)
                    val_losses.append(val_loss_i.item())

            avg_val_loss = sum(val_losses) / len(val_losses)
            val_loss_values.append(avg_val_loss)

            # Gesamtverlust ausgeben
            print(f'Epoch [{epoch + 1}/{num_epochs}], Total Loss: {avg_val_loss:.10f}')

            # Überprüfe, ob der Validierungsverlust verbessert wurde
            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                early_stop_counter = 0
            else:
                early_stop_counter += 1
                if early_stop_counter >= patience:
                    print(f'Early stopping at epoch {epoch + 1} with best validation loss: {best_val_loss:.10f}')
                    break

        # Pfad zum Speichern des Modells
        model_path = 'basic_model_turn.pth'

        # Modell und Optimizer speichern
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': avg_val_loss,
        }, model_path)

        return model, optimizer, epoch, avg_val_loss

def magnetic_field(df):
    mu_0 = 4 * 3.14159265358979323846e-7
    mu_r = 3000

    current = df['Current']
    turns = df['Turns']
    corediameter = df['Core Diameter']
    windowheight = df['Window Height']
    windowwidth = df['Window Width']
    number_airgaps = df['Air Gap Number']
    airgap_height = df['Air Gap Height']
    corearea = ((df['Core Diameter']/2)**2) * np.pi

    # length of core legs
    l_left_top = corediameter/2 + windowwidth + corediameter/2
    l_left_mid = corediameter/2 + windowheight + corediameter/2
    l_left_bot = corediameter/2 + windowwidth + corediameter/2

    l_left = l_left_top + l_left_mid + l_left_bot

    l_mid = corediameter/2 + windowheight + corediameter/2 - (number_airgaps * airgap_height)

    l_airgap = number_airgaps * airgap_height

    reluctance_left = l_left / (mu_0 * mu_r * corearea)

    reluctance_mid_core = l_mid / (mu_0 * mu_r * corearea)

    reluctance_mid_air = l_airgap / (mu_0 * corearea)

    reluctance_mid = reluctance_mid_core + reluctance_mid_air

    # Kernreluktanz links und recht sind gleich groß und parallel "geschaltet", weswegen die gesamte Reluktanz die Hälfte der linken beträgt
    reluctance_right_and_left = 0.5 * reluctance_left

    total_reluctance = reluctance_mid + reluctance_right_and_left

    magnetic_field = (current * turns) / (total_reluctance * corearea)

    return magnetic_field

# Lade oder erstelle DataFrame
df = load_or_create_df()
print(df.shape)

# DataFrame Filter
df = df[magnetic_field(df) < 0.3]
print(df.shape)

# DataFrame Filter
#df = df[df['Current'] < 1]
print(df.shape)

# DataFrame Filter
#df = df[df['Air Gap Number'] == 1]
#df = df[df['Turns'] * df['Current'] > 1]
df = df[df['P'] > 1e-6]
print(df.shape)

# DataFrame samplen
df = df.sample(n=150000, random_state=42)
print(df.shape)

# One-Hot-Encoding
df = pd.get_dummies(df, columns=['Conductor Arrangement'], prefix='Conductor_Arrangement')

# Vorbeugen von eventuellen Fehlern durch negativen Bereich
df['FluxIM'] = df['FluxIM'].abs()


# Inputs
features = df[["Frequency", "Current", "Turns", "Temperature", "Core Diameter", "Core Height", "Window Height",
               "Window Width", "Air Gap Number","Air Gap 1 Position", "Air Gap 2 Position","Air Gap 3 Position", "Air Gap Height",
               "Core Insulation Left Core", "Conductor Radius", "Conductor_Arrangement_Hexagonal", "Conductor_Arrangement_Square", "Conductor_Arrangement_SquareFullWidth"]].copy()

# Feature Engineering
features.loc[:, 'CoreArea'] = (features["Core Diameter"]/2 * features["Core Diameter"]/2 * 3.14159265358979323846)
features.loc[:, 'MFF'] = features["Turns"] * features["Current"]
features.loc[:, 'FrequencySquared'] = features['Frequency'] ** 2
features.loc[:, 'CurrentSquared'] = features['Current'] ** 2
features.loc[:, 'TurnsSquared'] = features['Turns'] ** 2
features.loc[:, 'Air Gap Height Squared'] = features['Air Gap Height'] ** 2
features.loc[:, 'CoreDiameterSquared'] = features['Core Diameter'] ** 2

# Outputs mit Windungsverluste (Auflistung der ausgeschlossenen Parameter)
target = df.drop(["Frequency", "Current", "Turns", "Temperature", "Core Diameter", "Core Height", "Window Height",
                    "Window Width", "Air Gap Height", "Core Insulation Left Core", "Inner Winding Insulation",
                    "Conductor Radius", "Fläche", "Air Gap Number","Air Gap 1 Position",
                    "Air Gap 2 Position","Air Gap 3 Position", "Conductor_Arrangement_Hexagonal", "Conductor_Arrangement_Square", "Conductor_Arrangement_SquareFullWidth",
                    "Core Part 1 Eddy Losses", "Core Part 1 Hyst Losses",
                    "Core Part 2 Eddy Losses", "Core Part 2 Hyst Losses",
                    "Core Part 3 Eddy Losses", "Core Part 3 Hyst Losses",
                    "Winding Window Split", "Core Part 1 Losses", "Core Part 2 Losses", "Core Part 3 Losses", "Volumen",
                    "Total Losses", "Total Core Losses", "Flux over Current RE", "Flux over Current IM", "S"], axis = 1)

# Outputs ohne Windungsverluste als Array
target_arr = ["P", "Q", "Winding Loss", "Total Core Eddy Losses", "Total Core Hyst Losses", "VoltageRE", "VoltageIM",
              "FluxRE", "FluxIM"]

# Datensatz in Trainings-, Validierungs- und Testdaten aufteilen
X_train, X_temp, y_train, y_temp = train_test_split(features, target, test_size=0.3, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42)

# Zur Kontrolle von Skallierung und Deskallierung
y_test_real = y_test
print(tabulate(y_test_real.head(9), headers='keys', tablefmt='pretty'))

# Spalten für Log-Skalierung und anschließender Min-Max-Skalierung
log_scaled_columns = ["P", "Q", "Winding Loss",
                      "Total Core Eddy Losses", "Total Core Hyst Losses",]

# Restliche Spalten NUR Min-Max-Skalierung
minmax_scaled_columns = [col for col in y_train.columns if col not in log_scaled_columns]

# Skalierurng der Eingangsdaten
scaler_features = MinMaxScaler()
X_train_scaled = scaler_features.fit_transform(X_train)
X_val_scaled = scaler_features.transform(X_val)
X_test_scaled = scaler_features.transform(X_test)

# Log-Skalierung für einige Ausgangswerte durchführen
y_train_log = np.log1p(y_train[log_scaled_columns])
y_val_log = np.log1p(y_val[log_scaled_columns])
y_test_log = np.log1p(y_test[log_scaled_columns])

# # Zusammenführen der unskalierten und log skalierten Daten
y_train = np.concatenate((y_train_log, y_train[minmax_scaled_columns]), axis=1)
y_val = np.concatenate((y_val_log, y_val[minmax_scaled_columns]), axis=1)
y_test = np.concatenate((y_test_log, y_test[minmax_scaled_columns]), axis=1)

# Min-Max-Skalierung der Ausgangsdaten
scaler_target_minmax = MinMaxScaler()
y_train_scaled = scaler_target_minmax.fit_transform(y_train)
y_val_scaled = scaler_target_minmax.transform(y_val)
y_test_scaled = scaler_target_minmax.transform(y_test)

# Konvertiere Daten in PyTorch-Tensoren
X_train_tensor = torch.tensor(X_train_scaled, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train_scaled, dtype=torch.float32)
X_val_tensor = torch.tensor(X_val_scaled, dtype=torch.float32)
y_val_tensor = torch.tensor(y_val_scaled, dtype=torch.float32)
X_test_tensor = torch.tensor(X_test_scaled, dtype=torch.float32)
y_test_tensor = torch.tensor(y_test_scaled, dtype=torch.float32)


# Modell initialisieren und trainieren
input_size = X_train.shape[1]
output_size = y_train.shape[1]

model, optimizer, epochs, loss = load_or_create_model(X_train_tensor, y_train_tensor, X_val_tensor, y_val_tensor,
                                                      input_size, output_size)