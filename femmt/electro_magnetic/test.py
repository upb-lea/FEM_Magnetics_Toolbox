import json
from self import self
import femmt as fmt
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class Database:
    """This is class to call data from the database."""

    def __init__(self):
        self.frequency = None
        self.temp = None



def database_to_pro(frequency: float, temp: int):
    # ---- input for Material data--------
    with open("materials_temp.pro", "w") as file:
        mat = str('N95')
        file.write(f'Material = %s; \n' % mat)
        file.write(f'Frequency = %d; \n' % frequency)
        file.write(f'Temperature = %d; \n' % temp)

    # -----Read the database-------
    with open('data.json') as database:
        m_data = json.load(database)
    freq_list = []
    # print(len(m_data["data"]))
    for i in range(len(m_data["data"])):
        freq_list.append(m_data["data"][i]["frequency"])
    # print(freq_list)
    n = len(freq_list)  # len of array

    # ------Remove Duplicate from freq array------
    def remove(arr, n):
        mp = {i: 0 for i in arr}
        for i in range(n):
            if mp[arr[i]] == 0:
                mp[arr[i]] = 1
                return mp

    freq_list = list(remove(freq_list, n))
    print(freq_list)

    # -----find nearby frequency---------
    def find_nearest(array, value):
        array = np.asarray(array)
        array.sort()
        idx = (np.abs(array - value)).argmin()
        if array[idx] > value:
            return array[idx - 1], array[idx]
        else:
            return array[idx], array[idx + 1]

    result = find_nearest(freq_list, frequency)
    print(result)

    f_l = result[0]
    f_h = result[1]

    # ------find nearby temperature------
    temp_list_l = []
    temp_list_h = []

    for i in range(len(m_data["data"])):
        if m_data["data"][i]["frequency"] == f_l:
            temp_list_l.append(m_data["data"][i]["temperature"])
    for i in range(len(m_data["data"])):
        if m_data["data"][i]["frequency"] == f_h:
            temp_list_h.append(m_data["data"][i]["temperature"])

    temp_list_l = find_nearest(temp_list_l, temp)
    temp_list_h = find_nearest(temp_list_h, temp)
    print(temp_list_l)
    print(temp_list_h)

    # -------get the data----------
    def getdata(variable, f, t_1, t_2):
        for k in range(len(m_data["data"])):
            if m_data["data"][k]["frequency"] == f and m_data["data"][k]["temperature"] == t_1:
                b_1 = m_data["data"][k]["b"]
                mu_real_1 = m_data["data"][k]["mu_real"]
                mu_imag_1 = m_data["data"][k]["mu_imag"]
                t_mu_imag_1 = interp1d(b_1, mu_imag_1)
                t_mu_real_1 = interp1d(b_1, mu_real_1)
            if m_data["data"][k]["frequency"] == f and m_data["data"][k]["temperature"] == t_2:
                b_2 = m_data["data"][k]["b"]
                mu_real_2 = m_data["data"][k]["mu_real"]
                mu_imag_2 = m_data["data"][k]["mu_imag"]
                t_mu_imag_2 = interp1d(b_2, mu_imag_2)
                t_mu_real_2 = interp1d(b_2, mu_real_2)

        # --------linear interpolation at constant freq-------------
        mu_i = []
        mu_r = []
        b_f = [0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]

        for j in range(len(b_f)):
            mu_r.append(
                t_mu_real_1(b_f[j]) + (t_mu_real_2(b_f[j]) - t_mu_real_1(b_f[j])) / (t_2 - t_1) * (variable - t_1))
            mu_i.append(
                t_mu_imag_1(b_f[j]) + (t_mu_imag_2(b_f[j]) - t_mu_imag_1(b_f[j])) / (t_2 - t_1) * (variable - t_1))
        return mu_r, mu_i

    # --------interpolated data at constant freq and nearby temp--------
    interpolate_temp_1 = getdata(temp, f_l, temp_list_l[0], temp_list_l[1])
    interpolate_temp_2 = getdata(temp, f_h, temp_list_h[0], temp_list_h[1])
    print(interpolate_temp_1)
    print(interpolate_temp_2)

    # ------linear interpolation at constant temp and nearby freq-----------------
    b_g = [0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
    f_mu_real_1 = interp1d(b_g, interpolate_temp_1[0])
    f_mu_imag_1 = interp1d(b_g, interpolate_temp_1[1])
    f_mu_real_2 = interp1d(b_g, interpolate_temp_2[0])
    f_mu_imag_2 = interp1d(b_g, interpolate_temp_2[1])
    mu_i_f = []
    mu_r_f = []
    for b in range(len(b_g)):
        mu_r_f.append(
            f_mu_real_1(b_g[b]) + (f_mu_real_2(b_g[b]) - f_mu_real_1(b_g[b])) / (f_h - f_l) * (frequency - f_l))
        mu_i_f.append(
            f_mu_imag_1(b_g[b]) + (f_mu_imag_2(b_g[b]) - f_mu_imag_1(b_g[b])) / (f_h - f_l) * (frequency - f_l))
    print(mu_r_f)
    print(mu_i_f)

    # ------------write the data to the .pro file---------
    with open("materials_temp.pro", "a") as file:
        file.write(f'B = %s; \n' % b_g)
        file.write(f'mu_real = %s; \n' % mu_r_f)
        file.write(f'mu_imag = %s; \n' % mu_i_f)


database_to_pro(frequency=150000, temp=40)
