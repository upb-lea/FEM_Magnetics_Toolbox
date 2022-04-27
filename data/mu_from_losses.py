from mu_extractor import *

k = 0.527
alpha = 1.507
beta = 2.382



# 100 kHz hyst loss
frequency = 100000

# Measurements
# b_amp = np.array([30, 40, 50, 60, 70, 80, 90, 100, 200, 300])  # mT
# p_hyst = np.array([4, 8, 15, 25, 40, 60, 80, 100,  700, 2000])  # kW/m^3
# mur =       np.array([2860, 3150, 3200, 3200, 2800, 2700, 2600, 2500, 2400, 2300])  # grad
# phi_mu =    np.array([ 1.4,    2,  4.5, None])  # grad from measurement
# b_amp =     np.array([  30,  100,  200,  300])  # mT
# mur =       np.array([2800, 3400, 3400, 3000])  #
# p_hyst =    np.array([   4,  100,  600, 2000])  # kW/m^3

# Datasheet
# 100 C
b_amp =     np.array([0.025, 0.05,  0.1,  0.2,  0.3])*1000
mur =       np.array([ 3000, 3000, 3000, 3000, 3000])
p_hyst =    np.array([    2,    9,   50,  400, 1000])  # 100 kHz


# 25 C
# b_amp =     np.array([   25,   50,  100,  200,  300])
# mur =       np.array([ 3000, 3000, 3000, 3000, 3000])
# p_hyst =    np.array([  3.5,   20,  100,  500, 1300])  # 100 kHz
# p_hyst =    np.array([    9,   50,  240, 1200, 3000])  # 200 kHz


# used for curve fitting:
# b_amp_extrapolated = np.concatenate((b_amp, max(b_amp)/len(b_amp) + np.linspace(max(b_amp), 350-max(b_amp)/len(b_amp), len(b_amp))))
print(f"{b_amp=}")
print(f"{b_amp/1000=}")

# plot losses
plt.plot(b_amp, p_hyst)
# plt.plot(b_amp, p_hyst_extrapolated)
# plt.plot(b_amp, p_hyst_extrapolated/b_amp)
plt.xlabel("B in mT")
plt.ylabel("P_hyst in kW/m^3")
plt.grid()
plt.show()


# mur_from_losses(b=b_amp, p=p_hyst, f=frequency, phi=phi_mu)
phis = phi_mu_from_losses(b=b_amp, p=p_hyst, f=frequency, mur=mur)
print(f"{phis=}")

# show plot without extrapolation
plt.plot(b_amp, phis, label="phi from p_hyst and mur")
# plt.plot(b_amp, phi_mu, label="phi direct")

# # extrapolate by linear curve fitting
# # Fit polynomial factors
# z_phi = np.polyfit(b_amp, phis, 1)
# z_phi = np.flip(z_phi, 0)
# phi_extrapolated = PolyCoefficients(x=b_amp, coeffs=z_phi)

# # plot extrapolated phi
# plt.plot(b_amp, phi_extrapolated)

# show both phi plots
plt.xlabel("B in mT")
plt.ylabel("phi_h in °")
plt.legend()
plt.grid()
plt.show()

# calc real and imaginary parts
mu_real = np.round(np.cos(np.deg2rad(phis)) * mur)
mu_imag = np.round(np.sin(np.deg2rad(phis)) * mur)
print(list(mu_real))
print(list(mu_imag))
plt.plot(b_amp, mu_real)
plt.plot(b_amp, mu_imag)
plt.grid()
plt.show()

# bring the other data into the right format
# print(b_amp_extrapolated/1000, phi_extrapolated)



