import numpy as np
from matplotlib import pyplot as plt
import femmt as fmt


phi = np.linspace(0, 2*np.pi, 100)
t =
i1 = 8*np.sin(phi) + 1*np.sin(2*phi+5)
i2 = 8*np.cos(phi) - 2*np.cos(2*phi+4)
print(np.mean(np.mean(i1)))
print(np.mean(np.mean(i2)))

i1 = i1 - np.mean(i1)
i2 = i2 - np.mean(i2)

res1 = fmt.fft(i1, mode='time')
print(res1)
fmt.fft(i2)
def wizard_n_ffts():
    # TODO: function that takes n ffts with a certain energy limit and createst the missing harmonic results!
    pass


plt.plot(phi, i1)
plt.plot(phi, i2)
plt.show()

print(np.mean(np.mean(i1)))
print(np.mean(np.mean(i2)))

i1.tofile('i1.csv', sep=',')
i2.tofile('i2.csv', sep=',')
