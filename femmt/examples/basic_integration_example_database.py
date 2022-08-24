import materialdatabase as mdb
import matplotlib.pyplot as plt
database = mdb.MaterialDatabase()

# -----Enter the freq and Temp-----------
mu_real = database.get_data_at_working_point(T=60, f=160000, material_name="N95")[0]
mu_imag = database.get_data_at_working_point(T=60, f=160000, material_name="N95")[1]
b = database.get_data_at_working_point(T=60, f=160000, material_name="N95")[2]
print(mu_imag)
print(b)
"""
plt.plot(b, mu_real)
plt.ylabel("Real part of permeability")
plt.xlabel('B in T')
plt.title("Real part of permeability")
plt.show()
"""

plt.plot(b, mu_imag)
plt.ylabel("Imaginary part of permeability")
plt.xlabel('B in T')
plt.title("Imaginary part of permeability")
plt.show()
