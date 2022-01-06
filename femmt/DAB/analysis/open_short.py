import numpy as np

dec = 2


# Open Circuit:

# long SP
# L11 = 830  # 2nd open: [L] = uH
# L22 = 76   # 1st open: [L] = uH

# short SP
L11 = 683  # 2nd open: [L] = uH
L22 = 67  # 1st open: [L] = uH



# Closed Circuit:

# Long SP
# Lk1 = 105  # [L] = uH

# Short SP
Lk1 = 84  # [L] = uH
# Lk2 = 8  # unused ... [L] = uH


# From Lk1 = (1 - k^2) * L11:
k = np.sqrt(1 - Lk1 / L11)

# Primary side concentrated Ls
Ls = np.round(Lk1, dec)

M = np.round(np.sqrt(k**2 * L11 * L22), dec)

n = np.round(k * np.sqrt(L11 / L22), dec)

Lh = np.round(k**2 * L11, dec)

print(f"Primary side concentrated "
      f"{Ls=}uH\n"
      f"{n=}\n"
      f"{Lh}uH\n"
      f"---\n"
      f"{k=}\n")

# Mprim = L11-Lk1 = 5
# Msec  = L22-Lk2 =