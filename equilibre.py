
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ==================================================
# 1. DONNEES DE L'EXERCICE DWSIM
# ==================================================

F = 60.0                 # Debit molaire total (mol/h)
zM = 0.40                # Fraction molaire methanol
zW = 0.60                # Fraction molaire eau
T = 80.0                 # Temperature (°C)
P = 760.0                # Pression (mmHg) = 1 atm

# Resultats de reference transmis depuis DWSIM
L_DWSIM = 35.4
V_DWSIM = 24.60
xM_DWSIM = 0.244
xW_DWSIM = 0.756
yM_DWSIM = 0.624
yW_DWSIM = 0.376
beta_DWSIM = 0.41

# Constantes K deduites des compositions DWSIM
K_M = yM_DWSIM / xM_DWSIM
K_W = yW_DWSIM / xW_DWSIM

# ==================================================
# 2. EQUATION DE RACHFORD-RICE
# ==================================================

def rachford_rice(beta):
    return (
        zM * (K_M - 1) / (1 + beta * (K_M - 1))
        + zW * (K_W - 1) / (1 + beta * (K_W - 1))
    )

beta = brentq(rachford_rice, 0, 1)

# Calcul des debits liquide et vapeur
V = beta * F
L = F - V

# Calcul des compositions
xM = zM / (1 + beta * (K_M - 1))
xW = zW / (1 + beta * (K_W - 1))

yM = K_M * xM
yW = K_W * xW

# ==================================================
# 3. AFFICHAGE DES RESULTATS
# ==================================================

print("===== FLASH EAU-METHANOL : RACHFORD-RICE =====")
print(f"Temperature : {T:.2f} °C")
print(f"Pression : {P/760:.2f} atm")
print(f"K methanol : {K_M:.5f}")
print(f"K eau : {K_W:.5f}")
print(f"Fraction vaporisee beta : {beta:.5f}")
print(f"Debit liquide : {L:.3f} mol/h")
print(f"Debit vapeur : {V:.3f} mol/h")
print(f"x methanol : {xM:.5f}")
print(f"x eau : {xW:.5f}")
print(f"y methanol : {yM:.5f}")
print(f"y eau : {yW:.5f}")

print("\n===== BILANS DE MATIERE =====")
print(f"xM + xW = {xM + xW:.6f}")
print(f"yM + yW = {yM + yW:.6f}")
print(f"Bilan methanol = {L*xM + V*yM:.4f} mol/h")
print(f"Bilan eau = {L*xW + V*yW:.4f} mol/h")

print("\n===== COMPARAISON DWSIM =====")
print(f"beta DWSIM : {beta_DWSIM:.3f}")
print(f"L DWSIM : {L_DWSIM:.2f} mol/h")
print(f"V DWSIM : {V_DWSIM:.2f} mol/h")
print(f"x methanol DWSIM : {xM_DWSIM:.3f}")
print(f"y methanol DWSIM : {yM_DWSIM:.3f}")

# ==================================================
# 4. DIAGRAMME T-x-y
# Modele ideal : loi de Raoult
# ==================================================

def Psat_M(temp):
    return 10**(8.08097 - 1582.271/(temp + 239.726))

def Psat_W(temp):
    return 10**(8.07131 - 1730.63/(temp + 233.426))

x_values = np.linspace(0, 1, 101)
T_bubble = []
y_values = []

for x in x_values:

    def equation_bulle(temp):
        return x * Psat_M(temp) + (1-x) * Psat_W(temp) - P

    Tb = brentq(equation_bulle, 50, 110)
    T_bubble.append(Tb)

    y = x * Psat_M(Tb) / P
    y_values.append(y)

# ==================================================
# 5. TRACAGE DU GRAPHIQUE
# ==================================================

plt.figure(figsize=(8, 6))

plt.plot(x_values, T_bubble, "b-", label="Courbe de bulle")
plt.plot(y_values, T_bubble, "r-", label="Courbe de rosee")

# Point liquide et point vapeur issus de DWSIM
plt.plot(xM_DWSIM, T, "bo", label="Liquide - DWSIM")
plt.plot(yM_DWSIM, T, "ro", label="Vapeur - DWSIM")

plt.xlabel("Fraction molaire du methanol (xM ou yM)")
plt.ylabel("Temperature (°C)")
plt.title("Diagramme T-x-y : eau-methanol")
plt.xlim(0, 1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()





