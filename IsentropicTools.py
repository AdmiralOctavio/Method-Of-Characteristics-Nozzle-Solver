import numpy as np
from scipy.optimize import brentq
from scipy.optimize import fsolve

def PM(M, g): 
    if M >= 1:
        return np.sqrt((g+1)/(g-1)) * np.arctan(np.sqrt((g-1)/(g+1) *(M**2-1))) - np.arctan(np.sqrt(M**2-1))
    else: return 0

def mu(M):
    if M < 1.0:
        return np.nan 
    return np.arcsin(1/M) 

def PMinv(nu_target, gamma, M_guess=1.1, M_min=1.0, M_max=10.0):
    def f(M, nu_target, gamma):
        return PM(M, gamma) - nu_target
    try:
        return brentq(f, M_min, M_max, args=(nu_target, gamma))
    except ValueError:
        print(f"Warning: Could not find M for nu={nu_target}. Returning M_max.")
        return M_max 

def Pressure(P0, g, M):
    return P0 / (1 + (g-1)/2 * M**2)**(g / (g-1))

def Temperature(T0, g, M):
    return T0 / (1 + (g-1)/2 * M**2)

def LocalSoS(g, Rs, T):
    return np.sqrt(g * Rs * T)

def AreaRatio(g, M):
    term1 = (2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)
    exponent = (g + 1) / (2 * (g - 1))
    return (1 / M) * term1**exponent

def AreaRatioInverse(Ar_target, g, M_guess=1.1, M_min=1.0, M_max=10.0):
    def f(M, Ar_target, g):
        return AreaRatio(g, M) - Ar_target
    try:
        if f(M_min, Ar_target, g) * f(M_max, Ar_target, g) > 0:
            return fsolve(f, M_guess, args=(Ar_target, g))[0]
        return brentq(f, M_min, M_max, args=(Ar_target, g))
    except ValueError:
        print(f"Warning: Could not find M for Area Ratio={Ar_target}. Returning M_max.")
        return M_max