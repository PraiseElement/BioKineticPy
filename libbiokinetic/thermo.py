import numpy as np
from scipy.stats import linregress

# Physical Constants
R = 8.314462618     # Gas Constant (J / mol K)
KB = 1.380649e-23   # Boltzmann Constant (J / K)
H = 6.62607015e-34  # Planck Constant (J s)

def solve_arrhenius(temps, rates):
    """
    Arrhenius Equation: k = A * exp(-Ea / RT)
    Linearization: ln(k) = ln(A) - (Ea/R) * (1/T)
    x = 1/T, y = ln(k)
    Slope = -Ea/R
    Intercept = ln(A)
    """
    x = 1.0 / temps
    y = np.log(rates)
    
    slope, intercept, r_val, p_val, std_err = linregress(x, y)
    
    ea = -slope * R  # J/mol
    a_factor = np.exp(intercept)
    
    return {
        "Ea": ea,
        "A": a_factor,
        "r_squared": r_val**2,
        "reg_x": x,
        "reg_y": y,
        "slope": slope,
        "intercept": intercept
    }

def solve_eyring(temps, rates):
    """
    Eyring Equation: k = (kB*T / h) * exp(-dH/RT) * exp(dS/R)
    Linearization: ln(k/T) = -dH/R * (1/T) + ln(kB/h) + dS/R
    x = 1/T, y = ln(k/T)
    Slope = -dH/R
    Intercept = ln(kB/h) + dS/R
    """
    x = 1.0 / temps
    y = np.log(rates / temps)
    
    slope, intercept, r_val, p_val, std_err = linregress(x, y)
    
    dh = -slope * R # J/mol
    
    ln_kb_h = np.log(KB / H)
    ds = (intercept - ln_kb_h) * R # J/mol K
    
    # Calculate dG at 25C (298.15 K)
    t_ref = 298.15
    dg_25 = dh - (t_ref * ds)
    
    return {
        "dH": dh,
        "dS": ds,
        "dG_25": dg_25,
        "r_squared": r_val**2,
        "reg_x": x,
        "reg_y": y,
        "slope": slope,
        "intercept": intercept
    }