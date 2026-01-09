import numpy as np

def michaelis_menten(s, vmax, km):
    return (vmax * s) / (km + s)

def competitive_inhibition(s, i, vmax, km, ki):
    km_app = km * (1 + (i / ki))
    return (vmax * s) / (km_app + s)

def uncompetitive_inhibition(s, i, vmax, km, ki):
    term = (1 + (i / ki))
    return (vmax * s) / (km + (s * term))

def non_competitive_inhibition(s, i, vmax, km, ki):
    vmax_app = vmax / (1 + (i / ki))
    return (vmax_app * s) / (km + s)

def mixed_inhibition(s, i, vmax, km, ki, alpha):
    km_app = km * (1 + (i / ki))
    vmax_app = vmax / (1 + (i / (alpha * ki)))
    return (vmax_app * s) / (km_app + s)

# --- ADVANCED MODELS ---

def hill_equation(s, vmax, khalf, n):
    """Sigmoidal Kinetics (Allostery)."""
    s_n = np.power(s, n)
    k_n = np.power(khalf, n)
    return (vmax * s_n) / (k_n + s_n)

def haldane_inhibition(s, i, vmax, km, ksi, ki):
    """
    Substrate Inhibition + Optional Competitive Inhibition.
    v = (Vmax * S) / (Km(1 + I/Ki) + S + S^2/Ksi)
    """
    # Term 1: Km increases due to external competitive inhibitor
    km_app = km * (1 + (i / ki))
    
    # Term 3: Substrate Inhibition term (S^2 / Ksi)
    # Note: Denominator = Km_app + S + S^2/Ksi
    denom = km_app + s + (np.square(s) / ksi)
    return (vmax * s) / denom

MODEL_REGISTRY = {
    "Michaelis-Menten": michaelis_menten,
    "Competitive": competitive_inhibition,
    "Uncompetitive": uncompetitive_inhibition,
    "Non-Competitive": non_competitive_inhibition,
    "Mixed": mixed_inhibition,
    "Hill (Allosteric)": hill_equation,
    "Substrate Inhibition": haldane_inhibition
}