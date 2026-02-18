import datetime

# --- EXTENDED THEORY DATABASE ---
MODEL_THEORY = {
    "Michaelis-Menten": (
        "THEORY AND MECHANISM:\n"
        "The Michaelis-Menten model describes the kinetics of enzymes exhibiting simple saturation. "
        "It assumes the formation of a specific Enzyme-Substrate (ES) complex that breaks down to form product. "
        "This model implies a hyperbolic relationship between velocity and substrate concentration.\n\n"
        "KEY PARAMETERS:\n"
        "- Vmax: The theoretical maximum velocity at infinite substrate concentration.\n"
        "- Km: The substrate concentration yielding half-maximal velocity. It is an inverse measure of affinity."
    ),
    "Competitive": (
        "THEORY AND MECHANISM:\n"
        "In Competitive Inhibition, the inhibitor resembles the substrate and binds reversibly to the active site "
        "of the free enzyme. This prevents the substrate from binding. However, because the binding is reversible, "
        "high concentrations of substrate can out-compete the inhibitor.\n\n"
        "DIAGNOSTIC SIGNATURE:\n"
        "- Vmax remains unchanged (infinite substrate eventually displaces all inhibitor).\n"
        "- The apparent Km increases (lower affinity)."
    ),
    "Non-Competitive": (
        "THEORY AND MECHANISM:\n"
        "In Non-Competitive Inhibition, the inhibitor binds to an allosteric site distinct from the active site. "
        "It can bind with equal affinity to both the free enzyme (E) and the Enzyme-Substrate complex (ES). "
        "Binding does not prevent substrate binding, but it prevents catalysis.\n\n"
        "DIAGNOSTIC SIGNATURE:\n"
        "- Vmax decreases (active enzyme concentration is effectively reduced).\n"
        "- Km remains unchanged."
    ),
    "Uncompetitive": (
        "THEORY AND MECHANISM:\n"
        "Uncompetitive Inhibitors bind ONLY to the Enzyme-Substrate (ES) complex, not the free enzyme. "
        "This often occurs in multi-substrate reactions or when substrate binding creates an allosteric site. "
        "Binding of the inhibitor traps the substrate on the enzyme.\n\n"
        "DIAGNOSTIC SIGNATURE:\n"
        "- Both Vmax and Km decrease by the same factor.\n"
        "- Lineweaver-Burk plots show parallel lines."
    ),
    "Mixed": (
        "THEORY AND MECHANISM:\n"
        "Mixed Inhibition is the general case where the inhibitor binds to both the free enzyme (E) and the "
        "Enzyme-Substrate complex (ES), but with unequal affinities. This is defined by the parameter Alpha.\n\n"
        "INTERPRETING ALPHA:\n"
        "- Alpha > 1: Inhibitor prefers Free Enzyme (Competitive-like).\n"
        "- Alpha < 1: Inhibitor prefers ES Complex (Uncompetitive-like).\n"
        "- Alpha = 1: Non-Competitive."
    ),
    "Hill (Allosteric)": (
        "THEORY AND MECHANISM:\n"
        "The Hill Equation describes enzymes with multiple binding sites that exhibit cooperativity (e.g., Hemoglobin). "
        "Binding of substrate to one site facilitates (or hinders) binding to subsequent sites.\n\n"
        "INTERPRETING 'n':\n"
        "- n > 1: Positive Cooperativity (Sigmoidal curve).\n"
        "- n < 1: Negative Cooperativity.\n"
        "- n = 1: No Cooperativity (Reduces to Michaelis-Menten)."
    ),
    "Substrate Inhibition": (
        "THEORY AND MECHANISM:\n"
        "Substrate Inhibition (Haldane Kinetics) occurs when a second substrate molecule binds to the ES complex, "
        "forming an inactive ESS complex. This creates a characteristic velocity peak followed by a decline at high "
        "substrate concentrations. This is often observed in enzymes with narrow active site tunnels."
    )
}

def generate_methods_report(best_result, dataset_size, c_unit, v_unit, weighting_mode, stats_res, derived_stats=None, ci_data=None):
    model_name = best_result['model']
    aic = best_result['aic']
    r2 = best_result['r_squared']
    params = best_result['parameters']
    
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    
    lines = []
    def add(text=""): lines.append(text)
    def section(title):
        add()
        add("=" * 70)
        add(f"  {title}")
        add("=" * 70)
        add()

    # --- HEADER ---
    add(f"BIOKINETICPY v1.0 - COMPREHENSIVE ANALYSIS REPORT")
    add(f"Generated: {timestamp}")
    add("-" * 70)
    
    # --- 1. OVERVIEW ---
    section("1. EXPERIMENTAL OVERVIEW")
    add(f"Dataset Size:           {dataset_size} data points")
    add(f"Concentration Unit:     {c_unit}")
    add(f"Velocity Unit:          {v_unit}")
    add(f"Regression Algorithm:   Levenberg-Marquardt (Damped Least Squares)")
    add(f"Weighting Strategy:     {weighting_mode}")
    add("Model Selection:        Akaike Information Criterion (AIC)")
    
    if derived_stats and derived_stats.get('et_val', 0) > 0:
        add(f"Enzyme Concentration:   {derived_stats['et_val']} {derived_stats['et_unit']}")

    # --- 2. MODEL SELECTION ---
    section("2. MODEL SELECTION & THEORY")
    add(f"SELECTED MODEL: {model_name.upper()}")
    add("-" * 30)
    add(MODEL_THEORY.get(model_name, "Standard kinetic model."))
    add()
    
    # --- 3. PARAMETERS ---
    section("3. KINETIC PARAMETERS")
    add("Values are reported with 95% Confidence Intervals (if computed).")
    add("-" * 70)
    
    for k, v in params.items():
        name = k.replace("vmax", "Vmax").replace("km", "Km").replace("ki", "Ki").replace("alpha", "Alpha").replace("khalf", "Khalf").replace("ksi", "Ksi")
        
        # Build CI String
        ci_str = ""
        if ci_data and k in ci_data:
            low, high = ci_data[k]
            # Calculate % error
            if v != 0:
                rel_err = (abs(high - low) / 2) / v * 100
                ci_str = f"  [95% CI: {low:.4f} - {high:.4f}] (±{rel_err:.1f}%)"
            else:
                ci_str = f"  [95% CI: {low:.4f} - {high:.4f}]"
        
        add(f"{name:<12}: {v:.4f} {ci_str}")
    
    add()
    if derived_stats and derived_stats.get('kcat', 0) > 0:
        add("DERIVED CONSTANTS (Based on [E]t)")
        add("-" * 30)
        add(f"kcat (Turnover Number):     {derived_stats['kcat']:.4f} s^-1")
        add("   -> Number of catalytic cycles per active site per second.")
        add(f"kcat/Km (Specificity):      {derived_stats['efficiency']:.2e} M^-1 s^-1")
        add("   -> The second-order rate constant for E + S reaction. A measure of catalytic efficiency.")

    # --- 4. STATISTICS ---
    section("4. STATISTICAL VALIDATION")
    add("Goodness-of-Fit Metrics:")
    add(f"- R-Squared (R²): {r2:.4f}")
    add("  (Indicates the proportion of variance explained by the model. >0.95 is ideal)")
    add(f"- AIC Score:      {aic:.2f}")
    add("  (Used to compare models. The lowest score indicates the most parsimonious fit.)")
    add()
    
    add("Residual Analysis (Post-Hoc Tests):")
    sw_p = stats_res.get('shapiro_p', 1.0)
    sw_status = "PASSED" if sw_p > 0.05 else "WARNING"
    add(f"1. Normality (Shapiro-Wilk): {sw_status} (p={sw_p:.3f})")
    if sw_p > 0.05:
        add("   -> The residuals appear to be normally distributed (Gaussian noise).")
    else:
        add("   -> The residuals deviate from normality. Check for outliers or systematic errors.")
        
    rz = abs(stats_res.get('runs_z', 0))
    run_status = "PASSED" if rz < 1.96 else "WARNING"
    add(f"2. Randomness (Runs Test):   {run_status} (Z-score={rz:.2f})")
    if rz < 1.96:
        add("   -> No systematic deviation detected. The model captures the trend well.")
    else:
        add("   -> Systematic deviation detected (e.g., residuals form a U-shape). The model may be incorrect.")

    add()
    add("=" * 70)
    add("Generated by BioKineticPy - Open Source Kinetic Analysis")
    
    return "\n".join(lines)

def generate_thermo_report(arr_res, eyr_res):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    lines = []
    def add(text=""): lines.append(text)
    
    add(f"BIOKINETICPY v1.0 - THERMODYNAMIC PROFILE")
    add(f"Generated: {timestamp}")
    add("=" * 70)
    add()
    
    add("This report analyzes the temperature dependence of the reaction rate using")
    add("Arrhenius and Transition State Theory (Eyring) models.")
    add()
    
    add("1. ARRHENIUS ANALYSIS")
    add("-" * 30)
    add(f"Activation Energy (Ea):  {arr_res['Ea']/1000:.2f} kJ/mol")
    add("   INTERPRETATION: This is the energy barrier that substrates must overcome")
    add("   to react. A higher Ea implies the reaction is more sensitive to temperature.")
    add()
    add(f"Pre-Exponential (A):     {arr_res['A']:.2e}")
    add("   INTERPRETATION: Represents the theoretical collision frequency.")
    add()
    add(f"Linear Fit R²:           {arr_res['r_squared']:.4f}")
    add()
    
    add("2. EYRING ANALYSIS (Transition State Theory)")
    add("-" * 30)
    add(f"Enthalpy of Activation (ΔH‡):  {eyr_res['dH']/1000:.2f} kJ/mol")
    add("   -> Reflects the bond breaking/making energy required to reach the transition state.")
    add()
    add(f"Entropy of Activation (ΔS‡):   {eyr_res['dS']:.2f} J/mol·K")
    add("   INTERPRETATION:")
    if eyr_res['dS'] < 0:
        add("   -> Negative Value: The transition state is MORE ordered than the ground state.")
        add("      This suggests a highly specific orientation is required (associative mechanism).")
    else:
        add("   -> Positive Value: The transition state is LESS ordered (dissociative mechanism).")
    add()
    add(f"Gibbs Free Energy (ΔG‡ 25°C):  {eyr_res['dG_25']/1000:.2f} kJ/mol")
    add("   -> The total energy barrier determining the rate at room temperature.")
    add()
    add("=" * 70)
    return "\n".join(lines)

def generate_simulation_report(model_name, params, c_unit, v_unit):
    """Generate a formatted simulation report consistent with the analysis reports."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    lines = []
    def add(text=""): lines.append(text)

    add(f"BIOKINETICPY v1.0 - SIMULATION REPORT")
    add(f"Generated: {timestamp}")
    add("=" * 70)
    add()

    add("1. SIMULATION OVERVIEW")
    add("-" * 30)
    add(f"Selected Model:         {model_name}")
    add(f"Concentration Unit:     {c_unit}")
    add(f"Velocity Unit:          {v_unit}")
    add()

    add("2. INPUT PARAMETERS")
    add("-" * 30)
    param_labels = {
        "vmax": "Vmax (Max Velocity)",
        "km": "Km (Michaelis Constant)",
        "ki": "Ki (Inhibition Constant)",
        "alpha": "Alpha (Binding Selectivity)",
        "khalf": "Khalf (Half-Saturation)",
        "n": "Hill Coefficient (n)",
        "ksi": "Ksi (Substrate Inhibition)",
        "i": "Inhibitor Concentration [I]",
    }
    for k, v in params.items():
        label = param_labels.get(k, k)
        if v != 0:
            add(f"  {label:<30}: {v:.4f}")
    add()

    add(MODEL_THEORY.get(model_name, "Standard kinetic model."))
    add()
    add("=" * 70)
    add("Generated by BioKineticPy - Open Source Kinetic Analysis")
    return "\n".join(lines)