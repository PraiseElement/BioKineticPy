import datetime

# ─────────────────────────────────────────────────────────────────────────────
#  MODEL THEORY DATABASE  (used by analysis and simulation reports)
# ─────────────────────────────────────────────────────────────────────────────
MODEL_THEORY = {
    "Michaelis-Menten": (
        "MECHANISM\n"
        "  The Michaelis-Menten model describes the kinetics of a simple, single-substrate\n"
        "  enzyme reaction: E + S ⇌ ES → E + P. It assumes rapid equilibrium or steady-state\n"
        "  ES formation and predicts a hyperbolic v vs [S] relationship.\n\n"
        "KEY PARAMETERS\n"
        "  Vmax  — Maximum velocity at saturating [S] (≡ kcat × [E]t)\n"
        "  Km    — [S] at half-maximal velocity; inverse measure of substrate affinity"
    ),
    "Competitive": (
        "MECHANISM\n"
        "  The inhibitor is structurally similar to the substrate and binds reversibly to\n"
        "  the active site of the free enzyme (E), forming an EI complex that cannot\n"
        "  bind S. High [S] can completely reverse inhibition.\n\n"
        "DIAGNOSTIC SIGNATURE\n"
        "  ✓ Vmax unchanged        (inhibitor fully outcompeted at high [S])\n"
        "  ✓ Km(apparent) ↑       (reduced apparent affinity for S)\n"
        "  ✓ Lineweaver-Burk      — lines converge on the y-axis\n"
        "  ✓ Dixon plot           — lines converge above the x-axis → x-intercept = −Ki"
    ),
    "Non-Competitive": (
        "MECHANISM\n"
        "  The inhibitor binds to an allosteric site on both the free enzyme (E) and the\n"
        "  ES complex with equal affinity. Substrate and inhibitor binding are independent;\n"
        "  the inhibitor reduces catalytic efficiency but not substrate binding.\n\n"
        "DIAGNOSTIC SIGNATURE\n"
        "  ✓ Vmax ↓               (active enzyme population reduced)\n"
        "  ✓ Km unchanged         (substrate affinity unaffected)\n"
        "  ✓ Lineweaver-Burk      — lines converge on the x-axis"
    ),
    "Uncompetitive": (
        "MECHANISM\n"
        "  The inhibitor binds ONLY to the Enzyme-Substrate (ES) complex, not the free\n"
        "  enzyme. Substrate binding creates or exposes the inhibitor site. Common in\n"
        "  multi-substrate reactions and enzymes with narrow active-site tunnels.\n\n"
        "DIAGNOSTIC SIGNATURE\n"
        "  ✓ Both Vmax and Km decrease by the same factor (ratio Vmax/Km preserved)\n"
        "  ✓ Lineweaver-Burk      — parallel lines (same slope, different intercepts)\n"
        "  ✓ Dixon plot           — lines converge below the x-axis"
    ),
    "Mixed": (
        "MECHANISM\n"
        "  The inhibitor binds to both the free enzyme (E) and the ES complex, but with\n"
        "  different affinities governed by the parameter Alpha (α).\n\n"
        "INTERPRETING ALPHA (α)\n"
        "  α > 1  — Inhibitor prefers free enzyme  (Competitive-like)\n"
        "  α < 1  — Inhibitor prefers ES complex   (Uncompetitive-like)\n"
        "  α = 1  — Equal affinity for E and ES    (Reduces to Non-Competitive)"
    ),
    "Hill (Allosteric)": (
        "MECHANISM\n"
        "  Describes cooperative enzymes with multiple substrate-binding sites (e.g.,\n"
        "  ATCase, Hemoglobin). Binding at one site alters affinity at neighbouring sites.\n\n"
        "INTERPRETING HILL COEFFICIENT (n)\n"
        "  n > 1  — Positive cooperativity  (sigmoidal curve, e.g., allosteric activators)\n"
        "  n < 1  — Negative cooperativity  (flatter than hyperbola)\n"
        "  n = 1  — No cooperativity        (reduces to Michaelis-Menten)"
    ),
    "Substrate Inhibition": (
        "MECHANISM\n"
        "  (Haldane Kinetics) A second substrate molecule binds the active-site cleft of\n"
        "  the ES complex, forming an inactive ESS ternary complex. Gives a characteristic\n"
        "  velocity peak followed by decline at high [S]. Common in narrow active-site\n"
        "  tunnels (e.g., acetylcholinesterase, tyrosinase).\n\n"
        "KEY PARAMETERS\n"
        "  Km    — Concentration for half-maximal activation\n"
        "  Ksi   — Inhibitory substrate dissociation constant (lower = stronger inhibition)"
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
#  RATE EQUATION STRINGS  (shown in the report header)
# ─────────────────────────────────────────────────────────────────────────────
MODEL_EQUATIONS = {
    "Michaelis-Menten":  "v = Vmax·[S] / (Km + [S])",
    "Competitive":       "v = Vmax·[S] / (Km·(1 + [I]/Ki) + [S])",
    "Non-Competitive":   "v = (Vmax/(1 + [I]/Ki)) · [S] / (Km + [S])",
    "Uncompetitive":     "v = Vmax·[S] / ((Km/(1+[I]/Ki)) + [S]·(1 + [I]/Ki))  ⟺  both Vmax & Km ÷ (1+[I]/Ki)",
    "Mixed":             "v = Vmax·[S] / (α·Km·(1+[I]/Ki) + [S]·(1 + [I]/(α·Ki)))",
    "Hill (Allosteric)": "v = Vmax·[S]ⁿ / (Khalf^n + [S]ⁿ)",
    "Substrate Inhibition": "v = Vmax·[S] / (Km + [S] + [S]²/Ksi)",
}

_DIVIDER  = "═" * 72
_SUBDIV   = "─" * 72
_THIN     = "·" * 72


def _section(title):
    return f"\n{_DIVIDER}\n  {title}\n{_DIVIDER}\n"


def generate_methods_report(best_result, dataset_size, c_unit, v_unit,
                             weighting_mode, stats_res,
                             derived_stats=None, ci_data=None,
                             all_results=None):
    """Generate a comprehensive kinetic analysis report (plain-text)."""
    model_name = best_result["model"]
    aicc       = best_result.get("aic", best_result.get("aicc", float("nan")))
    r2         = best_result["r_squared"]
    params     = best_result["parameters"]
    timestamp  = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")

    lines = []
    A = lines.append

    # ── HEADER ──────────────────────────────────────────────────────────────
    A(_DIVIDER)
    A("  BioKineticPy v1.1  ·  COMPREHENSIVE ANALYSIS REPORT")
    A(_DIVIDER)
    A(f"  Generated : {timestamp}")
    A(f"  Software  : BioKineticPy — Open Source Enzyme Kinetics Suite")
    A(f"  URL       : https://github.com/PraiseElement/BioKineticPy")
    A(_THIN)

    # ── 1. EXPERIMENTAL OVERVIEW ────────────────────────────────────────────
    A(_section("1.  EXPERIMENTAL OVERVIEW"))
    A(f"  Dataset size          :  {dataset_size} data points")
    A(f"  Concentration unit    :  {c_unit}")
    A(f"  Velocity unit         :  {v_unit}")
    A(f"  Regression algorithm  :  Levenberg-Marquardt (Damped Least Squares)")
    A(f"  Weighting strategy    :  {weighting_mode}")
    A(f"  Model selection       :  Corrected Akaike Information Criterion (AICc)")
    if derived_stats and derived_stats.get("et_val", 0) > 0:
        A(f"  Enzyme concentration  :  {derived_stats['et_val']} {derived_stats['et_unit']}")

    # ── 2. MODEL SELECTION RANKING ───────────────────────────────────────────
    if all_results:
        A(_section("2.  MODEL SELECTION RANKING  (all fitted models)"))
        A(f"  {'Rank':<5} {'Model':<22} {'AICc':>10} {'R²':>8}  {'ΔAIC':>8}")
        A("  " + _SUBDIV)
        best_aic = min(r["aic"] for r in all_results if r["aic"] is not None)
        ranked   = sorted(all_results, key=lambda r: r["aic"] if r["aic"] is not None else 1e9)
        for rank, res in enumerate(ranked, 1):
            marker = " ◀ BEST FIT" if res["model"] == model_name else ""
            delta  = res["aic"] - best_aic if res["aic"] is not None else float("nan")
            A(f"  {rank:<5} {res['model']:<22} {res['aic']:>10.2f} {res['r_squared']:>8.4f}  {delta:>8.2f}{marker}")
        A()
        A("  ΔAIC interpretation:")
        A("  0–2   → Model has substantial support; difficult to distinguish from best")
        A("  2–7   → Considerably less support")
        A("  > 10  → Essentially no support")
    else:
        A(_section("2.  SELECTED MODEL"))

    # ── 3. MODEL THEORY ──────────────────────────────────────────────────────
    A(_section("3.  SELECTED MODEL — THEORY & MECHANISM"))
    A(f"  Model name  :  {model_name.upper()}")
    eq = MODEL_EQUATIONS.get(model_name, "")
    if eq:
        A(f"  Rate law    :  {eq}")
    A()
    for line in MODEL_THEORY.get(model_name, "Standard kinetic model.").splitlines():
        A(f"  {line}")

    # ── 4. KINETIC PARAMETERS ────────────────────────────────────────────────
    A(_section("4.  KINETIC PARAMETERS"))
    nice = {
        "vmax": "Vmax  (Max velocity)",
        "km":   "Km    (Michaelis const.)",
        "ki":   "Ki    (Inhibition const.)",
        "alpha":"Alpha (Binding selectivity)",
        "khalf":"Khalf (Half-saturation)",
        "n":    "n     (Hill coefficient)",
        "ksi":  "Ksi   (Substrate inhib.)",
    }
    A(f"  {'Parameter':<28} {'Value':>12}   {'95 % CI':^26}  {'±%':>6}")
    A("  " + _SUBDIV)
    for k, v in params.items():
        label  = nice.get(k, k)
        ci_str = ""
        pct    = ""
        if ci_data and k in ci_data:
            lo, hi = ci_data[k]
            ci_str = f"[{lo:>10.4f} – {hi:<10.4f}]"
            if v:
                pct = f"{abs(hi - lo) / 2 / abs(v) * 100:>5.1f} %"
        A(f"  {label:<28} {v:>12.4f}   {ci_str:<26}  {pct}")

    if derived_stats and derived_stats.get("kcat", 0) > 0:
        A()
        A("  DERIVED CATALYTIC CONSTANTS  (requires [E]t input)")
        A("  " + _SUBDIV)
        A(f"  kcat   (turnover number)   :  {derived_stats['kcat']:.4f} s⁻¹")
        A(f"         Catalytic cycles per active site per second at saturating [S].")
        A(f"  kcat/Km (specificity const):  {derived_stats['efficiency']:.3e} M⁻¹ s⁻¹")
        A(f"         Apparent 2nd-order rate constant; upper limit set by diffusion (~10⁸–10⁹ M⁻¹s⁻¹).")

    # ── 5. GOODNESS OF FIT ──────────────────────────────────────────────────
    A(_section("5.  GOODNESS OF FIT"))
    r2_interp = (
        "Excellent — model explains > 99 % of variance" if r2 > 0.99 else
        "Good      — model explains > 97 % of variance" if r2 > 0.97 else
        "Acceptable— model explains > 95 % of variance" if r2 > 0.95 else
        "Poor      — consider alternative models or outlier removal"
    )
    A(f"  R²    :  {r2:.6f}   →  {r2_interp}")
    A(f"  AICc  :  {aicc:.3f}")
    A(f"         (Lower is better. Used internally to rank all 7 models.)")

    # ── 6. RESIDUAL DIAGNOSTICS ─────────────────────────────────────────────
    A(_section("6.  RESIDUAL DIAGNOSTICS"))
    sw_p   = stats_res.get("shapiro_p", 1.0)
    sw_ok  = sw_p > 0.05
    rz     = abs(stats_res.get("runs_z", 0))
    rz_ok  = rz < 1.96

    sw_sym = "✓  PASS" if sw_ok else "✗  WARNING"
    rz_sym = "✓  PASS" if rz_ok else "✗  WARNING"

    A("  Test 1 — Shapiro-Wilk  (residual normality)")
    A(f"    Result  :  {sw_sym}   (p = {sw_p:.4f}  |  threshold p > 0.05)")
    if sw_ok:
        A("    Meaning :  Residuals are consistent with Gaussian noise → fit is valid.")
    else:
        A("    Meaning :  Residuals deviate from normality. Possible causes:")
        A("                • Outliers in the dataset")
        A("                • Systematic error in experimental protocol")
        A("                • Wrong kinetic model selected")
        A("                • Heteroscedastic variance (try 1/v or 1/v² weighting)")
    A()
    A("  Test 2 — Runs Test  (residual randomness)")
    A(f"    Result  :  {rz_sym}   (|Z| = {rz:.3f}  |  threshold |Z| < 1.96)")
    if rz_ok:
        A("    Meaning :  No systematic pattern in residuals → model is appropriate.")
    else:
        A("    Meaning :  Residuals show systematic structure (e.g., U-shape or waves).")
        A("                • The fitted model may be mechanistically incorrect")
        A("                • Consider a model with an extra degree of freedom")

    # ── FOOTER ───────────────────────────────────────────────────────────────
    A()
    A(_THIN)
    A("  BioKineticPy v1.1  ·  github.com/PraiseElement/BioKineticPy")
    A("  Methods: Levenberg-Marquardt NLS · AICc model selection · Bootstrap 95% CI")
    A("  Author : Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DIVIDER)

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
def generate_thermo_report(arr_res, eyr_res):
    """Generate a comprehensive thermodynamic profile report (plain-text)."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    lines = []
    A = lines.append

    A(_DIVIDER)
    A("  BioKineticPy v1.1  ·  THERMODYNAMIC PROFILE REPORT")
    A(_DIVIDER)
    A(f"  Generated : {timestamp}")
    A(f"  Methods   : Arrhenius Analysis  ·  Eyring / Transition-State Theory")
    A(_THIN)

    # ── Arrhenius ────────────────────────────────────────────────────────────
    A(_section("1.  ARRHENIUS ANALYSIS"))
    ea_kj = arr_res["Ea"] / 1000
    A(f"  Activation Energy  Ea  :  {ea_kj:.3f} kJ/mol")
    A(f"  Pre-exponential    A   :  {arr_res['A']:.4e}  s⁻¹")
    A(f"  Linear fit         R²  :  {arr_res['r_squared']:.6f}")
    A()
    A("  INTERPRETATION")
    A(f"  Ea = {ea_kj:.1f} kJ/mol")
    if ea_kj < 30:
        A("  → Very low barrier. Reaction may be diffusion-limited or involves")
        A("    very weak bond rearrangements (e.g., isomerases).")
    elif ea_kj < 60:
        A("  → Typical for mesophilic enzymes. Temperature sensitivity is moderate.")
        A("    A 10 °C rise roughly doubles the rate (Q₁₀ ≈ 2).")
    elif ea_kj < 100:
        A("  → High barrier. Reaction involves significant bond breaking/forming.")
        A("    Strong temperature dependence; consider thermophilic variants.")
    else:
        A("  → Very high barrier. Unusual — verify data quality and temperature range.")
    A()
    A(f"  Arrhenius equation:  ln(k) = ln(A) − Ea / (R·T)")
    A(f"  where R = 8.314 J·mol⁻¹·K⁻¹ and T is in Kelvin.")

    # ── Eyring ───────────────────────────────────────────────────────────────
    A(_section("2.  EYRING / TRANSITION-STATE THEORY"))
    dH_kj = eyr_res["dH"] / 1000
    dS    = eyr_res["dS"]
    dG_kj = eyr_res["dG_25"] / 1000

    A(f"  Enthalpy of activation  ΔH‡   :  {dH_kj:.3f} kJ/mol")
    A(f"  Entropy of activation   ΔS‡   :  {dS:.3f} J/mol·K")
    A(f"  Gibbs free energy       ΔG‡   :  {dG_kj:.3f} kJ/mol  (at 25 °C / 298.15 K)")
    A()
    A("  INTERPRETATION OF ΔS‡")
    if dS < -20:
        A(f"  ΔS‡ = {dS:.1f} J/mol·K  →  Strongly negative (ordered transition state)")
        A("  The TS requires a high degree of geometric precision. Reactants, solvent,")
        A("  or co-factors must adopt a specific orientation before catalysis can occur.")
        A("  Hallmark of ASSOCIATIVE mechanisms (e.g., nucleophilic substitutions).")
    elif dS < 0:
        A(f"  ΔS‡ = {dS:.1f} J/mol·K  →  Mildly negative (slightly ordered TS)")
        A("  Moderate orientation requirement. Consistent with many bimolecular enzyme")
        A("  reactions.")
    else:
        A(f"  ΔS‡ = {dS:.1f} J/mol·K  →  Positive (disordered transition state)")
        A("  The TS is less ordered than the ground state. Solvent release, ring opening,")
        A("  or bond cleavage increases disorder. Hallmark of DISSOCIATIVE mechanisms.")
    A()
    A(f"  Eyring equation:  k = (kB·T/h) · exp(−ΔG‡/RT)")
    A(f"  where kB = 1.381×10⁻²³ J/K (Boltzmann), h = 6.626×10⁻³⁴ J·s (Planck)")
    A()
    A("  RELATIONSHIP BETWEEN Arrhenius and Eyring parameters:")
    A(f"  ΔH‡ ≈ Ea − RT  =  {ea_kj:.2f} − {8.314 * 298.15 / 1000:.2f}  =  {ea_kj - 8.314 * 298.15 / 1000:.2f} kJ/mol")

    # ── FOOTER ───────────────────────────────────────────────────────────────
    A()
    A(_THIN)
    A("  BioKineticPy v1.1  ·  github.com/PraiseElement/BioKineticPy")
    A("  Author: Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DIVIDER)

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
def generate_simulation_report(model_name, params, c_unit, v_unit):
    """Generate a formatted simulation report (plain-text)."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    lines = []
    A = lines.append

    A(_DIVIDER)
    A("  BioKineticPy v1.1  ·  SIMULATION REPORT")
    A(_DIVIDER)
    A(f"  Generated : {timestamp}")
    A(_THIN)

    A(_section("1.  SIMULATION OVERVIEW"))
    A(f"  Model             :  {model_name}")
    A(f"  Concentration unit:  {c_unit}")
    A(f"  Velocity unit     :  {v_unit}")
    eq = MODEL_EQUATIONS.get(model_name, "")
    if eq:
        A(f"  Rate law          :  {eq}")

    A(_section("2.  INPUT PARAMETERS"))
    nice = {
        "vmax":  "Vmax  (max velocity)",
        "km":    "Km    (Michaelis const.)",
        "ki":    "Ki    (inhibition const.)",
        "alpha": "Alpha (binding selectivity)",
        "khalf": "Khalf (half-saturation)",
        "n":     "n     (Hill coefficient)",
        "ksi":   "Ksi   (substrate inhibition)",
        "i":     "[I]   (inhibitor concentration)",
    }
    for k, v in params.items():
        if v != 0:
            A(f"  {nice.get(k, k):<30}:  {v:.4f}")

    A(_section("3.  KINETIC THEORY"))
    for line in MODEL_THEORY.get(model_name, "Standard kinetic model.").splitlines():
        A(f"  {line}")

    A()
    A(_THIN)
    A("  BioKineticPy v1.1  ·  github.com/PraiseElement/BioKineticPy")
    A("  Author: Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DIVIDER)

    return "\n".join(lines)