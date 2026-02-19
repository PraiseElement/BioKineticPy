"""
libbiokinetic/report.py
─────────────────────────────────────────────────────────────────────────────
Plain-text report generator for BioKineticPy v1.1.
All three report types (analysis, thermodynamics, simulation) are defined here.
Each section is written to be self-explanatory — the reader should not need a
textbook to understand the output.
"""
import datetime
import math

# ─────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
R_GAS  = 8.314          # J mol⁻¹ K⁻¹  (universal gas constant)
kB     = 1.380649e-23   # J K⁻¹         (Boltzmann constant)
h_P    = 6.62607015e-34 # J s            (Planck constant)
T_STD  = 298.15         # K              (25 °C standard temperature)

_W  = 72        # line width
_DV = "═" * _W  # heavy divider
_DH = "─" * _W  # light divider
_DT = "·" * _W  # dot divider

def _sec(title):
    return f"\n{_DV}\n  {title}\n{_DV}\n"

def _sub(title):
    return f"\n  ── {title} {'─' * max(0, _W - len(title) - 6)}\n"

def _wrap(text, indent=2, width=_W - 4):
    """Naive word-wrap for long strings."""
    words   = text.split()
    lines   = []
    current = " " * indent
    for w in words:
        if len(current) + len(w) + 1 <= width + indent:
            current += w + " "
        else:
            lines.append(current.rstrip())
            current = " " * indent + w + " "
    if current.strip():
        lines.append(current.rstrip())
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
#  MODEL THEORY  (mechanism + key diagnostic per model)
# ─────────────────────────────────────────────────────────────────────────────
MODEL_THEORY = {
    "Michaelis-Menten": {
        "equation": "v = Vmax · [S] / (Km + [S])",
        "assumptions": [
            "Single substrate, single active site",
            "Rapid-equilibrium or steady-state ES complex",
            "Negligible product inhibition (initial velocities only)",
            "[S] >> [E]t  (substrate greatly exceeds enzyme)",
        ],
        "mechanism": (
            "E + S  ⇌  ES  →  E + P\n"
            "  │              │\n"
            "  └─ kon/koff ─→ └─ kcat (rate-limiting step)\n\n"
            "  Km  =  (koff + kcat) / kon\n"
            "  When koff >> kcat:  Km ≈ Kd  (true dissociation constant)\n"
            "  When kcat >> koff:  Km >> Kd (Km overestimates substrate affinity)"
        ),
        "parameters": {
            "vmax": "Maximum catalytic velocity at saturating [S] ( = kcat × [E]t )",
            "km":   "Substrate concentration at half-maximal velocity (≈ substrate affinity)",
        },
        "diagnostic": "Hyperbolic v vs [S] curve. Linear Lineweaver-Burk (1/v vs 1/[S]) plot.",
        "biology": "Most single-substrate enzymes at sub-saturating [S] follow this model.",
    },
    "Competitive": {
        "equation": "v = Vmax · [S] / (Km·(1 + [I]/Ki) + [S])",
        "assumptions": [
            "Inhibitor and substrate compete for the SAME binding site",
            "Inhibitor binding is reversible",
            "Inhibitor does NOT bind the ES complex",
        ],
        "mechanism": (
            "E + S  ⇌  ES  →  E + P\n"
            "E + I  ⇌  EI   (dead-end complex — cannot bind S)\n\n"
            "  Km(app) = Km · α      where  α = 1 + [I]/Ki  ≥ 1\n"
            "  Vmax(app) = Vmax      (unchanged — excess S displaces I)\n\n"
            "  α is the competitive inhibition factor: at [I] = Ki, α = 2,\n"
            "  meaning Km doubles and you need twice as much substrate to\n"
            "  achieve half-maximal velocity."
        ),
        "parameters": {
            "vmax": "Unchanged — can always be recovered by raising [S]",
            "km":   "Apparent Km = Km · (1 + [I]/Ki)  →  increases with [I]",
            "ki":   "Inhibitor dissociation constant from the EI complex; lower Ki = tighter binding",
        },
        "diagnostic": "Vmax unchanged, Km increases. Lineweaver-Burk lines intersect on the y-axis.",
        "biology": "Statins (HMG-CoA reductase), methotrexate (DHFR), most competitive drug leads.",
    },
    "Non-Competitive": {
        "equation": "v = (Vmax / (1 + [I]/Ki)) · [S] / (Km + [S])",
        "assumptions": [
            "Inhibitor binds an ALLOSTERIC site (not the active site)",
            "Binds E and ES complex with equal affinity",
            "Does not affect substrate binding; blocks catalysis only",
        ],
        "mechanism": (
            "E + S  ⇌  ES  →  E + P\n"
            "│              │\n"
            "E + I ⇌ EI    ES + I ⇌ ESI  (both inactive)\n\n"
            "  Vmax(app) = Vmax / α     (decreases — fewer active complexes)\n"
            "  Km(app)   = Km           (unchanged — substrate affinity preserved)\n\n"
            "  The inhibitor effectively reduces the concentration of active\n"
            "  enzyme without competing with the substrate."
        ),
        "parameters": {
            "vmax": "Apparent Vmax = Vmax / (1 + [I]/Ki)  →  decreases with [I]",
            "km":   "Unchanged — substrate affinity unaffected by inhibitor",
            "ki":   "Inhibitor dissociation constant (equal for EI and ESI complexes)",
        },
        "diagnostic": "Vmax decreases, Km unchanged. LB lines cross on the x-axis (same slope, different intercepts).",
        "biology": "Heavy-metal poisoning; iodoacetamide on cysteine proteases; allosteric regulators.",
    },
    "Uncompetitive": {
        "equation": "v = (Vmax/(1+[I]/Ki)) · [S] / (Km/(1+[I]/Ki) + [S])",
        "assumptions": [
            "Inhibitor binds ONLY the ES complex, not free enzyme",
            "Substrate binding creates (or exposes) the inhibitor site",
            "Binding is reversible",
        ],
        "mechanism": (
            "E + S  ⇌  ES  →  E + P\n"
            "          │\n"
            "         ESI  (inactive;  I binds only after S)\n\n"
            "  α = 1 + [I]/Ki\n"
            "  Vmax(app) = Vmax / α     (decreases)\n"
            "  Km(app)   = Km / α       (ALSO decreases — counter-intuitive!)\n"
            "  Vmax/Km ratio = constant (unaffected by inhibitor)\n\n"
            "  PARADOX EXPLAINED: binding I to ES traps more S on the enzyme\n"
            "  (removes free ES from equilibrium), so apparent Km drops even\n"
            "  though catalysis is impaired."
        ),
        "parameters": {
            "vmax": "Decreases: Vmax / (1 + [I]/Ki)",
            "km":   "Also decreases: Km / (1 + [I]/Ki) — note: unusual direction",
            "ki":   "Inhibitor-ES dissociation constant; the only way to determine Ki here is from Dixon or LB slopes",
        },
        "diagnostic": "Both Vmax and Km decrease equally. LB lines PARALLEL (same slope). Dixon: lines converge below x-axis.",
        "biology": "Li⁺ on inositol monophosphatase; phenylalanine on alkaline phosphatase; drug targets in multi-substrate systems.",
    },
    "Mixed": {
        "equation": "v = Vmax · [S] / (α·Km·(1+[I]/(α·Ki)) + [S]·(1 + [I]/(α·Ki)))",
        "assumptions": [
            "Inhibitor binds both E and ES, but with DIFFERENT affinities",
            "α (alpha) captures the binding preference",
            "Includes Competitive (α→∞) and Uncompetitive (α→0) as limiting cases",
        ],
        "mechanism": (
            "E + S  ⇌   ES  →  E + P\n"
            "│               │\n"
            "EI ⇌ E + I     ESI (inhibitor binds both sides)\n\n"
            "  Ki(E)  = inhibitor dissociation constant from free enzyme\n"
            "  Ki(ES) = Ki / α  (inhibitor dissociation from ES complex)\n\n"
            "  α > 1 :  Ki(E) < Ki(ES) → prefers free enzyme  (competitive-like)\n"
            "  α < 1 :  Ki(E) > Ki(ES) → prefers ES complex   (uncompetitive-like)\n"
            "  α = 1 :  Ki(E) = Ki(ES) → pure non-competitive\n\n"
            "  Vmax(app ) = Vmax / (1 + [I]/(α·Ki))\n"
            "  Km(app)   = Km · α · (1 + [I]/Ki) / (1 + [I]/(α·Ki))"
        ),
        "parameters": {
            "vmax":  "Maximum velocity (Vmax from non-inhibited enzyme)",
            "km":    "Michaelis constant of uninhibited enzyme",
            "ki":    "Inhibitor dissociation constant from free enzyme",
            "alpha": "Binding selectivity (α).  >1 → competitive-like, <1 → uncompetitive-like",
        },
        "diagnostic": "Lines on LB plot converge in the second quadrant (α > 1) or third quadrant (α < 1).",
        "biology": "Many real drugs; ibuprofen on COX-2; calmodulin antagonists.",
    },
    "Hill (Allosteric)": {
        "equation": "v = Vmax · [S]ⁿ / (Khalf^n + [S]ⁿ)",
        "assumptions": [
            "Enzyme has two or more substrate-binding sites",
            "Sites communicate (positive or negative cooperativity)",
            "Does not distinguish between concerted (MWC) and sequential (KNF) mechanisms",
        ],
        "mechanism": (
            "  Empirical extension of Michaelis-Menten for cooperative enzymes:\n\n"
            "  log(v / (Vmax − v))  =  n · log[S]  −  n · log(Khalf)\n"
            "  (Hill plot — slope = n over the dynamic range)\n\n"
            "  n > 1  →  positive cooperativity: sigmoidal curve, switch-like response\n"
            "  n = 1  →  hyperbolic: reduces exactly to Michaelis-Menten\n"
            "  n < 1  →  negative cooperativity: sub-hyperbolic, flat curve\n\n"
            "  NOTE: n rarely equals the number of binding sites — it is a\n"
            "  phenomenological coefficient, not a stoichiometric one."
        ),
        "parameters": {
            "vmax":  "Maximum velocity at saturating substrate",
            "khalf": "[S] at half-maximal velocity ( ≡ Km only when n = 1 )",
            "n":     "Hill coefficient: extent of cooperativity",
        },
        "diagnostic": "Sigmoidal v vs [S] for n > 1. Hill plot linearises the relationship.",
        "biology": "Haemoglobin (n ≈ 2.8–3.0), ATCase, phosphofructokinase-1, Ca²⁺-calmodulin.",
    },
    "Substrate Inhibition": {
        "equation": "v = Vmax · [S] / (Km + [S] + [S]² / Ksi)",
        "assumptions": [
            "Second substrate molecule binds the active site when [S] is very high",
            "The resulting ESS complex is catalytically inactive",
            "Single-substrate reaction at lower [S] behaves as Michaelis-Menten",
        ],
        "mechanism": (
            "  E + S  ⇌  ES  →  E + P    (productive, low [S])\n"
            "  ES + S ⇌  ESS              (dead-end, high [S])\n\n"
            "  Velocity formula:\n"
            "    v  =  Vmax · [S]  /  (Km  +  [S]  +  [S]²/Ksi)\n\n"
            "  Optimal substrate concentration [S]opt:\n"
            "    [S]opt  =  √(Km · Ksi)\n\n"
            "  At [S] = [S]opt  the velocity is maximal. Beyond this point,\n"
            "  ESS accumulation causes velocity to fall — the enzyme is\n"
            "  paradoxically inhibited by its own substrate."
        ),
        "parameters": {
            "vmax": "Theoretical maximum velocity (at [S]opt, actual peak < Vmax)",
            "km":   "Michaelis constant for the productive ES step",
            "ksi":  "Substrate inhibition constant (lower Ksi = inhibition sets in at lower [S])",
        },
        "diagnostic": "Bell-shaped v vs [S] curve with a distinct velocity peak.",
        "biology": "Acetylcholinesterase, alcohol dehydrogenase, tyrosinase, phenol-oxidising enzymes.",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
#  PARAMETER NICE NAMES
# ─────────────────────────────────────────────────────────────────────────────
PARAM_LABELS = {
    "vmax":  "Vmax  (maximum velocity)",
    "km":    "Km    (Michaelis constant)",
    "ki":    "Ki    (inhibition constant)",
    "alpha": "Alpha (binding selectivity α)",
    "khalf": "Khalf (half-saturation conc.)",
    "n":     "n     (Hill coefficient)",
    "ksi":   "Ksi   (substrate inhibition const.)",
    "i":     "[I]   (inhibitor concentration)",
}


# ─────────────────────────────────────────────────────────────────────────────
#  ANALYSIS REPORT
# ─────────────────────────────────────────────────────────────────────────────
def generate_methods_report(best_result, dataset_size, c_unit, v_unit,
                             weighting_mode, stats_res,
                             derived_stats=None, ci_data=None,
                             all_results=None):
    """Generate a comprehensive, self-explanatory kinetic analysis report."""
    model_name = best_result["model"]
    aicc       = best_result.get("aic", best_result.get("aicc", float("nan")))
    r2         = best_result["r_squared"]
    params     = best_result["parameters"]
    timestamp  = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    theory     = MODEL_THEORY.get(model_name, {})

    lines = []
    A = lines.append

    # ── TITLE PAGE ─────────────────────────────────────────────────────────
    A(_DV)
    A("  BioKineticPy v1.1  ·  COMPREHENSIVE KINETIC ANALYSIS REPORT")
    A(_DV)
    A(f"  Generated  :  {timestamp}")
    A(f"  Software   :  BioKineticPy (open source, github.com/PraiseElement/BioKineticPy)")
    A(f"  Author     :  Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DT)
    A("  HOW TO READ THIS REPORT")
    A("  ─────────────────────────────────────────────────────────────────────")
    A("  1. Section 2 tells you WHICH model fits best and WHY it was chosen.")
    A("  2. Section 3 explains the biological mechanism behind that model.")
    A("  3. Section 4 gives the KINETIC PARAMETERS with confidence intervals.")
    A("  4. Section 5 evaluates HOW GOOD the fit is (R², AICc).")
    A("  5. Section 6 checks whether the fit is STATISTICALLY VALID.")
    A("  6. Section 7 provides PRACTICAL RECOMMENDATIONS based on the results.")
    A(_DT)

    # ── 1. EXPERIMENTAL OVERVIEW ────────────────────────────────────────────
    A(_sec("1.  EXPERIMENTAL OVERVIEW"))
    A(f"  Dataset size          :  {dataset_size} data points")
    A(f"  Concentration unit    :  {c_unit}")
    A(f"  Velocity unit         :  {v_unit}")
    A(f"  Regression algorithm  :  Levenberg-Marquardt (Damped Least Squares)")
    A(f"  Weighting strategy    :  {weighting_mode}")
    A(f"  Model selection       :  AICc (Corrected Akaike Information Criterion)")
    if derived_stats and derived_stats.get("et_val", 0) > 0:
        A(f"  Enzyme concentration  :  {derived_stats['et_val']} {derived_stats['et_unit']}")
    A("")
    A("  WEIGHTING STRATEGY EXPLAINED")
    A(_DH)
    wt = weighting_mode.lower()
    if "1/v²" in wt or "1/v2" in wt:
        A("  1/v² weighting: each point is weighted by the inverse square of its")
        A("  velocity. This is appropriate when the COEFFICIENT OF VARIATION (CV%)")
        A("  of repeated measurements is roughly constant across all velocities")
        A("  (multiplicative / log-normal error). Most common in enzyme kinetics.")
    elif "1/v" in wt:
        A("  1/v weighting: each point is weighted by 1/velocity. Intermediate between")
        A("  unweighted and 1/v² — use when variance scales with v (not v²).")
    else:
        A("  Unweighted (ordinary least squares): all points contribute equally.")
        A("  Appropriate when variance is roughly constant across all velocities.")
        A("  Less ideal if high-velocity points have larger absolute errors.")

    # ── 2. MODEL SELECTION RANKING ───────────────────────────────────────────
    A(_sec("2.  MODEL SELECTION RANKING"))
    A("  WHAT IS AICc?")
    A(_DH)
    A("  The Corrected Akaike Information Criterion penalises model complexity:")
    A("")
    A("     AICc  =  AIC  +  2k(k+1) / (n − k − 1)")
    A("     AIC   =  n·ln(RSS/n) + 2k")
    A("")
    A("  where:  n = number of data points        (more data = more power)")
    A("          k = number of free parameters    (more params = penalty)")
    A("          RSS = residual sum of squares     (smaller = better fit)")
    A("")
    A("  The correction term matters most when n is small (< 40 points).")
    A("  The model with the LOWEST AICc is the best-supported model.")
    A("")
    A("  ΔAIC INTERPRETATION GUIDE")
    A("  ─────────────────────────")
    A("  ΔAIC  0 – 2  :  Models are virtually indistinguishable. The simpler")
    A("                   model is usually preferred (Occam's razor).")
    A("  ΔAIC  2 – 7  :  Moderate evidence for the winning model.")
    A("  ΔAIC  7 – 10 :  Substantial evidence; lower-ranked model unlikely.")
    A("  ΔAIC  > 10   :  Decisive evidence — best model strongly preferred.")
    A("")

    if all_results:
        A(f"  {'Rank':<5} {'Model':<24} {'AICc':>10} {'R²':>8}  {'ΔAIC':>8}  Notes")
        A("  " + _DH)
        best_aic = min(r["aic"] for r in all_results if r["aic"] is not None)
        ranked   = sorted(all_results, key=lambda r: r["aic"] if r["aic"] is not None else 1e9)
        for rank, res in enumerate(ranked, 1):
            delta  = res["aic"] - best_aic if res["aic"] is not None else float("nan")
            note   = " ◀ SELECTED" if res["model"] == model_name else (
                     " (indistinct)"  if delta < 2 else
                     " (some support)" if delta < 7 else
                     " (little support)" if delta < 10 else
                     " (no support)")
            A(f"  {rank:<5} {res['model']:<24} {res['aic']:>10.2f} {res['r_squared']:>8.4f}  {delta:>8.2f}  {note}")
    else:
        A(f"  Selected model:  {model_name}  (AICc = {aicc:.2f},  R² = {r2:.4f})")

    # ── 3. MODEL MECHANISM ───────────────────────────────────────────────────
    A(_sec("3.  KINETIC MODEL — MECHANISM & THEORY"))
    A(f"  Model   :  {model_name.upper()}")
    eq = theory.get("equation", "")
    if eq:
        A(f"  Formula :  {eq}")
    A("")

    assump = theory.get("assumptions", [])
    if assump:
        A("  MODEL ASSUMPTIONS")
        A("  " + _DH)
        for a in assump:
            A(f"  • {a}")
        A("")

    mech = theory.get("mechanism", "")
    if mech:
        A("  MECHANISM (reaction scheme + mathematics)")
        A("  " + _DH)
        for line in mech.splitlines():
            A(f"  {line}")
        A("")

    diag = theory.get("diagnostic", "")
    if diag:
        A("  GRAPHICAL DIAGNOSTIC SIGNATURE")
        A("  " + _DH)
        A(f"  {diag}")

    bio = theory.get("biology", "")
    if bio:
        A("")
        A("  KNOWN BIOLOGICAL EXAMPLES")
        A("  " + _DH)
        A(f"  {bio}")

    # ── 4. KINETIC PARAMETERS ────────────────────────────────────────────────
    A(_sec("4.  KINETIC PARAMETERS"))
    A("  Values shown with 95% Bootstrap Confidence Intervals (if computed).")
    A("  The ± % column is the relative half-width of the CI — a measure of")
    A("  parameter precision. Values below 20% are generally acceptable.")
    A("")
    A(f"  {'Parameter':<28} {'Value':>12}   {'95% CI':^28}  ±%")
    A("  " + _DH)

    for k, v in params.items():
        label  = PARAM_LABELS.get(k, k)
        ci_str = ""
        pct    = ""
        if ci_data and k in ci_data:
            lo, hi = ci_data[k]
            ci_str = f"[{lo:>10.4f} – {hi:<10.4f}]"
            if v:
                pct = f"{abs(hi - lo) / 2 / abs(v) * 100:>5.1f}%"
        A(f"  {label:<28} {v:>12.4f}   {ci_str:<28}  {pct}")

    A("")
    A("  PARAMETER INTERPRETATIONS")
    A("  " + _DH)
    theory_params = theory.get("parameters", {})
    for k, v in params.items():
        label  = PARAM_LABELS.get(k, k).split()[0]
        interp = theory_params.get(k, "")
        if interp:
            A(f"  {label}")
            A(_wrap(interp, indent=6))
            A("")

    if derived_stats and derived_stats.get("kcat", 0) > 0:
        kcat = derived_stats["kcat"]
        eff  = derived_stats["efficiency"]
        A("")
        A("  DERIVED CATALYTIC CONSTANTS")
        A("  " + _DH)
        A("  These require [E]t (total enzyme concentration) to be specified.")
        A("")
        A(f"  kcat (turnover number)     :  {kcat:.4f} s⁻¹")
        A("")
        A("    Formula:  kcat = Vmax / [E]t")
        A("    Meaning:  Number of substrate molecules converted to product per")
        A("              active site per second at saturating [S]. This is the")
        A("              intrinsic speed of the enzyme regardless of [E].")
        A("")
        A("    Benchmark ranges:")
        A("    < 1 s⁻¹      — very slow (e.g., some isomerases, restriction enzymes)")
        A("    1–1000 s⁻¹   — typical range for most enzymes")
        A("    > 1000 s⁻¹   — fast (catalase ~4×10⁷, carbonic anhydrase ~10⁶)")
        A("")
        A(f"  kcat/Km (catalytic efficiency) :  {eff:.3e} M⁻¹ s⁻¹")
        A("")
        A("    Formula:  kcat/Km = kcat / Km  (units: M⁻¹ s⁻¹)")
        A("    Meaning:  The apparent second-order rate constant for the E + S")
        A("              reaction; captures both binding speed and catalytic rate.")
        A("              This is the most useful metric for comparing enzymes.")
        A("")
        A("    Benchmark ranges:")
        A("    < 10³ M⁻¹s⁻¹   — poor catalyst or non-optimal conditions")
        A("    10⁵–10⁷ M⁻¹s⁻¹ — typical enzyme (e.g., triosephosphate isomerase ~2×10⁸)")
        A("    > 10⁸ M⁻¹s⁻¹   — 'catalytically perfect' — diffusion-limited")

    # ── 5. GOODNESS OF FIT ──────────────────────────────────────────────────
    A(_sec("5.  GOODNESS OF FIT"))
    r2_pct   = r2 * 100
    r2_grade = (
        "EXCELLENT  (> 99% variance explained)" if r2 > 0.99 else
        "GOOD       (> 97% variance explained)" if r2 > 0.97 else
        "ACCEPTABLE (> 95% variance explained)" if r2 > 0.95 else
        "POOR       (< 95% — model may be wrong or data quality limited)"
    )
    A(f"  R²  =  {r2:.6f}   ({r2_pct:.3f}% of variance explained)   →  {r2_grade}")
    A("")
    A("  WHAT R² MEANS")
    A("  " + _DH)
    A("  R² = 1 − RSS / TSS")
    A("  where RSS = Σ(v_observed − v_predicted)²  (residual sum of squares)")
    A("        TSS = Σ(v_observed − mean(v))²      (total sum of squares)")
    A("")
    A("  R² = 1.0 → perfect fit (all points lie on the model curve)")
    A("  R² = 0.0 → model is no better than the mean of the data")
    A("")
    A("  IMPORTANT: R² alone is NOT sufficient to validate a model.")
    A("  A wrong model can still give a high R² if it has enough parameters.")
    A("  Always also check the residual diagnostics in Section 6.")
    A("")
    A(f"  AICc  =  {aicc:.3f}")
    A("  AICc is used for model COMPARISON (see Section 2), not absolute judgment.")

    # ── 6. RESIDUAL DIAGNOSTICS ─────────────────────────────────────────────
    A(_sec("6.  RESIDUAL DIAGNOSTICS"))
    A("  Residuals = v_observed − v_predicted for each data point.")
    A("  A valid kinetic fit should produce residuals that are:")
    A("  (a) Normally distributed (Gaussian measurement noise)")
    A("  (b) Randomly scattered around zero (no systematic pattern)")
    A("")

    sw_p   = stats_res.get("shapiro_p", 1.0)
    sw_ok  = sw_p > 0.05
    rz     = abs(stats_res.get("runs_z", 0))
    rz_ok  = rz < 1.96

    # ── Shapiro-Wilk ────────
    A("  TEST 1  —  SHAPIRO-WILK  (tests whether residuals are normally distributed)")
    A("  " + _DH)
    A("  Null hypothesis (H₀): residuals come from a Normal distribution.")
    A("  p > 0.05 → FAIL TO REJECT H₀ (residuals look Gaussian)  ← desired")
    A("  p ≤ 0.05 → REJECT H₀ (residuals deviate from normality)  ← investigate")
    A("")
    sw_sym = "✓  PASS" if sw_ok else "✗  WARNING"
    A(f"  Result  :  {sw_sym}   (p = {sw_p:.4f})")
    A("")
    if sw_ok:
        A("  Interpretation:")
        A("  The residuals are consistent with Gaussian (random) measurement noise.")
        A("  This confirms that the model captures the systematic trend in the data")
        A("  and leaves only random scatter behind.")
    else:
        A("  Interpretation:")
        A("  The residuals deviate significantly from normality. Likely causes:")
        A("")
        A("  1. OUTLIERS — one or more data points are far from the trend.")
        A("     Action: Plot residuals vs [S]. Identify and check outlier wells.")
        A("")
        A("  2. WRONG MODEL — the model shape is mechanistically incorrect.")
        A("     Action: Compare AICc values for alternative models (Section 2).")
        A("")
        A("  3. SYSTEMATIC EXPERIMENTAL ERROR — pipetting inaccuracy, evaporation,")
        A("     substrate depletion, or substrate depletion at high [S].")
        A("     Action: Repeat assays; check [S] vs product inhibition window.")
        A("")
        A("  4. HETEROSCEDASTIC ERROR — variance is not constant across [S].")
        A("     Action: Re-run with 1/v² weighting and compare AICc.")

    A("")
    # ── Runs Test ────────
    A("  TEST 2  —  WALD-WOLFOWITZ RUNS TEST  (tests for systematic patterns)")
    A("  " + _DH)
    A("  A 'run' is a consecutive sequence of residuals with the same sign (+ or −).")
    A("  Too few runs → U-shaped or systematic residuals (model underfits).")
    A("  Too many runs → residuals oscillate (model overparameterised or data noisy).")
    A("")
    A("  H₀: The residual signs are randomly ordered.")
    A("  |Z| < 1.96 → FAIL TO REJECT H₀ (random order)  ← desired")
    A("  |Z| ≥ 1.96 → REJECT H₀ at α = 0.05 (non-random)  ← investigate")
    A("")
    rz_sym = "✓  PASS" if rz_ok else "✗  WARNING"
    A(f"  Result  :  {rz_sym}   (|Z| = {rz:.3f})")
    A("")
    if rz_ok:
        A("  Interpretation:")
        A("  No systematic pattern detected. The model captures the underlying")
        A("  relationship well, and residuals appear randomly scattered.")
    else:
        A("  Interpretation:")
        A("  Systematic structure found in residuals. Likely causes:")
        A("")
        A("  1. MODEL UNDERFITTING — the model curve has a different SHAPE than")
        A("     the data (e.g., data is sigmoidal but Michaelis-Menten was used).")
        A("     Action: Try Hill (Allosteric) or Substrate Inhibition models.")
        A("")
        A("  2. WRONG MECHANISM — e.g., competitive inhibition data fitted with")
        A("     Michaelis-Menten shows U-shaped residuals at each inhibitor [I].")
        A("     Action: Check the AICc ranking for inhibition models.")
        A("")
        A("  3. MISSING VARIABLE — a third variable (pH, temperature, cofactor)")
        A("     is affecting velocity but is not captured in the model.")

    # ── 7. RECOMMENDATIONS ─────────────────────────────────────────────────
    A(_sec("7.  PRACTICAL RECOMMENDATIONS"))
    issues  = []
    if not sw_ok:
        issues.append("• Normality test failed — check for outliers and consider re-weighting.")
    if not rz_ok:
        issues.append("• Runs Test failed — the model shape may not match your data.")
    if r2 < 0.95:
        issues.append("• R² < 0.95 — consider alternative models or investigate data quality.")
    if ci_data:
        wide = [k for k, v_val in params.items()
                if ci_data.get(k) and v_val and
                abs(ci_data[k][1] - ci_data[k][0]) / 2 / abs(v_val) > 0.30]
        if wide:
            pnames = ", ".join(PARAM_LABELS.get(k, k).split()[0] for k in wide)
            issues.append(f"• Wide CIs on {pnames} — extended [S] range or more replicates recommended.")

    if not issues:
        A("  All statistical checks passed. The model is well-supported by the data.")
        A("")
        A("  Suggested next steps:")
        A("  • Compute 95% Bootstrap CIs if not already done (adds confidence intervals).")
        A("  • Repeat at a wider [S] range to confirm Vmax plateau.")
        if derived_stats and not derived_stats.get("kcat", 0):
            A("  • Enter [E]t in the sidebar to compute kcat and kcat/Km.")
    else:
        A("  The following issues require attention:\n")
        for iss in issues:
            A(f"  {iss}")
        A("")
        A("  General best-practice checklist:")
        A("  ── Use at least 6–8 [S] points spanning 0.2×Km to 10×Km.")
        A("  ── Include at least 2–3 inhibitor concentrations near Ki (if applicable).")
        A("  ── Take triplicate readings at each [S] for robust variance estimation.")
        A("  ── Verify initial velocity conditions (linear progress curve < 10% conversion).")

    # ── FOOTER ───────────────────────────────────────────────────────────────
    A("")
    A(_DT)
    A("  BioKineticPy v1.1  ·  github.com/PraiseElement/BioKineticPy")
    A("  Methods: Levenberg-Marquardt NLS · AICc model ranking · Bootstrap 95% CI")
    A("  Reference: Motulsky & Christopoulos, 'Fitting Models to Biological Data'")
    A("             Cornish-Bowden, 'Fundamentals of Enzyme Kinetics' (4th ed.)")
    A(_DV)

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
#  THERMODYNAMIC REPORT
# ─────────────────────────────────────────────────────────────────────────────
def generate_thermo_report(arr_res, eyr_res):
    """Generate a comprehensive thermodynamic profile report."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    lines = []
    A = lines.append

    A(_DV)
    A("  BioKineticPy v1.1  ·  THERMODYNAMIC PROFILE REPORT")
    A(_DV)
    A(f"  Generated  :  {timestamp}")
    A(f"  Software   :  BioKineticPy (github.com/PraiseElement/BioKineticPy)")
    A(_DT)
    A("  HOW TO READ THIS REPORT")
    A("  ─────────────────────────────────────────────────────────────────────")
    A("  Section 1 — Arrhenius analysis: activation energy and collision frequency.")
    A("  Section 2 — Eyring analysis: enthalpy, entropy, and Gibbs energy of activation.")
    A("  Section 3 — Comparative summary and practical interpretation.")
    A(_DT)

    ea_kj    = arr_res["Ea"] / 1000
    A_val    = arr_res["A"]
    r2_arr   = arr_res["r_squared"]
    dH_kj    = eyr_res["dH"] / 1000
    dS       = eyr_res["dS"]
    dG_kj    = eyr_res["dG_25"] / 1000

    # ── 1. Arrhenius ────────────────────────────────────────────────────────
    A(_sec("1.  ARRHENIUS ANALYSIS"))
    A("  THEORETICAL BACKGROUND")
    A("  " + _DH)
    A("  The Arrhenius equation describes how reaction rate depends on temperature:")
    A("")
    A("     k(T)  =  A · exp(−Ea / (R·T))")
    A("")
    A("     Linearised form used for fitting:")
    A("     ln(k)  =  ln(A)  −  Ea/R · (1/T)")
    A("")
    A("     where:  k   = observed rate constant at temperature T (Kelvin)")
    A("             A   = pre-exponential factor  (s⁻¹ for unimolecular reactions)")
    A("             Ea  = activation energy  (J mol⁻¹)")
    A("             R   = 8.314 J mol⁻¹ K⁻¹  (gas constant)")
    A("")
    A("  FITTING RESULTS")
    A("  " + _DH)
    A(f"  Activation Energy  Ea  :  {ea_kj:.3f} kJ/mol")
    A(f"  Pre-Exponential    A   :  {A_val:.4e} s⁻¹")
    A(f"  Linear fit quality R²  :  {r2_arr:.6f}")
    A("")
    A("  ACTIVATION ENERGY EXPLAINED")
    A("  " + _DH)
    A("  Ea is the minimum kinetic energy that reacting molecules must possess")
    A("  for the reaction to proceed. It governs temperature sensitivity:")
    A("")
    RT_25 = R_GAS * T_STD / 1000   # kJ/mol thermal energy at 25°C
    A(f"  Thermal energy at 25°C (RT):  {RT_25:.2f} kJ/mol")
    A(f"  Ea / RT  =  {ea_kj:.1f} / {RT_25:.2f}  =  {ea_kj/RT_25:.1f}  ×  thermal energy")
    A("")
    # Q10 estimate
    try:
        q10 = math.exp(ea_kj * 1000 * 10 / (R_GAS * T_STD * (T_STD + 10)))
        A(f"  Q₁₀ estimate (rate doubling per 10°C?):  {q10:.2f}")
        A(f"  → A 10°C rise multiplies the rate by approximately {q10:.1f}×.")
    except Exception:
        pass
    A("")
    if ea_kj < 20:
        A("  ▸ Very low Ea  (< 20 kJ/mol)")
        A("    Possibly diffusion-limited; the chemical step is not rate-limiting.")
        A("    Examples: proton transfers, diffusion-controlled radical reactions.")
    elif ea_kj < 50:
        A("  ▸ Low-moderate Ea  (20–50 kJ/mol)")
        A("    Common for enzyme-catalysed reactions involving proton transfer or")
        A("    conformational change steps. Weak temperature dependence.")
    elif ea_kj < 80:
        A("  ▸ Moderate Ea  (50–80 kJ/mol)")
        A("    Typical for covalent bond rearrangements in enzymatic mechanisms.")
        A("    Strong temperature dependence — Q10 ≈ 2–4.")
    elif ea_kj < 130:
        A("  ▸ High Ea  (80–130 kJ/mol)")
        A("    Significant bond breaking required. Reaction is strongly temperature-")
        A("    sensitive. Common in uncatalysed organic reactions.")
    else:
        A("  ▸ Very high Ea  (> 130 kJ/mol)")
        A("    Unusual for enzyme-catalysed reactions. Verify temperature range and")
        A("    data quality. Protein denaturation may be confounding the analysis.")

    A("")
    A("  PRE-EXPONENTIAL FACTOR A EXPLAINED")
    A("  " + _DH)
    A(f"  A  =  {A_val:.3e} s⁻¹")
    A("")
    A("  A is the theoretical maximum rate at infinite temperature (all collisions")
    A("  are productive). It incorporates collision frequency AND the probability")
    A("  that colliding molecules have the correct orientation.")
    A("")
    A("  Typical range for enzyme kcat:  10⁶ – 10¹³ s⁻¹")
    A("  (Eyring theory predicts maximum kT/h ≈ 6×10¹² s⁻¹ at 25°C)")

    # ── 2. Eyring ───────────────────────────────────────────────────────────
    A(_sec("2.  EYRING / TRANSITION-STATE THEORY  (TST)"))
    A("  THEORETICAL BACKGROUND")
    A("  " + _DH)
    A("  TST treats the reaction as passing through a high-energy TRANSITION STATE (‡):")
    A("")
    A("     Reactants  →  [TS‡]  →  Products")
    A("")
    A("  The Eyring equation:")
    A("")
    A("     k(T)  =  (kB·T / h) · exp(−ΔG‡ / RT)")
    A("            =  (kB·T / h) · exp(ΔS‡/R) · exp(−ΔH‡/(RT))")
    A("")
    A("     Linearised form used for fitting:")
    A("     ln(k/T)  =  ln(kB/h) + ΔS‡/R  −  ΔH‡/R · (1/T)")
    A("")
    A("     kB  =  1.381 × 10⁻²³ J K⁻¹  (Boltzmann constant)")
    A("     h   =  6.626 × 10⁻³⁴ J s    (Planck constant)")
    A("     kT/h at 25°C  =  6.25 × 10¹² s⁻¹")
    A("")
    A("  FITTING RESULTS")
    A("  " + _DH)
    A(f"  ΔH‡  (Enthalpy of activation)  :  {dH_kj:.3f} kJ/mol")
    A(f"  ΔS‡  (Entropy of activation)   :  {dS:.3f} J/mol·K")
    A(f"  ΔG‡  (Gibbs energy barrier)    :  {dG_kj:.3f} kJ/mol  (at 25°C, 298.15 K)")
    A("")
    A("  THERMODYNAMIC CYCLE RELATIONSHIP")
    A("  " + _DH)
    A("  ΔG‡  =  ΔH‡  −  T·ΔS‡")
    A(f"       =  {dH_kj:.2f}  −  298.15 × {dS/1000:.4f}")
    A(f"       =  {dH_kj:.2f}  −  {298.15 * dS/1000:.2f}")
    A(f"       =  {dG_kj:.3f} kJ/mol  ✓")
    A("")
    A("  ENTHALPY  ΔH‡  =  {:.3f} kJ/mol".format(dH_kj))
    A("  " + _DH)
    A("  ΔH‡ is the enthalpic component of the barrier — energy required to")
    A("  break bonds and reach the transition state geometry.")
    A("")
    A("  Relationship to Ea:")
    ea_from_dH = dH_kj + RT_25
    A(f"  Ea  ≈  ΔH‡ + RT  =  {dH_kj:.2f} + {RT_25:.2f}  =  {ea_from_dH:.2f} kJ/mol")
    A(f"  Reported Ea  =  {ea_kj:.2f} kJ/mol    Difference = {abs(ea_kj - ea_from_dH):.2f} kJ/mol")
    A("  (Small differences are expected depending on the fitting method.)")
    A("")
    A("  ENTROPY  ΔS‡  =  {:.3f} J/mol·K".format(dS))
    A("  " + _DH)
    A("  ΔS‡ is the entropic cost of forming the transition state geometry.")
    A("")
    if dS < -40:
        A("  ΔS‡ << 0  (strongly negative — {:.1f} J/mol·K)".format(dS))
        A("  ──────────────────────────────────────────────────────────────")
        A("  Interpretation: HIGHLY ORDERED transition state.")
        A("  Both substrate and enzyme must adopt a very precise geometry.")
        A("  Substrate, protein residues, and solvent (water) lose considerable")
        A("  freedom of motion as they lock into the TS conformation.")
        A("  ASSOCIATIVE MECHANISM: a new bond forms BEFORE the old one breaks.")
        A("  Examples: nucleophilic substitution (SN2-like), phosphoryl transfers.")
    elif dS < 0:
        A("  ΔS‡ < 0  (mildly negative — {:.1f} J/mol·K)".format(dS))
        A("  ──────────────────────────────────────────────────────────────")
        A("  Interpretation: Moderately ordered transition state.")
        A("  Some geometric constraint required, but not extreme.")
        A("  Consistent with bimolecular reactions or enzyme active-site closure.")
    else:
        A("  ΔS‡ > 0  (positive — {:.1f} J/mol·K)".format(dS))
        A("  ──────────────────────────────────────────────────────────────")
        A("  Interpretation: DISORDERED transition state (higher entropy than ground).")
        A("  Net release of structural constraint as TS forms — solvent dipoles")
        A("  are freed, ring opens, or a bond cleaves before a new one forms.")
        A("  DISSOCIATIVE MECHANISM: old bond breaks before new one forms.")
        A("  Examples: ring-opening reactions, carbocation intermediates.")
    A("")
    A("  GIBBS FREE ENERGY  ΔG‡  =  {:.3f} kJ/mol".format(dG_kj))
    A("  " + _DH)
    A("  ΔG‡ is the total energy barrier that determines the observable rate.")
    A("  The Eyring equation converts ΔG‡ directly to a rate constant:")
    try:
        k_25 = (kB * T_STD / h_P) * math.exp(-dG_kj * 1000 / (R_GAS * T_STD))
        A(f"  Predicted k at 25°C  =  {k_25:.3e} s⁻¹")
    except Exception:
        pass
    A("")
    A("  ΔG‡ scale for enzyme-catalysed reactions:")
    A("  40–60 kJ/mol   — fast enzymes (kcat > 100 s⁻¹)")
    A("  60–80 kJ/mol   — moderate enzymes (kcat 1–100 s⁻¹)")
    A("  80–100 kJ/mol  — slow enzymes (kcat < 1 s⁻¹)")

    # ── 3. Comparison ───────────────────────────────────────────────────────
    A(_sec("3.  COMPARATIVE SUMMARY & BIOLOGICAL CONTEXT"))
    A(f"  Ea   (Arrhenius)  =  {ea_kj:.2f} kJ/mol")
    A(f"  ΔH‡  (Eyring)     =  {dH_kj:.2f} kJ/mol")
    A(f"  ΔS‡  (Eyring)     =  {dS:.2f} J/mol·K")
    A(f"  ΔG‡  (Eyring)     =  {dG_kj:.2f} kJ/mol  (at 25°C)")
    A("")
    A("  BIOLOGICAL CONTEXT — typical mesophilic enzyme values:")
    A("  Ea      :  40–80 kJ/mol")
    A("  ΔH‡     :  similar to Ea (within RT ≈ 2.5 kJ/mol)")
    A("  ΔS‡     :  −80 to +20 J/mol·K (highly variable by mechanism)")
    A("  ΔG‡     :  50–90 kJ/mol")

    # ── FOOTER ─────────────────────────────────────────────────────────────
    A("")
    A(_DT)
    A("  BioKineticPy v1.1  ·  Levenberg-Marquardt NLS · AICc model selection")
    A("  Reference: Fersht, 'Structure and Mechanism in Protein Science'")
    A("             Atkins & de Paula, 'Physical Chemistry for the Life Sciences'")
    A("  Author: Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DV)

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
#  SIMULATION REPORT
# ─────────────────────────────────────────────────────────────────────────────
def generate_simulation_report(model_name, params, c_unit, v_unit):
    """Generate a comprehensive simulation/parameter-exploration report."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    theory    = MODEL_THEORY.get(model_name, {})
    lines = []
    A = lines.append

    A(_DV)
    A("  BioKineticPy v1.1  ·  SIMULATION REPORT")
    A(_DV)
    A(f"  Generated  :  {timestamp}")
    A(_DT)

    A(_sec("1.  SIMULATION OVERVIEW"))
    A(f"  Model             :  {model_name}")
    A(f"  Concentration unit:  {c_unit}")
    A(f"  Velocity unit     :  {v_unit}")
    eq = theory.get("equation", "")
    if eq:
        A(f"  Rate law          :  {eq}")

    A(_sec("2.  INPUT PARAMETERS"))
    for k, v in params.items():
        if v != 0:
            label  = PARAM_LABELS.get(k, k)
            interp = theory.get("parameters", {}).get(k, "")
            A(f"  {label:<30}:  {v:.4f}")
            if interp:
                A(_wrap(interp, indent=6))

    A(_sec("3.  KINETIC MODEL — THEORY & MECHANISM"))
    assump = theory.get("assumptions", [])
    if assump:
        A("  MODEL ASSUMPTIONS")
        for a in assump:
            A(f"  • {a}")
        A("")
    mech = theory.get("mechanism", "")
    if mech:
        A("  MECHANISM (reaction scheme + mathematics)")
        A("  " + _DH)
        for line in mech.splitlines():
            A(f"  {line}")
        A("")
    diag = theory.get("diagnostic", "")
    if diag:
        A("  GRAPHICAL SIGNATURE")
        A(f"  {diag}")
    bio = theory.get("biology", "")
    if bio:
        A("")
        A("  KNOWN BIOLOGICAL EXAMPLES")
        A(f"  {bio}")

    A("")
    A(_DT)
    A("  BioKineticPy v1.1  ·  github.com/PraiseElement/BioKineticPy")
    A("  Author: Chibuike Praise Okechukwu  ·  praizekene1@gmail.com")
    A(_DV)

    return "\n".join(lines)
