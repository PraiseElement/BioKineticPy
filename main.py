import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats 
import json
import random

from libbiokinetic.models import Dataset, KineticPoint
from libbiokinetic.solver import KineticSolver
from libbiokinetic.equations import MODEL_REGISTRY
from libbiokinetic.units import convert_velocity_to_standard, convert_param_to_user, CONC_TO_MOLAR
from libbiokinetic.report import generate_methods_report, generate_simulation_report, generate_thermo_report
from libbiokinetic.thermo import solve_arrhenius, solve_eyring

# ───────────────────────────────────────────────────────
#  AUTHENTIC EXAMPLE DATA GENERATORS
# ───────────────────────────────────────────────────────
def _noisy(values, cv=0.08):
    """Add realistic Gaussian noise (default 8 % coefficient of variation)."""
    arr = np.array(values, dtype=float)
    noise = np.random.normal(1.0, cv, size=arr.shape)
    return np.round(np.clip(arr * noise, 0.01, None), 2)

def _mm_v(s, vmax, km):
    return vmax * s / (km + s)

def _comp_v(s, i, vmax, km, ki):
    return vmax * s / (km * (1 + i / ki) + s)

def _ncomp_v(s, i, vmax, km, ki):
    return (vmax / (1 + i / ki)) * s / (km + s)

def _hill_v(s, vmax, khalf, n):
    return vmax * s**n / (khalf**n + s**n)

def _ksi_v(s, vmax, km, ksi):
    return vmax * s / (km + s + s**2 / ksi)

def make_example_matrix(mode="competitive", c_unit="μM"):
    """
    Generate an authentic, randomly-varied matrix-format kinetic dataset.
    mode: 'competitive' | 'noncompetitive' | 'michaelis'
    Returns a DataFrame and the inhibitor concentrations [I_A, I_B].
    """
    rng = np.random  # fresh seed each call
    # ── Pick biologically realistic parameters ──────────
    if mode == "competitive":
        vmax  = rng.uniform(80, 160)     # μM/min  – typical enzyme
        km    = rng.uniform(5, 20)        # μM
        ki    = rng.uniform(3, 15)        # μM
        i_a   = round(rng.uniform(0.5*ki, 1.2*ki), 1)
        i_b   = round(rng.uniform(1.5*ki, 3.0*ki), 1)
        s_arr = np.array([round(x, 1) for x in np.geomspace(km*0.1, km*8, 7)])
        v0  = _noisy(_comp_v(s_arr, 0,   vmax, km, ki))
        vA  = _noisy(_comp_v(s_arr, i_a, vmax, km, ki))
        vB  = _noisy(_comp_v(s_arr, i_b, vmax, km, ki))
        label = "Competitive"
    elif mode == "noncompetitive":
        vmax  = rng.uniform(80, 140)
        km    = rng.uniform(8, 25)
        ki    = rng.uniform(5, 20)
        i_a   = round(rng.uniform(0.5*ki, 1.2*ki), 1)
        i_b   = round(rng.uniform(1.5*ki, 2.5*ki), 1)
        s_arr = np.array([round(x, 1) for x in np.geomspace(km*0.1, km*8, 7)])
        v0  = _noisy(_ncomp_v(s_arr, 0,   vmax, km, ki))
        vA  = _noisy(_ncomp_v(s_arr, i_a, vmax, km, ki))
        vB  = _noisy(_ncomp_v(s_arr, i_b, vmax, km, ki))
        label = "Non-Competitive"
    else:  # michaelis
        vmax = rng.uniform(60, 150)
        km   = rng.uniform(5, 30)
        i_a, i_b = 0.0, 0.0
        s_arr = np.array([round(x, 1) for x in np.geomspace(km*0.1, km*10, 8)])
        v0  = _noisy(_mm_v(s_arr, vmax, km))
        vA  = np.array([None]*len(s_arr))
        vB  = np.array([None]*len(s_arr))
        label = "Michaelis-Menten"

    df = pd.DataFrame({
        "Include?": [True]*len(s_arr),
        f"Substrate [{c_unit}]": s_arr,
        "v (I=0)": v0,
        "v (I=A)": vA,
        "v (I=B)": vB,
        "v (I=C)": [None]*len(s_arr),
    })
    return df, i_a, i_b, label

def make_example_standard(mode="hill", c_unit="μM"):
    """
    Generate an authentic, randomly-varied standard-format kinetic dataset.
    mode: 'hill' | 'substrate_inhibition'
    """
    rng = np.random
    if mode == "hill":
        vmax  = rng.uniform(80, 160)
        khalf = rng.uniform(10, 30)
        n     = round(rng.uniform(1.5, 3.0), 2)
        s_arr = np.array([round(x, 1) for x in np.geomspace(khalf*0.05, khalf*10, 9)])
        v_arr = _noisy(_hill_v(s_arr, vmax, khalf, n))
        label = f"Hill (n ≈ {n:.1f})"
    else:
        vmax = rng.uniform(100, 200)
        km   = rng.uniform(3, 12)
        ksi  = rng.uniform(40, 120)
        s_arr = np.array([round(x, 1) for x in np.geomspace(km*0.2, ksi*3, 10)])
        v_arr = _noisy(_ksi_v(s_arr, vmax, km, ksi))
        label = f"Substrate Inhibition (Ksi ≈ {ksi:.0f} μM)"

    df = pd.DataFrame({
        "Include?": [True]*len(s_arr),
        f"S [{c_unit}]": s_arr,
        "v": v_arr,
        "I": [0]*len(s_arr),
    })
    return df, label

def make_example_thermo():
    """
    Generate authentic thermodynamic data using the Arrhenius equation + noise.
    Ea and A values span the typical range for mesophilic enzymes.
    """
    rng  = np.random
    R    = 8.314
    Ea   = rng.uniform(40_000, 80_000)   # J/mol  (40–80 kJ/mol typical)
    A    = rng.uniform(1e10, 1e13)        # s⁻¹ pre-exponential
    temps_c = np.array([20, 25, 30, 35, 40, 45, 50], dtype=float)
    temps_k = temps_c + 273.15
    k_true   = A * np.exp(-Ea / (R * temps_k))
    k_noisy  = _noisy(k_true, cv=0.05)
    return pd.DataFrame({"Temp (T)": temps_c, "Rate (k)": k_noisy}), Ea/1000


# ═══════════════════════════════════════════════════════
#  PAGE CONFIG & PREMIUM THEME
# ═══════════════════════════════════════════════════════
st.set_page_config(
    page_title="BioKineticPy v1.1 — Enzyme Kinetics Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- PREMIUM SCIENTIFIC CSS THEME ---
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;600;700&family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

:root {
    --bg-primary:    hsl(222, 28%, 5%);
    --bg-secondary:  hsl(222, 22%, 8%);
    --bg-card:       hsl(222, 20%, 10%);
    --bg-hover:      hsl(222, 20%, 13%);
    --accent-cyan:   #38bdf8;
    --accent-teal:   #2dd4bf;
    --accent-violet: #a78bfa;
    --accent-emerald:#34d399;
    --accent-amber:  #fbbf24;
    --accent-rose:   #fb7185;
    --text-primary:  #f0f6fc;
    --text-secondary:#94a3b8;
    --text-muted:    #64748b;
    --border:        hsl(222, 15%, 18%);
    --border-accent: rgba(56,189,248,0.25);
    --gradient-hero: linear-gradient(135deg, #38bdf8 0%, #a78bfa 55%, #fb7185 100%);
    --gradient-teal: linear-gradient(135deg, #2dd4bf 0%, #38bdf8 100%);
    --glow-cyan:     0 0 24px rgba(56,189,248,0.18);
    --glow-violet:   0 0 24px rgba(167,139,250,0.18);
    --radius-card:   14px;
    --radius-btn:    10px;
}

html, body, [data-testid="stAppViewContainer"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
    color: var(--text-primary);
    background: var(--bg-primary);
}
[data-testid="stAppViewContainer"] { background: var(--bg-primary); }
[data-testid="stHeader"] { background: transparent !important; border-bottom: none !important; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: var(--bg-secondary) !important;
    border-right: 1px solid var(--border) !important;
}
[data-testid="stSidebar"] * { color: var(--text-primary) !important; }
[data-testid="stSidebar"] hr { border-color: var(--border) !important; }

/* hide default radio dots, style as nav pills */
[data-testid="stSidebar"] .stRadio > div { gap: 6px !important; }
[data-testid="stSidebar"] .stRadio label {
    display: flex !important;
    align-items: center !important;
    width: 100% !important;
    padding: 10px 14px !important;
    border-radius: 10px !important;
    border: 1px solid transparent !important;
    background: transparent !important;
    cursor: pointer !important;
    transition: all 0.2s ease !important;
    font-weight: 500 !important;
    font-size: 0.9rem !important;
    color: var(--text-secondary) !important;
}
[data-testid="stSidebar"] .stRadio label:hover {
    background: rgba(56,189,248,0.07) !important;
    border-color: var(--border-accent) !important;
    color: var(--text-primary) !important;
}
[data-testid="stSidebar"] .stRadio label[data-baseweb="radio"]:has(input:checked),
[data-testid="stSidebar"] .stRadio div[aria-checked="true"] ~ label,
[data-testid="stSidebar"] .stRadio input:checked + div {
    background: rgba(56,189,248,0.1) !important;
    border-color: var(--border-accent) !important;
    color: var(--accent-cyan) !important;
}

/* ── Cards ── */
[data-testid="stExpander"] {
    background: var(--bg-card) !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius-card) !important;
}
div[data-testid="stMetric"] {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: var(--radius-card);
    padding: 18px 22px;
    transition: transform 0.2s ease, box-shadow 0.2s ease, border-color 0.2s ease;
}
div[data-testid="stMetric"]:hover {
    transform: translateY(-3px);
    box-shadow: var(--glow-cyan);
    border-color: var(--border-accent);
}
div[data-testid="stMetric"] label {
    color: var(--text-secondary) !important;
    font-family: 'Space Grotesk', sans-serif !important;
    font-weight: 500 !important;
    font-size: 0.72rem !important;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
div[data-testid="stMetric"] [data-testid="stMetricValue"] {
    color: var(--accent-cyan) !important;
    font-family: 'JetBrains Mono', monospace !important;
    font-weight: 600 !important;
    font-size: 1.55rem !important;
}

/* ── Tabs ── */
.stTabs [data-baseweb="tab-list"] {
    background: var(--bg-card);
    border-radius: var(--radius-card);
    padding: 5px;
    gap: 4px;
    border: 1px solid var(--border);
}
.stTabs [data-baseweb="tab"] {
    border-radius: 9px !important;
    color: var(--text-secondary) !important;
    font-weight: 500 !important;
    font-size: 0.88rem !important;
    transition: all 0.2s ease !important;
    padding: 8px 18px !important;
}
.stTabs [aria-selected="true"] {
    background: rgba(56,189,248,0.12) !important;
    color: var(--accent-cyan) !important;
    font-weight: 600 !important;
}
.stTabs [data-baseweb="tab-highlight"] { background-color: var(--accent-cyan) !important; }

/* ── Buttons ── */
.stButton > button {
    background: linear-gradient(135deg, #38bdf8 0%, #0ea5e9 100%) !important;
    color: #fff !important; border: none !important;
    border-radius: var(--radius-btn) !important;
    font-family: 'Space Grotesk', sans-serif !important;
    font-weight: 600 !important; font-size: 0.92rem !important;
    padding: 0.6rem 1.8rem !important; letter-spacing: 0.03em;
    transition: all 0.22s ease !important;
    box-shadow: 0 4px 18px rgba(56,189,248,0.3) !important;
}
.stButton > button:hover {
    transform: translateY(-2px) !important;
    box-shadow: 0 8px 28px rgba(56,189,248,0.42) !important;
    filter: brightness(1.08) !important;
}
.stDownloadButton > button {
    background: var(--bg-card) !important;
    color: var(--accent-cyan) !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius-btn) !important;
    font-weight: 500 !important; transition: all 0.22s ease !important;
}
.stDownloadButton > button:hover {
    border-color: var(--accent-cyan) !important;
    background: rgba(56,189,248,0.08) !important;
    box-shadow: var(--glow-cyan) !important;
}

/* ── Slider ── */
[data-testid="stSlider"] [data-baseweb="slider"] [data-testid="stThumbValue"] {
    background: var(--accent-cyan) !important; color: #000 !important;
    font-family: 'JetBrains Mono', monospace !important; font-size: 0.75rem !important;
}

/* ── Inputs ── */
[data-testid="stNumberInput"] input,
[data-testid="stTextInput"] input,
.stSelectbox [data-baseweb="select"],
[data-baseweb="input"] {
    background: var(--bg-card) !important;
    border-color: var(--border) !important;
    color: var(--text-primary) !important;
    border-radius: 8px !important;
    font-family: 'JetBrains Mono', monospace !important;
}
[data-testid="stNumberInput"] input:focus,
[data-testid="stTextInput"] input:focus {
    border-color: var(--accent-cyan) !important;
    box-shadow: 0 0 0 2px rgba(56,189,248,0.15) !important;
}

/* ── Data / Tables ── */
[data-testid="stDataFrame"], .stDataFrame {
    border: 1px solid var(--border) !important;
    border-radius: var(--radius-card) !important;
    overflow: hidden;
}

/* ── Alerts ── */
.stAlert { border-radius: 10px !important; border: 1px solid var(--border) !important; }

/* ── Custom Components ── */
.hero-title {
    font-family: 'Space Grotesk', sans-serif;
    background: var(--gradient-hero);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
    background-clip: text;
    font-weight: 700; font-size: 2.1rem; line-height: 1.15; margin-bottom: 0;
}
.hero-subtitle {
    color: var(--text-secondary); font-size: 0.92rem;
    font-weight: 400; margin-top: 6px; line-height: 1.6;
}
.section-badge {
    display: inline-flex; align-items: center; gap: 5px;
    background: rgba(56,189,248,0.08);
    color: var(--accent-cyan);
    border: 1px solid rgba(56,189,248,0.22);
    border-radius: 20px; padding: 4px 14px;
    font-size: 0.72rem; font-weight: 600;
    letter-spacing: 0.07em; text-transform: uppercase; margin-bottom: 10px;
}
.section-header {
    display: flex; align-items: center; gap: 10px;
    padding: 14px 18px; margin: 0 0 16px 0;
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-left: 3px solid var(--accent-cyan);
    border-radius: var(--radius-card);
}
.section-header h3 {
    font-family: 'Space Grotesk', sans-serif;
    font-size: 0.95rem; font-weight: 600;
    color: var(--text-primary); margin: 0;
}
.section-header p {
    font-size: 0.8rem; color: var(--text-secondary); margin: 2px 0 0 0;
}
.result-banner {
    padding: 20px 24px; margin: 16px 0;
    background: linear-gradient(135deg, rgba(56,189,248,0.08) 0%, rgba(167,139,250,0.05) 100%);
    border: 1px solid rgba(56,189,248,0.25);
    border-radius: var(--radius-card);
    display: flex; flex-wrap: wrap; gap: 24px; align-items: center;
}
.result-banner .rb-model {
    font-family: 'Space Grotesk', sans-serif;
    font-size: 1.3rem; font-weight: 700;
    background: var(--gradient-hero);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;
}
.result-banner .rb-stat {
    display: flex; flex-direction: column; gap: 2px;
}
.result-banner .rb-label {
    font-size: 0.68rem; text-transform: uppercase;
    letter-spacing: 0.08em; color: var(--text-muted); font-weight: 600;
}
.result-banner .rb-value {
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.05rem; font-weight: 500; color: var(--text-primary);
}
.chip { display:inline-flex; align-items:center; gap:5px;
    padding: 5px 12px; border-radius: 20px;
    font-size: 0.78rem; font-weight: 600; letter-spacing: 0.03em;
}
.chip-pass { background: rgba(52,211,153,0.12); color: #34d399; border: 1px solid rgba(52,211,153,0.3); }
.chip-warn { background: rgba(251,191,36,0.12); color: #fbbf24; border: 1px solid rgba(251,191,36,0.3); }
.chip-fail { background: rgba(251,113,133,0.12); color: #fb7185; border: 1px solid rgba(251,113,133,0.3); }
.param-card {
    background: var(--bg-card); border: 1px solid var(--border);
    border-radius: var(--radius-card); padding: 16px 20px; margin: 8px 0;
}
.param-card .pc-name {
    font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.08em;
    color: var(--text-muted); font-weight: 600; margin-bottom: 4px;
}
.param-card .pc-value {
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.35rem; font-weight: 500; color: var(--accent-cyan);
}
.eq-display {
    background: var(--bg-card); border: 1px solid var(--border-accent);
    border-left: 3px solid var(--accent-teal);
    border-radius: var(--radius-card); padding: 16px 20px; margin: 10px 0;
    font-family: 'JetBrains Mono', monospace; font-size: 0.95rem;
    color: var(--accent-teal); line-height: 1.8;
}
.eq-display .eq-label {
    font-size: 0.68rem; text-transform: uppercase; letter-spacing: 0.08em;
    color: var(--text-muted); font-family: 'Inter', sans-serif; margin-bottom: 6px;
}
.thermo-card {
    background: var(--bg-card); border: 1px solid var(--border);
    border-radius: var(--radius-card); padding: 20px; text-align: center;
}
.thermo-card .tc-label {
    font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.08em;
    color: var(--text-muted); font-weight: 600; margin-bottom: 8px;
}
.thermo-card .tc-value {
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.5rem; font-weight: 600; color: var(--accent-teal);
}
.thermo-card .tc-unit {
    font-size: 0.8rem; color: var(--text-secondary); margin-top: 2px;
}
.sidebar-logo {
    font-family: 'Space Grotesk', sans-serif;
    font-size: 1.5rem; font-weight: 700;
    background: var(--gradient-hero);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;
}
.version-badge {
    display: inline-block; background: rgba(56,189,248,0.1);
    color: var(--accent-cyan); border: 1px solid rgba(56,189,248,0.2);
    border-radius: 20px; padding: 2px 10px;
    font-size: 0.68rem; font-weight: 600; letter-spacing: 0.05em;
}

/* Plotly SVG */
.js-plotly-plot .plotly .main-svg { border-radius: 12px; }

/* Dividers & Scrollbar */
hr { border-color: var(--border) !important; }
::-webkit-scrollbar { width: 5px; }
::-webkit-scrollbar-track { background: var(--bg-primary); }
::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: var(--text-muted); }
</style>
""", unsafe_allow_html=True)

# ── Plotly Dark Theme Template ──
PLOTLY_DARK = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(13,17,28,0.7)",
    font=dict(family="'JetBrains Mono', monospace", color="#f0f6fc", size=12),
    title=dict(font=dict(family="'Space Grotesk', sans-serif", size=15, color="#94a3b8"), x=0.02, xanchor="left"),
    xaxis=dict(gridcolor="rgba(255,255,255,0.04)", zerolinecolor="rgba(56,189,248,0.2)", tickfont=dict(size=11)),
    yaxis=dict(gridcolor="rgba(255,255,255,0.04)", zerolinecolor="rgba(56,189,248,0.2)", tickfont=dict(size=11)),
    colorway=["#38bdf8", "#34d399", "#a78bfa", "#fbbf24", "#fb7185", "#f472b6", "#60a5fa", "#4ade80"],
    margin=dict(l=55, r=25, t=55, b=55),
    hoverlabel=dict(
        bgcolor="rgba(22,27,42,0.95)",
        bordercolor="rgba(56,189,248,0.3)",
        font=dict(family="'JetBrains Mono', monospace", size=12, color="#f0f6fc")
    ),
    legend=dict(
        bgcolor="rgba(13,17,28,0.7)",
        bordercolor="rgba(56,189,248,0.15)",
        borderwidth=1,
        font=dict(size=11)
    )
)

# ═══════════════════════════════════════════════════════
#  UTILS
# ═══════════════════════════════════════════════════════
def reset_editor_key():
    """Forces the data_editor to reload by changing its ID."""
    st.session_state.editor_key = random.randint(0, 100000)

if "editor_key" not in st.session_state:
    st.session_state.editor_key = 0

# ═══════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════
def generate_synthetic_curve(model_name, s_range, params):
    func = MODEL_REGISTRY[model_name]
    valid_args = func.__code__.co_varnames
    call_args = {k: v for k, v in params.items() if k in valid_args}
    if 'i' in valid_args and 'i' not in call_args: call_args['i'] = 0.0
    if 'ki' in valid_args and 'ki' not in call_args: call_args['ki'] = 1.0 
    return func(s=s_range, **call_args)

def generate_3d_surface(model_name, s_grid, i_grid, params):
    func = MODEL_REGISTRY[model_name]
    valid_args = func.__code__.co_varnames
    call_args = {k: v for k, v in params.items() if k in valid_args}
    call_args['s'] = s_grid
    if 'i' in valid_args: call_args['i'] = i_grid
    return func(**call_args)

def calculate_residuals(points, model_name, params):
    s_vals = np.array([p.substrate_conc for p in points])
    v_obs = np.array([p.velocity for p in points])
    v_pred = []
    func = MODEL_REGISTRY[model_name]
    valid_args = func.__code__.co_varnames
    for p in points:
        call_args = {k: v for k, v in params.items() if k in valid_args}
        call_args['s'] = p.substrate_conc
        if 'i' in valid_args: call_args['i'] = p.inhibitor_conc
        v_pred.append(func(**call_args))
    return s_vals, v_obs, np.array(v_pred), v_obs - np.array(v_pred)

def apply_plotly_theme(fig, height=500):
    """Apply the unified dark theme to a plotly figure."""
    fig.update_layout(**PLOTLY_DARK, height=height)
    return fig

# ═══════════════════════════════════════════════════════
#  SIDEBAR
# ═══════════════════════════════════════════════════════
with st.sidebar:
    # ── Logo & Brand ──
    st.markdown("""
    <div style="padding: 10px 4px 4px 4px;">
        <p class="sidebar-logo">🧬 BioKineticPy</p>
        <p style="color:var(--text-secondary);font-size:0.8rem;margin-top:2px;">
            Professional Enzyme Kinetics Suite
        </p>
        <span class="version-badge">v1.1</span>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("---")

    # ── Mode Navigation ──
    st.markdown('<span class="section-badge">🔬 Mode</span>', unsafe_allow_html=True)
    mode = st.radio(
        "mode",
        ["🔬  Analysis (Fit Data)", "📈  Simulation", "🌡️  Thermodynamics", "📖  Guide"],
        label_visibility="collapsed"
    )
    # Normalise to short names for downstream code
    if "Analysis" in mode:       mode = "Analysis (Fit Data)"
    elif "Simulation" in mode:   mode = "Predictive (Simulation)"
    elif "Thermodynamics" in mode: mode = "Thermodynamics"
    else:                          mode = "Guide"

    st.markdown("---")

    # ── Project Save/Load ──
    st.markdown('<span class="section-badge">💾 Project</span>', unsafe_allow_html=True)
    current_state = {
        "df_matrix":   st.session_state.get("df_matrix",   pd.DataFrame()).to_dict() if "df_matrix"   in st.session_state else None,
        "df_standard": st.session_state.get("df_standard", pd.DataFrame()).to_dict() if "df_standard" in st.session_state else None,
        "df_thermo":   st.session_state.get("df_thermo",   pd.DataFrame()).to_dict() if "df_thermo"   in st.session_state else None
    }
    st.download_button("💾  Save Project", json.dumps(current_state), "biokinetic_project.json", "application/json")
    uploaded_proj = st.file_uploader("Load Project (.json)", type=["json"], label_visibility="collapsed")
    if uploaded_proj:
        try:
            data = json.load(uploaded_proj)
            if data["df_matrix"]:   st.session_state.df_matrix   = pd.DataFrame.from_dict(data["df_matrix"])
            if data["df_standard"]: st.session_state.df_standard = pd.DataFrame.from_dict(data["df_standard"])
            if data.get("df_thermo"): st.session_state.df_thermo = pd.DataFrame.from_dict(data["df_thermo"])
            reset_editor_key()
            st.success("✅ Project loaded!")
        except Exception:
            st.error("Failed to parse project file.")

    if mode != "Thermodynamics":
        st.markdown("---")
        with st.expander("⚙️ Fit Settings", expanded=True):
            weighting_mode = st.selectbox("Weighting", ["None (Homoscedastic)", "1/v (Poisson)", "1/v² (Relative Error)"], label_visibility="collapsed")
        with st.expander("📏 Units", expanded=True):
            c_unit = st.selectbox("Concentration", ["mM", "μM", "nM", "M", "pM"], index=1)
            t_unit = st.selectbox("Time Base", ["min", "sec"], index=0)
            v_unit_display = f"{c_unit}/{t_unit}"
            st.markdown(f'<span class="section-badge" style="background:rgba(45,212,191,0.08);color:#2dd4bf;border-color:rgba(45,212,191,0.2);">v  →  {v_unit_display}</span>', unsafe_allow_html=True)
        with st.expander("🧫 Enzyme Properties"):
            use_et = st.checkbox("Define [E]ₜ (enzyme total)")
            et_val  = st.number_input("Enzyme Conc", 1.0) if use_et else 0
            et_unit = st.selectbox("Enzyme Unit", ["nM", "μM", "mM", "M"], 0) if use_et else "nM"

# ═══════════════════════════════════════════════════════
#  ANALYSIS MODE
# ═══════════════════════════════════════════════════════
if mode == "Analysis (Fit Data)":
    st.markdown('<p class="hero-title">🔬 Kinetic Analysis</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Automated multi-model fitting · AICc model selection · Residual diagnostics · Bootstrap CI</p>', unsafe_allow_html=True)
    st.markdown("")

    col1, col2 = st.columns([1.2, 2])
    
    with col1:
        st.markdown('<span class="section-badge">📥 Data Input</span>', unsafe_allow_html=True)
        if "input_format" not in st.session_state:
            st.session_state["input_format"] = "Matrix"
        input_format = st.radio("Layout", ["Matrix", "Standard"], horizontal=True, key="input_format")
        uploaded_file = st.file_uploader("Upload Data (CSV / Excel)", type=["csv", "xlsx"])

        # DEMO LOADERS
        with st.expander("📂 Load Example Data"):
            st.caption("Each click generates freshly randomised authentic data with realistic noise (~8% CV).")
            ex_r1c1, ex_r1c2, ex_r1c3 = st.columns(3)
            ex_r2c1, ex_r2c2, _ = st.columns(3)

            if ex_r1c1.button("Michaelis-Menten", use_container_width=True):
                df, i_a, i_b, lbl = make_example_matrix("michaelis", c_unit)
                st.session_state.df_matrix = df
                st.session_state.i_conc_A = float(i_a)
                st.session_state.i_conc_B = float(i_b)
                st.session_state.example_label = lbl
                reset_editor_key(); st.rerun()

            if ex_r1c2.button("Competitive", use_container_width=True):
                df, i_a, i_b, lbl = make_example_matrix("competitive", c_unit)
                st.session_state.df_matrix = df
                st.session_state.i_conc_A = float(i_a)
                st.session_state.i_conc_B = float(i_b)
                st.session_state.example_label = lbl
                reset_editor_key(); st.rerun()

            if ex_r1c3.button("Non-Competitive", use_container_width=True):
                df, i_a, i_b, lbl = make_example_matrix("noncompetitive", c_unit)
                st.session_state.df_matrix = df
                st.session_state.i_conc_A = float(i_a)
                st.session_state.i_conc_B = float(i_b)
                st.session_state.example_label = lbl
                reset_editor_key(); st.rerun()

            if ex_r2c1.button("Hill (Allosteric)", use_container_width=True):
                df, lbl = make_example_standard("hill", c_unit)
                st.session_state.df_standard = df
                st.session_state.example_label = lbl
                st.session_state["input_format"] = "Standard"
                reset_editor_key(); st.rerun()

            if ex_r2c2.button("Substrate Inhibition", use_container_width=True):
                df, lbl = make_example_standard("substrate_inhibition", c_unit)
                st.session_state.df_standard = df
                st.session_state.example_label = lbl
                st.session_state["input_format"] = "Standard"
                reset_editor_key(); st.rerun()

            if "example_label" in st.session_state:
                st.success(f"✅ Loaded: **{st.session_state.example_label}** — re-click to regenerate with new parameters")

        if input_format == "Matrix":
            if "df_matrix" not in st.session_state:
                st.session_state.df_matrix = pd.DataFrame({
                    "Include?": [True]*5,
                    f"Substrate [{c_unit}]": [1, 5, 10, 20, 50],
                    "v (I=0)": [10.5, 35.1, 55.0, 70.2, 90.5],
                    "v (I=A)": [4.8, 19.5, 32.9, 45.0, 60.1],
                    "v (I=B)": [None]*5, "v (I=C)": [None]*5
                })
            if uploaded_file:
                try: 
                    df = pd.read_csv(uploaded_file) if uploaded_file.name.endswith('.csv') else pd.read_excel(uploaded_file)
                    if "Include?" not in df.columns: df.insert(0, "Include?", True)
                    st.session_state.df_matrix = df.iloc[:, :6]
                    reset_editor_key()
                except Exception:
                    st.error("Could not parse file. Check format.")
            
            ic1, ic2, ic3 = st.columns(3)
            i_conc_A = ic1.number_input(f"[I] A ({c_unit})", value=5.0)
            i_conc_B = ic2.number_input(f"[I] B ({c_unit})", value=10.0)
            i_conc_C = ic3.number_input(f"[I] C ({c_unit})", value=20.0)
            
            edited_df = st.data_editor(
                st.session_state.df_matrix, 
                num_rows="dynamic", use_container_width=True,
                column_config={"Include?": st.column_config.CheckboxColumn(default=True)},
                key=f"editor_matrix_{st.session_state.editor_key}"
            )
        
        else: 
            if "df_standard" not in st.session_state:
                st.session_state.df_standard = pd.DataFrame({
                    "Include?": [True]*5, f"S [{c_unit}]": [1, 5, 10, 1, 5], f"v": [10, 35, 55, 4, 19], f"I": [0, 0, 0, 5, 5]
                })
            if uploaded_file:
                try:
                    df = pd.read_csv(uploaded_file) if uploaded_file.name.endswith('.csv') else pd.read_excel(uploaded_file)
                    if "Include?" not in df.columns: df.insert(0, "Include?", True)
                    st.session_state.df_standard = df.iloc[:, :4]
                    reset_editor_key()
                except Exception:
                    st.error("Could not parse file. Check format.")

            edited_df = st.data_editor(
                st.session_state.df_standard, 
                num_rows="dynamic", use_container_width=True,
                column_config={"Include?": st.column_config.CheckboxColumn(default=True)},
                key=f"editor_std_{st.session_state.editor_key}"
            )

        run_fit = st.button("🚀  Run Analysis", type="primary", use_container_width=True)

    with col2:
        if run_fit or st.session_state.get("fit_done", False):
            try:
                if run_fit:
                    st.session_state.fit_done = True
                    if "ci_data_user" in st.session_state: del st.session_state.ci_data_user

                points = []
                c_factor = CONC_TO_MOLAR[c_unit]
                def norm_s(val): return float(val) * c_factor
                def norm_v(val): return convert_velocity_to_standard(float(val), c_unit, t_unit)
                
                if input_format == "Matrix":
                    cols = edited_df.columns
                    inh_map = {2: 0.0, 3: i_conc_A, 4: i_conc_B, 5: i_conc_C}
                    for _, row in edited_df.iterrows():
                        if not row["Include?"] or pd.isna(row[cols[1]]): continue
                        s_val = norm_s(row[cols[1]])
                        for c_idx, i_val in inh_map.items():
                            if c_idx < len(cols) and not pd.isna(row[cols[c_idx]]):
                                points.append(KineticPoint(substrate_conc=s_val, velocity=norm_v(row[cols[c_idx]]), inhibitor_conc=norm_s(i_val)))
                else:
                    cols = edited_df.columns
                    for _, row in edited_df.iterrows():
                        if not row["Include?"] or pd.isna(row[cols[1]]): continue
                        points.append(KineticPoint(substrate_conc=norm_s(row[cols[1]]), velocity=norm_v(row[cols[2]]), inhibitor_conc=norm_s(row.get(cols[3], 0))))

                if len(points) < 3: st.error("❌ Insufficient data — need at least 3 points."); st.stop()

                dataset = Dataset(points=points)
                solver = KineticSolver(dataset)
                results = solver.run_model_competition(weighting_mode)
                best_model = results[0]
                ambiguity_msg = solver.check_ambiguity(results)
                
                if run_fit:
                    n_pts = len(points); n_params = len(best_model['parameters'])
                    aicc = best_model['aic'] + (2*n_params*(n_params+1))/(max(n_pts-n_params-1,1))
                    st.markdown(f"""
                    <div class="result-banner">
                        <div>
                            <div style="font-size:0.68rem;text-transform:uppercase;letter-spacing:0.08em;color:var(--text-muted);font-weight:600;">Best Model</div>
                            <div class="rb-model">{best_model['model']}</div>
                        </div>
                        <div class="rb-stat"><div class="rb-label">AICc</div><div class="rb-value">{aicc:.2f}</div></div>
                        <div class="rb-stat"><div class="rb-label">R²</div><div class="rb-value">{best_model['r_squared']:.4f}</div></div>
                        <div class="rb-stat"><div class="rb-label">Data Points</div><div class="rb-value">{n_pts}</div></div>
                    </div>
                    """, unsafe_allow_html=True)
                    if ambiguity_msg: st.warning(ambiguity_msg)
                
                disp_params = {}
                for k, v in best_model['parameters'].items():
                    if k in ['vmax']: disp_params[k] = convert_param_to_user(v, 'velocity', c_unit, t_unit)
                    elif k in ['km', 'ki', 'khalf', 'ksi']: disp_params[k] = convert_param_to_user(v, 'concentration', c_unit, t_unit)
                    else: disp_params[k] = v 
                
                s_vals, v_obs, v_pred, residuals = calculate_residuals(points, best_model['model'], best_model['parameters'])
                stat_tests = solver.diagnose_residuals(residuals)
                derived_stats = {}
                if use_et and et_val > 0:
                    vmax_m_s = best_model['parameters']['vmax']
                    et_m = et_val * CONC_TO_MOLAR[et_unit]
                    if et_m > 0: derived_stats = {"kcat": vmax_m_s/et_m, "efficiency": (vmax_m_s/et_m)/best_model['parameters'].get('km', 1e-9)}

                if st.button("📊  Compute 95% Confidence Intervals"):
                    with st.spinner("Running Monte Carlo bootstrap…"):
                        ci_raw = solver.bootstrap_uncertainty(best_model['model'], best_model['parameters'], 50, weighting_mode)
                        ci_converted = {}
                        for k, (l, h) in ci_raw.items():
                            if k in ['vmax']: ci_converted[k] = (convert_param_to_user(l, 'velocity', c_unit, t_unit), convert_param_to_user(h, 'velocity', c_unit, t_unit))
                            elif k in ['km', 'ki', 'khalf', 'ksi']: ci_converted[k] = (convert_param_to_user(l, 'concentration', c_unit, t_unit), convert_param_to_user(h, 'concentration', c_unit, t_unit))
                            else: ci_converted[k] = (l, h)
                        st.session_state.ci_data_user = ci_converted
                        st.success("✅ Confidence Intervals computed.")

                report = generate_methods_report(
                    {"model": best_model['model'], "aic": best_model['aic'], "r_squared": best_model['r_squared'], "parameters": disp_params},
                    len(points), c_unit, v_unit_display, weighting_mode, stat_tests, derived_stats, 
                    st.session_state.get("ci_data_user", None)
                )
                st.download_button("📄  Download Full Report", report, "biokinetic_report.txt")

                user_s = np.array([convert_param_to_user(p.substrate_conc, 'concentration', c_unit, t_unit) for p in points])
                user_v = np.array([convert_param_to_user(p.velocity, 'velocity', c_unit, t_unit) for p in points])
                user_i = np.array([convert_param_to_user(p.inhibitor_conc, 'concentration', c_unit, t_unit) for p in points])
                unique_i = sorted(list(set(user_i)))

                tab1, tab2, tab3, tab4, tab5 = st.tabs(["📈 Fit", "🩺 Diagnostics", "📐 Linearizations", "🧪 Inhibition", "🌐 3D Landscape"])
                
                with tab1:
                    st.markdown("**Fitted Parameters**")
                    param_table = []
                    ci_info = st.session_state.get("ci_data_user", {})
                    for k, v in disp_params.items():
                        row = {"Parameter": k.upper(), "Value": v}
                        if ci_info and k in ci_info:
                            row["Lower 95%"] = ci_info[k][0]; row["Upper 95%"] = ci_info[k][1]
                        param_table.append(row)
                    st.dataframe(pd.DataFrame(param_table), hide_index=True, use_container_width=True, column_config={"Value": st.column_config.NumberColumn(format="%.4f")})
                    if use_et and derived_stats: st.metric("kcat (Turnover)", f"{derived_stats.get('kcat', 0):.2f} s⁻¹")

                    fig = go.Figure()
                    s_smooth = np.linspace(0, max(user_s)*1.1, 100)
                    s_molar = s_smooth * c_factor
                    for i_val in unique_i:
                        idx = [j for j, val in enumerate(user_i) if val == i_val]
                        fig.add_trace(go.Scatter(x=user_s[idx], y=user_v[idx], mode='markers', name=f"Data [I]={i_val:.1f}", marker=dict(size=9, line=dict(width=1, color='rgba(255,255,255,0.3)'))))
                        p_calc = best_model['parameters'].copy(); p_calc['i'] = i_val * c_factor
                        v_m = generate_synthetic_curve(best_model['model'], s_molar, p_calc)
                        v_u = [convert_param_to_user(v, 'velocity', c_unit, t_unit) for v in v_m]
                        fig.add_trace(go.Scatter(x=s_smooth, y=v_u, mode='lines', name=f"Fit [I]={i_val:.1f}", line=dict(width=2.5)))
                    fig.update_layout(xaxis_title=f"[S] ({c_unit})", yaxis_title=f"Velocity ({v_unit_display})", title=f"Kinetic Fit — {best_model['model']}")
                    st.plotly_chart(apply_plotly_theme(fig, 480), use_container_width=True)

                with tab2:
                    sw_p = stat_tests.get('shapiro_p', 0); rz = stat_tests.get('runs_z', 0)
                    r2_val = best_model['r_squared']
                    norm_chip  = 'chip-pass' if sw_p > 0.05 else 'chip-fail'
                    norm_txt   = f'✅ Normal (p={sw_p:.3f})' if sw_p > 0.05 else f'⚠️ Non-Normal (p={sw_p:.3f})'
                    rand_chip  = 'chip-pass' if abs(rz) < 1.96 else 'chip-warn'
                    rand_txt   = f'✅ Random (Z={rz:.2f})' if abs(rz) < 1.96 else f'⚠️ Systematic (Z={rz:.2f})'
                    r2_chip    = 'chip-pass' if r2_val > 0.95 else ('chip-warn' if r2_val > 0.80 else 'chip-fail')
                    r2_txt     = f'R² = {r2_val:.4f}'
                    st.markdown(f"""
                    <div style="display:flex;gap:12px;flex-wrap:wrap;margin:12px 0 20px 0;">
                        <div>
                            <div style="font-size:0.68rem;color:var(--text-muted);text-transform:uppercase;letter-spacing:0.07em;margin-bottom:5px;">Normality (Shapiro-Wilk)</div>
                            <span class="chip {norm_chip}">{norm_txt}</span>
                        </div>
                        <div>
                            <div style="font-size:0.68rem;color:var(--text-muted);text-transform:uppercase;letter-spacing:0.07em;margin-bottom:5px;">Randomness (Runs Test)</div>
                            <span class="chip {rand_chip}">{rand_txt}</span>
                        </div>
                        <div>
                            <div style="font-size:0.68rem;color:var(--text-muted);text-transform:uppercase;letter-spacing:0.07em;margin-bottom:5px;">Goodness of Fit</div>
                            <span class="chip {r2_chip}">{r2_txt}</span>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    fig_diag = make_subplots(rows=1, cols=2, subplot_titles=("Residuals vs [S]", "Q-Q Plot"))
                    fig_diag.add_trace(go.Scatter(x=s_vals, y=residuals, mode='markers', marker=dict(color='#f85149', size=8, line=dict(width=1, color='rgba(255,255,255,0.2)'))), row=1, col=1)
                    fig_diag.add_hline(y=0, line_dash="dash", line_color="#58a6ff", row=1, col=1)
                    res_sorted = np.sort(residuals); theo = stats.norm.ppf(np.linspace(0.01, 0.99, len(res_sorted)))
                    fig_diag.add_trace(go.Scatter(x=theo, y=res_sorted, mode='markers', marker=dict(color='#bc8cff', size=8)), row=1, col=2)
                    # add reference line for QQ
                    qq_min, qq_max = min(theo), max(theo)
                    fig_diag.add_trace(go.Scatter(x=[qq_min, qq_max], y=[qq_min * np.std(residuals) + np.mean(residuals), qq_max * np.std(residuals) + np.mean(residuals)], mode='lines', line=dict(dash='dash', color='#58a6ff'), showlegend=False), row=1, col=2)
                    st.plotly_chart(apply_plotly_theme(fig_diag, 400), use_container_width=True)

                with tab3:
                    lt1, lt2, lt3 = st.tabs(["📐 Lineweaver-Burk", "📐 Hanes-Woolf", "📐 Eadie-Hofstee"])

                    with lt1:
                        fig_lb = go.Figure()
                        for i_val in unique_i:
                            idx = [j for j, val in enumerate(user_i) if val == i_val]
                            valid = [j for j in idx if user_s[j] > 0 and user_v[j] > 0]
                            fig_lb.add_trace(go.Scatter(x=1/user_s[valid], y=1/user_v[valid], mode='markers', name=f"[I]={i_val:.1f}"))
                            p_calc = best_model['parameters'].copy(); p_calc['i'] = i_val * c_factor
                            v_m = generate_synthetic_curve(best_model['model'], s_molar, p_calc)
                            v_u = np.array([convert_param_to_user(v, 'velocity', c_unit, t_unit) for v in v_m])
                            mask = (s_smooth > 1e-9) & (v_u > 1e-9)
                            fig_lb.add_trace(go.Scatter(x=1.0/s_smooth[mask], y=1.0/v_u[mask], mode='lines', name=f"Fit [I]={i_val:.1f}"))
                        fig_lb.update_layout(title="Lineweaver-Burk", xaxis_title=f"1/[S]  (1/{c_unit})", yaxis_title=f"1/v  (1/({v_unit_display}))")
                        st.plotly_chart(apply_plotly_theme(fig_lb, 480), use_container_width=True)
                        st.caption("Double-reciprocal plot. Lines that converge on the y-axis → competitive inhibition. Parallel lines → uncompetitive.")

                    with lt2:
                        fig_hw = go.Figure()
                        for i_val in unique_i:
                            idx = [j for j, val in enumerate(user_i) if val == i_val]
                            valid = [j for j in idx if user_v[j] > 0]
                            fig_hw.add_trace(go.Scatter(x=user_s[valid], y=user_s[valid]/user_v[valid], mode='markers', name=f"[I]={i_val:.1f}"))
                            p_calc = best_model['parameters'].copy(); p_calc['i'] = i_val * c_factor
                            v_m = generate_synthetic_curve(best_model['model'], s_molar, p_calc)
                            v_u = np.array([convert_param_to_user(v, 'velocity', c_unit, t_unit) for v in v_m])
                            mask = (v_u > 1e-9)
                            fig_hw.add_trace(go.Scatter(x=s_smooth[mask], y=s_smooth[mask]/v_u[mask], mode='lines', name=f"Fit [I]={i_val:.1f}"))
                        fig_hw.update_layout(title="Hanes-Woolf", xaxis_title=f"[S]  ({c_unit})", yaxis_title=f"[S]/v  ({c_unit}/({v_unit_display}))")
                        st.plotly_chart(apply_plotly_theme(fig_hw, 480), use_container_width=True)
                        st.caption("More robust to error at low [S] than Lineweaver-Burk. Slope = 1/Vmax; x-intercept = −Km.")

                    with lt3:
                        fig_eh = go.Figure()
                        for i_val in unique_i:
                            idx = [j for j, val in enumerate(user_i) if val == i_val]
                            valid = [j for j in idx if user_s[j] > 0]
                            fig_eh.add_trace(go.Scatter(x=user_v[valid]/user_s[valid], y=user_v[valid], mode='markers', name=f"[I]={i_val:.1f}"))
                            p_calc = best_model['parameters'].copy(); p_calc['i'] = i_val * c_factor
                            v_m = generate_synthetic_curve(best_model['model'], s_molar, p_calc)
                            v_u = np.array([convert_param_to_user(v, 'velocity', c_unit, t_unit) for v in v_m])
                            mask = (s_smooth > 1e-9)
                            fig_eh.add_trace(go.Scatter(x=v_u[mask]/s_smooth[mask], y=v_u[mask], mode='lines', name=f"Fit [I]={i_val:.1f}"))
                        fig_eh.update_layout(title="Eadie-Hofstee", xaxis_title=f"v/[S]  ({v_unit_display}/{c_unit})", yaxis_title=f"v  ({v_unit_display})")
                        st.plotly_chart(apply_plotly_theme(fig_eh, 480), use_container_width=True)
                        st.caption("Best for visually distinguishing inhibition types. Slope = −Km; y-intercept = Vmax; concave curves indicate cooperativity.")

                with tab4:
                    it1, it2 = st.tabs(["🧪 Dixon Plot", "🧪 Cornish-Bowden"])

                    with it1:
                        fig_dixon = go.Figure()
                        for s_val in sorted(list(set([round(s, 6) for s in user_s]))):
                            idx = [j for j, s in enumerate(user_s) if abs(s - s_val) < 1e-9 and user_v[j] > 0]
                            if len(idx) > 1:
                                fig_dixon.add_trace(go.Scatter(x=user_i[idx], y=1/user_v[idx], mode='markers+lines', name=f"[S]={s_val:.1f}"))
                        fig_dixon.update_layout(title="Dixon Plot", xaxis_title=f"[I]  ({c_unit})", yaxis_title=f"1/v  (1/({v_unit_display}))")
                        st.plotly_chart(apply_plotly_theme(fig_dixon, 480), use_container_width=True)
                        st.caption("1/v vs [I] at each fixed [S]. Lines converging above the x-axis → competitive inhibition. Converging below → uncompetitive. The x-intercept of the intersection gives −Ki.")

                    with it2:
                        fig_cb = go.Figure()
                        for s_val in sorted(list(set([round(s, 6) for s in user_s]))):
                            idx = [j for j, s in enumerate(user_s) if abs(s - s_val) < 1e-9 and user_v[j] > 0]
                            if len(idx) > 1:
                                fig_cb.add_trace(go.Scatter(x=user_i[idx], y=s_val/user_v[idx], mode='markers+lines', name=f"[S]={s_val:.1f}"))
                        fig_cb.update_layout(title="Cornish-Bowden", xaxis_title=f"[I]  ({c_unit})", yaxis_title=f"[S]/v  ({c_unit}/({v_unit_display}))")
                        st.plotly_chart(apply_plotly_theme(fig_cb, 480), use_container_width=True)
                        st.caption("[S]/v vs [I] at each fixed [S]. More robust than Dixon — unaffected by competitive inhibition artefacts. Parallel lines → competitive. Converging lines → non-competitive/mixed.")

                with tab5:
                    st.markdown("**Interactive Velocity Surface**")
                    s_range_3d = np.linspace(0, max(user_s)*1.2, 30)
                    i_range_3d = np.linspace(0, max(user_i)*1.2, 30)
                    S_MESH, I_MESH = np.meshgrid(s_range_3d, i_range_3d)
                    S_MOLAR = S_MESH * c_factor; I_MOLAR = I_MESH * c_factor
                    params_3d = best_model['parameters'].copy()
                    Z_MOLAR = generate_3d_surface(best_model['model'], S_MOLAR, I_MOLAR, params_3d)
                    Z_USER = convert_param_to_user(Z_MOLAR, 'velocity', c_unit, t_unit)
                    fig_3d = go.Figure(data=[go.Surface(z=Z_USER, x=s_range_3d, y=i_range_3d, colorscale='Viridis', opacity=0.85)])
                    fig_3d.add_trace(go.Scatter3d(x=user_s, y=user_i, z=user_v, mode='markers', marker=dict(size=5, color='#f85149', line=dict(width=1, color='white')), name='Experimental'))
                    fig_3d.update_layout(title=f"3D Kinetic Landscape — {best_model['model']}", scene=dict(xaxis_title=f"[S] ({c_unit})", yaxis_title=f"[I] ({c_unit})", zaxis_title=f"v ({v_unit_display})"))
                    st.plotly_chart(apply_plotly_theme(fig_3d, 700), use_container_width=True)

            except Exception as e: st.error(f"⚠️ Analysis Failed: {str(e)}")

# ═══════════════════════════════════════════════════════
#  PREDICTIVE MODE
# ═══════════════════════════════════════════════════════
elif mode == "Predictive (Simulation)":
    st.markdown('<p class="hero-title">📈 Kinetic Simulator</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Explore enzyme mechanisms · Tune parameters interactively · Visualise velocity surfaces</p>', unsafe_allow_html=True)
    st.markdown("")

    col_input, col_graph = st.columns([1, 3])
    with col_input:
        st.markdown('<span class="section-badge">🔧 Parameters</span>', unsafe_allow_html=True)
        model_select = st.selectbox("🔧 Mechanism", list(MODEL_REGISTRY.keys()))

        # ── Kinetic Equation Display ──
        EQ_MAP = {
            "Michaelis-Menten":      "v = Vmax · [S] / (Km + [S])",
            "Competitive":           "v = Vmax · [S] / (Km·(1 + [I]/Ki) + [S])",
            "Uncompetitive":         "v = Vmax · [S] / (Km + [S]·(1 + [I]/Ki))",
            "Non-Competitive":       "v = (Vmax/(1+[I]/Ki)) · [S] / (Km + [S])",
            "Mixed":                 "v = Vmax’ · [S] / (Km’ + [S])  [α-mixed]",
            "Hill (Allosteric)":     "v = Vmax · [S]ⁿ / (Khalfⁿ + [S]ⁿ)",
            "Substrate Inhibition":  "v = Vmax · [S] / (Km + [S] + [S]²/Ksi)",
        }
        st.markdown(f"""
        <div class="eq-display">
            <div class="eq-label">Rate Equation</div>
            {EQ_MAP.get(model_select, "v = f([S])")}
        </div>
        """, unsafe_allow_html=True)
        vmax = st.slider(f"Vmax ({v_unit_display})", 0.0, 200.0, 100.0)
        km = 0; ki = 0; alpha = 1; ksi = 0; n = 1; i_conc = 0; khalf = 0
        
        if "Hill" in model_select:
            khalf = st.slider(f"Khalf ({c_unit})", 0.1, 50.0, 10.0); n = st.slider("Hill Coeff (n)", 0.1, 5.0, 1.0)
        elif "Substrate Inhibition" in model_select:
            km = st.slider(f"Km ({c_unit})", 0.1, 50.0, 10.0); ksi = st.slider(f"Ksi ({c_unit})", 0.1, 200.0, 50.0)
            st.markdown("---")
            st.markdown("**External Inhibitor**")
            i_conc = st.number_input(f"[I] ({c_unit})", 0.0, 100.0, 0.0); ki = st.slider(f"Ki ({c_unit})", 0.1, 50.0, 5.0)
        else:
            km = st.slider(f"Km ({c_unit})", 0.1, 50.0, 10.0)
            if model_select != "Michaelis-Menten": i_conc = st.number_input(f"[I] ({c_unit})", 0.0, 100.0, 10.0); ki = st.slider(f"Ki ({c_unit})", 0.1, 50.0, 5.0)
            if model_select == "Mixed": alpha = st.slider("Alpha", 0.1, 10.0, 1.0)
        st.markdown("---")
        sim_params = {"vmax": vmax, "km": km, "ki": ki, "alpha": alpha, "i": i_conc, "khalf": khalf, "n": n, "ksi": ksi}
        st.download_button("📄  Download Parameters", generate_simulation_report(model_select, sim_params, c_unit, v_unit_display), "biokineticpy_simulation.txt", use_container_width=True)

    with col_graph:
        tab_2d, tab_3d = st.tabs(["📈 2D Plot", "🌐 3D Landscape"])
        c_fac = CONC_TO_MOLAR[c_unit]; v_fac = c_fac / (60 if t_unit == 'min' else 1)
        params_molar = {"vmax": vmax * v_fac, "km": km * c_fac, "ki": ki * c_fac, "alpha": alpha, "i": i_conc * c_fac, "khalf": khalf * c_fac, "n": n, "ksi": ksi * c_fac}
        
        with tab_2d:
            s_sim = np.linspace(0, 100, 200)
            v_sim_user = generate_synthetic_curve(model_select, s_sim * c_fac, params_molar) / v_fac 
            fig = go.Figure()
            if "Substrate Inhibition" in model_select:
                 if i_conc > 0: fig.add_trace(go.Scatter(x=s_sim, y=generate_synthetic_curve(model_select, s_sim * c_fac, {**params_molar, 'i': 0.0}) / v_fac, mode='lines', name="No External Inhibitor", line=dict(color='#6e7681', dash='dash')))
                 else: fig.add_trace(go.Scatter(x=s_sim, y=(vmax*s_sim)/(km+s_sim), mode='lines', name="No Substrate Inhibition", line=dict(color='#6e7681', dash='dash')))
            elif "Hill" in model_select: fig.add_trace(go.Scatter(x=s_sim, y=(vmax*s_sim)/(khalf+s_sim), mode='lines', name="Hyperbolic (n=1)", line=dict(color='#6e7681', dash='dash')))
            else: fig.add_trace(go.Scatter(x=s_sim, y=generate_synthetic_curve(model_select, s_sim * c_fac, {**params_molar, 'i': 0.0}) / v_fac, mode='lines', name="No Inhibitor", line=dict(color='#6e7681', dash='dash')))
            fig.add_trace(go.Scatter(x=s_sim, y=v_sim_user, mode='lines', name=f"{model_select}", line=dict(width=3)))
            fig.update_layout(xaxis_title=f"[S] ({c_unit})", yaxis_title=f"v ({v_unit_display})", title=f"Simulation — {model_select}")
            st.plotly_chart(apply_plotly_theme(fig, 550), use_container_width=True)
            
        with tab_3d:
            s_range_3d = np.linspace(0, 100, 30); i_range_3d = np.linspace(0, 100, 30)
            S_MESH, I_MESH = np.meshgrid(s_range_3d, i_range_3d)
            Z_MOLAR = generate_3d_surface(model_select, S_MESH * c_fac, I_MESH * c_fac, params_molar)
            Z_USER = Z_MOLAR / v_fac
            fig_3d = go.Figure(data=[go.Surface(z=Z_USER, x=s_range_3d, y=i_range_3d, colorscale='Viridis', opacity=0.85)])
            fig_3d.update_layout(title=f"Simulation Surface — {model_select}", scene=dict(xaxis_title=f"[S] ({c_unit})", yaxis_title=f"[I] ({c_unit})", zaxis_title=f"v ({v_unit_display})"))
            st.plotly_chart(apply_plotly_theme(fig_3d, 700), use_container_width=True)

# ═══════════════════════════════════════════════════════
#  THERMODYNAMICS MODE
# ═══════════════════════════════════════════════════════
elif mode == "Thermodynamics":
    st.markdown('<p class="hero-title">🔥 Thermodynamic Analysis</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Arrhenius activation energy · Eyring transition-state theory · ΔH‡, ΔS‡, ΔG‡ profiling</p>', unsafe_allow_html=True)
    st.markdown("")

    col1, col2 = st.columns([1, 2])
    with col1:
        st.markdown('<span class="section-badge">📥 Data Entry</span>', unsafe_allow_html=True)
        with st.expander("📂 Load Example"):
            st.caption("Generates Arrhenius-based data with realistic noise. Each click uses new random parameters.")
            if st.button("Load Thermo Example", use_container_width=True):
                df_th, ea_kj = make_example_thermo()
                st.session_state.df_thermo = df_th
                st.session_state.thermo_ea_hint = ea_kj
                reset_editor_key(); st.rerun()
            if "thermo_ea_hint" in st.session_state:
                st.info(f"💡 True Ea ≈ **{st.session_state.thermo_ea_hint:.1f} kJ/mol** — see how close the fit gets!")

        if "df_thermo" not in st.session_state:
            st.session_state.df_thermo = pd.DataFrame({
                "Temp (T)": [25, 30, 35, 40, 45],
                "Rate (k)": [10.5, 15.2, 22.1, 31.5, 42.0]
            })
        temp_unit = st.radio("Temperature Unit", ["Celsius (°C)", "Kelvin (K)"], horizontal=True)
        edited_df = st.data_editor(st.session_state.df_thermo, num_rows="dynamic", use_container_width=True, key=f"editor_thermo_{st.session_state.editor_key}")
        run_thermo = st.button("🚀  Calculate", type="primary", use_container_width=True)

    with col2:
        if run_thermo:
            try:
                raw_t = edited_df.iloc[:, 0].values.astype(float)
                raw_k = edited_df.iloc[:, 1].values.astype(float)
                temps_k = raw_t + 273.15 if "Celsius" in temp_unit else raw_t
                
                arr = solve_arrhenius(temps_k, raw_k)
                eyr = solve_eyring(temps_k, raw_k)
                st.success("✅ Thermodynamic analysis complete!")
                
                report_text = generate_thermo_report(arr, eyr)
                st.download_button("📄  Download Report", report_text, "biokinetic_thermo.txt")
                
                t1, t2 = st.tabs(["🌡️ Arrhenius", "⚗️ Eyring"])
                with t1:
                    ea_kj = arr['Ea']/1000
                    st.markdown(f"""
                    <div style="display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:12px 0 20px 0;">
                        <div class="thermo-card">
                            <div class="tc-label">Activation Energy</div>
                            <div class="tc-value">{ea_kj:.2f}</div>
                            <div class="tc-unit">kJ / mol</div>
                        </div>
                        <div class="thermo-card">
                            <div class="tc-label">Pre-Exponential (A)</div>
                            <div class="tc-value">{arr['A']:.2e}</div>
                            <div class="tc-unit">s⁻¹</div>
                        </div>
                        <div class="thermo-card">
                            <div class="tc-label">Linear Fit R²</div>
                            <div class="tc-value">{arr['r_squared']:.4f}</div>
                            <div class="tc-unit">goodness of fit</div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=arr['reg_x'], y=arr['reg_y'], mode='markers', name="Data", marker=dict(size=10, line=dict(width=1, color='rgba(255,255,255,0.3)'))))
                    x_l = np.linspace(min(arr['reg_x']), max(arr['reg_x']), 10)
                    fig.add_trace(go.Scatter(x=x_l, y=arr['slope']*x_l + arr['intercept'], mode='lines', name="Linear Fit", line=dict(width=2.5)))
                    fig.update_layout(title="Arrhenius Plot", xaxis_title="1/T (K⁻¹)", yaxis_title="ln(k)")
                    st.plotly_chart(apply_plotly_theme(fig, 480), use_container_width=True)
                
                with t2:
                    entropy_note = "Associative — ordered TS" if eyr['dS'] < 0 else "Dissociative — disordered TS"
                    st.markdown(f"""
                    <div style="display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:12px 0 20px 0;">
                        <div class="thermo-card">
                            <div class="tc-label">ΔH‡ Enthalpy</div>
                            <div class="tc-value">{eyr['dH']/1000:.2f}</div>
                            <div class="tc-unit">kJ / mol</div>
                        </div>
                        <div class="thermo-card">
                            <div class="tc-label">ΔS‡ Entropy</div>
                            <div class="tc-value" style="color:{'#fb7185' if eyr['dS']<0 else '#fbbf24'}">{eyr['dS']:.2f}</div>
                            <div class="tc-unit">J / mol·K &mdash; {entropy_note}</div>
                        </div>
                        <div class="thermo-card">
                            <div class="tc-label">ΔG‡ at 25°C</div>
                            <div class="tc-value">{eyr['dG_25']/1000:.2f}</div>
                            <div class="tc-unit">kJ / mol</div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=eyr['reg_x'], y=eyr['reg_y'], mode='markers', name="Data", marker=dict(size=10, line=dict(width=1, color='rgba(255,255,255,0.3)'))))
                    x_l = np.linspace(min(eyr['reg_x']), max(eyr['reg_x']), 10)
                    fig.add_trace(go.Scatter(x=x_l, y=eyr['slope']*x_l + eyr['intercept'], mode='lines', name="Linear Fit", line=dict(width=2.5)))
                    fig.update_layout(title="Eyring Plot", xaxis_title="1/T (K⁻¹)", yaxis_title="ln(k/T)")
                    st.plotly_chart(apply_plotly_theme(fig, 480), use_container_width=True)

            except Exception as e: st.error(f"⚠️ Error: {e}")

# ═══════════════════════════════════════════════════════
#  GUIDE MODE
# ═══════════════════════════════════════════════════════
elif mode == "Guide":
    st.markdown('<p class="hero-title">📖 User Guide</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Everything you need to use BioKineticPy — from zero to publication-ready results</p>', unsafe_allow_html=True)
    st.markdown("")

    # ── Quick Start ──
    st.markdown('<span class="section-badge">🚀 Quick Start</span>', unsafe_allow_html=True)
    c1, c2, c3 = st.columns(3)
    c1.markdown("""<div class="param-card" style="border-left:3px solid #38bdf8;">
<div class="pc-name">Step 1 — Enter Data</div>
<p style="color:var(--text-secondary);font-size:0.87rem;margin-top:8px;line-height:1.6;">
Go to <strong>🔬 Analysis</strong>. Type your <strong>[S]</strong> and <strong>v</strong> values into
the table, or click <em>Load Example Data</em> to try a built-in dataset instantly.
</p></div>""", unsafe_allow_html=True)
    c2.markdown("""<div class="param-card" style="border-left:3px solid #a78bfa;">
<div class="pc-name">Step 2 — Run Analysis</div>
<p style="color:var(--text-secondary);font-size:0.87rem;margin-top:8px;line-height:1.6;">
Click <strong>🚀 Run Analysis</strong>. All 7 kinetic models are tested automatically.
The best model is selected by <strong>AICc</strong> and shown instantly with fitted parameters.
</p></div>""", unsafe_allow_html=True)
    c3.markdown("""<div class="param-card" style="border-left:3px solid #34d399;">
<div class="pc-name">Step 3 — Download</div>
<p style="color:var(--text-secondary);font-size:0.87rem;margin-top:8px;line-height:1.6;">
Review the charts and diagnostics, then click <strong>📄 Download Full Report</strong> for a
publication-ready Methods text with all parameters and statistical tests.
</p></div>""", unsafe_allow_html=True)

    st.markdown("---")

    # ── Mode Explanations ──
    st.markdown('<span class="section-badge">🔬 The Three Modes</span>', unsafe_allow_html=True)
    tab_am, tab_sm, tab_tm = st.tabs(["🔬 Analysis", "📈 Simulation", "🌡️ Thermodynamics"])

    with tab_am:
        cola, colb = st.columns(2)
        cola.markdown("""
**What you need:**
- Substrate concentrations **[S]**
- Measured velocities **v**
- Optionally inhibitor concentrations **[I]**

**Data layouts:**
- **Matrix** — each column is a different [I]; good for side-by-side inhibitor experiments  
- **Standard** — each row is one (S, v, I) triple; good for long-form spreadsheet data

**CSV / Excel upload:** include a header row. App accepts `.csv` and `.xlsx`.
""")
        colb.markdown("""
**Output tabs:**

| Tab | Contents |
|-----|----------|
| 📈 Fit | Best-fit curve over your data |
| 🩺 Diagnostics | Normality & randomness chips |
| 📐 Linearizations | LB, HW, Eadie-Hofstee plots |
| 🧪 Inhibition | Dixon & Cornish-Bowden plots |
| 🌐 3D Landscape | Velocity surface over [S] and [I] |

**Confidence Intervals:** click *Compute 95% CI* to run a 50-iteration Monte Carlo bootstrap for parameter uncertainty.
""")

    with tab_sm:
        cola, colb = st.columns(2)
        cola.markdown("""
**How to use:**
1. Choose a **mechanism** from the dropdown
2. Drag the sliders to set Vmax, Km, Ki, etc.
3. The 2D/3D plots update in real time
4. The **rate equation** card updates too — useful for teaching

**When useful?**
- Teaching enzyme kinetics concepts  
- Designing [S] ranges before an experiment  
- Verifying whether estimated parameters are biologically reasonable
""")
        colb.markdown("""
**Slider Parameters:**

| Symbol | Meaning |
|--------|---------|
| Vmax | Rate at saturating [S] |
| Km | [S] at half-Vmax |
| Ki | [I] that halves activity |
| Khalf | Like Km for allosteric enzymes |
| n | Hill coefficient (>1 = cooperativity) |
| Ksi | [S] at which substrate self-inhibits |
| α | Mixed inhibition preference |
""")

    with tab_tm:
        cola, colb = st.columns(2)
        cola.markdown("""
**What you need:**
- Temperature values (°C or K)  
- Rate constants **k** at each temperature  
- Minimum 4 points spanning a meaningful range (e.g. 20–50 °C)

**Arrhenius output:**
- **Ea** — activation energy (kJ/mol): energy barrier to catalysis  
- **A** — pre-exponential: theoretical collision frequency  

**Eyring output:**
- **ΔH‡** — enthalpy of activation  
- **ΔS‡** — entropy of activation (negative = ordered TS; positive = disordered)  
- **ΔG‡** — Gibbs energy barrier at 25 °C  
""")
        colb.markdown("""
**Interpreting the plots:**

*Arrhenius Plot:* ln(k) vs 1/T(K⁻¹)  
→ slope = **−Ea/R**. Straight line confirms the model holds.

*Eyring Plot:* ln(k/T) vs 1/T(K⁻¹)  
→ slope = **−ΔH‡/R**. Deviations suggest mechanism changes.

> **Tip:** a curvature in either plot at high temperature often signals  
> protein unfolding — not a simple thermodynamic transition.
""")

    st.markdown("---")

    # ── Kinetic Models ──
    st.markdown('<span class="section-badge">⚗️ Kinetic Models Explained</span>', unsafe_allow_html=True)
    st.markdown("Expand any model below for its equation, mechanism, diagnostic signature, and a real-world example.")

    models_info = [
        ("Michaelis-Menten", "#38bdf8",
         "v = Vmax·[S] / (Km + [S])",
         "The foundational model for a single-substrate enzyme with no inhibitor. Produces a **hyperbolic** curve that plateaus at Vmax. Assumes: rapid-equilibrium or steady-state ES complex, single substrate, no product inhibition.",
         "Vmax **unchanged**, Km **unchanged**. This is the baseline.",
         "Most simple one-substrate enzymes under standard conditions."),
        ("Competitive", "#34d399",
         "v = Vmax·[S] / (Km·(1+[I]/Ki) + [S])",
         "The inhibitor **competes** with substrate for the active site. Because binding is reversible, adding more substrate displaces the inhibitor.",
         "Vmax **unchanged** · Km **increases**. LB lines converge on the **y-axis**.",
         "Methotrexate vs. dihydrofolate reductase; statins vs. HMG-CoA reductase."),
        ("Non-Competitive", "#a78bfa",
         "v = (Vmax/(1+[I]/Ki))·[S] / (Km + [S])",
         "Inhibitor binds an allosteric site on both free enzyme and ES complex with equal affinity. Does not prevent substrate binding but blocks catalysis.",
         "Vmax **decreases** · Km **unchanged**. LB lines converge on the **x-axis**.",
         "Heavy metal inhibition; iodoacetamide on cysteine-protease active sites."),
        ("Uncompetitive", "#fbbf24",
         "v = Vmax·[S] / (Km + [S]·(1+[I]/Ki))",
         "Inhibitor binds **only** to the ES complex — it cannot bind free enzyme. Paradoxically, both Vmax and Km decrease by the same factor.",
         "Both Vmax and Km **decrease** equally. LB lines are **parallel** (same slope).",
         "Phenylalanine inhibition of alkaline phosphatase; Li⁺ on inositol monophosphatase."),
        ("Mixed", "#fb7185",
         "v = Vmax'·[S] / (Km' + [S])  [α controls E vs ES preference]",
         "General case: inhibitor binds both free enzyme and ES complex with **different** affinities (α). α > 1 → prefers free enzyme (competitive-like). α < 1 → prefers ES (uncompetitive-like). α = 1 → pure non-competitive.",
         "Both Vmax and Km change. LB lines intersect at a **non-axis** point.",
         "Many real-world inhibitors; the α parameter distinguishes the binding preference."),
        ("Hill (Allosteric)", "#f472b6",
         "v = Vmax·[S]ⁿ / (Khalfⁿ + [S]ⁿ)",
         "Describes enzymes with multiple cooperative binding sites. n>1: positive cooperativity (sigmoidal S-curve). n<1: negative cooperativity. n=1: reduces to Michaelis-Menten.",
         "Sigmoidal curve shape. n is the **cooperativity index**.",
         "Haemoglobin (n≈2.8), aspartate transcarbamoylase (ATCase), phosphofructokinase."),
        ("Substrate Inhibition", "#60a5fa",
         "v = Vmax·[S] / (Km + [S] + [S]²/Ksi)",
         "At high [S], a second substrate molecule binds the ES complex forming an inactive ESS complex. Velocity rises to a peak then falls — a bell-shaped curve.",
         "**Bell-shaped curve** with a velocity maximum followed by decline at high [S].",
         "Alcohol dehydrogenase at high ethanol; phenol hydroxylase; some CYP enzymes."),
    ]

    for name, color, eq, theory, signature, example in models_info:
        with st.expander(f"**{name}**"):
            mc1, mc2 = st.columns([1.1, 1])
            mc1.markdown(f"""<div class="eq-display" style="border-left-color:{color};">
<div class="eq-label">Rate Equation</div>{eq}</div>

**Mechanism:** {theory}""", unsafe_allow_html=True)
            mc2.markdown(f"""**Diagnostic signature:**  
{signature}

**Real example:**  
*{example}*""")

    st.markdown("---")

    # ── Interpreting Results ──
    st.markdown('<span class="section-badge">📊 Interpreting Results</span>', unsafe_allow_html=True)
    ir1, ir2 = st.columns(2)
    ir1.markdown("""
**AICc — Model Selection**

The AICc (corrected Akaike Information Criterion) rewards fit quality and penalises complexity.
The model with the **lowest AICc** wins.

| ΔAIC vs best | Interpretation |
|-------------|----------------|
| < 2 | Models virtually indistinguishable |
| 4 – 7 | Modest evidence for the winner |
| > 10 | Strong evidence — clear winner |

An "ambiguity" warning appears when ΔAIC < 2 between the top models.
""")
    ir2.markdown("""**R² — Goodness of Fit**""")
    ir2.markdown("""<div style="display:flex;flex-direction:column;gap:8px;margin-top:4px;">
  <div style="display:flex;align-items:center;gap:10px;"><span class="chip chip-pass">R² > 0.95</span>
    <span style="color:var(--text-secondary);font-size:0.87rem;">Excellent — model captures the data well</span></div>
  <div style="display:flex;align-items:center;gap:10px;"><span class="chip chip-warn">0.80 – 0.95</span>
    <span style="color:var(--text-secondary);font-size:0.87rem;">Acceptable — inspect residuals for patterns</span></div>
  <div style="display:flex;align-items:center;gap:10px;"><span class="chip chip-fail">R² < 0.80</span>
    <span style="color:var(--text-secondary);font-size:0.87rem;">Poor — try another model or remove outliers</span></div>
</div>""", unsafe_allow_html=True)

    st.markdown("")
    dg1, dg2 = st.columns(2)
    dg1.markdown("""
**Shapiro-Wilk Test (Normality)**  
Tests whether residuals follow a Gaussian distribution — the assumption behind least-squares fitting.

- ✅ **p > 0.05** — residuals are normal; model is appropriate  
- ⚠️ **p ≤ 0.05** — non-normal residuals; possible outliers, wrong model, or unequal variance → try weighted fitting (1/v or 1/v²)
""")
    dg2.markdown("""
**Runs Test (Randomness)**  
Tests whether residuals are randomly scattered above/below zero, or show a systematic U-shape or wave.

- ✅ **|Z| < 1.96** — random residuals; model fits the trend  
- ⚠️ **|Z| ≥ 1.96** — systematic deviation; the model may be **mechanistically incorrect** — consider switching to a different model
""")

    st.markdown("---")

    # ── Data Format ──
    st.markdown('<span class="section-badge">📋 Data Format Guide</span>', unsafe_allow_html=True)
    df1, df2 = st.columns(2)
    df1.markdown("""**Matrix Format** — multiple inhibitor columns side by side

| Include? | Substrate [μM] | v (I=0) | v (I=5) | v (I=10) |
|----------|---------------|---------|---------|----------|
| ✅ | 1 | 10.5 | 4.8 | 2.5 |
| ✅ | 5 | 35.1 | 19.5 | 12.1 |
| ✅ | 20 | 70.2 | 45.0 | 30.5 |

Set the [I] values using the **[I] A / B / C** inputs above the table. Leave a column blank if not measured.
""")
    df2.markdown("""**Standard Format** — one row per observation

| Include? | S [μM] | v | I |
|----------|--------|---|---|
| ✅ | 1 | 10.5 | 0 |
| ✅ | 5 | 35.1 | 0 |
| ✅ | 1 | 4.8 | 5 |
| ✅ | 5 | 19.5 | 5 |

Set **I = 0** for no inhibitor. Mix any combination of [S] and [I] rows. Use the **Include?** checkbox to exclude outliers without deleting rows.
""")

    st.markdown("---")

    # ── Glossary ──
    st.markdown('<span class="section-badge">📚 Glossary</span>', unsafe_allow_html=True)
    with st.expander("📚 Scientific Glossary (click to expand)"):
        glossary = [
            ("Vmax", "Maximum reaction velocity at saturating substrate. Units: concentration/time (e.g. μM/min)."),
            ("Km", "Michaelis constant — [S] at half-Vmax. Inverse measure of substrate affinity; lower Km = tighter binding."),
            ("Ki", "Inhibition constant — [I] producing 50% reduction in apparent affinity or Vmax. Lower = more potent inhibitor."),
            ("kcat", "Turnover number — substrate molecules converted per active site per second. Requires [E]total to be defined."),
            ("kcat/Km", "Catalytic efficiency — second-order rate constant for E+S→P. 'Gold standard' of enzyme performance. Theoretical max ~10⁹ M⁻¹s⁻¹."),
            ("AICc", "Akaike Information Criterion (small-sample corrected). Lower = better. Penalises extra parameters to prevent overfitting."),
            ("R²", "Coefficient of determination. Fraction of velocity variance explained by the model. 1.0 = perfect fit."),
            ("Residual", "Difference between observed and predicted velocity: residual = v_obs − v_pred. Should be random and normally distributed."),
            ("Levenberg-Marquardt", "The non-linear least-squares optimisation algorithm used for fitting. Combines gradient descent and Newton's method."),
            ("Bootstrap CI", "Uncertainty estimation by fitting the model to many simulated noisy datasets. Gives 95% confidence bounds for each parameter."),
            ("Lineweaver-Burk", "Double-reciprocal plot: 1/v vs 1/[S]. Historically used to estimate Vmax and Km; sensitive to error at low [S]."),
            ("Eadie-Hofstee", "Plot of v vs v/[S]. More balanced error distribution. Good for distinguishing inhibition types visually."),
            ("Hanes-Woolf", "Plot of [S]/v vs [S]. Least distorted by measurement error; best for Km estimation at high [S]."),
            ("Dixon Plot", "1/v vs [I] for each [S]. Intersection gives Ki and reveals inhibition type."),
            ("Cornish-Bowden", "[S]/v vs [I]. More robust than Dixon for identifying inhibition mechanism."),
            ("Arrhenius Ea", "Activation energy from Arrhenius equation: slope of ln(k) vs 1/T × (−R). Higher Ea = more temperature-sensitive reaction."),
            ("ΔH‡", "Enthalpy of activation (Eyring). Bond-breaking/forming energy at the transition state."),
            ("ΔS‡", "Entropy of activation. Negative = ordered transition state (associative). Positive = disordered (dissociative)."),
            ("ΔG‡", "Gibbs free energy of activation = ΔH‡ − TΔS‡. The total kinetic barrier at a given temperature."),
        ]
        for term, defn in glossary:
            st.markdown(f"""<div style="display:flex;gap:14px;padding:9px 0;border-bottom:1px solid var(--border);">
<div style="min-width:155px;font-family:'JetBrains Mono',monospace;font-size:0.83rem;
color:var(--accent-cyan);font-weight:500;padding-top:2px;">{term}</div>
<div style="color:var(--text-secondary);font-size:0.87rem;line-height:1.6;">{defn}</div>
</div>""", unsafe_allow_html=True)

# ── Footer ──
st.markdown("---")
st.markdown("""
<div style="text-align:center; padding: 16px 0 8px 0;">
    <span class="version-badge" style="margin-right:8px;">🧬 BioKineticPy v1.1</span>
    <span style="color:var(--text-muted);font-size:0.8rem;">Open Source Enzyme Kinetics &amp; Thermodynamics Suite</span>
    <br><br>
    <span style="color:var(--text-muted);font-size:0.75rem;">
        Built by <strong style="color:var(--text-secondary);">Chibuike Praise Okechukwu</strong>
        &nbsp;·&nbsp; Levenberg-Marquardt regression
        &nbsp;·&nbsp; AICc model selection
        &nbsp;·&nbsp; Bootstrap CI
    </span>
</div>
""", unsafe_allow_html=True)