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

# ═══════════════════════════════════════════════════════
#  PAGE CONFIG & PREMIUM THEME
# ═══════════════════════════════════════════════════════
st.set_page_config(
    page_title="BioKineticPy v1.0 — Enzyme Kinetics Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- PREMIUM CSS THEME ---
st.markdown("""
<style>
/* ── Google Font ── */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');

/* ── Root Variables ── */
:root {
    --bg-primary: #0d1117;
    --bg-secondary: #161b22;
    --bg-card: #1c2333;
    --bg-hover: #242d3d;
    --accent-cyan: #58a6ff;
    --accent-emerald: #3fb950;
    --accent-purple: #bc8cff;
    --accent-amber: #d29922;
    --accent-rose: #f85149;
    --text-primary: #e6edf3;
    --text-secondary: #8b949e;
    --text-muted: #6e7681;
    --border: #30363d;
    --gradient-hero: linear-gradient(135deg, #58a6ff 0%, #bc8cff 50%, #f778ba 100%);
    --glass-bg: rgba(22, 27, 34, 0.75);
    --glass-border: rgba(88, 166, 255, 0.15);
}

/* ── Global ── */
html, body, [data-testid="stAppViewContainer"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
    color: var(--text-primary);
}

[data-testid="stAppViewContainer"] {
    background: var(--bg-primary);
}

[data-testid="stHeader"] {
    background: transparent !important;
}

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: var(--bg-secondary) !important;
    border-right: 1px solid var(--border) !important;
}
[data-testid="stSidebar"] [data-testid="stMarkdown"] p,
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] .stRadio label span {
    color: var(--text-primary) !important;
}
[data-testid="stSidebar"] hr {
    border-color: var(--border) !important;
}

/* ── Cards / Containers ── */
[data-testid="stExpander"],
[data-testid="stForm"] {
    background: var(--bg-card) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
}

div[data-testid="stMetric"] {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 16px 20px;
    transition: transform 0.2s ease, box-shadow 0.2s ease;
}
div[data-testid="stMetric"]:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 24px rgba(88,166,255,0.08);
}
div[data-testid="stMetric"] label {
    color: var(--text-secondary) !important;
    font-weight: 500 !important;
    font-size: 0.78rem !important;
    text-transform: uppercase;
    letter-spacing: 0.06em;
}
div[data-testid="stMetric"] [data-testid="stMetricValue"] {
    color: var(--accent-cyan) !important;
    font-weight: 700 !important;
}

/* ── Tabs ── */
.stTabs [data-baseweb="tab-list"] {
    background: var(--bg-card);
    border-radius: 12px;
    padding: 4px;
    gap: 4px;
    border: 1px solid var(--border);
}
.stTabs [data-baseweb="tab"] {
    border-radius: 8px !important;
    color: var(--text-secondary) !important;
    font-weight: 500 !important;
    transition: all 0.2s ease !important;
    padding: 8px 16px !important;
}
.stTabs [aria-selected="true"] {
    background: rgba(88,166,255,0.12) !important;
    color: var(--accent-cyan) !important;
}
.stTabs [data-baseweb="tab-highlight"] {
    background-color: var(--accent-cyan) !important;
}

/* ── Buttons ── */
.stButton > button {
    background: linear-gradient(135deg, #58a6ff 0%, #388bfd 100%) !important;
    color: #fff !important;
    border: none !important;
    border-radius: 10px !important;
    font-weight: 600 !important;
    padding: 0.55rem 1.6rem !important;
    letter-spacing: 0.02em;
    transition: all 0.25s ease !important;
    box-shadow: 0 4px 14px rgba(88,166,255,0.25) !important;
}
.stButton > button:hover {
    transform: translateY(-1px) !important;
    box-shadow: 0 6px 20px rgba(88,166,255,0.35) !important;
    filter: brightness(1.08) !important;
}

/* ── Download Button ── */
.stDownloadButton > button {
    background: var(--bg-card) !important;
    color: var(--accent-cyan) !important;
    border: 1px solid var(--border) !important;
    border-radius: 10px !important;
    font-weight: 500 !important;
    transition: all 0.25s ease !important;
}
.stDownloadButton > button:hover {
    border-color: var(--accent-cyan) !important;
    background: rgba(88,166,255,0.08) !important;
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
}

/* ── Data Editor / Tables ── */
[data-testid="stDataFrame"],
.stDataFrame {
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    overflow: hidden;
}

/* ── Alerts ── */
.stAlert {
    border-radius: 10px !important;
    border: 1px solid var(--border) !important;
}

/* ── Plotly overrides ── */
.js-plotly-plot .plotly .main-svg {
    border-radius: 12px;
}

/* ── Hero header ── */
.hero-title {
    background: var(--gradient-hero);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    font-weight: 800;
    font-size: 2rem;
    line-height: 1.15;
    margin-bottom: 0;
}
.hero-subtitle {
    color: var(--text-secondary);
    font-size: 0.95rem;
    font-weight: 400;
    margin-top: 4px;
}
.section-badge {
    display: inline-block;
    background: rgba(88,166,255,0.1);
    color: var(--accent-cyan);
    border: 1px solid rgba(88,166,255,0.2);
    border-radius: 20px;
    padding: 4px 14px;
    font-size: 0.75rem;
    font-weight: 600;
    letter-spacing: 0.05em;
    text-transform: uppercase;
    margin-bottom: 8px;
}

/* ── Dividers ── */
hr {
    border-color: var(--border) !important;
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 6px; }
::-webkit-scrollbar-track { background: var(--bg-primary); }
::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: var(--text-muted); }
</style>
""", unsafe_allow_html=True)

# ── Plotly Dark Theme Template ──
PLOTLY_DARK = dict(
    paper_bgcolor="rgba(28,35,51,0)",
    plot_bgcolor="rgba(28,35,51,0.6)",
    font=dict(family="Inter, sans-serif", color="#e6edf3", size=12),
    xaxis=dict(gridcolor="rgba(48,54,61,0.6)", zerolinecolor="#30363d"),
    yaxis=dict(gridcolor="rgba(48,54,61,0.6)", zerolinecolor="#30363d"),
    colorway=["#58a6ff", "#3fb950", "#bc8cff", "#d29922", "#f85149", "#f778ba", "#79c0ff", "#56d364"],
    margin=dict(l=50, r=20, t=50, b=50),
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
    st.markdown('<p class="hero-title" style="font-size:1.6rem;">🧬 BioKineticPy</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Professional Enzyme Kinetics Suite</p>', unsafe_allow_html=True)
    st.markdown("---")
    
    # PROJECT SAVE/LOAD
    st.markdown('<span class="section-badge">💾 Project</span>', unsafe_allow_html=True)
    current_state = {
        "df_matrix": st.session_state.get("df_matrix", pd.DataFrame()).to_dict() if "df_matrix" in st.session_state else None,
        "df_standard": st.session_state.get("df_standard", pd.DataFrame()).to_dict() if "df_standard" in st.session_state else None,
        "df_thermo": st.session_state.get("df_thermo", pd.DataFrame()).to_dict() if "df_thermo" in st.session_state else None
    }
    st.download_button("💾  Save Project", json.dumps(current_state), "biokinetic_project.json", "application/json")
    uploaded_proj = st.file_uploader("Load Project", type=["json"], label_visibility="collapsed")
    if uploaded_proj:
        try:
            data = json.load(uploaded_proj)
            if data["df_matrix"]: st.session_state.df_matrix = pd.DataFrame.from_dict(data["df_matrix"])
            if data["df_standard"]: st.session_state.df_standard = pd.DataFrame.from_dict(data["df_standard"])
            if data.get("df_thermo"): st.session_state.df_thermo = pd.DataFrame.from_dict(data["df_thermo"])
            reset_editor_key()
            st.success("✅ Project loaded successfully!")
        except Exception:
            st.error("Failed to load project file.")

    st.markdown("---")
    mode = st.radio("🔬 Operation Mode", ["Analysis (Fit Data)", "Predictive (Simulation)", "Thermodynamics"])
    
    if mode != "Thermodynamics":
        st.markdown("---")
        st.markdown('<span class="section-badge">⚙️ Settings</span>', unsafe_allow_html=True)
        weighting_mode = st.selectbox("Weighting Method", ["None (Homoscedastic)", "1/v (Poisson)", "1/v² (Relative Error)"])
        st.markdown('<span class="section-badge">📏 Units</span>', unsafe_allow_html=True)
        c_unit = st.selectbox("Concentration Unit", ["mM", "μM", "nM", "M", "pM"], index=1)
        t_unit = st.selectbox("Time Base", ["min", "sec"], index=0)
        v_unit_display = f"{c_unit}/{t_unit}"
        st.info(f"Velocity → **{v_unit_display}**")
        with st.expander("🧫 Enzyme Properties"):
            use_et = st.checkbox("Define [E]total")
            et_val = st.number_input("Enzyme Conc", 1.0) if use_et else 0
            et_unit = st.selectbox("Enzyme Unit", ["nM", "μM", "mM", "M"], 0) if use_et else "nM"

# ═══════════════════════════════════════════════════════
#  ANALYSIS MODE
# ═══════════════════════════════════════════════════════
if mode == "Analysis (Fit Data)":
    st.markdown('<p class="hero-title">🔬 Kinetic Analysis</p>', unsafe_allow_html=True)
    st.markdown('<p class="hero-subtitle">Automated model selection & statistical validation</p>', unsafe_allow_html=True)
    st.markdown("")

    col1, col2 = st.columns([1.2, 2])
    
    with col1:
        st.markdown('<span class="section-badge">📥 Data Input</span>', unsafe_allow_html=True)
        input_format = st.radio("Layout", ["Matrix", "Standard"], horizontal=True)
        uploaded_file = st.file_uploader("Upload Data (CSV / Excel)", type=["csv", "xlsx"])

        # DEMO LOADERS
        with st.expander("📂 Load Example Data"):
            c_ex1, c_ex2 = st.columns(2)
            if c_ex1.button("Competitive", use_container_width=True):
                st.session_state.df_matrix = pd.DataFrame({
                    "Include?": [True]*5,
                    f"Substrate [{c_unit}]": [1, 5, 10, 20, 50],
                    "v (I=0)": [10.5, 35.1, 55.0, 70.2, 90.5],
                    "v (I=5)": [4.8, 19.5, 32.9, 45.0, 60.1],
                    "v (I=10)": [2.5, 12.1, 21.0, 30.5, 42.0], 
                    "v (I=C)": [None]*5
                })
                reset_editor_key()
                st.rerun()
            if c_ex2.button("Substrate Inhib", use_container_width=True):
                st.session_state.df_standard = pd.DataFrame({
                    "Include?": [True]*6,
                    f"S [{c_unit}]": [1, 5, 10, 50, 100, 200],
                    "v": [10, 45, 80, 60, 40, 25],
                    "I": [0]*6
                })
                reset_editor_key()
                st.rerun()

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
                    st.success(f"✅ Best Fit: **{best_model['model']}**  ·  AIC: {best_model['aic']:.1f}  ·  R²: {best_model['r_squared']:.4f}")
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
                    c1, c2, c3 = st.columns(3)
                    sw_p = stat_tests.get('shapiro_p', 0); rz = stat_tests.get('runs_z', 0)
                    c1.metric("Normality (SW p)", f"{sw_p:.3f}", delta="Normal" if sw_p > 0.05 else "Non-Normal", delta_color="normal")
                    c2.metric("Randomness (Runs Z)", f"{rz:.2f}", delta="Random" if abs(rz) < 1.96 else "Systematic", delta_color="inverse")
                    c3.metric("R²", f"{best_model['r_squared']:.4f}", delta="Good" if best_model['r_squared'] > 0.95 else "Low", delta_color="normal")
                    
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
                    cols = st.columns(2)
                    # Lineweaver-Burk
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
                    fig_lb.update_layout(title="Lineweaver-Burk", xaxis_title="1/[S]", yaxis_title="1/v")
                    cols[0].plotly_chart(apply_plotly_theme(fig_lb, 400), use_container_width=True)
                    
                    # Hanes-Woolf
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
                    fig_hw.update_layout(title="Hanes-Woolf", xaxis_title="[S]", yaxis_title="[S]/v")
                    cols[1].plotly_chart(apply_plotly_theme(fig_hw, 400), use_container_width=True)
                    
                    # Eadie-Hofstee
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
                    fig_eh.update_layout(title="Eadie-Hofstee", xaxis_title="v/[S]", yaxis_title="v")
                    st.plotly_chart(apply_plotly_theme(fig_eh, 400), use_container_width=True)

                with tab4:
                    cols = st.columns(2)
                    # Dixon
                    fig_dixon = go.Figure()
                    for s_val in sorted(list(set([round(s, 6) for s in user_s]))):
                        idx = [j for j, s in enumerate(user_s) if abs(s - s_val) < 1e-9 and user_v[j] > 0]
                        if len(idx) > 1:
                            fig_dixon.add_trace(go.Scatter(x=user_i[idx], y=1/user_v[idx], mode='markers+lines', name=f"[S]={s_val:.1f}"))
                    fig_dixon.update_layout(title="Dixon Plot", xaxis_title="[I]", yaxis_title="1/v")
                    cols[0].plotly_chart(apply_plotly_theme(fig_dixon, 400), use_container_width=True)
                    
                    # Cornish-Bowden
                    fig_cb = go.Figure()
                    for s_val in sorted(list(set([round(s, 6) for s in user_s]))):
                        idx = [j for j, s in enumerate(user_s) if abs(s - s_val) < 1e-9 and user_v[j] > 0]
                        if len(idx) > 1:
                            fig_cb.add_trace(go.Scatter(x=user_i[idx], y=s_val/user_v[idx], mode='markers+lines', name=f"[S]={s_val:.1f}"))
                    fig_cb.update_layout(title="Cornish-Bowden", xaxis_title="[I]", yaxis_title="[S]/v")
                    cols[1].plotly_chart(apply_plotly_theme(fig_cb, 400), use_container_width=True)

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
    st.markdown('<p class="hero-subtitle">Explore enzyme mechanisms with interactive parameter tuning</p>', unsafe_allow_html=True)
    st.markdown("")

    col_input, col_graph = st.columns([1, 3])
    with col_input:
        st.markdown('<span class="section-badge">🔧 Parameters</span>', unsafe_allow_html=True)
        model_select = st.selectbox("Mechanism", list(MODEL_REGISTRY.keys()))
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
    st.markdown('<p class="hero-subtitle">Arrhenius & Eyring transition-state profiling</p>', unsafe_allow_html=True)
    st.markdown("")

    col1, col2 = st.columns([1, 2])
    with col1:
        st.markdown('<span class="section-badge">📥 Data Entry</span>', unsafe_allow_html=True)
        with st.expander("📂 Load Example"):
            if st.button("Load Thermo Data", use_container_width=True):
                st.session_state.df_thermo = pd.DataFrame({
                    "Temp (T)": [25, 30, 37, 42, 50],
                    "Rate (k)": [12.5, 18.2, 35.1, 52.5, 90.0]
                })
                reset_editor_key()
                st.rerun()

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
                    c = st.columns(3)
                    c[0].metric("Ea", f"{arr['Ea']/1000:.2f} kJ/mol")
                    c[1].metric("A (Pre-exp)", f"{arr['A']:.2e}")
                    c[2].metric("R²", f"{arr['r_squared']:.4f}")
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=arr['reg_x'], y=arr['reg_y'], mode='markers', name="Data", marker=dict(size=10, line=dict(width=1, color='rgba(255,255,255,0.3)'))))
                    x_l = np.linspace(min(arr['reg_x']), max(arr['reg_x']), 10)
                    fig.add_trace(go.Scatter(x=x_l, y=arr['slope']*x_l + arr['intercept'], mode='lines', name="Linear Fit", line=dict(width=2.5)))
                    fig.update_layout(title="Arrhenius Plot", xaxis_title="1/T (K⁻¹)", yaxis_title="ln(k)")
                    st.plotly_chart(apply_plotly_theme(fig, 480), use_container_width=True)
                
                with t2:
                    c = st.columns(3)
                    c[0].metric("ΔH‡", f"{eyr['dH']/1000:.2f} kJ/mol")
                    c[1].metric("ΔS‡", f"{eyr['dS']:.2f} J/mol·K")
                    c[2].metric("ΔG‡ (25°C)", f"{eyr['dG_25']/1000:.2f} kJ/mol")
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=eyr['reg_x'], y=eyr['reg_y'], mode='markers', name="Data", marker=dict(size=10, line=dict(width=1, color='rgba(255,255,255,0.3)'))))
                    x_l = np.linspace(min(eyr['reg_x']), max(eyr['reg_x']), 10)
                    fig.add_trace(go.Scatter(x=x_l, y=eyr['slope']*x_l + eyr['intercept'], mode='lines', name="Linear Fit", line=dict(width=2.5)))
                    fig.update_layout(title="Eyring Plot", xaxis_title="1/T (K⁻¹)", yaxis_title="ln(k/T)")
                    st.plotly_chart(apply_plotly_theme(fig, 480), use_container_width=True)

            except Exception as e: st.error(f"⚠️ Error: {e}")

# ── Footer ──
st.markdown("---")
st.markdown(
    '<p style="text-align:center; color:var(--text-muted); font-size:0.8rem;">'
    '🧬 BioKineticPy v1.0 — Open Source Enzyme Kinetics &amp; Thermodynamics Suite '
    '· Built by <strong>Chibuike Praise Okechukwu</strong></p>',
    unsafe_allow_html=True
)