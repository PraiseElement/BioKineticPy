# BioKineticPy v1.1 🧬

### _A Professional Computational Workbench for Enzyme Kinetics & Thermodynamics_

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Framework](https://img.shields.io/badge/built%20with-Streamlit-red.svg)
![Version](https://img.shields.io/badge/version-1.1-brightgreen.svg)

---

## 📖 Project Overview

**BioKineticPy** is a specialized data analysis suite designed for biochemists, enzymologists, and biophysicists. Unlike general-purpose graphing software, BioKineticPy is built specifically for the non-linear regression challenges inherent in enzyme kinetics.

It moves beyond simple Michaelis-Menten fitting to handle complex scenarios like **mixed inhibition**, **allostery (Hill equation)**, and **substrate inhibition (Haldane kinetics)**. Furthermore, it bridges the gap between kinetics and thermodynamics, offering automated **Arrhenius** and **Eyring** analyses to determine activation energy ($E_a$), enthalpy ($\Delta H^\ddagger$), and entropy ($\Delta S^\ddagger$).

## ✨ Key Features

### 🔬 Advanced Kinetic Analysis

- **Automated Model Selection:** Fits data against 7 mechanistic models simultaneously, ranked by **AICc** (corrected Akaike Information Criterion).
- **Global Matrix Fitting:** Analyse full inhibition datasets (multiple [I] concentrations) in a single global fit.
- **Statistical Rigor:**
  - **Bootstrap Resampling (Monte Carlo):** 95% Confidence Intervals for $V_{max}$, $K_m$, and $K_i$.
  - **Residual Diagnostics:** Shapiro-Wilk (normality) and Runs Test (randomness) with pass/warn/fail chip badges.
  - **Weighted Regression:** Supports $1/v$ and $1/v^2$ weighting for heteroscedastic error.

### 🔥 Thermodynamic Profiling

- **Arrhenius Analysis:** Activation energy $E_a$ and pre-exponential factor $A$.
- **Eyring / Transition-State Theory:** $\Delta H^\ddagger$, $\Delta S^\ddagger$, $\Delta G^\ddagger$ with colour-coded entropy interpretation.
- **Unit Intelligence:** Auto-converts between Celsius/Kelvin and handles energy units.

### 📊 Visualization Suite

- **Premium Plotly Dark Theme:** Consistent styling using Space Grotesk and JetBrains Mono fonts, glass hover labels, and teal zero lines.
- **3D Kinetic Landscapes:** Interactive surface plots of velocity over [S] and [I] space.
- **Diagnostic Linearizations:** Lineweaver-Burk, Hanes-Woolf, Eadie-Hofstee, Dixon, and Cornish-Bowden plots.

### 📖 In-App User Guide _(v1.1 new)_

Built-in Guide mode accessible from the sidebar, covering:

- 3-step Quick Start for new users
- All three modes explained in plain language
- All 7 kinetic models with equations, mechanisms, and real examples
- Result interpretation (AICc, R², diagnostic tests)
- Data format reference (Matrix vs Standard)
- 19-term Scientific Glossary

---

## 🚀 Installation & Usage

### Prerequisites

- Python 3.8 or higher
- pip (Python Package Installer)

### Setup Guide

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/PraiseElement/BioKineticPy.git
   cd BioKineticPy
   ```

2. **Create a Virtual Environment (Recommended):**

   ```bash
   python -m venv venv
   # Windows:
   venv\Scripts\activate
   # Mac/Linux:
   source venv/bin/activate
   ```

3. **Install Dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

4. **Launch the Application:**
   ```bash
   streamlit run main.py
   ```

---

## 📚 Theoretical Models Supported

1. **Michaelis-Menten** — `v = Vmax·[S] / (Km + [S])`
2. **Competitive Inhibition** — `v = Vmax·[S] / (Km·(1+[I]/Ki) + [S])`
3. **Non-Competitive Inhibition** — `v = (Vmax/(1+[I]/Ki))·[S] / (Km + [S])`
4. **Uncompetitive Inhibition** — `v = Vmax·[S] / (Km + [S]·(1+[I]/Ki))`
5. **Mixed Inhibition** — general case with α parameter controlling E vs ES preference
6. **Hill Equation** — `v = Vmax·[S]ⁿ / (Khalfⁿ + [S]ⁿ)` — cooperative/allosteric
7. **Substrate Inhibition (Haldane)** — `v = Vmax·[S] / (Km + [S] + [S]²/Ksi)`

---

## 🛠 Changelog

### v1.1 — UI/UX Overhaul & In-App Guide _(2026-02-19)_

#### 🎨 Premium Design System

| #   | Enhancement                                                                                                                   |
| --- | ----------------------------------------------------------------------------------------------------------------------------- |
| 1   | **New HSL-based deep dark palette** — richer backgrounds (`hsl(222,28%,5%)`) with layered glassmorphism cards                 |
| 2   | **Typography upgrade** — Space Grotesk (headings), JetBrains Mono (parameters/numbers), Inter (body) loaded from Google Fonts |
| 3   | **CSS variable system** — full design token set: `--accent-cyan`, `--accent-violet`, `--glow-cyan`, `--gradient-hero`, etc.   |
| 4   | **Micro-animations** — hover lift on cards, glow shadow on primary button, smooth transitions throughout                      |
| 5   | **Custom scrollbar** matching the dark theme                                                                                  |

#### 🧭 Sidebar Navigation

| #   | Enhancement                                                                                                                |
| --- | -------------------------------------------------------------------------------------------------------------------------- |
| 6   | **Branded sidebar logo** — gradient text with vertical accent bar and version badge pill                                   |
| 7   | **Nav-pill mode selector** — radio buttons now rendered as styled pill buttons with hover highlight and active teal border |
| 8   | **Grouped settings expanders** — ⚙️ Fit Settings · 📏 Units · 🧫 Enzyme Properties (replaces flat layout)                  |
| 9   | **📖 Guide mode** added as a fourth navigation option                                                                      |

#### 🔬 Analysis Mode

| #   | Enhancement                                                                                                                          |
| --- | ------------------------------------------------------------------------------------------------------------------------------------ |
| 10  | **Result banner** — full-width gradient card showing model name (gradient text), AICc, R², and data-point count                      |
| 11  | **AICc substituted for AIC** — corrected information criterion displayed throughout                                                  |
| 12  | **Chip badges for diagnostics** — Normality, Randomness, and R² now show colour-coded `chip-pass` / `chip-warn` / `chip-fail` badges |
| 13  | **Enhanced hero subtitle** — lists all capabilities inline                                                                           |

#### 📈 Simulation Mode

| #   | Enhancement                                                                                                     |
| --- | --------------------------------------------------------------------------------------------------------------- |
| 14  | **Live kinetic equation display** — `eq-display` card updates automatically when the mechanism dropdown changes |
| 15  | **Mechanism label prefix** added to selectbox for clarity                                                       |

#### 🔥 Thermodynamics Mode

| #   | Enhancement                                                                                                                                    |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| 16  | **Arrhenius thermo-cards** — 3-column styled grid (Ea, A, R²) replacing plain `st.metric`                                                      |
| 17  | **Eyring thermo-cards** — 3-column grid (ΔH‡, ΔS‡, ΔG‡) with colour-coded entropy (rose = ordered TS, amber = disordered) and mechanistic note |

#### 📊 Plotly Theme

| #   | Enhancement                                                           |
| --- | --------------------------------------------------------------------- |
| 18  | **JetBrains Mono** for tick labels; Space Grotesk for chart titles    |
| 19  | **Glass hover labels** — dark background, cyan border, monospace font |
| 20  | **Legend styled as glass card**                                       |
| 21  | **Subtle gridlines** and teal zero lines                              |
| 22  | Chart titles repositioned to **top-left**                             |

#### 📖 In-App Guide (new)

| #   | Section                                                                                                |
| --- | ------------------------------------------------------------------------------------------------------ |
| 23  | **Quick Start** — 3-step colour-coded cards for first-time users                                       |
| 24  | **Mode explanations** — tabbed guide for Analysis, Simulation, Thermodynamics with input/output tables |
| 25  | **7 Model cards** — each with rate equation, mechanism, diagnostic signature, and real enzyme example  |
| 26  | **Result interpretation** — AICc ΔAIC table, R² chip key, Shapiro-Wilk and Runs Test explanations      |
| 27  | **Data format guide** — side-by-side Matrix vs Standard annotated tables                               |
| 28  | **Scientific Glossary** — 19-term expandable reference (Vmax, Km, kcat, AICc, ΔG‡, etc.)               |

#### 🦶 Footer

| #   | Enhancement                                                                                                           |
| --- | --------------------------------------------------------------------------------------------------------------------- |
| 29  | **Redesigned footer** — version badge, method credits (Levenberg-Marquardt · AICc · Bootstrap CI), author attribution |

---

### v1.0.1 — Bug Fixes & Infrastructure

| #   | File                        | Fix                                                              |
| --- | --------------------------- | ---------------------------------------------------------------- |
| 1   | `main.py`                   | Bare `except: pass` → `except Exception` throughout              |
| 2   | `main.py`                   | `width="stretch"` → `use_container_width=True` on data editors   |
| 3   | `README.md`                 | Wrong launch command `app.py` → `main.py`                        |
| 4   | `.gitignore`                | Added `!requirements.txt` exception                              |
| 5   | `libbiokinetic/__init__.py` | Created missing `__init__.py` to expose public API               |
| 6   | `libbiokinetic/models.py`   | Added `@classmethod` to Pydantic v2 validators                   |
| 7   | `libbiokinetic/report.py`   | Expanded `generate_simulation_report()` from stub to full report |

---

## 🤝 Support & Contact

This project is maintained by **Chibuike Praise Okechukwu**.

- **Email:** [praizekene1@gmail.com](mailto:praizekene1@gmail.com)
- **LinkedIn:** [linkedin.com/in/chukwubuikem-okechukwu](https://linkedin.com/in/chukwubuikem-okechukwu)
- **GitHub:** [github.com/PraiseElement](https://github.com/PraiseElement)

---

**Disclaimer:** This software is intended for research and educational purposes. While it uses rigorous statistical methods, experimental validation is always recommended.
