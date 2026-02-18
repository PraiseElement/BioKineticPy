# BioKineticPy v1.0 🧬

### _A Professional Computational Workbench for Enzyme Kinetics & Thermodynamics_

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Framework](https://img.shields.io/badge/built%20with-Streamlit-red.svg)

---

## 📖 Project Overview

**BioKineticPy** is a specialized data analysis suite designed for biochemists, enzymologists, and biophysicists. Unlike general-purpose graphing software, BioKineticPy is built specifically for the non-linear regression challenges inherent in enzyme kinetics.

It moves beyond simple Michaelis-Menten fitting to handle complex scenarios like **mixed inhibition**, **allostery (Hill equation)**, and **substrate inhibition (Haldane kinetics)**. Furthermore, it bridges the gap between kinetics and thermodynamics, offering automated **Arrhenius** and **Eyring** analyses to determine activation energy ($E_a$), enthalpy ($\Delta H^\ddagger$), and entropy ($\Delta S^\ddagger$).

## ✨ Key Features

### 🔬 Advanced Kinetic Analysis

- **Automated Model Selection:** The system fits your data against 7 distinct mechanistic models simultaneously and ranks them using the **Akaike Information Criterion (AIC)**. It tells you _which_ model fits best, not just _how_ it fits.
- **Global Matrix Fitting:** Analyze entire datasets (e.g., substrate titration curves at 3 different inhibitor concentrations) in a single global fit, reducing parameter error.
- **Statistical Rigor:**
  - **Bootstrap Resampling (Monte Carlo):** Generates true 95% Confidence Intervals for $V_{max}$, $K_m$, and $K_i$.
  - **Residual Diagnostics:** Automated Shapiro-Wilk (normality) and Runs Test (randomness) to validate fit quality.
  - **Weighted Regression:** Supports $1/v$ and $1/v^2$ weighting to account for heteroscedastic experimental error.

### 🔥 Thermodynamic Profiling

- **Temperature Dependence:** Input rate constants ($k_{cat}$ or $V_{max}$) vs. Temperature.
- **Transition State Theory:** Automatically calculates Gibbs Free Energy of Activation ($\Delta G^\ddagger$) at standard conditions.
- **Unit Intelligence:** Auto-converts between Celsius/Kelvin and handles energy units (Joules).

### 📊 Visualization Suite

- **3D Kinetic Landscapes:** Interactive surface plots to visualize velocity as a function of both Substrate and Inhibitor simultaneously.
- **Diagnostic Linearizations:** Auto-generates Lineweaver-Burk, Hanes-Woolf, Eadie-Hofstee, Dixon, and Cornish-Bowden plots for visual verification.

---

## 🚀 Installation & Usage

### Prerequisites

- Python 3.8 or higher
- pip (Python Package Installer)

### Setup Guide

1.  **Clone the Repository:**

    ```bash
    git clone https://github.com/PraiseElement/BioKineticPy.git
    cd BioKineticPy
    ```

2.  **Create a Virtual Environment (Recommended):**

    ```bash
    python -m venv venv
    # Windows:
    venv\Scripts\activate
    # Mac/Linux:
    source venv/bin/activate
    ```

3.  **Install Dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

4.  **Launch the Application:**
    ```bash
    streamlit run main.py
    ```

---

## 📚 Theoretical Models Supported

1.  **Michaelis-Menten:** The fundamental model of enzyme kinetics.
2.  **Competitive Inhibition:** Inhibitor binds to the active site ($K_m$ increases).
3.  **Non-Competitive Inhibition:** Inhibitor binds allosterically ($V_{max}$ decreases).
4.  **Uncompetitive Inhibition:** Inhibitor binds only the ES complex ($V_{max}$ & $K_m$ decrease).
5.  **Mixed Inhibition:** A general case of non-competitive inhibition ($\alpha \neq 1$).
6.  **Hill Equation:** For cooperative / allosteric enzymes (Sigmoidal curves).
7.  **Substrate Inhibition (Haldane):** When excess substrate acts as an inhibitor.

---

## 🛠 Corrections & Updates (v1.0.1)

The following corrections and updates were applied to the codebase:

### 🐛 Bug Fixes

| #   | File         | Issue                                                                                                                         | Fix                                                          |
| --- | ------------ | ----------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| 1   | `main.py`    | **Bare `except: pass` clauses** silently swallowed all errors (including `SystemExit`, `KeyboardInterrupt`), hiding real bugs | Changed to `except Exception` throughout                     |
| 2   | `main.py`    | **`width="stretch"` is not a valid Streamlit parameter** for `st.data_editor` — data tables may not render at full width      | Replaced with `use_container_width=True` on all data editors |
| 3   | `README.md`  | **Wrong launch command** — `streamlit run app.py` referenced a non-existent file                                              | Corrected to `streamlit run main.py`                         |
| 4   | `README.md`  | **Nested markdown link** in clone URL rendered incorrectly                                                                    | Simplified to plain URL                                      |
| 5   | `.gitignore` | **`*.txt` glob excluded `requirements.txt`** from Git — new clones would be missing the dependency file                       | Added `!requirements.txt` exception                          |

### 🏗 Infrastructure Improvements

| #   | File                        | Change                                                                                                                                                                                            |
| --- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 6   | `libbiokinetic/__init__.py` | **Created missing `__init__.py`** — Python did not recognize `libbiokinetic` as a proper package; this could break imports with certain tooling or editable installs. Now exports the public API. |
| 7   | `libbiokinetic/main.py`     | **Converted absolute imports to relative imports** (`from libbiokinetic.models` → `from .models`) for consistency with the rest of the package                                                    |
| 8   | `libbiokinetic/models.py`   | **Added `@classmethod` decorator** to the Pydantic `@field_validator` to suppress deprecation warnings in Pydantic v2                                                                             |
| 9   | `libbiokinetic/report.py`   | **Expanded `generate_simulation_report()`** from a one-line stub into a full formatted report matching the style of the analysis and thermodynamic reports                                        |

### 🎨 UI Overhaul

| #   | Enhancement                                                                                              |
| --- | -------------------------------------------------------------------------------------------------------- |
| 10  | **Premium dark theme** with custom CSS — deep navy backgrounds, glassmorphic cards, and polished borders |
| 11  | **Gradient hero titles** using `Inter` web font for all section headers                                  |
| 12  | **Unified Plotly dark template** applied to all charts for visual consistency                            |
| 13  | **Pill-shaped tab styling** with accent highlight for the active tab                                     |
| 14  | **Gradient primary buttons** with hover lift animation and glow shadow                                   |
| 15  | **Styled metric cards** with hover elevation and uppercase labels                                        |
| 16  | **Q-Q plot reference line** added to diagnostics for better visual validation                            |
| 17  | **R² metric card** added to the Diagnostics tab alongside Normality and Randomness                       |
| 18  | **Branded footer** with author attribution at the bottom of every page                                   |
| 19  | **Section badges** (pill-shaped labels) for visual grouping in sidebar and content areas                 |
| 20  | **Custom scrollbar** matching the dark theme                                                             |

---

## 🤝 Support & Contact

This project is maintained by **Chibuike Praise Okechukwu**.

If you encounter errors, have questions about the kinetic models, to support, or to advise me, or need advice on interpreting your data, please reach out.

- **Email:** [praizekene1@gmail.com](mailto:praizekene1@gmail.com)
- **LinkedIn:** [linkedin.com/in/chukwubuikem-okechukwu](https://linkedin.com/in/chukwubuikem-okechukwu)
- **GitHub:** [github.com/PraiseElement](https://github.com/PraiseElement)

---

**Disclaimer:** This software is intended for research and educational purposes. While it uses rigorous statistical methods, experimental validation is always recommended.
