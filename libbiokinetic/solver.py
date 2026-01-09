import numpy as np
from lmfit import Model, Parameters
from scipy import stats
from .equations import MODEL_REGISTRY
from .models import Dataset

class KineticSolver:
    def __init__(self, data: Dataset):
        self.s, self.v, self.i = data.to_arrays()

    def _estimate_initial_guesses(self):
        v_max_guess = np.max(self.v)
        half_v = v_max_guess / 2
        idx = (np.abs(self.v - half_v)).argmin()
        km_guess = self.s[idx] if self.s[idx] > 0 else 0.1
        return v_max_guess, km_guess

    def _calculate_weights(self, mode):
        """
        Generates weights for regression.
        lmfit 'weights' are 1/sigma. 
        So for 1/v^2 weighting (relative error), we pass 1/v.
        """
        # Avoid division by zero
        safe_v = np.maximum(np.abs(self.v), 1e-9)
        
        if mode == "None (Homoscedastic)":
            return None
        elif mode == "1/v (Poisson)":
            # Variance proportional to v -> sigma proportional to sqrt(v) -> weight = 1/sqrt(v)
            return 1.0 / np.sqrt(safe_v)
        elif mode == "1/v² (Relative Error)":
            # Variance proportional to v^2 -> sigma proportional to v -> weight = 1/v
            return 1.0 / safe_v
        return None

    def fit_model(self, model_name: str, weighting_mode: str = "None (Homoscedastic)"):
        equation = MODEL_REGISTRY[model_name]
        arg_names = equation.__code__.co_varnames[:equation.__code__.co_argcount]
        indep_vars = ['s']
        if 'i' in arg_names: indep_vars.append('i')

        gmodel = Model(equation, independent_vars=indep_vars)
        est_vmax, est_km = self._estimate_initial_guesses()
        params = Parameters()
        
        has_inhibitor = np.max(self.i) > 0

        # Parameter Definitions
        if 'vmax' in arg_names: params.add('vmax', value=est_vmax, min=1e-9)
        if 'km' in arg_names: params.add('km', value=est_km, min=1e-9)
        if 'ki' in arg_names:
            if has_inhibitor: params.add('ki', value=est_km, min=1e-9)
            else: params.add('ki', value=1e9, vary=False)
        if 'alpha' in arg_names: params.add('alpha', value=1.0, min=0.01, max=100.0)
        if 'khalf' in arg_names: params.add('khalf', value=est_km, min=1e-9)
        if 'n' in arg_names: params.add('n', value=1.0, min=0.1, max=10.0)
        if 'ksi' in arg_names: params.add('ksi', value=est_km * 10, min=1e-9)

        try:
            fit_kwargs = {'s': self.s}
            if 'i' in indep_vars: fit_kwargs['i'] = self.i
            
            # Apply Weights
            weights = self._calculate_weights(weighting_mode)
            if weights is not None:
                fit_kwargs['weights'] = weights
            
            result = gmodel.fit(self.v, params, **fit_kwargs)
            
            param_errors = {name: (p.stderr if p.stderr else 0.0) for name, p in result.params.items()}

            return {
                "model": model_name,
                "aic": result.aic,
                "r_squared": result.rsquared,
                "parameters": result.best_values,
                "errors": param_errors,
                "fit_object": result
            }
        except Exception:
            return {"model": model_name, "aic": 999999, "r_squared": 0, "parameters": {}, "errors": {}}

    def run_model_competition(self, weighting_mode="None (Homoscedastic)"):
        results = []
        has_inhibitor = np.max(self.i) > 0
        core_models = ["Michaelis-Menten", "Hill (Allosteric)", "Substrate Inhibition"]
        if has_inhibitor:
            core_models.extend(["Competitive", "Uncompetitive", "Non-Competitive", "Mixed"])
        
        for m in core_models:
            if not has_inhibitor and m in ["Competitive", "Uncompetitive", "Non-Competitive", "Mixed"]: continue
            results.append(self.fit_model(m, weighting_mode))
            
        results.sort(key=lambda x: x['aic'])
        return results

    def check_ambiguity(self, sorted_results):
        if len(sorted_results) < 2: return None
        delta_aic = abs(sorted_results[1]['aic'] - sorted_results[0]['aic'])
        if delta_aic < 2.0:
            return f"⚠️ **Ambiguity Detected:** '{sorted_results[0]['model']}' and '{sorted_results[1]['model']}' are statistically indistinguishable (ΔAIC = {delta_aic:.2f})."
        return None

    def bootstrap_uncertainty(self, model_name, best_params, n_iter=50, weighting_mode="None (Homoscedastic)"):
        equation = MODEL_REGISTRY[model_name]
        arg_names = equation.__code__.co_varnames[:equation.__code__.co_argcount]
        
        call_args = {k: v for k, v in best_params.items() if k in arg_names}
        call_args['s'] = self.s
        if 'i' in arg_names: call_args['i'] = self.i
        v_pred = equation(**call_args)
        residuals = self.v - v_pred
        
        collected_params = {k: [] for k in best_params.keys()}
        indep_vars = ['s']; 
        if 'i' in arg_names: indep_vars.append('i')
        gmodel = Model(equation, independent_vars=indep_vars)
        weights = self._calculate_weights(weighting_mode)

        for _ in range(n_iter):
            res_sample = np.random.choice(residuals, size=len(residuals), replace=True)
            v_synthetic = v_pred + res_sample
            
            params = Parameters()
            for k, v in best_params.items():
                if k == 'ki' and np.max(self.i) == 0: params.add(k, value=1e9, vary=False)
                else: 
                    params.add(k, value=v, min=1e-9)
                    if k == 'alpha': params[k].max = 100.0
                    if k == 'n': params[k].max = 10.0

            fit_kwargs = {'s': self.s}
            if 'i' in indep_vars: fit_kwargs['i'] = self.i
            if weights is not None: fit_kwargs['weights'] = weights
            
            try:
                res = gmodel.fit(v_synthetic, params, **fit_kwargs, method='leastsq')
                for k in collected_params: collected_params[k].append(res.best_values[k])
            except: continue
                
        ci_results = {}
        for k, values in collected_params.items():
            if len(values) > 5:
                ci_results[k] = (np.percentile(values, 2.5), np.percentile(values, 97.5))
            else:
                ci_results[k] = (best_params[k], best_params[k])
        return ci_results

    def diagnose_residuals(self, residuals):
        """
        Performs statistical tests on residuals.
        1. Shapiro-Wilk (Normality): p < 0.05 means NOT normal.
        2. Runs Test (Randomness): Z-score > 1.96 means systematic deviation.
        """
        # Shapiro-Wilk
        if len(residuals) >= 3:
            shapiro_stat, shapiro_p = stats.shapiro(residuals)
        else:
            shapiro_stat, shapiro_p = 0, 1.0

        # Runs Test (Manual Implementation)
        signs = np.sign(residuals)
        # Remove zeros
        signs = signs[signs != 0]
        n = len(signs)
        if n > 1:
            n_pos = np.sum(signs > 0)
            n_neg = np.sum(signs < 0)
            runs = 1 + np.sum(signs[1:] != signs[:-1])
            
            # Expected runs
            mu = (2 * n_pos * n_neg) / n + 1
            var = (mu - 1) * (mu - 2) / (n - 1)
            
            if var > 0: z_score = (runs - mu) / np.sqrt(var)
            else: z_score = 0
        else:
            z_score = 0
            
        return {
            "shapiro_p": shapiro_p,
            "runs_z": z_score
        }