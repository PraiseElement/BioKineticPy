import unittest
import numpy as np
import sys
import os

# Add parent directory to path so we can import the app modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from libbiokinetic.models import Dataset, KineticPoint
from libbiokinetic.solver import KineticSolver
from libbiokinetic.equations import MODEL_REGISTRY

class TestKineticSolver(unittest.TestCase):

    def test_michaelis_menten_fit(self):
        """
        Test if the solver correctly identifies Vmax=100 and Km=10
        from synthetic data with no noise.
        """
        # 1. Generate Synthetic Data
        true_vmax = 100.0
        true_km = 10.0
        s_values = [1, 2, 5, 10, 20, 50, 100]
        points = []
        
        mm_func = MODEL_REGISTRY["Michaelis-Menten"]
        
        for s in s_values:
            # v = (Vmax * S) / (Km + S)
            v = mm_func(s, true_vmax, true_km)
            points.append(KineticPoint(substrate_conc=s, velocity=v, inhibitor_conc=0.0))
            
        # 2. Run Solver
        dataset = Dataset(points=points)
        solver = KineticSolver(dataset)
        result = solver.fit_model("Michaelis-Menten")
        
        # 3. Assertions (Check if result is close to truth)
        # We use assertAlmostEqual because floating point math is rarely exact
        self.assertAlmostEqual(result['parameters']['vmax'], true_vmax, places=4)
        self.assertAlmostEqual(result['parameters']['km'], true_km, places=4)
        self.assertGreater(result['r_squared'], 0.999)
        print(f"\n✅ Michaelis-Menten Test Passed: Vmax={result['parameters']['vmax']:.2f}, Km={result['parameters']['km']:.2f}")

    def test_competitive_inhibition(self):
        """
        Test if solver detects Competitive Inhibition correctly.
        True Params: Vmax=100, Km=10, Ki=5
        """
        true_vmax = 100.0
        true_km = 10.0
        true_ki = 5.0
        s_values = [1, 5, 10, 50]
        i_values = [0, 5, 10]
        points = []
        
        comp_func = MODEL_REGISTRY["Competitive"]
        
        for i in i_values:
            for s in s_values:
                # v = (Vmax * S) / (Km(1 + I/Ki) + S)
                v = comp_func(s, i, true_vmax, true_km, true_ki)
                points.append(KineticPoint(substrate_conc=s, velocity=v, inhibitor_conc=i))
                
        dataset = Dataset(points=points)
        solver = KineticSolver(dataset)
        result = solver.fit_model("Competitive")
        
        self.assertAlmostEqual(result['parameters']['ki'], true_ki, places=3)
        print(f"✅ Competitive Inhibition Test Passed: Ki={result['parameters']['ki']:.2f} (Expected 5.0)")

    def test_haldane_substrate_inhibition(self):
        """
        Test Substrate Inhibition (Haldane).
        True Params: Vmax=100, Km=10, Ksi=50
        """
        true_vmax = 100.0
        true_km = 10.0
        true_ksi = 50.0
        s_values = [1, 5, 10, 20, 50, 100, 200] # Range goes high to show inhibition
        points = []
        
        haldane_func = MODEL_REGISTRY["Substrate Inhibition"]
        
        for s in s_values:
            # v = (Vmax * S) / (Km + S + S^2/Ksi)
            # Note: Our updated haldane signature includes 'i' and 'ki', need to pass defaults
            v = haldane_func(s=s, i=0, vmax=true_vmax, km=true_km, ksi=true_ksi, ki=1e9)
            points.append(KineticPoint(substrate_conc=s, velocity=v, inhibitor_conc=0.0))
            
        dataset = Dataset(points=points)
        solver = KineticSolver(dataset)
        result = solver.fit_model("Substrate Inhibition")
        
        self.assertAlmostEqual(result['parameters']['ksi'], true_ksi, places=2)
        print(f"✅ Substrate Inhibition Test Passed: Ksi={result['parameters']['ksi']:.2f} (Expected 50.0)")

if __name__ == '__main__':
    unittest.main()