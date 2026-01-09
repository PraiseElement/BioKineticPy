from libbiokinetic.models import Dataset, KineticPoint
from libbiokinetic.solver import KineticSolver

# 1. Simulate data (Competitive Inhibition: Vmax=100, Km=10, Ki=5)
# Using "Dirty Data" approach - no perfect lines
raw_data = [
    KineticPoint(substrate_conc=1, velocity=8.5, inhibitor_conc=0),
    KineticPoint(substrate_conc=5, velocity=33.1, inhibitor_conc=0),
    KineticPoint(substrate_conc=10, velocity=49.8, inhibitor_conc=0),
    KineticPoint(substrate_conc=1, velocity=4.8, inhibitor_conc=5), # Inhibited
    KineticPoint(substrate_conc=5, velocity=19.5, inhibitor_conc=5),
    KineticPoint(substrate_conc=10, velocity=32.9, inhibitor_conc=5),
]

dataset = Dataset(points=raw_data)

# 2. Run Solver
solver = KineticSolver(dataset)
results = solver.run_model_competition()

# 3. Print Report
print(f"Best Fitting Model: {results[0]['model']}")
print(f"AIC: {results[0]['aic']:.2f}")
print("Parameters:")
for k, v in results[0]['parameters'].items():
    print(f"  {k}: {v:.4f}")