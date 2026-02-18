from pydantic import BaseModel, Field, field_validator
from typing import List, Optional
import numpy as np

class KineticPoint(BaseModel):
    substrate_conc: float
    velocity: float
    inhibitor_conc: float = 0.0

    @field_validator('substrate_conc', 'velocity')
    @classmethod
    def must_be_non_negative(cls, v):
        if v < 0:
            raise ValueError("Concentration and velocity must be non-negative.")
        return v

class Dataset(BaseModel):
    """Represents a full experimental dataset ready for fitting."""
    points: List[KineticPoint]
    
    # Helper to export as numpy arrays for the solver
    def to_arrays(self):
        s = np.array([p.substrate_conc for p in self.points])
        v = np.array([p.velocity for p in self.points])
        i = np.array([p.inhibitor_conc for p in self.points])
        return s, v, i