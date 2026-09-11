import math
from typing import Literal
from pydantic import BaseModel, Field, model_validator


class SimulationConfig(BaseModel):
    dataset_id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_-]+$")
    engine: Literal["openmm", "gromacs"] = "openmm"
    name: str = Field(default="My simulation", min_length=1, max_length=100)
    duration_ps: float = Field(default=10, gt=0, le=1000, allow_inf_nan=False)
    temperature_k: float = Field(default=300, ge=50, le=500, allow_inf_nan=False)
    timestep_fs: float = Field(default=2, ge=0.1, le=2, allow_inf_nan=False)
    report_interval: int = Field(default=50, ge=1, le=100000)
    friction_ps: float = Field(default=1, ge=0.01, le=100, allow_inf_nan=False)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)
    solvent: Literal["implicit", "explicit"] = "implicit"
    minimize: bool = True
    equilibration_steps: int = Field(default=100, ge=0, le=100_000)
    padding_nm: float = Field(default=1, ge=1, le=3, allow_inf_nan=False)

    @model_validator(mode="after")
    def check_limits(self):
        steps = round(self.duration_ps * 1000 / self.timestep_fs)
        if steps < 1 or steps > 2_000_000:
            raise ValueError("Choose a duration and timestep giving 1 to 2,000,000 production steps.")
        if math.ceil(steps / self.report_interval) + 1 > 10_000:
            raise ValueError("Increase the frame interval: this local version stores at most 10,000 frames.")
        if self.engine == "gromacs" and self.solvent != "explicit":
            raise ValueError("The GROMACS workflow supports explicit TIP3P solvent only.")
        return self


class MeasurementRequest(BaseModel):
    kind: Literal["distance", "angle", "dihedral", "hbond"]
    atoms: list[int] = Field(min_length=2, max_length=4)
