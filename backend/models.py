import math
from typing import Annotated, Literal
from pydantic import BaseModel, Field, model_validator

from . import config

LigandOverrides = dict[Annotated[str, Field(min_length=1, max_length=100)], Annotated[str, Field(min_length=1, max_length=10000)]]


class LiveMeasurement(BaseModel):
    id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_-]+$")
    kind: Literal["distance", "angle", "dihedral", "hbond"]
    atoms: list[Annotated[int, Field(strict=True, ge=0)]] = Field(min_length=2, max_length=4)
    label: str = Field(min_length=1, max_length=160)
    color: str | None = Field(default=None, pattern=r"^#[0-9a-fA-F]{6}$")
    trackDuringRun: bool | None = None

    @model_validator(mode="after")
    def check_selection(self):
        expected = {"distance": 2, "angle": 3, "dihedral": 4, "hbond": 3}[self.kind]
        if len(self.atoms) != expected or len(set(self.atoms)) != expected:
            raise ValueError(f"Select {expected} distinct atoms for a {self.kind} measurement.")
        return self


class RemapMeasurementsRequest(BaseModel):
    source_dataset_id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_-]+$")
    measurements: list[LiveMeasurement] = Field(default_factory=list, max_length=100)


class SimulationConfig(BaseModel):
    dataset_id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_-]+$")
    engine: Literal["openmm", "gromacs"] = "openmm"
    name: str = Field(default="My simulation", min_length=1, max_length=100)
    duration_ps: float = Field(default=10, gt=0, le=config.MAX_DURATION_PS, allow_inf_nan=False)
    temperature_k: float = Field(default=300, ge=50, le=500, allow_inf_nan=False)
    timestep_fs: float = Field(default=2, ge=0.1, le=2, allow_inf_nan=False)
    report_interval: int = Field(default=50, ge=1, le=100000)
    friction_ps: float = Field(default=1, ge=0.01, le=100, allow_inf_nan=False)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)
    solvent: Literal["implicit", "explicit"] = "implicit"
    minimize: bool = True
    equilibration_steps: int = Field(default=100, ge=0, le=100_000)
    padding_nm: float = Field(default=1, ge=1, le=3, allow_inf_nan=False)
    measurements: list[LiveMeasurement] = Field(default_factory=list, max_length=12)

    @model_validator(mode="after")
    def check_limits(self):
        if len({measurement.id for measurement in self.measurements}) != len(self.measurements):
            raise ValueError("Live measurement IDs must be unique.")
        steps = round(self.duration_ps * 1000 / self.timestep_fs)
        if steps < 1 or steps > config.MAX_STEPS:
            raise ValueError(f"Choose a duration and timestep giving 1 to {config.MAX_STEPS:,} production steps.")
        if math.ceil(steps / self.report_interval) + 1 > config.MAX_FRAMES:
            raise ValueError(f"Increase the frame interval: this local version stores at most {config.MAX_FRAMES:,} frames.")
        if self.engine == "gromacs" and self.solvent != "explicit":
            raise ValueError("The GROMACS workflow supports explicit TIP3P solvent only.")
        return self


class MeasurementRequest(BaseModel):
    kind: Literal["distance", "angle", "dihedral", "hbond"]
    atoms: list[int] = Field(min_length=2, max_length=4)


class StructureFetchRequest(BaseModel):
    provider: Literal["pdb", "pubchem"] = "pdb"
    identifier: str = Field(min_length=1, max_length=200)
    name: str | None = Field(default=None, max_length=100)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)


class SmilesRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=10000)
    name: str | None = Field(default=None, max_length=100)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)


class PreparationRequest(BaseModel):
    dataset_id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_-]+$")
    name: str = Field(default="Protein preparation", min_length=1, max_length=100)
    ph: float = Field(default=7, ge=0, le=14, allow_inf_nan=False)
    add_missing_atoms: bool = True
    build_missing_residues: bool = False
    optimize_sidechains: bool = True
    remove_waters: bool = True
    remove_heterogens: bool = False
    ligand_overrides: LigandOverrides = Field(default_factory=dict, max_length=100)
    ligand_actions: dict[str, Literal["repair", "remove"]] = Field(default_factory=dict, max_length=100)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)


class SolvationRequest(BaseModel):
    padding_nm: float = Field(default=1, ge=1, le=3, allow_inf_nan=False)
    ph: float = Field(default=7, ge=0, le=14, allow_inf_nan=False)
    seed: int = Field(default=2026, ge=1, le=2_147_483_646)


class InspectionRequest(BaseModel):
    ph: float = Field(default=7, ge=0, le=14, allow_inf_nan=False)
    ligand_overrides: LigandOverrides = Field(default_factory=dict, max_length=100)
    ligand_actions: dict[str, Literal["repair", "remove"]] = Field(default_factory=dict, max_length=100)
