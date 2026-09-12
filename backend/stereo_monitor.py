"""Check residue stereochemistry after native relaxation and in saved MD frames."""
import numpy as np


class StereoMonitor:
    def __init__(self, topology, positions_nm, template_source):
        from openmm import unit
        from pdbfixer import PDBFixer
        from .preparation_worker import stereochemistry_report
        templates = PDBFixer(filename=str(template_source)).templates
        report = stereochemistry_report(topology, np.asarray(positions_nm) * unit.nanometer, templates)
        self.centers = report['centers']
        self.indices = np.asarray([c['atom_indices'] for c in self.centers], dtype=int).reshape(-1, 4)
        self.expected = np.asarray([c['template_signed_volume_nm3'] for c in self.centers])

    def check(self, xyz_nm, box_nm=None, stage='Dynamics'):
        xyz = np.asarray(xyz_nm, dtype=float)
        points = xyz[self.indices]
        vectors = points[:, 1:] - points[:, :1]
        if box_nm is not None:
            box = np.asarray(box_nm, dtype=float)
            if box.shape != (3, 3) or not np.isfinite(box).all() or np.linalg.det(box) <= 0:
                raise ValueError('A finite periodic box is required for the stereochemistry check.')
            fractional = vectors @ np.linalg.inv(box)
            # Bonds are local; use their nearest periodic images, never whole-box
            # wrapped atom positions when determining handedness.
            vectors = (fractional - np.rint(fractional)) @ box
        with np.errstate(invalid='ignore'):
            volumes = np.linalg.det(vectors)
        invalid = ~np.isfinite(volumes) | (np.abs(volumes) < 1e-4) | (volumes * self.expected <= 0)
        violations = [{**self.centers[i], 'signed_volume_nm3': float(volumes[i]) if np.isfinite(volumes[i]) else None}
                      for i in np.flatnonzero(invalid)]
        return {'stage': stage, 'checked_centers': len(self.centers), 'violations': violations,
                'passed': not violations,
                'method': 'Local minimum-image signed stereocenter volumes versus residue templates; no coordinate changes.'}


def require_valid(report):
    if report['violations']:
        first = report['violations'][0]
        raise ValueError(f"{report['stage']} inverted or flattened residue stereochemistry at "
                         f"{first['chain']}:{first['resid']}{first.get('insertion_code', '')} "
                         f"{first['residue']} {first['center']}. Dynamics stopped; "
                         "the prepared geometry needs repair. The coordinates were retained for inspection.")
