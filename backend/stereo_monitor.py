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

    def protection_force(self, xyz_nm, strength=1e6):
        """A temporary flat-bottom chirality restraint over the checked centers.

        Gradient minimization from a strained start can drive a standard-residue
        stereocenter through a planar intermediate into its mirror image. This
        force penalizes only the wrong-sign or near-planar side of each center's
        signed volume, so the minimizer can still relieve clashes but cannot
        invert or flatten a center. It is added for minimization only and never
        becomes part of the dynamics force field. The signed-volume expression
        matches ``check`` exactly; signs and magnitudes come from the supplied
        (already validated) reference coordinates.
        """
        import openmm as mm
        xyz = np.asarray(xyz_nm, dtype=float)
        # Normalize the flat-bottom penalty by the reference signed volume so the
        # force constant k is a real energy scale (signed volumes are tiny nm^3
        # numbers, so an unnormalized k*(Δv)^2 barrier is negligible against a
        # steric clash). The penalty is zero while the center keeps its sign and
        # at least a quarter of its reference magnitude.
        force = mm.CustomCompoundBondForce(4,
            '0.5*k*(min(0,s*v-vmin)/vref)^2; '
            'v=(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))'
            '-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))'
            '+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))')
        for name in ('k', 's', 'vmin', 'vref'):
            force.addPerBondParameter(name)
        protected = []
        for center, idx in zip(self.centers, self.indices):
            p = xyz[idx]
            volume = float(np.linalg.det([p[1] - p[0], p[2] - p[0], p[3] - p[0]]))
            if not np.isfinite(volume) or abs(volume) < 1e-4:
                continue
            vmin = 0.25 * abs(volume)
            force.addBond([int(i) for i in idx], [float(strength), float(np.sign(volume)), vmin, abs(volume)])
            protected.append({'chain': center['chain'], 'resid': center['resid'],
                              'residue': center['residue'], 'center': center['center'],
                              'reference_signed_volume_nm3': volume, 'minimum_oriented_volume_nm3': vmin})
        return force, protected

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
