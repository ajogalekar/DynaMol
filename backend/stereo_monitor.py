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

    def chirality_guard_force(self, xyz_nm, k=5000.0, vmin_fraction=0.4):
        """A permanent, gentle flat-bottom chirality restraint for the MD force field.

        The AMBER force field has no explicit term that keeps a Cα (or Cβ)
        stereocenter left-handed: chirality is maintained only by the local
        bonded geometry. Where a strained prepared geometry (a residue next to a
        low-confidence rebuilt loop) distorts that geometry, the force field's
        own local minimum can be the *inverted* (D) center -- free minimization
        then drives it there, and dynamics inverts it in the first steps. This
        restraint is added to the system before minimization and kept for all
        dynamics. It is a flat-bottom well that is exactly zero while a center
        keeps its correct sign and at least ``vmin_fraction`` of its ideal
        (template) signed-volume magnitude, and rises only as the center
        approaches planarity -- so it never biases the many healthy centers and
        only resists racemization of the few strained ones. ``k`` is gentle
        enough to be integrator-safe yet, because the penalty is normalized by
        the tiny reference volume, easily holds a center against the large
        distorting force a strained rebuilt loop exerts on it.
        The signed-volume expression matches ``check`` exactly. Returns
        ``(force, guarded)``.
        """
        import openmm as mm
        xyz = np.asarray(xyz_nm, dtype=float)
        force = mm.CustomCompoundBondForce(4,
            '0.5*k_chi*(min(0,s*v-vmin)/vref)^2; '
            'v=(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))'
            '-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))'
            '+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))')
        force.addGlobalParameter('k_chi', float(k))
        for name in ('s', 'vmin', 'vref'):
            force.addPerBondParameter(name)
        guarded = []
        for center, idx, ideal in zip(self.centers, self.indices, self.expected):
            p = xyz[idx]
            volume = float(np.linalg.det([p[1] - p[0], p[2] - p[0], p[3] - p[0]]))
            ideal = float(ideal)
            if not np.isfinite(volume) or abs(volume) < 1e-4 or not np.isfinite(ideal) or abs(ideal) < 1e-4:
                continue
            if volume * ideal <= 0:  # already inverted vs template; the before-check owns this
                continue
            vref = abs(ideal)
            force.addBond([int(i) for i in idx], [float(np.sign(volume)), vmin_fraction * vref, vref])
            guarded.append({'chain': center['chain'], 'resid': center['resid'],
                            'residue': center['residue'], 'center': center['center'],
                            'prepared_signed_volume_nm3': volume, 'ideal_signed_volume_nm3': ideal,
                            'floor_oriented_volume_nm3': vmin_fraction * vref})
        return force, guarded

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
