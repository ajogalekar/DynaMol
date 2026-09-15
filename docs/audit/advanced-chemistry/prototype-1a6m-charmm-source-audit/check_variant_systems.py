"""Ensure the two isolated published-model definitions differ only as declared."""
from pathlib import Path
import json
import hashlib
import openmm
from openmm import XmlSerializer, NonbondedForce, unit

HERE = Path(__file__).resolve().parent
paths = [HERE / ('published-horse-' + v) / 'system.xml' for v in ['baseline', 'updated-charges']]
systems = [XmlSerializer.deserialize(p.read_text()) for p in paths]
assert systems[0].getNumParticles() == systems[1].getNumParticles() == 27386
assert systems[0].getNumForces() == systems[1].getNumForces()
for a, b in zip(systems[0].getForces(), systems[1].getForces()):
    assert type(a) is type(b)
    if not isinstance(a, NonbondedForce):
        assert XmlSerializer.serialize(a) == XmlSerializer.serialize(b)
nb = [next(f for f in s.getForces() if isinstance(f, NonbondedForce)) for s in systems]
charges, changed_particles = [], []
for i in range(27386):
    a, b = [f.getParticleParameters(i) for f in nb]
    assert a[1:] == b[1:]
    assert systems[0].getParticleMass(i) == systems[1].getParticleMass(i)
    qa, qb = [v[0].value_in_unit(unit.elementary_charge) for v in [a, b]]
    charges.append((qa, qb))
    if qa != qb:
        changed_particles.append(i)
declared = json.loads((HERE / 'published-horse-updated-charges/audit.json').read_text())['charge_modifications']
assert changed_particles == [r['index_1based'] - 1 for r in declared]
assert nb[0].getNumExceptions() == nb[1].getNumExceptions()
changed_exceptions = 0
max_unscaled_14_error = 0.
for k in range(nb[0].getNumExceptions()):
    a, b = [f.getExceptionParameters(k) for f in nb]
    assert a[:2] == b[:2] and a[3:] == b[3:]
    qa, qb = [p[2].value_in_unit(unit.elementary_charge**2) for p in [a, b]]
    if qa == 0:
        assert qb == 0
    else:
        i, j = a[:2]
        errors = [abs(qa - charges[i][0] * charges[j][0]), abs(qb - charges[i][1] * charges[j][1])]
        max_unscaled_14_error = max(max_unscaled_14_error, *errors)
        assert max(errors) < 1e-12
    if qa != qb:
        changed_exceptions += 1
report = {'openmm_version': openmm.version.version, 'particles': 27386,
          'changed_particle_charges': len(changed_particles), 'changed_exception_charge_products': changed_exceptions,
          'nonbonded_exceptions': nb[0].getNumExceptions(),
          'maximum_unscaled_nonzero_14_charge_product_error_e2': max_unscaled_14_error,
          'unchanged': ['All other force XML definitions, including zero-amplitude site terms', 'Particle masses and Lennard-Jones terms', 'Exception pairs and Lennard-Jones terms', 'Excluded charge products'],
          'system_sha256': [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths],
          'scope': 'Two OpenMM import definitions are internally consistent with declared charge-only changes. No independent native-engine energy/force parity, no Context, and no model-accuracy validation.'}
(HERE / 'variant-mechanics-comparison.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
