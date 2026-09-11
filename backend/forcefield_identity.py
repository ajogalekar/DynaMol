"""Keep native per-atom ligand parameters attached to their recorded identities.

OpenMM's normal residue matching intentionally ignores atom names. A graph
automorphism can therefore exchange distinct native Amber improper assignments.
This adapter changes only the match for explicitly recorded ligand residues;
protein, modified-residue and solvent parameter generators remain unchanged.
"""
from __future__ import annotations

from openmm import app


class IdentityForceField(app.ForceField):
    def __init__(self, *files, ligand_atom_maps=()):
        super().__init__(*files)
        self._ligand_identities = {}
        for record in ligand_atom_maps:
            try:
                key = tuple(record['residue_key'])
                template_name = record['template_name']
                atoms = record['atoms']
                names = [atom['prepared_name'] for atom in atoms]
                template_names = [atom['template_atom_name'] for atom in atoms]
                indices = [atom['native_index'] for atom in atoms]
            except (KeyError, TypeError):
                raise ValueError('The prepared ligand atom identity map is malformed. Prepare the complex again.') from None
            if len(key) != 4 or not all(isinstance(value, str) for value in key):
                raise ValueError('Invalid prepared ligand residue identity. Prepare the complex again.')
            if key in self._ligand_identities:
                raise ValueError('Duplicate prepared ligand residue identity. Prepare the complex again.')
            if not isinstance(template_name, str):
                raise ValueError('The recorded ligand template name is invalid.')
            template = self._templates.get(template_name)
            if template is None:
                raise ValueError('The recorded ligand template is missing from its parameter bundle.')
            if (not all(isinstance(name, str) and name for name in names + template_names)
                    or not all(type(index) is int for index in indices)
                    or len(atoms) != len(template.atoms) or len(set(names)) != len(names)
                    or sorted(indices) != list(range(len(template.atoms)))
                    or set(template_names) != {atom.name for atom in template.atoms}
                    or len(set(template_names)) != len(template_names)):
                raise ValueError('The prepared ligand atom identity map is incomplete or ambiguous.')
            if template.externalBonds:
                raise ValueError('Native noncovalent ligand templates cannot contain external bonds.')
            template_indices = {atom.name: i for i, atom in enumerate(template.atoms)}
            if any(template_indices[atom['template_atom_name']] != atom['native_index'] for atom in atoms):
                raise ValueError('The ligand atom identity map conflicts with its native parameter indices.')
            self._ligand_identities[key] = (template, {
                atom['prepared_name']: template_indices[atom['template_atom_name']] for atom in atoms
            })

    @staticmethod
    def _residue_key(residue):
        if not hasattr(residue, 'chain'):
            return None
        return (residue.chain.id, residue.id, (residue.insertionCode or '').strip(), residue.name)

    def _identity_match(self, residue, bonded_to_atom):
        identity = self._ligand_identities.get(self._residue_key(residue))
        if identity is None:
            return None
        template, by_name = identity
        atoms = list(residue.atoms())
        names = [atom.name for atom in atoms]
        error = f'Prepared ligand {residue.name} {residue.chain.id}:{residue.id} lost its recorded atom identity'
        if len(names) != len(by_name) or len(set(names)) != len(names) or set(names) != set(by_name):
            raise ValueError(error + ': atom names or counts changed. Prepare the complex again.')
        matches = [by_name[atom.name] for atom in atoms]
        if any(atom.element != template.atoms[index].element for atom, index in zip(atoms, matches)):
            raise ValueError(error + ': atom elements changed. Prepare the complex again.')
        mapped = {atom.index: index for atom, index in zip(atoms, matches)}
        bonds = set()
        for atom in atoms:
            for neighbor in bonded_to_atom[atom.index]:
                if neighbor not in mapped:
                    raise ValueError(error + ': unexpected external/covalent bond.')
                bonds.add(tuple(sorted((mapped[atom.index], mapped[neighbor]))))
        if bonds != {tuple(sorted(bond)) for bond in template.bonds}:
            raise ValueError(error + ': ligand connectivity changed. Prepare the complex again.')
        return [template, matches]

    # OpenMM's public template matcher can choose a template but cannot return
    # an atom map: it invokes graph matching again. Keep this private hook small
    # and cover it with native-energy, atom-order and Modeller/reload regressions.
    def _getResidueTemplateMatches(self, res, bondedToAtom, templateSignatures=None,
                                   ignoreExternalBonds=False, ignoreExtraParticles=False):
        matched = self._identity_match(res, bondedToAtom)
        if matched is not None:
            return matched
        matched = super()._getResidueTemplateMatches(
            res, bondedToAtom, templateSignatures, ignoreExternalBonds, ignoreExtraParticles)
        if matched[0] is not None and matched[0].name.startswith('DML_'):
            raise ValueError('A DynaMol ligand has no matching recorded atom identity. Prepare the complex again.')
        return matched

    def _matchAllResiduesToTemplates(self, data, topology, residueTemplates, ignoreExternalBonds,
                                     ignoreExtraParticles=False, recordParameters=True):
        # An explicit template selection must not bypass the recorded atom map.
        selected = dict(residueTemplates)
        for residue, name in residueTemplates.items():
            identity = self._ligand_identities.get(self._residue_key(residue))
            if identity is not None:
                if name != identity[0].name:
                    raise ValueError('An explicit template conflicts with the prepared ligand atom identity.')
                del selected[residue]
            elif name.startswith('DML_'):
                raise ValueError('An explicit DynaMol ligand template requires its recorded atom identity.')
        return super()._matchAllResiduesToTemplates(
            data, topology, selected, ignoreExternalBonds, ignoreExtraParticles, recordParameters)
