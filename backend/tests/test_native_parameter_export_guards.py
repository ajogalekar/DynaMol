"""Native unsupported-map guards must fail before exporting partial evidence."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[2]/'docs/audit/advanced-chemistry/export_native_parameter_manifest.py'
spec = spec_from_file_location('native_manifest_export_guard_test', PATH)
module = module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize('prefix', ['', 'CHARMM_'])
def test_present_native_cmap_cannot_be_exported_as_ordinary_terms(prefix):
    raw = {prefix+'CMAP_COUNT': [1, 1], prefix+'CMAP_RESOLUTION': [24],
           prefix+'CMAP_INDEX': [1, 2, 3, 4, 5, 1], prefix+'CMAP_PARAMETER_01': [0.0]*576}
    with pytest.raises(ValueError, match='unsupported CMAP'):
        module.reject_unsupported_native_terms(raw)


@pytest.mark.parametrize('prefix', ['', 'CHARMM_'])
def test_explicitly_empty_cmap_count_has_no_missing_force_term(prefix):
    module.reject_unsupported_native_terms({prefix+'CMAP_COUNT': [0, 0], prefix+'CMAP_INDEX': []})


@pytest.mark.parametrize('prefix', ['', 'CHARMM_'])
@pytest.mark.parametrize('count', [None, [0, 0]])
def test_orphan_or_false_zero_map_is_rejected(prefix, count):
    raw = {prefix+'CMAP_PARAMETER_01': [0.0]*576}
    if count is not None:
        raw[prefix+'CMAP_COUNT'] = count
    with pytest.raises(ValueError, match='CMAP'):
        module.reject_unsupported_native_terms(raw)


@pytest.mark.parametrize('count', [[-1, 0], [1.0, 1], [False, 0], [0], '0 0'])
def test_malformed_map_count_is_not_empty_evidence(count):
    with pytest.raises(ValueError, match='Malformed'):
        module.reject_unsupported_native_terms({'CMAP_COUNT': count})
