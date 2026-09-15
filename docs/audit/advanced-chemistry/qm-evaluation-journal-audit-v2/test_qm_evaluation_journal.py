"""Synthetic callback bookkeeping only; no QM, molecular model or native imports.

Run with any Python 3 stdlib:
    python test_qm_evaluation_journal.py -v
"""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "qm_evaluation_journal", Path(__file__).with_name("qm_evaluation_journal.py"))
journal = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(journal)


class MockArray:
    """Exercise the ndarray-like callback interface without importing numpy."""
    def __init__(self, value):
        self.value = value

    def tolist(self):
        return copy.deepcopy(self.value)


class JournalTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.request = {
            "schema_version": 1,
            "atom_ids": ["SYNTHETIC:O", "SYNTHETIC:H1", "SYNTHETIC:H2"],
            "elements": ["O", "H", "H"], "charge": 0, "spin": 0,
            "coords_bohr": [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            "method": "RHF", "basis": "MOCK-BASIS-NOT-EVALUATED",
            "density_fit": True, "auxbasis": "MOCK-AUX-NOT-EVALUATED",
            "scf_conv_tol": 1e-10, "scf_conv_tol_grad": 1e-7,
            "optimization": {"enabled": True, "freeze_atom_ids": ["SYNTHETIC:O"],
                             "convergence_energy": 1e-6},
        }
        self.binding = journal.make_binding(
            self.request, source_input_sha256="a" * 64, worker_source_sha256="b" * 64)
        self.path = self.root / "synthetic-journal"
        self.writer = journal.EvaluationJournal(self.path, self.binding, create=True)
        self.args = {
            "cycle": 1, "coords": MockArray([[0, 0, 0], [1.1, 0, 0], [0, 1.2, 0]]),
            "gradient": MockArray([[.1, -.2, .3], [-.1, .2, -.3], [0, 0, 0]]),
            "energy": -1.234567890123, "units": dict(journal.UNITS),
            "observed_identity": copy.deepcopy(self.binding["identity"]),
            "requested_method_fingerprint_sha256": self.binding["requested_method_fingerprint_sha256"],
            "scf_converged": True,
        }

    def append(self, **changes):
        return self.writer.append(**{**self.args, **changes})

    def rewrite_record(self, change, *, rehash=True):
        path = self.path / "evaluation-000001.json"
        envelope = json.loads(path.read_bytes())
        change(envelope["payload"])
        if rehash:
            envelope["payload_sha256"] = journal.fingerprint(envelope["payload"])
        path.write_bytes(journal._canonical(envelope) + b"\n")

    def test_full_same_evaluation_is_preserved_with_positive_gradient_and_hashes(self):
        publication = self.append()
        raw = Path(publication["path"]).read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), publication["file_sha256"])
        [record] = self.writer.records()
        self.assertEqual(record["coords_bohr"], self.args["coords"].tolist())
        self.assertEqual(record["gradient_hartree_per_bohr"], self.args["gradient"].tolist())
        self.assertEqual(record["energy_hartree"], self.args["energy"])
        self.assertEqual(record["identity"], self.binding["identity"])
        self.assertEqual(record["units"], journal.UNITS)
        self.assertEqual(record["manifest_sha256"], journal.fingerprint(self.binding))
        for key, value in journal.STATUS.items():
            self.assertEqual(record[key], value)
        self.assertEqual(record["status"], "UNCONVERGED")
        self.assertFalse(record["checkpoint_reuse_authorized"])

    def test_callback_flat_coordinates_are_reshaped_without_conversion(self):
        self.append(coords=MockArray([0, 0, 0, 1.1, 0, 0, 0, 1.2, 0]))
        self.assertEqual(self.writer.records()[0]["coords_bohr"], self.args["coords"].tolist())

    def test_every_noncoordinate_requested_setting_changes_method_fingerprint(self):
        for name, value in {
            "basis": "DIFFERENT", "auxbasis": "DIFFERENT", "spin": 2,
            "functional": "NEW", "unknown_future_setting": {"value": 1},
            "optimization": {"enabled": True, "convergence_energy": 1e-5},
        }.items():
            with self.subTest(name=name):
                request = {**self.request, name: value}
                changed = journal.make_binding(request, source_input_sha256="a" * 64,
                                               worker_source_sha256="b" * 64)
                self.assertNotEqual(changed["requested_method_fingerprint_sha256"],
                                    self.binding["requested_method_fingerprint_sha256"])

    def test_geometry_changes_resolved_request_digest_but_not_method_digest(self):
        request = copy.deepcopy(self.request)
        request["coords_bohr"][1][0] += .1
        changed = journal.make_binding(request, source_input_sha256="a" * 64,
                                       worker_source_sha256="b" * 64)
        self.assertNotEqual(changed["resolved_request_sha256"], self.binding["resolved_request_sha256"])
        self.assertEqual(changed["requested_method_fingerprint_sha256"],
                         self.binding["requested_method_fingerprint_sha256"])

    def test_raw_source_and_worker_hashes_are_separate_binding_fields(self):
        self.assertEqual(self.binding["source_input_file_sha256"], "a" * 64)
        self.assertEqual(self.binding["worker_source_file_sha256"], "b" * 64)
        with self.assertRaises(ValueError):
            journal.make_binding(self.request, source_input_sha256="bad", worker_source_sha256="b" * 64)

    def test_changed_binding_cannot_reopen_a_run(self):
        for key in ("source_input_file_sha256", "worker_source_file_sha256", "resolved_request_sha256"):
            with self.subTest(key=key):
                binding = {**self.binding, key: "c" * 64}
                with self.assertRaises(ValueError):
                    journal.EvaluationJournal(self.path, binding)

    def test_binding_and_arrays_are_snapshotted_not_retained_mutable_references(self):
        self.append()
        self.binding["identity"]["atom_ids"][0] = "MUTATED"
        self.args["coords"].value[0][0] = 99
        record = self.writer.records()[0]
        self.assertEqual(record["identity"]["atom_ids"][0], "SYNTHETIC:O")
        self.assertEqual(record["coords_bohr"][0][0], 0)

    def test_mismatched_identity_order_elements_or_state_refused_before_writing(self):
        for key, value in {
            "atom_ids": ["SYNTHETIC:O", "SYNTHETIC:H2", "SYNTHETIC:H1"],
            "elements": ["H", "O", "H"], "charge": 1, "spin": 2,
        }.items():
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.append(observed_identity={**self.args["observed_identity"], key: value})
        self.assertEqual(self.writer.records(), [])

    def test_duplicate_and_missing_identity_atoms_refused(self):
        for ids in (["O", "H", "H"], ["SYNTHETIC:O"], []):
            with self.subTest(ids=ids), self.assertRaises(ValueError):
                self.append(observed_identity={**self.args["observed_identity"], "atom_ids": ids})

    def test_method_mismatch_refused(self):
        with self.assertRaisesRegex(ValueError, "method fingerprint"):
            self.append(requested_method_fingerprint_sha256="c" * 64)

    def test_scf_flag_must_be_true_not_truthy(self):
        for value in (False, None, 0, 1, "true"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.append(scf_converged=value)
        self.assertEqual(self.writer.records(), [])

    def test_wrong_units_missing_units_and_force_units_refused(self):
        for units in ({}, None, {**journal.UNITS, "coordinates": "angstrom"},
                      {**journal.UNITS, "gradient": "force_hartree/bohr"},
                      {**journal.UNITS, "energy": "kcal/mol"},
                      {**journal.UNITS, "additional": "ignored"}):
            with self.subTest(units=units), self.assertRaises(ValueError):
                self.append(units=units)
        self.assertEqual(self.writer.records(), [])

    def test_nonfinite_boolean_string_and_complex_energy_refused(self):
        for value in (float("nan"), float("inf"), float("-inf"), True, "1.0", 1j):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.append(energy=value)

    def test_nonfinite_and_invalid_coordinate_or_gradient_components_refused(self):
        for field in ("coords", "gradient"):
            for value in (float("nan"), float("inf"), True, "0", 1j):
                with self.subTest(field=field, value=value):
                    array = self.args[field].tolist()
                    array[1][2] = value
                    with self.assertRaises(ValueError):
                        self.append(**{field: array})

    def test_malformed_shapes_and_flat_gradient_refused(self):
        for field in ("coords", "gradient"):
            for array in ([], [[0, 0, 0]], [[0, 0]] * 3, [[0, 0, 0, 0]] * 3, None):
                with self.subTest(field=field, array=array), self.assertRaises(ValueError):
                    self.append(**{field: array})
        with self.assertRaises(ValueError):
            self.append(gradient=[0] * 9)

    def test_cycle_bounds_refused(self):
        for value in (0, -1, True, 1.0, "1", 1000000):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.append(cycle=value)

    def test_create_refuses_existing_directory_and_cycle_commit_never_overwrites(self):
        with self.assertRaises(FileExistsError):
            journal.EvaluationJournal(self.path, self.binding, create=True)
        publication = self.append()
        before = Path(publication["path"]).read_bytes()
        with self.assertRaises(FileExistsError):
            self.append(energy=999)
        self.assertEqual(Path(publication["path"]).read_bytes(), before)
        self.assertEqual(len(self.writer.records()), 1)
        self.assertEqual(len(list(self.path.glob(".evaluation-*.tmp"))), 1)

    def test_partial_write_interruption_preserves_evidence_and_existing_records(self):
        self.append()

        def partial_then_interrupt(fd, raw):
            os.write(fd, raw[:31])
            raise KeyboardInterrupt("synthetic interrupted callback write")

        with patch.object(journal, "_write_all", partial_then_interrupt):
            with self.assertRaises(KeyboardInterrupt):
                self.append(cycle=2)
        [record] = journal.EvaluationJournal(self.path, self.binding).records()
        self.assertEqual(record["cycle"], 1)
        [partial] = list(self.path.glob(".evaluation-000002.*.tmp"))
        self.assertEqual(partial.stat().st_size, 31)
        self.assertFalse((self.path / "evaluation-000002.json").exists())

    def test_file_fsync_failure_leaves_no_committed_record(self):
        with patch.object(journal.os, "fsync", side_effect=OSError("synthetic fsync failure")):
            with self.assertRaises(OSError):
                self.append()
        self.assertEqual(self.writer.records(), [])
        self.assertEqual(len(list(self.path.glob(".evaluation-*.tmp"))), 1)

    def test_interrupt_before_atomic_link_leaves_no_committed_record(self):
        with patch.object(journal.os, "link", side_effect=KeyboardInterrupt("before commit")):
            with self.assertRaises(KeyboardInterrupt):
                self.append()
        self.assertEqual(self.writer.records(), [])
        self.assertEqual(len(list(self.path.glob(".evaluation-*.tmp"))), 1)

    def test_failure_after_commit_recovers_complete_record_and_preserves_temporary(self):
        for error in (OSError("directory sync failed"), KeyboardInterrupt("after commit")):
            with self.subTest(error=type(error).__name__):
                cycle = len(self.writer.records()) + 1
                with patch.object(journal, "_sync_directory", side_effect=error):
                    with self.assertRaises(type(error)):
                        self.append(cycle=cycle)
                self.assertEqual(self.writer.records()[-1]["cycle"], cycle)
                self.assertTrue(list(self.path.glob(f".evaluation-{cycle:06d}.*.tmp")))

    def test_link_commit_is_only_after_file_fsync(self):
        original_sync, original_link = journal.os.fsync, journal.os.link
        order = []

        def sync(fd):
            order.append("sync")
            return original_sync(fd)

        def link(source, target):
            order.append("link")
            return original_link(source, target)

        with patch.object(journal.os, "fsync", sync), patch.object(journal.os, "link", link):
            self.append()
        self.assertEqual(order, ["sync", "link", "sync"])

    def test_short_writes_are_completed(self):
        original_write = journal.os.write
        with patch.object(journal.os, "write", lambda fd, raw: original_write(fd, raw[:17])):
            self.append()
        self.assertEqual(len(self.writer.records()), 1)

    def test_zero_progress_write_refused(self):
        with patch.object(journal.os, "write", return_value=0):
            with self.assertRaises(OSError):
                self.append()
        self.assertEqual(self.writer.records(), [])

    def test_no_latest_pointer_and_uncommitted_complete_temporary_is_ignored(self):
        self.append(cycle=3)
        self.append(cycle=1)
        self.assertEqual([x["cycle"] for x in self.writer.records()], [1, 3])
        self.assertFalse((self.path / "latest.json").exists())
        (self.path / ".evaluation-000099.valid.tmp").write_bytes(
            (self.path / "evaluation-000003.json").read_bytes())
        self.assertEqual([x["cycle"] for x in self.writer.records()], [1, 3])

    def test_tampered_payload_or_manifest_refused(self):
        self.append()
        self.rewrite_record(lambda r: r.update(energy_hartree=1), rehash=False)
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.writer.records()
        manifest = self.path / "manifest.json"
        envelope = json.loads(manifest.read_bytes())
        envelope["payload"]["source_input_file_sha256"] = "d" * 64
        envelope["payload_sha256"] = journal.fingerprint(envelope["payload"])
        manifest.write_bytes(journal._canonical(envelope))
        with self.assertRaisesRegex(ValueError, "manifest binding"):
            self.append(cycle=2)

    def test_rehashed_wrong_status_schema_units_or_binding_refused(self):
        self.append()
        original = (self.path / "evaluation-000001.json").read_bytes()
        for changes in ({"accepted": True}, {"accepted": 0}, {"status": "CONVERGED"},
                        {"physical_acceptance": True}, {"checkpoint_reuse_authorized": True},
                        {"scf_converged": False}, {"manifest_sha256": "d" * 64},
                        {"requested_method_fingerprint_sha256": "e" * 64},
                        {"schema_version": 99}, {"cycle": 2}, {"extra_field": 1},
                        {"units": {**journal.UNITS, "coordinates": "angstrom"}}):
            with self.subTest(changes=changes):
                (self.path / "evaluation-000001.json").write_bytes(original)
                self.rewrite_record(lambda r: r.update(changes))
                with self.assertRaises(ValueError):
                    self.writer.records()

    def test_rehashed_changed_geometry_or_gradient_refused_by_component_hash(self):
        self.append()
        for field in ("coords_bohr", "gradient_hartree_per_bohr"):
            with self.subTest(field=field):
                original = (self.path / "evaluation-000001.json").read_bytes()
                self.rewrite_record(lambda r: r[field][0].__setitem__(0, 9.99))
                with self.assertRaisesRegex(ValueError, "geometry or gradient hash"):
                    self.writer.records()
                (self.path / "evaluation-000001.json").write_bytes(original)

    def _assert_rehashed_identity_type_refused(self, key, value):
        self.append()

        def change(record):
            record["identity"][key] = value
            record["geometry_sha256"] = journal.fingerprint({
                "identity": record["identity"], "unit": "bohr", "coords": record["coords_bohr"]})

        # Both geometry and envelope hashes are consistent with the malformed
        # state: this must fail the strict identity check, not a hash comparison.
        self.rewrite_record(change)
        with self.assertRaisesRegex(ValueError, "explicit integers"):
            self.writer.records()

    def test_rehashed_boolean_charge_is_refused_even_when_equal_to_integer_zero(self):
        self._assert_rehashed_identity_type_refused("charge", False)

    def test_rehashed_float_spin_is_refused_even_when_equal_to_integer_zero(self):
        self._assert_rehashed_identity_type_refused("spin", 0.0)

    def test_rehashed_boolean_record_schema_is_refused_even_when_equal_to_integer_one(self):
        self.append()
        self.rewrite_record(lambda r: r.update(schema_version=True))
        with self.assertRaisesRegex(ValueError, "Evaluation binding"):
            self.writer.records()

    def test_boolean_manifest_schema_is_refused_even_with_matching_expected_binding(self):
        binding = {**self.binding, "schema_version": True}
        with self.assertRaisesRegex(ValueError, "Invalid journal binding"):
            journal.EvaluationJournal(self.root / "invalid-schema-journal", binding, create=True)

    def test_corrupt_final_record_is_not_silently_skipped(self):
        self.append()
        (self.path / "evaluation-000002.json").write_bytes(b"{partial")
        with self.assertRaises(ValueError):
            self.writer.records()

    def test_regular_file_and_filename_binding_required(self):
        self.append()
        path = self.path / "evaluation-000001.json"
        moved = self.root / "moved.json"
        path.rename(moved)
        path.symlink_to(moved)
        with self.assertRaisesRegex(ValueError, "nonsymlink"):
            self.writer.records()
        path.unlink()
        moved.rename(self.path / "evaluation-bad.json")
        with self.assertRaisesRegex(ValueError, "Malformed"):
            self.writer.records()

    def test_fresh_reopen_verifies_records_without_any_checkpoint(self):
        self.append()
        (self.path / "scf.chk").write_bytes(b"synthetic mixed-geometry checkpoint never read")
        reopened = journal.EvaluationJournal(self.path, self.binding)
        [record] = reopened.records()
        self.assertEqual(record["scf_checkpoint_same_geometry_validity"], "NOT_ESTABLISHED")
        self.assertFalse(record["accepted"])
        self.assertFalse(record["optimization_converged"])


if __name__ == "__main__":
    unittest.main()
