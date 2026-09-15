"""Regression checks for rounded native charges without relaxing state identity."""
import tempfile
import unittest
from pathlib import Path

from prepare_adduct import validate_native_charge_total


class ChargePrecisionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.native=self.root/'capped-native';self.native.mkdir()
        self.sqm=('Mulliken Charge\n1 C -0.242\n2 H 0.241\n'
                  'Total Mulliken Charge = 0.000\nCalculation Completed\n')
        (self.native/'sqm.out').write_text(self.sqm)
        for filename in ['ANTECHAMBER_AM1BCC_PRE.AC','ANTECHAMBER_AM1BCC.AC']:
            (self.native/filename).write_text('ATOM 1 C1 MOL 1 0 0 0 -0.300000 c3\n'
                                             'ATOM 2 H1 MOL 1 1 0 0 0.299000 h1\n')
        (self.native/'charged.mol2').write_text('@<TRIPOS>ATOM\n'
            '1 C1 0 0 0 c3 1 MOL -0.300000\n2 H1 1 0 0 h1 1 MOL 0.299000\n')

    def tearDown(self):
        self.temp.cleanup()

    def test_evidenced_rounding_admitted_without_changing_charges(self):
        charges=[-.3,.299]
        result=validate_native_charge_total(self.root,charges,0)
        self.assertEqual(charges,[-.3,.299])
        self.assertFalse(result['charges_modified_by_admission'])
        self.assertAlmostEqual(result['input_error_e'],-.001)

    def test_wrong_declared_state_rejected(self):
        with self.assertRaises(ValueError):validate_native_charge_total(self.root,[-.3,.299],1)

    def test_redistributed_candidate_cannot_hide_behind_same_total(self):
        with self.assertRaises(ValueError):validate_native_charge_total(self.root,[-.2,.199],0)

    def test_incomplete_quantum_output_rejected(self):
        (self.native/'sqm.out').write_text(self.sqm.replace('Calculation Completed','Interrupted'))
        with self.assertRaises(ValueError):validate_native_charge_total(self.root,[-.3,.299],0)

    def test_bcc_charge_loss_rejected(self):
        p=self.native/'ANTECHAMBER_AM1BCC.AC'
        p.write_text(p.read_text().replace('0.299000','0.289000'))
        with self.assertRaises(ValueError):validate_native_charge_total(self.root,[-.3,.299],0)

    def test_nonfinite_charge_rejected(self):
        with self.assertRaises(ValueError):validate_native_charge_total(self.root,[float('nan'),0],0)


if __name__=='__main__':unittest.main()
