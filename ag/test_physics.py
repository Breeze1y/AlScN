"""Physical invariants and solver checks, runnable with unittest."""
import unittest
from dataclasses import replace
import numpy as np
from kmc_model import Config, Model, EPS0, KB


class PhysicsTests(unittest.TestCase):
    def test_layered_polarization_against_analytic(self):
        c = Config(nx=5, disorder_eV=0, voltage=1.2)
        m = Model(c)
        h = c.spacing_nm*1e-9
        thickness = c.layer_cells*h
        eps = EPS0*np.array([c.eps_aln, c.eps_alscn, c.eps_aln])
        pol = np.array([0, c.polarization_C_m2, 0])
        d = (c.voltage + np.sum(pol*thickness/eps))/np.sum(thickness/eps)
        e = (d-pol)/eps
        expected = []
        for k in range(3):
            expected.extend(c.voltage-np.sum(e[:k]*thickness)-e[k]*(np.arange(c.layer_cells)+.5)*h)
        # High-epsilon metal rows approximate ideal surfaces to O(1e-5).
        np.testing.assert_allclose(m.phi.reshape(m.shape)[1:-1, 2], expected, atol=2e-5)
        self.assertLess(m.field_residual, 1e-10)

    def test_short_circuit_compensation_and_reversal(self):
        c = Config(nx=4, voltage=0, disorder_eV=0)
        a = Model(c); b = Model(replace(c, polarization_C_m2=-c.polarization_C_m2))
        np.testing.assert_allclose(a.phi, -b.phi, atol=1e-12)
        phi = a.phi.reshape(a.shape)[:, 1]
        self.assertLess(phi[1], phi[10])  # Outer and middle fields oppose.
        self.assertGreater(phi[25], phi[35])

    def test_hop_detailed_balance(self):
        m = Model(Config(nx=4))
        for dg in [-.2, -.03, 0, .07, .2]:
            forward = m.rate(.5+max(dg, 0))
            reverse = m.rate(.5+max(-dg, 0))
            self.assertAlmostEqual(float(np.log(forward/reverse)), -dg/(KB*m.c.temperature_K), places=11)

    def test_mass_and_no_preloaded_filament(self):
        m = Model(Config(nx=8, layer_cells=3, max_events=1500, max_time_s=.1))
        self.assertFalse(m.connected)
        self.assertEqual(np.sum(m.state[m.interior]), 0)
        result = m.run()
        self.assertEqual(result['ag_balance_error'], 0)
        self.assertEqual(result['counts']['inject']-result['counts']['absorb'], m.net_source)

    def test_topology_does_not_confuse_ions_with_metal(self):
        m = Model(Config(nx=5, layer_cells=2))
        m.state.reshape(m.shape)[1:-1, 2] = 1
        m.update_structure(); self.assertFalse(m.connected)
        m.state.reshape(m.shape)[1:-1, 2] = 2
        m.update_structure(); self.assertTrue(m.connected)

if __name__ == '__main__': unittest.main()
