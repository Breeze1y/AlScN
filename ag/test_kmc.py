"""Restoration checks plus the original physics, reaction and exact-cache tests."""
from pathlib import Path
from dataclasses import replace
import hashlib
import json
import unittest
from tempfile import TemporaryDirectory
import numpy as np

from kmc_model import Config, Model
from run_kmc import restore_filament,read_forming,hashes,save_dynamic,load_dynamic

ROOT=Path(__file__).resolve().parent


class RestorationTests(unittest.TestCase):
    def test_frozen_source_files_restored_byte_for_byte(self):
        expected=json.loads((ROOT/'evidence/T25_forming_V1_20261109/manifest.json').read_text())['source_and_configuration_sha256']
        for name,digest in hashes().items():
            if name in expected:self.assertEqual(digest,expected[name])

    def test_infinite_reservoir_inject_return_and_fixed_boundary(self):
        m=Model(Config(nx=8,layer_cells=4,seed=17,record_events=True))
        self.assertEqual(m.shape,(14,8));self.assertEqual(np.count_nonzero(m.state[m.interior]),0)
        a=m.c.nx+3
        m.pending_event=dict(kind='inject',a=a,b=-1,time_s=.001)
        m.advance(.001,10)
        self.assertEqual(m.net_source,1);self.assertEqual(m.state[a],1)
        self.assertTrue(np.all(m.state[:m.c.nx]==2))
        m.pending_event=dict(kind='absorb',a=a,b=-1,time_s=.002)
        m.advance(.002,10)
        self.assertEqual(m.net_source,0);self.assertEqual(m.state[a],0)
        self.assertEqual(m.result('check')['ag_balance_error'],0)
        self.assertEqual(set(m.counts),{'inject','absorb','hop','reduce','oxidize','migrate'})

    def test_saved_filament_geometry_and_inventory_recovered_exactly(self):
        path=ROOT/'evidence/F03_cycle_initial_states_20261109/forming_seed_17.npz'
        m,audit=restore_filament(path)
        with np.load(path) as archive:np.testing.assert_array_equal(m.state.reshape(m.shape),archive['state'])
        self.assertEqual(m.shape,(62,150));self.assertTrue(m.connected)
        self.assertEqual(int(np.sum(m.state[m.interior]==2)),1704)
        self.assertEqual(int(np.sum(m.state==1)),243);self.assertEqual(m.net_source,1947)
        self.assertTrue(audit['exact_state_match']);self.assertFalse(audit['rng_continuation_available'])

    def test_restored_current_matches_original_frozen_network(self):
        m,audit=restore_filament(ROOT/'evidence/F03_cycle_initial_states_20261109/forming_seed_17.npz')
        self.assertAlmostEqual(m.first_bridge['command_voltage_V'],1.305,places=12)
        self.assertAlmostEqual(m.operating['current_A'],.001,places=12)
        self.assertAlmostEqual(m.operating['device_voltage_V'],.03433688148564609,places=12)

    def test_fresh_default_keeps_50hz_full_loop_and_no_preloaded_ag(self):
        settings=json.loads((ROOT/'configs/current_kmc.json').read_text())
        curve=read_forming(ROOT/settings['protocol']['source'],'full')
        self.assertEqual(settings['protocol']['dwell_s'],.02)
        self.assertEqual(len(curve['points']),202);self.assertEqual(curve['points'][-1]['voltage_V'],0)
        self.assertEqual(curve['points'][100]['voltage_V'],1.5);self.assertEqual(curve['points'][101]['voltage_V'],1.5)
        m=Model(Config(**settings['model']))
        self.assertEqual(np.count_nonzero(m.state[m.interior]),0)
        self.assertTrue(np.all(m.state[:150]==2))

    def test_finite_electrode_stage_cannot_be_imported_as_legacy_shape(self):
        # Rejected finite-electrode schema needs no multi-MB historical run.
        with TemporaryDirectory() as folder:
            path=Path(folder)/'finite_electrode.npz'
            np.savez(path,state=np.zeros((4,4),dtype=np.uint8))
            with self.assertRaisesRegex(ValueError,'archived infinite-source'):
                restore_filament(path)

    def test_reference_morphology_import_keeps_fast_protocol_timing(self):
        path=ROOT/'results/morphology_reproduction_20261112/vacancy_penalty_seed_17_filament.npz'
        m,audit=restore_filament(path)
        self.assertEqual(audit['historical_dwell_s'],.0001)
        self.assertTrue(m.connected)
        self.assertEqual(int(np.sum(m.state[m.interior]==2)),2020)
        self.assertEqual(int(np.sum(m.state==1)),8)
        self.assertAlmostEqual(m.operating['device_voltage_V'],.9736083787755132,places=12)


class DynamicTests(unittest.TestCase):
    def model(self):
        from dynamic_kmc import DynamicCachedModel
        from forming_iv import ElectronicConfig
        cfg=Config(nx=4,layer_cells=1,voltage=1.5,source_barrier_eV=.45,hop_aln_eV=.25,hop_alscn_eV=.25,
                   reaction_barrier_eV=.3,pt_nucleation_eV=.02,max_events=20000,seed=17,disorder_seed=811)
        return DynamicCachedModel(cfg,ElectronicConfig(series_resistance_ohm=0))

    def test_prebridge_path_matches_original_then_continues(self):
        from fast_kmc import CachedFormingModel
        m=self.model();old=CachedFormingModel(replace(m.c),m.electronics)
        old.advance(.001,20000);self.assertTrue(old.connected)
        m.advance(old.time,old.steps)
        np.testing.assert_array_equal(m.state,old.state)
        self.assertEqual(m.counts,old.counts)
        self.assertAlmostEqual(m.first_bridge['time_s'],old.first_bridge['time_s'],places=14)
        first=m.steps;m.advance(.001,20000)
        self.assertEqual(m.time,.001);self.assertGreater(m.steps,first)
        self.assertEqual(m.result('check')['ag_balance_error'],0)

    def test_network_potential_agrees_with_original_conductance(self):
        from dynamic_kmc import dynamic_network
        from forming_iv import electronic_network
        m=self.model();s=m.state.reshape(m.shape).copy();s[:,2]=2
        old=electronic_network(s,m.c.spacing_nm,m.electronics)
        new=dynamic_network(s,m.c.spacing_nm,m.electronics)
        self.assertAlmostEqual(new['conductance_S'],old['conductance_S'],places=12)
        self.assertTrue(np.any((new['metal_phi_at_1V'][1:-1,2]>0)&(new['metal_phi_at_1V'][1:-1,2]<1)))

    def test_bridged_field_uses_finite_metal_voltage_profile(self):
        m=self.model();m.state.reshape(m.shape)[:,2]=2
        m.net_source=int(np.count_nonzero(m.state[m.interior]));m.update_structure()
        metal=m.state==2;expected=m.network['metal_phi_at_1V'].ravel()*m.c.voltage
        np.testing.assert_allclose(m.phi[metal],expected[metal],atol=1e-12)
        self.assertLess(m.field_residual,1e-8)
        inner=metal&m.interior
        self.assertGreater(float(m.phi[inner].max()),0)

    def test_metal_redox_uses_its_own_electron_potential_at_contact(self):
        m=self.model();m.state.reshape(m.shape)[:,1]=2;m.update_structure()
        i=m.c.nx+1
        local=m.electron_reservoir_potential(np.array([i]))[0]
        self.assertAlmostEqual(local,m.phi[i],places=12)
        donor=m.valid[i]&(m.state[m.safe[i]]==2)&m.cathode[m.safe[i]]
        arithmetic=float(np.mean(m.phi[m.safe[i][donor]]))
        self.assertGreater(abs(arithmetic-local),.01)
        m._build_fixed();sites=m._catalogue_fixed[3];ox=m._catalogue_fixed[5]
        k=int(np.flatnonzero(sites==i)[0])
        dg=m.c.metal_standard_eV-m.binding(np.array([i]))[0]-m.site[i]
        expected=m.rate(np.array([m.c.reaction_barrier_eV+max(-dg,0)]))[0]
        self.assertAlmostEqual(ox[k],expected,places=10)

    def test_virtual_reduction_node_weights_contact_resistance(self):
        m=self.model();m.state.reshape(m.shape)[:,1]=2;m.state.reshape(m.shape)[1,2]=1;m.update_structure()
        i=m.c.nx+2
        neighbors=m.safe[i];donor=m.valid[i]&(m.state[neighbors]==2)&m.cathode[neighbors]
        phi=m.network['metal_phi_at_1V'].ravel()*m.c.voltage
        terminal=(neighbors<m.c.nx)|(neighbors>=m.n-m.c.nx)
        resistance=m.electronics.metal_resistivity_ohm_m/(m.electronics.effective_depth_nm*1e-9)+terminal*m.electronics.contact_resistance_ohm
        w=1/resistance[donor];expected=float(np.sum(w*phi[neighbors[donor]])/w.sum())
        actual=m.electron_reservoir_potential(np.array([i]))[0]
        self.assertAlmostEqual(actual,expected,places=12)
        self.assertGreater(abs(actual-float(phi[neighbors[donor]].mean())),.01)

    def test_local_redox_forward_reverse_ratio_survives_contact_correction(self):
        from kmc_model import KB
        m=self.model();m.state.reshape(m.shape)[:,1]=2;m.state.reshape(m.shape)[1,2]=1;m.update_structure();m._build_fixed()
        sites=m._catalogue_fixed[3];red=m._catalogue_fixed[4];ox=m._catalogue_fixed[5]
        dg=m.c.metal_standard_eV-m.binding(sites)-m.site[sites]-m.phi[sites]+m.electron_reservoir_potential(sites)
        np.testing.assert_allclose(np.log(red/ox),-dg/(KB*m.c.temperature_K),rtol=1e-12,atol=1e-12)

    def test_return_to_zero_continues_time_and_keeps_inventory(self):
        m=self.model();m.advance(.0001,20000);self.assertIsNotNone(m.first_bridge)
        old=m.steps;m.apply_command(0);m.advance(.0003,20000)
        self.assertEqual(m.time,.0003);self.assertGreater(m.steps,old)
        self.assertEqual(m.result('check')['ag_balance_error'],0)
        self.assertEqual(m.operating['current_A'],0)

    def test_dynamic_checkpoint_preserves_random_path(self):
        m=self.model();m.advance(.0002,20000)
        with TemporaryDirectory() as folder:
            folder=Path(folder)
            save_dynamic(folder,m,dict(phase='waveform',index=0,deadline_s=.0004),[],[m.state.reshape(m.shape).copy()],['test'],'test',{})
            n,_,_,_,_=load_dynamic(folder)
        m.advance(.0004,20000);n.advance(.0004,20000)
        np.testing.assert_array_equal(m.state,n.state);self.assertEqual(m.counts,n.counts)
        self.assertEqual(m.pending_event,n.pending_event)
        self.assertEqual(m.topology_log,n.topology_log)

    def test_stability_requires_observed_hold_and_no_hidden_breaks(self):
        from dynamic_kmc import stability
        rows=[dict(time_s=i*.02,network_conductance_S=.01,metal=1000,ions=10,connected=True) for i in range(6)]
        self.assertFalse(stability(rows[:3],[])['assessed'])
        self.assertTrue(stability(rows,[])['filament_stable'])
        self.assertFalse(stability(rows,[dict(time_s=.05,connected=False)])['filament_stable'])
        rows[-1]['network_conductance_S']=.008
        self.assertFalse(stability(rows,[])['filament_stable'])

    def test_ion_stationarity_is_separate_from_metal_stability(self):
        from dynamic_kmc import stability
        rows=[dict(time_s=i*.02,network_conductance_S=.01,metal=1000,ions=100-i*10,connected=True) for i in range(6)]
        out=stability(rows,[])
        self.assertTrue(out['filament_stable']);self.assertFalse(out['ions_stationary'])


# Maintained invariants recovered from the 2026-11-09 source snapshot: test_physics.py
"""Physical invariants and solver checks, runnable with unittest."""
import unittest
from dataclasses import replace
import numpy as np
from kmc_model import Config, Model, EPS0, KB


class PhysicsTests(unittest.TestCase):
    def test_vacancy_mask_reproducibility_and_zero_penalty_control(self):
        cfg = Config(nx=8, layer_cells=4, voltage=1.5, seed=18, max_events=500,
                     max_time_s=.05, record_events=True)
        annotated = replace(cfg, vacancy_affected_fraction_by_layer=(0, .5, 0), vacancy_seed=42)
        baseline, control = Model(cfg), Model(annotated)
        other_seed = Model(replace(annotated, seed=19, vacancy_ag_ion_penalty_eV=.05))
        np.testing.assert_array_equal(control.vacancy_mask, other_seed.vacancy_mask)
        mask = control.vacancy_mask
        self.assertTrue(mask.any())
        self.assertTrue(np.all(control.interior[mask] & (control.zone[mask] == 1)))
        np.testing.assert_array_equal(baseline.site, control.site)
        baseline.run()
        control.run()
        self.assertGreater(len(control.event_log), 0)
        self.assertEqual(baseline.event_log, control.event_log)
        self.assertEqual(baseline.rng.bit_generator.state, control.rng.bit_generator.state)
        np.testing.assert_array_equal(baseline.state, control.state)

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


# Maintained invariants recovered from the 2026-11-09 source snapshot: test_reactions.py
"""Surface admissibility and contact-energy invariants of the local model."""
from dataclasses import replace
import unittest
import numpy as np
from kmc_model import Config,Model,KB


def rate(model,kind,a,b=-1):
    return sum(float(r) for k,aa,bb,rr in model.catalogue() if k==kind
               for x,y,r in zip(aa,bb,rr) if x==a and y==b)


class ReactionTests(unittest.TestCase):
    def test_vacancy_penalty_enters_all_ionic_forward_reverse_pairs(self):
        penalty = .05
        cfg = Config(nx=4, layer_cells=2, voltage=0, polarization_C_m2=0, disorder_eV=0,
                     vacancy_affected_fraction_by_layer=(1, 0, 1), vacancy_ag_ion_penalty_eV=penalty)
        model = Model(cfg)
        kT = KB * cfg.temperature_K
        source = cfg.nx + 1
        injection = rate(model, 'inject', source)
        model.state[source] = 1
        self.assertAlmostEqual(np.log(injection / rate(model, 'absorb', source)), -penalty / kT)
        model.state[source] = 0

        marked, unmarked = 2 * cfg.nx + 1, 3 * cfg.nx + 1
        model.state[marked] = 1
        out_of_defect = rate(model, 'hop', marked, unmarked)
        model.state[marked], model.state[unmarked] = 0, 1
        self.assertAlmostEqual(np.log(out_of_defect / rate(model, 'hop', unmarked, marked)), penalty / kT)
        model.state[unmarked] = 0

        contact = model.n - 2 * cfg.nx + 1
        model.state[contact] = 1
        reduction = rate(model, 'reduce', contact)
        model.state[contact] = 2
        model.update_structure()
        delta_g = cfg.metal_standard_eV - cfg.pt_adhesion_eV - penalty
        self.assertAlmostEqual(np.log(reduction / rate(model, 'oxidize', contact)), -delta_g / kT)

    def test_buried_redox_forbidden_in_both_directions(self):
        m=Model(Config(nx=5,layer_cells=2,voltage=0,polarization_C_m2=0))
        m.state.reshape(m.shape)[3:-1,1:4]=2;m.update_structure()
        a=5*m.c.nx+2
        self.assertEqual(rate(m,'oxidize',a),0)
        m.state[a]=1;m.update_structure()
        self.assertEqual(rate(m,'reduce',a),0)
        m.c=replace(m.c,surface_redox_only=False)
        self.assertGreater(rate(m,'reduce',a),0)

    def test_exposed_redox_pair_with_distinct_adhesion(self):
        m=Model(Config(nx=5,layer_cells=2,voltage=.3,polarization_C_m2=0,pt_adhesion_eV=.09))
        a=m.n-2*m.c.nx+2;m.state[a]=1
        phi=m.phi.copy()
        dg=m.c.metal_standard_eV-m.binding(np.array([a]))[0]-m.site[a]-phi[a]
        f=rate(m,'reduce',a)
        m.state[a]=2;m.update_structure();m.phi=phi
        r=rate(m,'oxidize',a)
        self.assertGreater(f,0);self.assertGreater(r,0)
        self.assertAlmostEqual(np.log(f/r),-dg/(KB*m.c.temperature_K),places=10)

    def test_migration_preserves_contact_energy_balance(self):
        m=Model(Config(nx=5,layer_cells=2,voltage=0,polarization_C_m2=0,pt_adhesion_eV=.10))
        a=m.n-2*m.c.nx+2;b=a-m.c.nx
        m.state[a]=2;m.update_structure();phi=m.phi.copy()
        dg=float(m.binding(np.array([a]))[0]-(m.binding(np.array([b]))[0]-m.c.bond_eV))
        f=rate(m,'migrate',a,b)
        m.state[a]=0;m.state[b]=2;m.update_structure();m.phi=phi
        r=rate(m,'migrate',b,a)
        self.assertAlmostEqual(np.log(f/r),-dg/(KB*m.c.temperature_K),places=10)


# Maintained invariants recovered from the 2026-11-09 source snapshot: test_fast_kmc.py
"""Caching must preserve full random trajectories, not merely means."""
from dataclasses import replace
import unittest
import numpy as np
from kmc_model import Config, Model
from forming_iv import FormingModel, ElectronicConfig
from fast_kmc import CachedModel,CachedFormingModel


class CacheTests(unittest.TestCase):
    def compare(self,a,b):
        np.testing.assert_array_equal(a.state,b.state)
        np.testing.assert_array_equal(a.phi,b.phi)
        self.assertEqual(a.time,b.time)
        self.assertEqual(a.event_log,b.event_log)
        self.assertEqual(a.counts,b.counts)
        self.assertEqual(a.pending_event,b.pending_event)
        self.assertEqual(a.rng.bit_generator.state,b.rng.bit_generator.state)

    def test_changed_bias_and_moving_ions_keep_exact_path(self):
        cfg=Config(nx=8,layer_cells=4,max_events=4000,max_time_s=.1,seed=18,record_events=True,
                   vacancy_affected_fraction_by_layer=(.1,.5,.1),vacancy_ag_ion_penalty_eV=.05)
        a,b=Model(cfg),CachedModel(cfg)
        for voltage,time in [(0,.005),(.8,.02),(1.5,.04),(.3,.08)]:
            for m in [a,b]:m.set_voltage(voltage);m.advance(time)
            self.compare(a,b)
            np.testing.assert_array_equal(a.vacancy_ion_residence_s,b.vacancy_ion_residence_s)

    def test_circuit_feedback_and_neutral_migration_keep_exact_path(self):
        cfg=Config(nx=5,layer_cells=2,max_events=2000,max_time_s=.05,seed=7,
                   record_events=True,migration_eV=.5,voltage=.8)
        e=ElectronicConfig(leakage_per_column_S=1e-4,compliance_A=1e-5)
        a,b=FormingModel(cfg,e),CachedFormingModel(cfg,e)
        for voltage,time in [(.8,.01),(1.2,.03),(2.,.05)]:
            for m in [a,b]:m.apply_command(voltage);m.advance(time)
            self.compare(a,b)
            self.assertEqual(a.iv,b.iv)

    def test_all_channels_match_on_mixed_states_and_after_ion_motion(self):
        cfg=Config(nx=8,layer_cells=4,voltage=1.2,migration_eV=.5)
        a,b=Model(cfg),CachedModel(cfg)
        rng=np.random.default_rng(912)
        state=rng.choice([0,1,2],size=(12,8),p=[.65,.2,.15]).astype(np.uint8)
        for m in [a,b]:m.state.reshape(m.shape)[1:-1]=state;m.update_structure()
        def same_catalogue():
            aa,bb=a.catalogue(),b.catalogue()
            self.assertEqual([r[0] for r in aa],[r[0] for r in bb])
            self.assertIn('migrate',[r[0] for r in aa])
            for ga,gb in zip(aa,bb):
                for x,y in zip(ga[1:],gb[1:]):np.testing.assert_array_equal(x,y)
        same_catalogue()
        ion=int(np.flatnonzero(a.state==1)[0]);host=int(np.flatnonzero(a.interior&(a.state==0))[0])
        for m in [a,b]:m.state[ion]=0;m.state[host]=1
        same_catalogue()
        for m in [a,b]:m.c=replace(m.c,temperature_K=310)
        same_catalogue()

if __name__=='__main__':unittest.main()
