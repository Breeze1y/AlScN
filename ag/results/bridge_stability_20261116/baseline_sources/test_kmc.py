"""Restoration checks plus the original physics, reaction and exact-cache tests."""
from pathlib import Path
from dataclasses import replace
import hashlib
import json
import sys
import types
import unittest
import zipfile
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
        with self.assertRaisesRegex(ValueError,'archived infinite-source'):
            restore_filament(ROOT/'results/current_kmc_forming_exactbatch_uponly_50Hz_20261112/checkpoint.npz')

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


def load_tests(loader,standard_tests,pattern):
    # Keep one maintained test entry. The source snapshot is an immutable fixture.
    suite=unittest.TestSuite([standard_tests])
    with zipfile.ZipFile(ROOT/'results/forming_refit_20261109/source_snapshot.zip') as archive:
        for filename in ('test_physics.py','test_reactions.py','test_fast_kmc.py'):
            name='_restored_'+filename[:-3]
            module=types.ModuleType(name);module.__file__=str(ROOT/filename);sys.modules[name]=module
            exec(compile(archive.read(filename),str(ROOT/filename),'exec'),module.__dict__)
            suite.addTests(loader.loadTestsFromModule(module))
    return suite


if __name__=='__main__':unittest.main()
