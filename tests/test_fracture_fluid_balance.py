import copy
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest
spec=importlib.util.spec_from_file_location('balance',Path(__file__).resolve().parents[1]/'scripts/fracture_fluid_balance.py')
balance=importlib.util.module_from_spec(spec);spec.loader.exec_module(balance)
def rows():
    result=[]
    for step in range(2):
        for phase in balance.PHASES:
            r={k:0.0 for k in balance.REQUIRED-balance.TEXT}
            r.update(schema=1,step=step,rank=0,ranks=1,phase=phase,control='rate_well',start_s=step*10,end_s=(step+1)*10,dt_s=10,fracture_dt_s=10,ncells=2,applicable=1,storage_enabled=1,inner_converged=1,volume_begin_m3=step+1,volume_old_on_current_m3=step+1,volume_end_m3=step+2,well_inflow_m3_s=.3,leakoff_m3_s=.2,storage_model_m3_s=.1,storage_physical_m3_s=.1,flow_tolerance_m3_s=1e-6)
            result.append(r)
    return result
class BalanceTests(unittest.TestCase):
    def test_closed_budget(self):
        report=balance.audit(rows())['accepted_end']
        self.assertAlmostEqual(report['inflow_m3'],6)
        self.assertAlmostEqual(report['leakoff_m3'],4)
        self.assertAlmostEqual(report['cumulative_defect_m3'],0)
    def test_duplicate_owner_is_not_double_counted(self):
        data=rows();data.append(copy.deepcopy(data[0]))
        with self.assertRaisesRegex(ValueError,'Duplicate'):balance.audit(data)
    def test_missing_phase(self):
        with self.assertRaisesRegex(ValueError,'Missing'):balance.audit(rows()[1:])
    def test_checkpoint_jump(self):
        data=rows();data[2]['volume_begin_m3']+=.1;data[3]['volume_begin_m3']+=.1
        with self.assertRaisesRegex(ValueError,'Checkpoint'):balance.audit(data)
    def test_time_gap(self):
        data=rows();data[2]['start_s']+=1;data[3]['start_s']+=1
        with self.assertRaisesRegex(ValueError,'Gap'):balance.audit(data)
    def test_source_residual_mismatch(self):
        data=rows();data[0]['well_inflow_m3_s']+=.1
        with self.assertRaisesRegex(ValueError,'Rate columns'):balance.audit(data)
    def test_post_flow_changes_are_not_hidden(self):
        data=rows()
        for r in data:
            if r['phase']=='flow_converged':
                r['volume_end_m3']-=.1;r['storage_physical_m3_s']-=.01;r['storage_model_m3_s']-=.01;r['physical_defect_m3_s']=.01;r['equation_defect_m3_s']=.01;r['assembled_defect_m3_s']=.01
        report=balance.audit(data)
        self.assertAlmostEqual(report['post_flow_solve']['summed_volume_change_m3'],.2)
        self.assertAlmostEqual(report['flow_converged']['cumulative_defect_m3'],.2)
    def test_nonfinite_csv_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case_fluid_balance_rank0.csv';data=rows();data[0]['well_inflow_m3_s']=float('nan')
            with p.open('w',newline='') as f:
                writer=csv.DictWriter(f,sorted(balance.REQUIRED));writer.writeheader();writer.writerows(data)
            with self.assertRaisesRegex(ValueError,'Non-finite'):balance.load([p],1)
if __name__=='__main__':unittest.main()
