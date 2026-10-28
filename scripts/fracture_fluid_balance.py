#!/usr/bin/env python3
"""Audit accepted-step legacy fracture volume balances; never retune tolerances."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

PHASES = ('flow_converged', 'accepted_end')
TEXT = {'phase', 'control'}
REQUIRED = set(('schema step rank ranks phase start_s end_s dt_s fracture_dt_s ncells control '
                'applicable storage_enabled inner_converged volume_begin_m3 volume_old_on_current_m3 '
                'volume_end_m3 well_inflow_m3_s leakoff_m3_s storage_model_m3_s storage_physical_m3_s '
                'equation_defect_m3_s physical_defect_m3_s assembled_defect_m3_s cell_defect_linf_m3_s '
                'input_well_reservoir_water_rate_m3_s input_perf_rate_m3_s well_pressure_Pa '
                'flow_tolerance_m3_s').split())

def near(a, b, absolute=1e-9, relative=1e-10):
    return abs(a-b) <= absolute + relative*max(abs(a), abs(b))

def load(paths, ranks):
    groups = {}
    for path in paths:
        identity = re.sub(r'_fluid_balance_rank\d+\.csv$', '', path.name)
        rows = groups.setdefault(identity, [])
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream)
            if set(reader.fieldnames or []) != REQUIRED:
                raise ValueError(f'Unexpected balance schema: {path}')
            for raw in reader:
                row = {k: v if k in TEXT else float(v) for k,v in raw.items()}
                if any(not math.isfinite(v) for k,v in row.items() if k not in TEXT):
                    raise ValueError('Non-finite balance value')
                if row['schema'] != 1 or row['ranks'] != ranks or not 0 <= row['rank'] < ranks:
                    raise ValueError('Wrong schema or MPI rank metadata')
                if row['phase'] not in PHASES or row['dt_s'] <= 0 or row['ncells'] < 1 or row['flow_tolerance_m3_s'] <= 0:
                    raise ValueError('Invalid phase, timestep or cell count')
                if not near(row['end_s']-row['start_s'], row['dt_s']):
                    raise ValueError('Timestep start/end disagree')
                if row['storage_enabled'] and not near(row['fracture_dt_s'], row['dt_s']):
                    raise ValueError('Fracture timestep differs from accepted timestep')
                rows.append(row)
    if not groups or any(not rows for rows in groups.values()):
        raise ValueError('No balance rows')
    return groups

def audit(rows):
    by_key = {}
    for row in rows:
        key = (row['step'], row['phase'])
        if key in by_key:
            raise ValueError('Duplicate accepted phase (possibly duplicate MPI fracture ownership)')
        by_key[key] = row
    steps = sorted({row['step'] for row in rows})
    for step in steps:
        if any((step,phase) not in by_key for phase in PHASES):
            raise ValueError('Missing accepted-step phase')
        a,b = [by_key[step,phase] for phase in PHASES]
        if any(not near(a[k], b[k]) for k in ['start_s','end_s','volume_begin_m3']):
            raise ValueError('Two phases do not describe the same step')
    output = {}
    for phase in PHASES:
        chain = [by_key[step,phase] for step in steps]
        for previous, current in zip(chain, chain[1:]):
            if not near(previous['end_s'], current['start_s'], absolute=1e-7):
                raise ValueError('Gap or overlap in accepted-step times')
        # The next checkpoint is the accepted END state, not the pre-end solve.
        for step, next_step in zip(steps,steps[1:]):
            if not near(by_key[step,'accepted_end']['volume_end_m3'], by_key[next_step,phase]['volume_begin_m3'], absolute=1e-8):
                raise ValueError('Checkpoint volume is not continuous across accepted steps')
        for row in chain:
            physical = (row['volume_end_m3']-row['volume_begin_m3'])/row['dt_s']
            if not near(physical,row['storage_physical_m3_s']):
                raise ValueError('Physical storage column disagrees with volumes')
            for kind in ['physical','equation']:
                storage = row['storage_physical_m3_s' if kind=='physical' else 'storage_model_m3_s']
                if not near(row['well_inflow_m3_s']-row['leakoff_m3_s']-storage,row[f'{kind}_defect_m3_s']):
                    raise ValueError('Rate columns disagree with reported balance defect')
        inflow = sum(r['well_inflow_m3_s']*r['dt_s'] for r in chain)
        leakoff = sum(r['leakoff_m3_s']*r['dt_s'] for r in chain)
        delta = sum(r['volume_end_m3']-r['volume_begin_m3'] for r in chain)
        defect = inflow-leakoff-delta
        output[phase] = dict(steps=len(chain),start_s=chain[0]['start_s'],end_s=chain[-1]['end_s'],
            inflow_m3=inflow,leakoff_m3=leakoff,summed_step_volume_change_m3=delta,
            first_volume_m3=chain[0]['volume_begin_m3'],last_volume_m3=chain[-1]['volume_end_m3'],
            cumulative_defect_m3=defect,
            sum_abs_step_defects_m3=sum(abs(r['physical_defect_m3_s'])*r['dt_s'] for r in chain),
            integrated_equation_defect_m3=sum(r['equation_defect_m3_s']*r['dt_s'] for r in chain),
            independent_and_assembled_defects_agree_within_solver_tolerance=all(abs(r['equation_defect_m3_s']-r['assembled_defect_m3_s']) <= r['flow_tolerance_m3_s'] for r in chain),
            relative_defect_to_inflow=abs(defect)/max(abs(inflow),1e-30),
            relative_defect_to_storage_change=abs(defect)/max(abs(delta),1e-30),
            max_abs_physical_defect_m3_s=max(abs(r['physical_defect_m3_s']) for r in chain),
            max_abs_equation_defect_m3_s=max(abs(r['equation_defect_m3_s']) for r in chain),
            max_abs_domain_volume_change_m3=max(abs(r['volume_old_on_current_m3']-r['volume_begin_m3']) for r in chain),
            max_independent_vs_assembled_defect_m3_s=max(abs(r['equation_defect_m3_s']-r['assembled_defect_m3_s']) for r in chain),
            max_cell_residual_over_solver_tolerance=max(r['cell_defect_linf_m3_s']/r['flow_tolerance_m3_s'] for r in chain),
            all_applicable=all(r['applicable']==1 for r in chain),
            all_inner_converged=all(r['inner_converged']==1 for r in chain),
            all_cell_residuals_within_solver_tolerance=all(r['cell_defect_linf_m3_s'] <= r['flow_tolerance_m3_s']*1.000001+1e-10 for r in chain))
    end_changes = [by_key[step,'accepted_end']['volume_end_m3']-by_key[step,'flow_converged']['volume_end_m3'] for step in steps]
    output['post_flow_solve'] = dict(max_abs_volume_change_m3=max(map(abs,end_changes)),
        summed_volume_change_m3=sum(end_changes))
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--log',type=Path,required=True)
    p.add_argument('--ranks',type=int,required=True)
    p.add_argument('--report-dir',type=Path,required=True)
    args=p.parse_args()
    paths=sorted(args.run_dir.glob('*_fluid_balance_rank*.csv'))
    log=args.log.read_text(errors='replace')
    if 'End of simulation' not in log or [int(n) for n in re.findall(r'Number of MPI processes:\s*(\d+)',log)] != [args.ranks]:
        raise ValueError('Missing completed simulation or mismatched MPI footer')
    groups=load(paths,args.ranks)
    report=dict(schema=1,status='measured',ranks=args.ranks,
        fractures={name:audit(rows) for name,rows in groups.items()},
        interpretation='Volumetric backward-Euler budget at both accepted-step phases. Metrics are measurements; no conservation tolerance has been invented. Compare cell residuals with the configured solver tolerance, and report cumulative defect relative to both throughput and storage.',
        analyzed_inputs={str(path.resolve()):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths+[args.log]})
    args.report_dir.mkdir(parents=True,exist_ok=False)
    (args.report_dir/'fluid-balance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['fractures'],indent=2))

if __name__=='__main__':
    main()
