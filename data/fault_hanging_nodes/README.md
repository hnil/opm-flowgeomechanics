# SIMPLE_MECH with a fault

`SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ` with a vertical fault between rows J=8 and J=9 and a
10 m (half-layer) throw, so every node on the fault on one side lies in the middle of an edge
on the other side (hanging nodes). `FAULT_GRID.INC` is made by `make_fault_grid.py 8 10`.

    flow_energy_geomech --parameter-file=../myparamsfrac_seq.ini SIMPLE_MECH_FAULT_SEQ.DATA \
      --fracture-param-file=../fracture_simple_bhp_seq_final.json \
      --linear-solver-mech=../amgmechsolverhypre_geosxcpu_new_nv.json \
      --enable-vtk-output=true --vtk-write-geo-mech-nodal=true --output-dir=<scratch>

The parameter file sets `EdgeConformal=true`. 600 of the 624 nodes on the fault plane hang.
The schedule is the 15-day one of the original deck: it ends at propagation onset.
