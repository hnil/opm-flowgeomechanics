## Geomechanical responses

### Volumetric (3D) fields/vectors

These vectors are automatically output to the restart file when the simulation run enables geomechanical analysis through the `MECH` keyword in the `RUNSPEC` section and have one value for each active cell.

| Vector | Dimension | Description |
| ------ | --------- | ----------- |
| `STRESSuv` | Pressure | Component $uv$ of the symmetric stress tensor, $\sigma_{uv}$, $u,v\in\{X,Y,Z\}$ |
| `FRCSTRuv` | Pressure | Component $uv$ of the fracture stress tensor, $u,v\in\{X,Y,Z\}$ |
| `LINSTRuv` | Pressure | Component $uv$ of the mathematical stress tensor, $u,v\in\{X,Y,Z\}$ |
| `DELSTRuv` | Pressure | Component $uv$ of the stress tensor change, $u,v\in\{X,Y,Z\}$ |
| `STRAINuv` | | Component $uv$ of the strain tensor, $u,v\in\{X,Y,Z\}$ |
| `DISPu` | Length | Component $u$ of the displacment field |
| `MECHPOTF` | Pressure | Total mechanical potential |
| `PRESPOTF` | Pressure | Potential due to pressure differences |
| `TEMPPOTF` | Pressure | Potential due to temperature differences |

### Block level summary vectors

| Vector | Dimension | Description |
| ------ | --------- | ----------- |
| `BSTRSSXX` | Pressure | Component $XX$ of the symmetric stress tensor, $\sigma_{XX}$ |
| `BSTRSSYY` | Pressure | Component $YY$ of the symmetric stress tensor, $\sigma_{YY}$ |
| `BSTRSSZZ` | Pressure | Component $ZZ$ of the symmetric stress tensor, $\sigma_{ZZ}$ |
| `BSTRSSXY` | Pressure | Component $XY$ of the symmetric stress tensor, $\sigma_{XY}$ |
| `BSTRSSXZ` | Pressure | Component $XZ$ of the symmetric stress tensor, $\sigma_{XZ}$ |
| `BSTRSSYZ` | Pressure | Component $YZ$ of the symmetric stress tensor, $\sigma_{YZ}$ |
