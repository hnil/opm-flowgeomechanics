## Summary vectors for filtrate injection

These summary vectors show how injected filtrate particles develop over time and partition between the well bore and fractures.

The naming convention is the same as for other summary vectors, with names starting with `C` corresponding to the connection level and names starting with `W` corresponding to the well level. Furthermore, if the connection block is defaulted in a connection level vector, as in
```
CFCFFVIR
  I1  * /
/
```
then additional vectors will be allocated for fracture connections as fracturing develops.

| Vector | Dimension | Description |
| ------ | --------- | ----------- |
| `[CW]FCFFVIR` | Volume/Time | Volumetric filtrate injection rate attributable to fracture |
| `[CW]FCFFVIT` | Volume | Total filtrate injection volume attributable to fracture |
| `[CW]FCWFVIR` | Volume/Time | Volumetric filtrate injection rate attributable to well |
| `[CW]FCWFVIT` | Volume | Total filtrate injection volume attributable to well |
| `[CW]FCFVIR` | Volume/Time | Volumetric filtrate injection rate attributable to both fracture and well (`[CW]FCFFVIR + [CW]FCWFVIR`) |
| `[CW]FCFVIT` | Volume | Total filtrate injection volume attributable to both fracture and well (`[CW]FCFFVIT + [CW]FCWFVIT`) |
| `[CW]INJFVR` | Volume/Time | Filtrate injection rate |
| `[CW]INJFVT` | Volume | Total injected filtrate volume |
| `WINJFC` | | Filtrate injection concentration |
| `CFCSKIN` | | Skin factor due to filter cake |
| `CFCWIDTH` | Length | Filter cake thickness |
| `CFCPERM` | Permeability | Filter cake permeability |
| `CFCPORO` | | Filter cake porosity |
| `CFCRAD` | Length | Filter cake radius |
| `CFCAOF` | Area | Filter cake area of flow |

## Examples

### Show volumetric filtrate injection rate attributable to fracture in well I1
```
WFCFFVIR
  'I1' /
```

### Show total injected filtrate volume attributable to fracture in all wells matching the name pattern 'I-*'
```
WFCFFVIT
  'I-*' /
```

### Show volumetric filtrate injection attributable to the well in connection (11,22,33) of well I1
```
CFCWFVIR
  'I1'  11 22 33 /
/
```

### Show total injected filtrate volume in connection (10,10,10) of well I1 and all connections, including those created by fracturing, in well I2
```
CFCFVIT
  'I1'  10 10 10 /
  'I2'  * /
/
```
