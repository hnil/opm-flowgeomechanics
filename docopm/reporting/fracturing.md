## Summary vectors for fracturing process

These summary vectors show how fractures develop over time.

The naming convention is the same as for other summary vectors, with names starting with `C` corresponding to the connection level and names starting with `W` corresponding to the well level. Furthermore, if the connection block is defaulted in a connection level vector, as in
```
CFRWI
  'I1' * /
/
```
then additional vectors will be allocated for fracture connections as fracturing develops.

| Vector | Dimension | Description |
| ------ | --------- | ----------- |
| `[CW]WIRFRAC` | Surface condition liquid flow rate | Water injection rate going into fracture |
| `[CW]WITFRAC` | Liquid surface volume | Total injected water volume going into fracture |
| `CFRWI` | Transmissibility | Connection transmissibility factor due to fracturing |
| `CFRAREA` | Area | Total fracture area |
| `CFRHEIGH` | Length | Fracture half-thickness |
| `CFRLENGT` | Length | Fracture half-length |
| `CFRVOLUM` | Volume | Total fracture volume |
| `CFRFVOLU` | Volume | Total filter cake volume on fracture wall |
| `CFCFFRAC` | | Flow rate fraction affected by filtrate/filter cake |
| `CFRAVG` | Length | Average fracture thickness |
| `CFRAVGFW` | Length | Average filter cake thickness |
| `CFCFRATE` | Volume/Time | Connection flow rate attributable to fracture |
| `CFRPMAX` | Pressure | Maximum pressure across fracture |
| `CFRPMIN` | Pressure | Mimimum pressure across fracture |
| `CFRPAVG` | Pressure | Average pressure across fracture |
| `CFRPSTD` | Pressure | Pressure standard deviation across fracture |
| `CFRIRMAX` | Surface condition liquid flow rate | Maximum flow rate across fracture |
| `CFRIRMIN` | Surface condition liquid flow rate | Mimimum flow rate across fracture |
| `CFRIRAVG` | Surface condition liquid flow rate | Average flow rate across fracture |
| `CFRIRSTD` | Surface condition liquid flow rate | Flow rate standard deviation across fracture |
| `CFRWDMAX` | Length | Maximum width across fracture |
| `CFRWDMIN` | Length | Mimimum width across fracture |
| `CFRWDAVG` | Length | Average width across fracture |
| `CFRWDSTD` | Length | Width standard deviation across fracture |
| `CFRINJPR` | Pressure | Injection pressure |
| `CFRINJBH` | Pressure | Injection bottom-hole pressure |
| `CFRINJRA` | Volume/Time | Injection rate |

## Examples

### Show surface condition liquid flow rate going into fractures from well I1
```
WWIRFRAC
  'I1' /
```

### Show total injected volume, at surface conditions, going into fracture from connection (8,8,16) of well I1.
```
CWITFRAC
  'I1'  8  8  16 /
/
```

### Show connection transmissibility factors due to fracturing in all connections, including fracture connections, from well I2
```
CFRWI
  I2 * /
/
```
