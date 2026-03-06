# CSTRESS

The CSTRESS item is used to set the critcal stress in mechanics models.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's critical stress value. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*0.123E4`.

## Sections

- GRID

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

The input values are treated as stresses in the SI units of Pascal, regardless of the model's unit system.

## Example

Set the rock's critical stress to 550,000 Pascal (0.55 MPa, 5.5 bar) in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
CSTRESS
  1331*0.55E6 /
```
