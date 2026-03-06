# STRESSEQUILNUM

Assigns stress equilibration region for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one positive integer for each grid block in the current input box, specifying the block's stress equilibration region. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `123*4`.

## Sections

- REGIONS

## Data

This is a *data* keyword containing a flat array of `INT` values.

## Examples

Suppose the model has 100-by-100-by-50 input cells, for a total of 500,000 grid blocks. The following statement places all grid blocks in stress equilibration region 1.
```
STRESSEQUILNUM
  500000*1 /
```
