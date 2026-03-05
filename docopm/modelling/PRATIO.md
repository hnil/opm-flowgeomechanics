# PRATIO

Assigns Poisson's ratio for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's Poisson ratio. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*0.25`.

## Sections

- GRID

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

The Poisson ratio is a pure scalar with no associated unit of measurement.

## Notes

Set the Poisson ratio to `0.2` in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
PRATIO
  1331*0.2 /
```
