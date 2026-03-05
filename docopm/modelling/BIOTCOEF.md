# BIOTCOEF

Assigns the Biot coefficient for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's Biot coefficient. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*1.23`.

## Sections

- GRID

## Dependencies

This keyword is mutually exclusive with the `POELCOEF` keyword, as the two keywords are alternative ways of defining the same underlying physical property. In particular, the poro-elastic coefficient $A_p$ is defined in terms of the Biot coefficient $\beta$ through the Poisson ratio $\nu$ as
$$
A_p = \frac{1-2\nu}{1-\nu}\beta.
$$

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

The Biot coefficient is a pure scalar with no associated unit of measurement.

## Notes

Set the Biot coefficient to `0.123e-4` in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
BIOTCOEF
  1331*0.123E-4 /
```
