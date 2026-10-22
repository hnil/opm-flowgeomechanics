# YMODULE

Assigns Young's modulus for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's Young modulus. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*6.78`.

## Sections

- GRID

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

Young's modulus is entered in units of giga-Pascal (GPa) in all unit systems.

## Notes

Set Young's modulus to `6.78 GPa` in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
YMODULE
  1331*6.78 /
```
