# THERMEXR

Assigns thermal expansion ratio for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's thermal expansion ratio. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*1.23E-4`.

## Sections

- GRID

## Dependencies

This keyword is mutually exclusive with the `THELCOEF` keyword, as the two keywords are alternative ways of defining the same underlying physical property, namely the rock's thermoelastic constant $A_T$.

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

The values should be input in units of reciprocal absolute temperature, meaning reciprocal Kelvin for the METRIC, LAB, and PVT-M unit conventions and as reciprocal degree Rankine for the FIELD unit conventions.

## Example

Set the thermal expansion ratio to 1.23 per Kelvin in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
THERMEXR
  1331*1.23 /
```
