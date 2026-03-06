# THELCOEF

Assigns thermal expansion coefficient for cells in simulation models that include geo-mechanical effects.

The keyword should be followed by one real number for each grid block in the current input box, specifying the rock's thermal expansion coefficient. Grid blocks are ordered with the X index cycling the fastest, followed by the Y index and finally the Z index. Repeat counts may be used for repeated values, for example `1000*1.23`.

## Sections

- GRID

## Dependencies

This keyword is mutually exclusive with the `THELMEXR` keyword, as the two keywords are alternative ways of defining the same underlying physical property, namely the rock's thermoelastic constant $A_T$.

## Data

This is a *data* keyword containing a flat array of `DOUBLE` values.

## Units

The values should be input in units of pressure per absolute temperature, meaning for example bars per Kelvin for the METRIC unit conventions and as psi per degree Rankine for the FIELD unit conventions.

## Example

Set the thermal expansion coefficient to 1.23 bars per Kelvin in the current input box of (`10..20`, `10..20`, `10..20`)--a total of 1331 grid blocks
```
THELCOEF
  1331*1.23 /
```
