# STREQUIL

Defines equilibration parameters for computing the initial state of mechanical quantities--notably stresses--in simulation models that include geo-mechanical effects.

The keyword should followed by exactly one record for each stress equilibration region defined by the STRESSEQUILNUM keyword in the REGIONS section. There is no special termination character for the `STREQUIL` block.

## Format

Records following `STREQUIL` are given in fixed field order. Leading optional fields may be explicitly defaulted using the notation `1*` in which case they will take their default values. Trailing optional fields may be omitted in which case they will take their default values.

Field order:

- `DATUM_DEPTH` (DOUBLE, default 0.0): Depth of input stress reference point.
- `DATUM_POSX` (DOUBLE, default 0.0): $X$ position of stress reference point.
- `DATUM_POSY` (DOUBLE, default 0.0): $Y$ position of stress reference point.
- `STRESSXX` (DOUBLE): Stress tensor $XX$ component, $\sigma_{xx}$, at input reference point.
- `STRESSXXGRAD` (DOUBLE): Derivative of stress tensor $XX$ component $\sigma_{xx}$ with respect to depth at input reference point.
- `STRESSYY` (DOUBLE): Stress tensor $YY$ component, $\sigma_{yy}$, at input reference point.
- `STRESSYYGRAD` (DOUBLE): Derivative of stress tensor $YY$ component $\sigma_{yy}$ with respect to depth at input reference point.
- `STRESSZZ` (DOUBLE): Stress tensor $ZZ$ component, $\sigma_{zz}$, at input reference point.
- `STRESSZZGRAD` (DOUBLE): Derivative of stress tensor $ZZ$ component $\sigma_{zz}$ with respect to depth at input reference point.
- `STRESSXY` (DOUBLE, default 0.0): Stress tensor $XY$ component, $\sigma_{xy}$, at input reference point.
- `STRESSXYGRAD` (DOUBLE, default 0.0): Derivative of stress tensor $XY$ component $\sigma_{xy}$ with respect to depth at input reference point.
- `STRESSXZ` (DOUBLE, default 0.0): Stress tensor $XZ$ component, $\sigma_{xz}$, at input reference point.
- `STRESSXZGRAD` (DOUBLE, default 0.0): Derivative of stress tensor $XZ$ component $\sigma_{xz}$ with respect to depth at input reference point.
- `STRESSYZ` (DOUBLE, 0.0): Stress tensor $YZ$ component, $\sigma_{yz}$, at input reference point.
- `STRESSYZGRAD` (DOUBLE, 0.0): Derivative of stress tensor $YZ$ component $\sigma_{yz}$ with respect to depth at input reference point.

Each record must end with a `/`.

## Sections

- SOLUTION

## Items

| # | Name | Type | Default | Dimension | Description |
|---|------|------|---------|-----------|-------------|
| 1 | `DATUM_DEPTH` | DOUBLE | 0 | Length | Depth of stress reference point. |
| 2 | `DATUM_POSX` | DOUBLE | 0 | Length | $X$ position of stress reference point. |
| 3 | `DATUM_POSY` | DOUBLE | 0 | Length | $Y$ position of stress reference point. |
| 4 | `STRESSXX` | DOUBLE |  | Pressure | Stress tensor $XX$ component, $\sigma_{xx}$, at input reference point. |
| 5 | `STRESSXXGRAD` | DOUBLE |  | Pressure/Length | Derivative of stress tensor $XX$ component $\sigma_{xx}$ with respect to depth at input reference point. |
| 6 | `STRESSYY` | DOUBLE |  | Pressure | Stress tensor $YY$ component, $\sigma_{yy}$, at input reference point. |
| 7 | `STRESSYYGRAD` | DOUBLE |  | Pressure/Length | Derivative of stress tensor $YY$ component $\sigma_{yy}$ with respect to depth at input reference point. |
| 8 | `STRESSZZ` | DOUBLE |  | Pressure | Stress tensor $ZZ$ component, $\sigma_{zz}$, at input reference point. |
| 9 | `STRESSZZGRAD` | DOUBLE |  | Pressure/Length | Derivative of stress tensor $ZZ$ component $\sigma_{zz}$ with respect to depth at input reference point. |
| 10 | `STRESSXY` | DOUBLE | 0 | Pressure | Stress tensor $XY$ component, $\sigma_{xy}$, at input reference point. |
| 11 | `STRESSXYGRAD` | DOUBLE | 0 | Pressure/Length | Derivative of stress tensor $XY$ component $\sigma_{xy}$ with respect to depth at input reference point. |
| 10 | `STRESSXZ` | DOUBLE | 0 | Pressure | Stress tensor $XZ$ component, $\sigma_{xz}$, at input reference point. |
| 11 | `STRESSXZGRAD` | DOUBLE | 0 | Pressure/Length | Derivative of stress tensor $XZ$ component $\sigma_{xz}$ with respect to depth at input reference point. |
| 12 | `STRESSYZ` | DOUBLE | 0 | Pressure | Stress tensor $YZ$ component, $\sigma_{yz}$, at input reference point. |
| 13 | `STRESSYZGRAD` | DOUBLE | 0 | Pressure/Length | Derivative of stress tensor $YZ$ component $\sigma_{yz}$ with respect to depth at input reference point. |

## Units

Numeric pressure/stress values use the simulator's active pressure unit (for example METRIC → bars). Numeric length values use the simulator's active length unit (for example, METRIC → metres). Default numeric values shown above assume metres.

## Example

### Define a zero-stress point at the origin with a diagonal stress tensor increase with respect to depth.

Single stress equilibration region.
```
STREQIL
  1*  1*  1*  0.0  0.167  0.0  0.167  0.0  0.22 /
```

### Define a zero-stress point at (1000, 1000, 500) with gradually increasing stress as a function of depth

Two stress equilibration regions.
```
STREQUIL
  1000.0  1000.0  500.0 -- Input stress reference point location.
  0.0  0.1              -- XX component and its derivative with respect to depth.
  0.0  0.1              -- YY component and its derivative with respect to depth.
  0.0  0.15             -- ZZ component and its derivative with respect to depth.
  0.0  0.05             -- XY component and its derivative with respect to depth.
  0.0  0.01             -- XZ component and its derivative with respect to depth.
  0.0  0.04             -- YZ component and its derivative with respect to depth.
/ -- First equilibration region.
  1000.0  1000.0  500.0 -- Input stress reference point location.
  0.0  0.2              -- XX component and its derivative with respect to depth.
  0.0  0.2              -- YY component and its derivative with respect to depth.
  0.0  0.25             -- ZZ component and its derivative with respect to depth.
  0.0  0.15             -- XY component and its derivative with respect to depth.
  0.0  0.11             -- XZ component and its derivative with respect to depth.
  0.0  0.14             -- YZ component and its derivative with respect to depth.
/ -- Second equilibration region.
```