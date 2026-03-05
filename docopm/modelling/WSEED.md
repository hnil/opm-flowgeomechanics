# WSEED

Defines a fracture seed point associated with a well–reservoir connection. A seed point locates an initial fracture at a specific Cartesian cell given by indices (I, J, K) for a named well connection. Each seed record ends with a slash (`/`) and the `WSEED` block is terminated by a line containing only `/`. Seed records that reference non-existent well–reservoir connections are ignored; the simulator may emit a warning when this occurs.

This keyword applies to standard (single-segment) wells only; multi-segment well seed points are not supported.

## Sections

- SCHEDULE

## Format

Records following `WSEED` are given in fixed field order. Fields shown as optional may be omitted and will take their default values.

Field order (one record per line, terminated by `/`):

- `WELL` (STRING): well name identifying the connection.
- `I` `J` `K` (INT): grid cell indices of the seed point (model indexing; see "Indexing & coordinates" below).
- `NORMAL_X` `NORMAL_Y` `NORMAL_Z` (DOUBLE): components of the fracture-plane normal vector (direction only).
- `SIZE_Z` (DOUBLE, default 0.3): seed length along the normal direction (length units).
- `SIZE_H` (DOUBLE, default 0.3): in-plane seed half-length orthogonal to the normal (length units).
- `WIDTH` (DOUBLE, default 1.0E-4): initial fracture width/thickness (length units).

Each record must end with a `/`. The `WSEED` block is closed by a line containing only `/`.

## Items

| # | Name | Type | Default | Dimension | Description |
|---:|------|------:|--------:|-----------|-------------|
| 1 | `WELL` | STRING | — | — | Well name for the seed point (string literal). |
| 2 | `I` | INT | — | — | Grid index I of the cell containing the seed. |
| 3 | `J` | INT | — | — | Grid index J of the cell containing the seed. |
| 4 | `K` | INT | — | — | Grid index K of the cell containing the seed. |
| 5 | `NORMAL_X` | DOUBLE | 1.0 | Length | X component of the fracture-plane normal (direction). |
| 6 | `NORMAL_Y` | DOUBLE | 0.0 | Length | Y component of the fracture-plane normal (direction). |
| 7 | `NORMAL_Z` | DOUBLE | 0.0 | Length | Z component of the fracture-plane normal (direction). |
| 8 | `SIZE_Z` | DOUBLE | 0.3 | Length | Seed length parallel to the normal vector. |
| 9 | `SIZE_H` | DOUBLE | 0.3 | Length | Seed length in the fracture plane (orthogonal to normal). |
| 10 | `WIDTH` | DOUBLE | 1.0E-4 | Length | Initial fracture aperture/width. |

## Indexing & coordinates

- I/J/K follow the model's grid indexing convention (replace below with your model's convention if different). In ECLIPSE-style grids, indices are 1-based and `K = 1` denotes the top layer with `K` increasing downward.
- The seed is located at the centroid of the cell identified by `(I,J,K)`.

## Normal vector semantics

- `NORMAL_X/Y/Z` specify the direction of the fracture-plane normal. Only the direction (unit vector) is used; components are treated as direction cosines and the vector will be normalized internally. The magnitude provided is ignored, but the relative magnitude between the components matters.

## Units

Numeric length values use the simulator's active length unit (for example, METRIC → metres). Default numeric values shown above assume metres.

## Behavior for invalid references

Seed records that reference well–reservoir connections that do not exist will be ignored. A warning may be emitted; check the simulator log for details.

## Examples

### Single seed (minimal fields)

Seed a fracture for well `I1` at cell (11, 22, 33). The normal is (1,1,1) and the seed length along the normal is 0.25 (length units). Other numeric fields use defaults.

```
WSEED
 'I1'   11  22  33   1.0  1.0  1.0   0.25  /
/
```

### Two seeds (explicit fields)

First record uses defaults for `SIZE_H` and `WIDTH`. Second record specifies in-plane seed and width explicitly.

```
WSEED
 'I2'   10  11  12   1.0  1.0  1.0                  /
 'I2'   10  11  20   1.0  1.0  0.0   0.5  0.5  1.0E-3 /
/
```