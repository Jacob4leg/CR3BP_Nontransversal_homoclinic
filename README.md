# CR3BP: non-transversal homoclinic orbit to a vertical Lyapunov orbit

Source code of the computer-assisted proof in the paper

> *Non-transversal homoclinic connection to a periodic orbit in the Circular
> Restricted Three-Body Problem — a computer-assisted proof.*

The program proves Theorem 6 of the paper (proof in Section 5). For the
Sun–Jupiter mass parameter and the energy interval

```
mu = 0.00095388114032796904
E  = 3.0268216642156945825 + 1e-8 * [-1,1]
```

it shows that
- for every `e` in `E` there is a vertical Lyapunov orbit near `L_1` on the
  level set `{C_mu = e}` of the Jacobi integral, given by a fixed point `v(e)`
  of `F = P_y^4` on the section `Pi_y = {y = 0}`;
- for some `e` in `E` this orbit has a non-transversal homoclinic connection
  to itself.

All rigorous computations use interval arithmetic and the
[CAPD](https://github.com/CAPDGroup/CAPD) library.

## Requirements

- CAPD (built from source, see the CAPD README),
- a C++17 compiler (tested with GCC 16.2),
- GNU Make.

## Build

If `capd-config` is in `PATH`:

```
make
```

Otherwise, set `CAPDBINDIR` to the directory that contains `capd-config`:

```
make CAPDBINDIR=/path/to/CAPD/build/bin
```

## Run

```
./main
```

The program uses all hardware threads. On an AMD Ryzen 7 7735U (16 threads)
the full run takes about 7 minutes. If the proof succeeds, the last line is

```
The existence of a homoclinic orbit proved successfully
```

If some verification fails, the program prints the reason for each failed piece
and does not print this line.

### Input: `init_data.txt`

Approximate (non-rigorous) points `v(e)`, i.e. points on vertical Lyapunov
orbits at energy `e`. One record per line:

```
tau e x y z dx dy dz
```

`(x, y, z, dx, dy, dz)` approximates `v(e)`. `(tau, e)` is a point on the
boundary of `Omega = (tau~, e~) + (-1e-10, 1e-10) x (-1e-8, 1e-8)` and
marks the piece of `Omega` on which this `v(e)` is used. Two consecutive
records give one piece `[tau_l, tau_r] x [e_l, e_r]`, and the hull of their
`v(e)` is the initial box for Step 1.

The file has 4 blocks separated by empty lines: the bottom and top walls of `Omega` 
(`e` fixed, 301 records each) and the left and right walls (`tau` fixed, 8001 records each).
Consecutive records give the 300 + 300 + 8000 + 8000 = 16 600 boundary
pieces. The right wall also gives the 2000 strips that cover `Omega`.

If `init_data.txt` does not exist, `main` generates it first
(`NonRigorousSolver::save_init_data_to_file`). To regenerate an existing file,
set `generate_file = 1` in `src/main.cpp`. The rigorous part does not depend
on the accuracy of these points: they only define the boxes in which the
proof is attempted.

### Output: `output.txt`

One line per boundary piece, in the order of `init_data.txt`:

```
x_left x_right y_left y_right
```

This box encloses `H([0,1], piece)`: it is the hull of `B(piece - (tau~, e~))` 
and the enclosure of `G(piece)`, i.e. of the `(dx, dz)` components of `J P_y P_x(w)` (Step 3).

## Structure

| File | Contents |
| --- | --- |
| `src/main.cpp` | Reads `init_data.txt`, builds the 16 600 boundary pieces and 2000 strips, checks that the Step 1 Jacobian is regular on the hull of all their Step 1 boxes (so all of them give the same `v(e)`), runs them in parallel, writes `output.txt`. Computes the matrices `J` and `L`. |
| `include/ComputeImageTask.h` | One piece of `Omega`: Step 1 (`prove_periodic_orbit`), Step 2 (`get_proper_bound_for_z`, `cone_conditions`; cone conditions are checked on the strips), Step 3 (enclosure of `H` and the check `0 ∉ H`). `is_jacobian_regular`: the hull check called from `main`. |
| `include/NonRigorousSolver.h` | Non-rigorous part: approximate `v(e)`, `init_data.txt`, the change of basis `L`, the matrix `J^{-1}`. |
| `include/cr3bp.h` | Vector field, Jacobi integral, Poincaré maps on `{x = 0}`, `{y = 0}`, `{z = 0}`. |
| `include/Constants.h` | `mu`, the approximate points `ṽ`, `w̃`, `(τ̃, ẽ)`, cone parameters `alpha_h`, `alpha_v`, `m_h`, `m_v`, size and subdivision of `Omega`, `tau_max`. |
| `include/MatrixTools.h` | Small matrix helpers. |

Notation in the code vs. the paper:
- `pm_y` crosses `Pi_y` only from `y > 0` to `y < 0`, so `(P_y)^2` in the code
  is `F = P_y^4` in the paper;
- the code works in coordinates `(x, y, z, dx, dy, dz)`; the paper drops `y = 0`
  and writes points of `Pi_y` as `(x, z, dx, dy, dz)`.
