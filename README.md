# 3D Spherical FDTD Simulation

FDTD simulator code for lightning-based passive radar. This is a three-dimensional
finite-difference time-domain (FDTD) solver in spherical coordinates. It propagates the
low-frequency field of a lightning return stroke over the sea surface and past perfectly
conducting (PEC) boxes that stand in for ships.

The solver models a thin slab of atmosphere on a spherical Earth. With the default
settings the slab is 50 km long (east-west), 7 km wide (north-south) and 5 km high. Cells
are about 30 m, one twentieth of the wavelength at 500 kHz, which gives roughly
66 million cells. The model includes:

- a lightning return-stroke current source (MTLE: modified transmission line model with
  exponential decay)
- CPML (convolutional perfectly matched layer) absorbing layers on the four side walls
  and the top
- a surface impedance boundary condition (SIBC) for the seawater surface, in a
  single-frequency form or a dispersive form fitted with Vector Fitting
- grid-aligned PEC boxes for a target and for clutter
- probes that record the radial electric field `E_r` along the propagation path

## Repository layout

| Path | Contents |
|---|---|
| `spherical_fdtd_cpml_allwalls_SIBC.m` | The solver. A single function; the physical settings are variables near the top of the file. |
| `setup_clutter_and_sweep.m` | Script that builds a fixed random clutter field and an 11-position target sweep and saves them to `clutter_sweep.mat`. |
| `clutter_sweep_driver.m` | Script that loads `clutter_sweep.mat` and runs the solver for a background case, a clutter-only case and each target position. |
| `visualize_boxes.m` | Function that draws a plan view of a box list: footprints, CPML walls, keep-out margins and the target path. |
| `vfit3/` | Third-party Vector Fitting package. The solver calls `vectfit3.m`; the other files are the package's own examples, data and papers. |
| `Narrow Slice Probe Results/` | Six saved probe records from July 2026. See [Example results](#example-results). |

## Requirements

- **MATLAB R2021a or newer.** The solver uses `Name=Value` argument syntax, which older
  releases do not accept.
- **Signal Processing Toolbox**, for the plots made at the end of every run (`butter`,
  `zp2sos`, `filtfilt`, `spectrogram`, `pow2db`).
- **The `vfit3` folder on the MATLAB path.** The solver fits the seawater surface
  impedance with `vectfit3` at the start of every run, whichever SIBC setting is chosen.
- **Memory.** On the default grid the solver keeps roughly 7 GB of arrays for the whole
  run, and each field update creates about a dozen full-size temporary arrays, roughly
  6 GB more. These figures are estimated from the array sizes and have not been measured.

## Quick start

Run these from the repository root.

```matlab
addpath(pwd, fullfile(pwd, 'vfit3'));

% One run with no scatterers
cfg = struct('box_list', zeros(0, 6));
spherical_fdtd_cpml_allwalls_SIBC(cfg);
```

The solver must be given a struct that has a `box_list` field, and `vfit3` must already
be on the path. Both points are explained under [Known issues](#known-issues).

To add one PEC box and save the probe record:

```matlab
cfg = struct();
cfg.box_list  = [1 2 115 122 1503 1504];
cfg.save_path = 'single_box.mat';
spherical_fdtd_cpml_allwalls_SIBC(cfg);
```

On the default grid that row is a box 30 m high, about 210 m long in theta and 30 m long
in phi, sitting on the sea surface at the 45 km mark. It is what the solver's built-in
single-box settings work out to.

At start-up the solver prints the grid size, time step and boundary settings. During the
run it prints a progress line with the elapsed time and an estimated time remaining about
every 2% of the steps. The default run is 9,128 steps of 54.8 ns, which is 0.50 ms of
simulated time.

## What a run produces

All output is written relative to MATLAB's current folder.

- **`probe_data<timestamp>.mat`**, written by every run. It holds `probe_record` (one row
  per time step, one column per probe, `E_r` in V/m), `del_t` (the time step in seconds),
  `probe_ranges_km`, `sibc_enabled` and the six `cpml_ct_*` wall thicknesses. Row `n` of
  `probe_record` is at time `(n-1)*del_t`.
- **The file named in `cfg.save_path`**, if one is given. It holds `probe_record`, `del_t`
  and `box_list`, saved in `-v7.3` format. It is written as soon as time stepping ends,
  before the plots.
- **Five figures:** each probe in dB against time; each probe raw and band-pass filtered
  to 100-300 kHz; a spectrogram of each filtered probe; the channel-base current; and the
  periodogram of that current.

## Configuring a run

### The `ext_cfg` argument

| Field | Required | Meaning |
|---|---|---|
| `box_list` | Yes | N-by-6 array with one PEC box per row, `[i1 i2 j1 j2 k1 k2]`. Use `zeros(0,6)` for no boxes. |
| `save_path` | No | File to save `probe_record`, `del_t` and `box_list` to. |
| `viz_enable` | No | `true` shows three slices of log10 of the `E_r` magnitude every 75 steps. Default `false`. |

A box row gives node index ranges in r, theta and phi. The box fills the cells between
those nodes, so it is `i2-i1` cells high, `j2-j1` cells long in theta and `k2-k1` cells
long in phi. `i1 = 1` puts its base on the sea surface. Each time step the solver sets the
electric field to zero inside the box and the tangential electric field to zero on its
faces. It warns, without stopping, when a box is degenerate, lies outside the grid or
covers the source column.

### Settings inside the solver file

Everything else is a variable near the top of `spherical_fdtd_cpml_allwalls_SIBC.m`.

| Variable | Default | Meaning |
|---|---|---|
| `f_max` | `500e3` Hz | Highest frequency resolved. Sets the cell size `ds` to one twentieth of the wavelength, 29.98 m. |
| `start_point` | `[0, 0]` | Latitude and longitude, in degrees, of the centre of the grid's east edge. |
| `dist_max` | `50e3` m | Length of the grid along the surface, westward from `start_point`. |
| `width_max` | `7e3` m | North-south width, centred on the start latitude. |
| `height_max` | `5e3` m | Height above the surface. |
| `below_ground` | `0` m | Depth of grid below the surface. Leave at 0: the sea is represented by the bottom boundary. |
| `cfl_safety` | `0.95` | Fraction of the CFL limit used for the time step. |
| `cpml_ct_west`, `_east`, `_north`, `_south`, `_top` | `10` cells | CPML thickness on each wall. `0` turns that wall's CPML off. |
| `cpml_m`, `cpml_kappa_max`, `cpml_alpha_max` | `3`, `2.07`, `2*pi*f_max*eps0` | CPML grading order, maximum kappa and maximum alpha. |
| `sibc_enabled` | `true` | `false` makes the sea surface a perfect conductor. |
| `sibc_mode` | `'conv'` | `'conv'` is the dispersive SIBC. `'mono'` is the single-frequency SIBC. |
| `sibc_f0` | `200e3` Hz | Design frequency of the `'mono'` SIBC. |
| `sibc_sigma` | `4.0` S/m | Seawater conductivity used by the SIBC. |
| `vf_Np`, `vf_flo`, `vf_fhi` | `6`, `1e4` Hz, `1e6` Hz | Number of poles and frequency band of the surface impedance fit. |
| `test_source_enabled` | `false` | `true` replaces the lightning source with a differentiated-Gaussian point source at the centre of the grid, peaked at `test_f0` (200 kHz). |
| `v_return`, `lambda_decay` | `(7/12)*c`, `1e3` m | MTLE return-stroke speed and current decay height. |
| `channel_base_current` | 10 kA double exponential | Local function at the end of the file. Decay time constant 50 µs, rise time constant 0.3 µs. |
| `probe_ranges_km` | `[1 5 10 15 20 25 30 35 40 45]` | Probe positions. See [Ranges](#ranges-are-measured-from-the-east-edge). |
| `viz_every`, `viz_save_pngs`, `viz_png_dir` | `75`, `false`, `'./fdtd_frames'` | Live-plot interval in steps, and optional saving of each frame as a PNG. |

The number of time steps is `ceil(3 * (dist_max/c) / del_t)`, three transit times across
the grid.

## Geometry and conventions

- **Coordinates.** r is radius, theta is colatitude (0 at the north pole, pi/2 at the
  equator) and phi increases westward.
- **Indices.** `i` runs along r with `i = 1` at the sea surface. `j` runs along theta with
  `j = 1` at the north edge. `k` runs along phi with `k = 1` at the east edge.
- **Default grid.** `Nr = 168`, `Nthe = 235`, `Nphi = 1671`. The solver prints these at
  start-up.
- **Source.** A vertical channel from the surface to about 1 km (34 cells) on the centre
  line `j = round(Nthe/2)`. Its phi index is the east CPML thickness plus the far-field
  distance of the 1 km channel at the shortest wavelength, which comes to `k = 127`.
- **Probes.** `E_r` in the first cell above the surface (`i = 1`) on the centre line, at
  the phi node nearest each entry of `probe_ranges_km`.

### Ranges are measured from the east edge

Probe ranges, the built-in box range and the sweep's `tgt_range_km` are all measured
along the surface from the east edge of the grid (`k = 1`), not from the source. With the
default settings the source column is 3.77 km west of that edge, so:

| Label | phi index `k` | Distance from the source |
|---|---|---|
| Source | 127 | 0 |
| 1 km probe | 34 | 2.79 km east of the source, on the opposite side from the other probes |
| 5 km probe | 168 | 1.23 km |
| 10 km probe | 335 | 6.23 km |
| 45 km probe | 1503 | 41.22 km |

A probe that lies inside a PEC box records exactly zero.

## Target-in-clutter sweep

The sweep keeps one clutter field fixed and moves a single target across theta, so that
later processing can separate the target from the clutter.

```matlab
addpath(pwd, fullfile(pwd, 'vfit3'));
setup_clutter_and_sweep      % writes clutter_sweep.mat and plots the layout
clutter_sweep_driver         % 13 solver runs; results go to clutter_runs/
```

**`setup_clutter_and_sweep`** builds, with its default values:

- a target nominally 30 m high, 200 m long in theta and 30 m long in phi (1 x 7 x 1
  cells) at the 45 km mark, at 11 positions across theta with centres
  `j = 20, 40, 59, 79, 98, 118, 138, 157, 177, 196, 216`
- 20 clutter boxes of the same size, each placed at random with its long side along theta
  or phi, using the fixed seed `rng(98457)` so the field is reproducible
- keep-out rules for the clutter: at least 15 cells from every side wall (10 CPML cells
  plus 5 cells for a ship's 147 m far-field distance), a gap of at least one cell to any
  other box and to the strip the target sweeps through, and nothing at a phi index below
  `phi_min = 235`

It saves the clutter list, the target positions and the grid constants to
`clutter_sweep.mat` and plots the layout with `visualize_boxes`.

**`clutter_sweep_driver`** loads that file and calls the solver for an empty background,
for the clutter alone, and for the clutter plus the target at each position. The switches
`run_positions`, `do_background` and `do_clutteronly` near the top choose which runs to
do; setting `run_positions = 6` tries only the centre position. Each run is a complete
simulation. The script starts with `clear`, so it empties the workspace. Results are
saved to:

| File in `clutter_runs/` | Run |
|---|---|
| `run_background.mat` | No boxes |
| `run_clutteronly.mat` | Clutter only |
| `run_posPP_jJJJJ.mat` | Clutter plus the target at position `PP`, theta centre index `JJJJ` |

Each file holds `probe_record`, `del_t` and `box_list`.

Things to know about the sweep:

- **The setup script has its own copy of the grid.** `Nr`, `Nthe`, `Nphi`, `del_the`,
  `del_phi`, `cpml_cells` and `phi_min` are typed into `setup_clutter_and_sweep.m` and
  currently match the solver defaults. If you change `f_max`, `dist_max`, `width_max`,
  `height_max` or the CPML thickness in the solver, change them there too.
- **The centre position covers the 45 km probe.** At position 6 the target sits on the
  45 km probe, which records zero for that run.
- **Figures and files pile up.** Every solver call opens five figures and writes a
  `probe_data<timestamp>.mat` file. The driver closes figures only after the background
  run, after the clutter-only run and after the last position, so a full sweep ends with
  more than 60 figures open.
- **The post-processing script is not here.** The driver's last message points to
  `integrate_clutter_suite.m`, which is not in this repository.

## Example results

`Narrow Slice Probe Results/` holds six `probe_data*.mat` files written between 16 and
21 July 2026. Each has a 14,604-by-16 `probe_record` for probes labelled 1, 5, 10, ...,
75 km, with `del_t` = 54.8 ns, the SIBC enabled and 10-cell CPML on every wall except the
bottom. The step count matches an 80 km grid, not the current 50 km default. In five of
the files one probe is zero throughout (35, 45, 55, 65 or 75 km), as happens when a PEC
box sits on a probe; the 17 July file has no such probe. The files do not record the grid
width or the box that was used.

```matlab
S = load(fullfile('Narrow Slice Probe Results', 'probe_data17_Jul_202609_56_42.mat'));
t_us = (0:size(S.probe_record, 1) - 1) * S.del_t * 1e6;
plot(t_us, S.probe_record(:, S.probe_ranges_km == 20));
xlabel('Time (\mus)'); ylabel('E_r at the 20 km probe (V/m)');
```

## Known issues

1. **The solver stops unless it is given `box_list`.** The line
   `box_enabled = ~isempty(ext_cfg.box_list);` (line 537) runs before the check that an
   argument was supplied. Calling the solver with no argument, or with a struct that has
   no `box_list` field, ends in an error there. The built-in single-box settings just
   below that line are therefore never used; pass the box explicitly as shown in
   [Quick start](#quick-start).
2. **The path to `vfit3` is hard-coded.** Line 382 adds
   `D:\GitHub Files\DAIRS Analysis\Simulation Suite\FDTD_Simulation\vfit3` to the path.
   Where that folder does not exist MATLAB prints a warning and carries on, and
   `vectfit3` is found only if `vfit3` was added to the path beforehand.
3. **Vector Fitting draws into figure 1.** The solver passes `'spy', 0` to `vectfit3`,
   but the option that turns its plot off is `spy2`, which defaults to 1. Each of the six
   fitting passes therefore plots into figure 1 and replaces whatever was there, for
   example a `visualize_boxes` layout.
4. **Range labels are not distances from the source.** See
   [Ranges](#ranges-are-measured-from-the-east-edge).
5. **Some comments describe older versions.** The solver header says to run
   `spherical_fdtd_cpml_allwalls`, says the bottom is PEC and says a call with no argument
   works. The driver header calls the file `run_sweep_driver.m` and asks for
   `spherical_fdtd_cpml_allwalls_SIBC_v2`. The setup header refers to
   `example_ext_cfg.m`. Those names are not files in this repository.
6. **Unused code remains in the helper files.** `setup_clutter_and_sweep.m` defines
   `phys2box` with a function `deal_box` that does not exist, and builds an `opts` struct
   it never uses. `visualize_boxes.m` documents a `.guard_cells` option that the code does
   not read. None of this runs.

## Third-party code

`vfit3/` is the Vector Fitting package (Fast Relaxed Vector Fitting, version 1.0, 2008) by
B. Gustavsen, SINTEF Energy Research. The header of `vectfit3.m` restricts the program to
non-commercial use and asks that scientific work using it cite:

1. B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by
   Vector Fitting," IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.
2. B. Gustavsen, "Improving the pole relocating properties of vector fitting," IEEE Trans.
   Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.
3. D. Deschrijver, M. Mrozowski, T. Dhaene and D. De Zutter, "Macromodeling of Multiport
   Systems Using a Fast Implementation of the Vector Fitting Method," IEEE Microwave and
   Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.

`vfit3/vectfit3_userguide.pdf` documents the function and its options.

## References

Background for the methods named in the code comments:

- **CPML.** J. A. Roden and S. D. Gedney, "Convolution PML (CPML): An efficient FDTD
  implementation of the CFS-PML for arbitrary media," Microwave and Optical Technology
  Letters, vol. 27, no. 5, pp. 334-339, 2000.
- **SIBC in FDTD.** J. H. Beggs, R. J. Luebbers, K. S. Yee and K. S. Kunz,
  "Finite-difference time-domain implementation of surface impedance boundary
  conditions," IEEE Trans. Antennas and Propagation, vol. 40, no. 1, pp. 49-56, 1992.
- **MTLE return-stroke model.** C. A. Nucci, G. Diendorfer, M. A. Uman, F. Rachidi,
  M. Ianoz and C. Mazzetti, "Lightning return stroke current models with specified
  channel-base current: A review and comparison," Journal of Geophysical Research,
  vol. 95, no. D12, pp. 20395-20408, 1990.
