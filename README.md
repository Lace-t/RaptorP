# RaptorP

**RaptorP** (Raptor-PLUTO) is a polarized radiative transfer code based on [RAPTOR](https://github.com/jordydavelaar/raptor) (Bronzwaer et al. 2018, A&A 613, A2; Davelaar et al. 2018, CompAC 5, 1), adapted for performing **special-relativistic radiative transfer** on PLUTO MHD simulation data in **Cartesian coordinates**.

The code solves the polarized radiative transfer equation along geodesics in Minkowski spacetime (Cartesian Kerr-Schild metric as a special case), computing Stokes **I, Q, U, V** emission, absorption, and Faraday conversion/rotation through user-defined electron distribution functions.

---

## 1. Overview

RaptorP reads PLUTO (https://plutocode.ph.unito.it/) RMHD simulation output in `.dbl` format with a corresponding `grid.out` file, and produces synthetic polarized intensity maps and spectra. It supports:

- **Cartesian coordinates** 
- **Uniform and non-uniform PLUTO grids** — automatically detects grid uniformity and uses binary search or direct index computation accordingly
- **Multiple electron distribution functions** — thermal, kappa, power-law, and hybrid models informed by PIC simulations
- **Nonthermal electron physics** from state-of-the-art particle-in-cell (PIC) studies (turbulence and magnetic reconnection)
- **Polarized radiative transfer** — full Stokes I, Q, U, V with thermal and nonthermal emission, absorption, and Faraday effects

---

## 2. Dependencies

- HDF5 (parallel `h5cc` compiler)
- GSL (GNU Scientific Library)
- OpenMP

Compilation uses the `makefile`. Set environment variables as needed for your system.

---

## 3. `definitions.h` — Compile-time Parameters

This file controls the physics, metric, and numerical configuration. The key parameter groups are described below.

### 3.1 Distribution Function Choices — `DF`

The macro `DF` selects the electron distribution function used for computing emission and absorption coefficients. Five options are available:

| Macro           | Value | Description |
|:----------------|:------|:------------|
| `TH`            | 1     | **Thermal (Maxwell-Jüttner) distribution** |
| `KAPPA`         | 0     | **Kappa distribution** (single fixed κ) |
| `POWER`         | 2     | **Power-law distribution** (fixed spectral index and Lorentz factor range) |
| `VAR_KAPPA`     | 3     | **Variable kappa** — thermal + kappa hybrid, with nonthermal fraction determined by PIC-based efficiency |
| `VAR_POWER`     | 4     | **Variable power-law** — thermal + power-law hybrid, with nonthermal fraction determined by PIC-based efficiency |

The active distribution is set via:
```c
#define DF (TH)   // change to KAPPA, POWER, VAR_KAPPA, or VAR_POWER
```

#### Thermal Distribution (`TH`)

The relativistic Maxwell-Jüttner distribution for electrons at temperature $`\theta_e = k_B T_e / (m_e c^2)`$:

$f(\gamma) \propto \gamma^2 \beta \, e^{-\gamma/\theta_e}$


Emission and absorption coefficients use the fitting functions from Leung et al. (2011) and Pandya et al. (2016).

#### Kappa Distribution (`KAPPA`)

The κ-distribution describes a thermal core with a power-law tail:

$f(\gamma) \propto \frac{\gamma(\gamma^2 - 1)^{1/2}}{(1 + \gamma/\kappa w)^{\kappa+1}}$


where $`w = (\kappa-3)\theta_e/\kappa`$. The κ-index is set by `kappa_const` (default 4.0). Emission and absorption use the fitting functions from Pandya et al. (2016) and Davelaar et al. (2019).

#### Power-Law Distribution (`POWER`)

A pure power-law distribution between $`\gamma_{\rm min}`$ and $`\gamma_{\rm max}`$:

$f(\gamma) \propto \gamma^{-p}, \quad \gamma_{\rm min} \leq \gamma \leq \gamma_{\rm max}$


where $`p`$ = `Power` (default 2.3), $`\gamma_{\rm min}`$ = `Gamma_min` (default 1e3), $`\gamma_{\rm max}`$ = `Gamma_max` (default 1e6).

#### Variable Kappa (`VAR_KAPPA`)

Hybrid model: a fraction $`1-\mathcal{E}`$ of electrons follow a thermal distribution and a fraction $`\mathcal{E}`$ follow a κ-distribution. The κ index and nonthermal efficiency $`\mathcal{E}`$ are computed from the local plasma parameters $`\beta`$ and $`\sigma`$ using PIC simulation fits.

#### Variable Power-Law (`VAR_POWER`)

Hybrid model: thermal + power-law, with:

- $`\mathcal{E}`$ = nonthermal efficiency (fraction of electrons accelerated)
- $`\eta`$ = partition fraction of nonthermal electrons (including energy conservation)
- `Power` = power-law index $`p`$
- `Gamma_min`, `Gamma_max` = Lorentz factor range
- `Power_type` = `cooling` or `nth` — controls whether a cooling break is included (break Lorentz factor $`\gamma_{\rm br} \propto B^{-2}`$)

[warning: the $\gamma_{br}$ should be below 10 to guarantee the validity of the coefficients]

### 3.2 Nonthermal Mechanism — `Epsilon`

The macro `Epsilon` selects which PIC-based fitting function to use for the nonthermal efficiency $`\mathcal{E}`$, the κ index or power-law index $`p`$:

| Mechanism       | Value | Reference |
|:----------------|:------|:----------|
| `turbulence`    | 1     | Meringolo et al. 2023, ApJ 944, 122 — [DOI: 10.3847/1538-4357/acaefe](https://doi.org/10.3847/1538-4357/acaefe) |
| `reconnection`  | 0     | Ball et al. 2018, ApJ 862, 80 — [DOI: 10.3847/1538-4357/aac820](https://doi.org/10.3847/1538-4357/aac820) |

The fitting functions provide $`\mathcal{E}(\beta, \sigma)`$ and $`\kappa(\beta, \sigma)`$ (or $`p(\beta, \sigma)`$) in the trans-relativistic regime. For the **turbulence** model, these come from 2D PIC simulations of relativistic turbulence (Meringolo+2023); for the **reconnection** model, from PIC simulations of trans-relativistic magnetic reconnection (Ball+2018).

**Turbulence model (Meringolo+2023)** — applicable for $0.1 < \sigma < 10$, $10^{-4} < \beta < 2$:

Efficiency:
$\mathcal{E}(\sigma, \beta) = e_0 + \frac{e_1}{\sqrt{\sigma}} + e_2 \sigma^{1/10} \tanh\!\big(e_3 \beta \sigma^{1/10}\big)$
with $e_0 = 1$, $e_1 = -0.23$, $e_2 = 0.5$, $e_3 = -10.18$.

Kappa index:
$\kappa(\sigma, \beta) = a_0 + \frac{a_1}{\sqrt{\sigma}} + a_2 \sigma^{-0.6} \tanh\!\big(a_3 \beta \sigma^{0.1}\big)$
with $a_0 = 2.8$, $a_1 = 0.2$, $a_2 = 1.6$, $a_3 = 2.25$. Clipped to $3.1 \leq \kappa \leq 7.5$.

Note: The power-law index `Power` is user-specified (default 2.3) and not computed from PIC fits for the turbulence case with `VAR_POWER`.

**Reconnection model (Ball+2018)** — applicable for $0.01 < \sigma < 7.2$, $10^{-4} < \beta < 2.5$:

Efficiency:
$\mathcal{E}(\sigma, \beta) = e_0 + \frac{e_1}{4.2\sigma^{0.55} + 1} + e_2 \sigma^{0.07} \tanh\!\big(e_3 \beta \sigma^{0.13}\big)$
with $e_0 = 1$, $e_1 = -1$, $e_2 = 0.64$, $e_3 = -68$.

Kappa index:
$\kappa(\sigma, \beta) = a_0 + \frac{a_1}{\sqrt{\sigma}} + a_2 \sigma^{-0.19} \tanh\!\big(a_3 \beta \sigma^{0.26}\big)$
with $a_0 = 2.8$, $a_1 = 0.7$, $a_2 = 3.7$, $a_3 = 23.4$. Clipped to $3.1 \leq \kappa \leq 7.5$.

Power-law index for reconnection (`VAR_POWER`):
$p(\sigma, \beta) = a_0 + \frac{a_1}{\sqrt{\sigma}} + a_2 \sigma^{-0.19} \tanh\!\big(a_3 \beta \sigma^{0.26}\big)$
Clipped to $2.1 \leq p \leq 4.6$.

### 3.3 Other Important Parameters

| Parameter          | Description |
|:-------------------|:------------|
| `IMGFILE` | whether generate 2D data for plotting image(defualt:1) |
| `SPECFILE` | whether save all stokes parameters for every frequency(defualt:1) |
| `POL` | whether calculate polarization(defualt:1) |
| `num_pixels_1d`    | Number of pixels per dimension (default 20) |
| `num_frequencies`  | Number of frequencies (default 1) |
| `FREQFILE`/`FREQLOG` | Frequency sampling mode: file-based or logarithmic |
| `rcam`             | Camera distance from the origin (in $R_g$) |
| `int_method`       | Integration method: `RK2`, `RK4`, `VER`, `RK45` |
| `max_steps`        | Maximum number of integration steps |
| `RT_OUTER_CUTOFF`  | Outer boundary for radiative transfer |
| `EPS_RTE`          | Limitation for adaptive optical depth stepping |

---

## 4. `model_jet.in` — Run-time Parameters

This input file specifies the model and camera parameters without requiring recompilation.

### Example
```
MBH             (Msun)   6.5e9
DISTANCE        (kpc)    16.8e3
M_UNIT          (g)      1.e26
R_LOW           (-)      1
R_HIGH          (-)      5

INCLINATION     (deg)    90
azmuth          (deg)    0
IMG_WIDTH       (pixels) 500
IMG_HEIGHT      (pixels) 600
SHIFT           (Rg)     0
CAM_SIZE_X      (Rg)     10
CAM_SIZE_Y      (Rg)     60
FREQS_PER_DEC   (-)      1
FREQ            (Hz)     230e9
STEPSIZE        (-)      0.02
STEPSIZE_MAX    (-)      0.02
STEPSIZE_MIN    (-)      0.01
```

### Parameter Descriptions

#### Black Hole & Physical Parameters

| Parameter     | Unit    | Description |
|:--------------|:--------|:------------|
| `MBH`         | $M_\odot$ | Black hole mass in solar masses |
| `DISTANCE`    | kpc     | Distance to the source |
| `M_UNIT`      | g       | Mass unit used to set density and magnetic field scaling: $`\rho_{\rm unit} = M_{\rm unit} / L_{\rm unit}^3`$, $`B_{\rm unit} = c \sqrt{4\pi \rho_{\rm unit}}`$ |
| `R_LOW`       | —       | Low plasma-β temperature ratio parameter (Moscibrodzka et al. 2016) |
| `R_HIGH`      | —       | High plasma-β temperature ratio parameter (Moscibrodzka et al. 2016) |

The electron–proton temperature ratio is given by:

$$
\frac{T_p}{T_e} = R_{\rm high} \frac{\beta^2}{1+\beta^2} + R_{\rm low} \frac{1}{1+\beta^2}
$$
The electron temperature is then $`\theta_e = \frac{m_p}{m_e} \frac{u}{\rho} / (1 + T_p/T_e)`$ where $u/\rho$ is the specific internal energy.

#### Camera Parameters

| Parameter      | Unit    | Description |
|:---------------|:--------|:------------|
| `INCLINATION`  | deg     | Observer inclination angle |
| `azmuth`       | deg     | Observer azimuthal angle |
| `IMG_WIDTH`    | pixels  | Number of pixels in the x-direction |
| `IMG_HEIGHT`   | pixels  | Number of pixels in the y-direction |
| `CAM_SIZE_X`   | $R_g$   | Field of view in the x-direction |
| `CAM_SIZE_Y`   | $R_g$   | Field of view in the y-direction |
| `SHIFT`        | $R_g$   | Center offset in the x3-direction (applied after grid read-in) |

#### Numerical Parameters

| Parameter       | Unit  | Description |
|:----------------|:------|:------------|
| `FREQS_PER_DEC`  | —     | Number of frequency points per logarithmic decade (for multi-frequency runs) |
| `FREQ`          | Hz    | Observation frequency |
| `STEPSIZE`      | —     | Base integration step size (in $R_g$) |
| `STEPSIZE_MAX`  | —     | Maximum allowed step size |
| `STEPSIZE_MIN`  | —     | Minimum allowed step size |

---

## 5. `model.c` — Model Initialization and Data Handling

### 5.1 `init_model()`

The function `init_model()` in `model.c` selects the data reading routine. Four options are available within `init_model()`:

```c
void init_model() {
    set_units(M_UNIT);

    fprintf(stderr, "\nStarting read in of PLUTO RMHD data...\n");

    init_rmhd_data(RMHD_FILE);                    // (1) Standard 3D data
    //init_axis_data(RMHD_FILE);                   // (2) z-axisymmetric (mirror across z=0)
    //init_trace_data(RMHD_FILE);                  // (3) Includes tracer to mask ambient
    //init_axis_trace_data(RMHD_FILE);             // (4) z-axisymmetric + tracer
}
```

| Function                  | File       | Description |
|:--------------------------|:-----------|:------------|
| `init_rmhd_data()`        | `model.c`  | Standard 3D PLUTO data read. Reads all 8 primitive variables (ρ, u, v¹, v², v³, B¹, B², B³) from a `.dbl` file. Internal energy is converted: $u_{\rm int} = u / (\gamma_{\rm ad} - 1)$ with $\gamma_{\rm ad} = 4/3$. |
| `init_axis_data()`        | `pluto.c`  | Same as above, but mirrors the data across the z=0 plane (doubling N3). Velocity and magnetic field components are sign-flipped appropriately for symmetry. |
| `init_trace_data()`       | `pluto.c`  | Includes an additional tracer variable from the `.dbl` file. The tracer is used to mask the ambient (inactive) region to zero density and internal energy. |
| `init_axis_trace_data()`  | `pluto.c`  | Combines z-axisymmetry with tracer masking. |

### 5.2 `grid.out` Path

The `grid.out` file contains the cell interface positions for the PLUTO grid. **All four init functions read `grid.out` from a hard-coded relative path**:

Adjust these paths as needed for your directory structure. The `grid.out` format follows the standard PLUTO output: comments starting with `#`, then lines specifying the number of cells N1, N2, N3, followed by each cell's (index, left_edge, right_edge).

### 5.3 Unit System

The physical units are set in `set_units()`:

$$
\begin{aligned}
L_{\rm unit} &= G M_{\rm BH} / c^2 \quad [\text{cm}] \\
T_{\rm unit} &= L_{\rm unit} / c \quad [\text{s}] \\
\rho_{\rm unit} &= M_{\rm unit} / L_{\rm unit}^3 \quad [\text{g/cm}^3] \\
U_{\rm unit} &= \rho_{\rm unit} c^2 \quad [\text{erg/cm}^3] \\
B_{\rm unit} &= c \sqrt{4\pi \rho_{\rm unit}} \quad [\text{Gauss}] \\
n_{e,\rm unit} &= \rho_{\rm unit} / (m_p + m_e) \quad [\text{cm}^{-3}]
\end{aligned}
$$

## 6. How to Run

### 6.1 Compilation

```bash
cd /path/to/RaptorP
make
```

This produces the executable `RAPTOR` in the current directory.

### 6.2 Running the Code

```bash
./RAPTOR model_jet.in <path/to/data.dbl> <output_index>
```

**Arguments:**
1. `model_jet.in` — Parameter file (see Section 4)
2. `<path/to/data.dbl>` — Path to the PLUTO binary data file (e.g., `../df/data.0096.dbl`)
3. `<output_index>` — Integer index appended to output filenames

### 6.3 Example — `run.sh`

```bash
#!/bin/bash
#export OMP_NUM_THREADS=1
./RAPTOR model_jet.in ../df/data.0096.dbl 40
```

This runs the code with:
- Parameter file: `model_jet.in`
- Data file: `../df/data.0096.dbl`
- Output index: 40

The `model_jet.in` file should be in the same directory as the executable. The output consists of HDF5 files containing the Stokes parameter maps and integrated spectra.

### 6.4 Visualization

The repository includes several Python scripts for analysis:

- `rapplot.py` — General-purpose plotting library (adapted from RAPTOR)
- `plot.py` — Quick plotting script

## 7. Citation

If you use RaptorP in your research, please cite the original RAPTOR papers and the application paper:

**RaptorP application:**
> Hu, X., et al. 2025, ApJ, 995, 76
> [https://ui.adsabs.harvard.edu/abs/2025ApJ...995...76H/abstract](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...76H/abstract)

**RAPTOR code papers:**

- Bronzwaer, Davelaar, et al. 2018, A&A, 613, A2
- Davelaar, Bronzwaer, et al. 2018, CompAC, 5, 1
- Bronzwaer, Younsi, Davelaar, et al. 2020, A&A, 641, A126
- Davelaar & Haiman 2022, PRD, 105, 10, 103010

**PIC model references:**
- **Turbulence**: Meringolo, C., Cruz-Osorio, A., Rezzolla, L., & Servidio, S. 2023, ApJ, 944, 122 — [DOI: 10.3847/1538-4357/acaefe](https://doi.org/10.3847/1538-4357/acaefe)
- **Reconnection**: Ball, D., Sironi, L., & Özel, F. 2018, ApJ, 862, 80 — [DOI: 10.3847/1538-4357/aac820](https://doi.org/10.3847/1538-4357/aac820)

---

## 8. License

RaptorP is a derivative of RAPTOR (Radboud Polarized Integrator), which is released under the GNU General Public License. The original RAPTOR code is copyright 2014–2021 Black Hole Cam (ERC Synergy Grant), authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi. PLUTO-specific modifications and additions by Xufan Hu (2025).
