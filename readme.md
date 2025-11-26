# RaptorP

**RaptorP** is a Special Relativity Radiation Transfer (SRRT) extension of the General Relativity Radiation Transfer (GRRT) code [RAPTOR](https://github.com/jordydavelaar/raptor). It adopts (Relativistic)-Magnetohydrodynamics ((R)MHD) simulated data from open-sourced numerical code [PLUTO](https://plutocode.ph.unito.it/) .

## Overview

**RaptorP** inherits most of beneficial functions from  [RAPTOR](https://github.com/jordydavelaar/raptor) . It takes a full form of radiation transfer formula, including synchrotron emissivity, absorption, and Faraday rotation. For electron distribution function (eDF), it supports Maxwell-J\"uttner(thermal), $\kappa$(hybrid), and broken power-law(non-thermal) distributions.

The other features are listed below:

(1) Reading and processing both uniform and non-uniform **Cartesian** PLUTO grids.

(2) Supporting a tracer to exclude the ambient medium. From one-sided jet simulation data, it can generate a symmetric two-sided jet.

(3) Estimating the nonthermal electron fraction from PIC-informed  prescriptions, e.g. magnetic reconnection ([Ball et al. 2018](https://iopscience.iop.org/article/10.3847/1538-4357/aac820)) , turbulence ([Meringolo et al. 2023](https://iopscience.iop.org/article/10.3847/1538-4357/acaefe)); and proposing a self-consistent method for evaluating 'w' parameter in the $\kappa$ distribution (see the Appendix.B of the code paper for details).

(4) Using a rectangle camera.

(5) Shifting camera position to image locally.

(6) Recording light paths.

(7) The plotting package `rapplot.py` is reconstructed with Python class. It is easier to use. 

## Quick Start 

### Dependence

gsl hdf5 gnu mpi

### Running RaptorP

1. open `model.c` 

	you should modify the path to `grid.out` in `init_rmnd_data()`

2. make

3. run RaptorP

	```shell
	./RAPTOR model.in path2/data.00xx.dbl xx
	```

## Defining radiation transfer:`definitions.h`

open `definitions.h`, and you can set up output files, eDF, etc.

### output

`IMGFILE`  : output radiation image (.h5)

`SPECFILE` : output the Stokes parameters for each frequency

```
freq  I Q U V
```

### frequeny

`num_frequencies` : number of frequencies

(1) `FREQFILE` : use `frequencies.txt` to list multi-frequencies for calculation.

(2) `FREQLOG` : if you only consider one frequency, you just need to list it in `model.in`. 

`FREQS` : way to input frequency , `FREQFILE` or `FREQLOG`

#### electron distribution function(eDF)

`DF` controls the eDF. The choices are:

1. `TH` : Maxwell-J\"uttner distribution

	$\frac{dN}{d\gamma}=\frac{\gamma^2v}{c\Theta_eK_2(1/\Theta_e)}e^{-\gamma/\Theta_e}$

	with $K_2$ the modified Bessel function of the second kind, and $\Theta_e=\frac{m_ec^2}{kT}$ 

3. `KAPPA` : $\kappa$ distribution ([Xiao 2006](https://iopscience.iop.org/article/10.1088/0741-3335/48/2/003))

$\frac{d N}{d\gamma} \propto \gamma\sqrt{\gamma^2-1}(1+\frac{\gamma-1}{\kappa w})^{-(\kappa+1)}$

 Specifying $\kappa$ by setting `kappa_const` , and  $w$  is determined automatically

3. `POWER` :  broke powerlaw

	$\frac{dN}{d\gamma}\propto \gamma^{-\alpha}\ \ \ \  \gamma_{min}<\gamma<\gamma_{max}$

	`Power` : paw-law spectrum index

	`Gamma_min` : $\gamma_{min}$

	`Gamma_max` : $\gamma_{max}$

3. `VAR_KAPPA` :  Estimating the nonthermal energy fraction through PIC simulations

	(1) `reconnection` : magnetic reconnection  ([Ball et al. 2018](https://iopscience.iop.org/article/10.3847/1538-4357/aac820)) 

	(2)  `turbulence` : turbulence ([Meringolo et al. 2023](https://iopscience.iop.org/article/10.3847/1538-4357/acaefe))

	`Epsilon`:  mirco-physical mechanism, `reconnection` or `reconnection`

	$w$ is determined automatically

### marking radiation transfer zone

`LIGHTPATH` : record lightpaths in file `lightpath.txt`

`RT_OUTER_CUTOFF` :  Outer boundary of rad transfer computation. In RaptorP, the RT radius  is defined as  $\sqrt{x^2+y^2}$.

`rcam` : the location of the camera, usually far larger than the simulation scale

## Configuring the simulation model: `model.in`

`MBH` :  Black hole mass in $M_\odot$, used for setting the length unit ($r_g=\frac{GM}{c^2}$)

`DISTANCE` : distance to the source in kpc

`M_UNIT`  :  Mass scaling in gram

`R_HIGH`  :  parameter in $R_\beta$ model

`R_LOW`  :   parameter in $R_\beta$ model



`INCLINATION`  :  Observer inclination in degrees

`IMG_WIDTH`  :  Amount of pixels in the x direction

`IMG_HEIGHT`  :  Amount of pixels in the y direction

`SHIFT` : Shifting the image center along the z-axis. So you can image locally. 

`CAM_SIZE_X`  :  The FOV in the x direction in $r_g$

`CAM_SIZE_Y`  :  The FOV in the y direction in $r_g$

`FREQ`  :  Starting frequencies of the frequency array in Hz

`STEPSIZE`  :  Step size scaling, suggested value $\le0.01$

<span style='color:red;'>[WARNING] both the `STEPSIZE	` parameter and  `src/gr_integrator.c/stepsize()` algorithm can affect the ray-trace sensitively. You may need to adjust them according to your grid.</span>

`STEPSIZE_MAX` : the maximum step size with `RT_OUTER_CUTOFF`

`STEPSIZE_MIN` : the minimum step size with `RT_OUTER_CUTOFF`

## Augmenting simulation data

In `model.c/init_model()`, you can choice the other initial function for following purposes:

(1) ` init_axis_data()` : Generate data that is mirror-symmetrical relative to the XY plane. As an example, from one-sided jet simulation data, it can generate a symmetric two-sided jet.

(2) ` init_trace_data()`: use a tracer (from PLUTO) to exclude the ambient medium.

(3) `init_axis_trace_data()` : perform (1)&(2) simultaneously

## Plotting

`rapplot.py` provides a portable approach to plot . Please check comments in it carefully. `plot.py` is a simple example.

## Code paper

[Quasi-Periodic Polarized Emissions from Kink Structure in Magnetized Relativistic Jets]()

## Contact me


Xufan Hu （胡旭凡） tie@mail.ustc.edu.cn 


