# ElSeisHom3D

This repository contains the code that comes with the publication

> Slob, E. and M. Mulder, 2015,  
> Seismoelectromagnetic homogeneous space Green's functions:  
> Geophysics, 81(4), F27-F40;
> DOI: [10.1190/geo2015-0337.1](https://doi.org/10.1190/geo2015-0337.1).

The original code was published with the article in a static archive by the
SEG, see
[software.seg.org/2016/0003](https://software.seg.org/2016/0003/index.html).
The version in the SEG archive is (c) 2016 by the Society of Exploration
Geophysicists, for more info consult the file
[software.seg.org/disclaimer.txt](http://software.seg.org/disclaimer.txt).

We[^1] release our source code here under the CC0-1.0 license, see the file
`LICENSE`. Please cite the above article if you use the code in your research.


# Manual

ElSeisHom3D stands for computing the Green's functions of electromagnetic
fields and elastic waves that are coupled in a homogeneous porous medium. The
fields are the electric field, the magnetic field, the particle velocity, the
filtration velocity, the acoustic pressure in the fluid, and bulk stress. The
sources that generate these fields are the electric current source, the
magnetic current source, the force source acting on the bulk, the force source
acting on the fluid, the acoustic volume injection rate source acting on the
fluid, and the deformation rate source acting on the bulk. The impulse
responses that connect the fields to these sources are the Green's functions.

The expressions on which these codes are based are given in Slob and Mulder
(2015). The codes are given in wavenumber-frequency domain and in space-time
domain. The space-time results are obtained by coding the space-frequency
domain expressions given in Slob and Mulder (2015), the results of which are
then transformed to space-time domain using a Fast Inverse Fourier
Transformation routine. All codes are given in Matlab and the algorithms can be
directly transferred to any other computer code language.

## Folder structure

The folder is called `3Dhom`. Inside this folder two folders are present. One
is called `3Dhomkf`, which contains the codes in wavenumber-frequency domain.
These codes are simply meant to demonstrate correctness of the obtained
coefficients and fields generated by a single component from each source. This
is achieved by demonstrating that the given fields satisfy the basic equations
through showing the numerical precision with which these are satisfied. These
codes serve to validate the Green's functions expressions. The second folder is
called `3Dhomxt`, which contains the codes to compute the space-time domain
impulse responses.

## Folder content `3Dhomkf`

This folder contains Matlab scripts in the main six files. For each source the six basic equations are tested for the Green's functions for a single source component.

The scripts are (`File name :` source type):

- `ElSeishom3Delcurrent.m       :` electric current source
- `ElSeishom3Dmagcurrent.m      :` magnetic current source
- `ElSeishom3Dforcebulk.m       :` force source acting on the bulk
- `ElSeishom3Dforcefluid.m      :` force source acting on the fluid
- `ElSeishom3Dvolumeinjection.m :` volume injection rate acting on fluid
- `ElSeishom3Ddefrate.m         :` deformation rate acting on the bulk

In each of these script files first the parameters are defined and the
parameters are described in the paper. The first is frequency, which is a
single number in these codes, and then all medium parameters are computed. Each
file starts in line 67 with computing the spherical wavenumbers and two
coefficients that are used in all Green's functions.

Then the wavenumber vector for kx, ky, and kz is defined, with 2D arrays
containing all horizontal wavenumbers kx and ky and a loop to run over all
kz-values is coded. In the code and in the paper subscript notation is used so
you will find k1=kx, k2=ky, and k3=kz. The radial wavenumber is computed
because all coefficients only depend on that parameter. All wavenumber tensor
components, scalar Green’s functions, and coefficients are computed for all
fields generated by a specific source having the x-component. The relations to
the equations in the paper are given in the codes. Once these are computed the
fields are computed and inserted in the expression of the six basic equations.
For each equation a plot is generated to show the result, which should always
come out as zero. Each plot shows the results of the components of each
equation and contains nine subplots, each as documented in the codes. This
demonstrates the precision with which these basic equations are satisfied.

## Folder content `3Dhomxt`

This folder contains the Matlab scripts in six main files. For each source the Green's functions are computed for the six fields. In these codes all field vector/tensor components are computed for all source vector/tensor components. In the paper four function files are present, each of which contains gradients to the Green's functions from a single gradient to four gradients.

The main scripts are (`File name :` source type):

- `ElSeishom3Delcurrentxt.m       :` electric current source
- `ElSeishom3Dmagcurrentxt.m      :` magnetic current source
- `ElSeishom3Dforcebulkxt.m       :` force source acting on the bulk
- `ElSeishom3Dforcefluidxt.m      :` force source acting on the fluid
- `ElSeishom3Dvolumeinjectionxt.m :` volume injection rate acting on fluid
- `ElSeishom3Ddefratext.m         :` deformation rate acting on the bulk

The function scripts are (`File name :` number of gradients):

- `greeni.m    :` one
- `greenij.m   :` two
- `greenijk.m  :` three
- `greenijkl.m :` four

The scripts to compute the figure plots in the paper are (`File name :`
frequency band):

- `ElSeishom3Delcurrentxtseparatedwaves.m   :` seismic
- `ElSeishom3Dmagcurrentxtseparatedwaves.m  :` seismic
- `ElSeishom3DelcurrentxtseparatedwavesM.m  :` ultrasonic
- `ElSeishom3DmagcurrentxtseparatedwavesM.m :` utrasonic


## The main scripts

In each of these script files first the parameters are defined and the
parameters are described in the paper. The first is frequency, which is a range
of numbers that allow for transforming all results to time domain with an
ifft-routine. Then all medium parameters are computed. For time domain results
we use a Ricker wavelet as source-time function and it is computed as the
source signature in the frequency domain just before the spherical wavenumbers
are computed. The structure is similar to the script files for the
wavenumber-frequency computations. Now a three-dimensional matrix is defined
containing information on all x1-, x2-, and frequency-components. On this grid
all Green's functions are computed. Depending on the field-source combinations
several functions are called that compute one or more gradients of the scalar
Green's functions and these are all referred to the equation numbers in the
paper. Then the tensor-components of the Green's functions for each of the
fields are computed. Some of these are multiplied with source signature and
transformed back to time domain. Of these the script plot one field component
to show the result. These are not shown in the paper.

## The functions computing the gradients on the Green’s functions

These function files compute one or more gradients as explained in the function- script files and reference is made to the equations in the paper.

## The script-files used to make the plots in the paper

These files have the same structure as the main script files and compute for
several fields the contributions from the fast- and slow-P waves, the S-wave,
and the electromagnetic field as separate arrays. The contributions from these
four wave types are shown as snapshots in time for a receiver at a certain
vertical distance to the source for all computed horizontal distances in a
plane. This is done for the particle velocity and the filtration velocity
generated by electric and magnetic current source. This is done at seismic
frequencies with a Ricker wavelet having a center frequency of 40 Hz and at
ultrasonic frequencies with a Ricker wavelet having a center frequency of 400
kHz.

## The medium parameters

The medium parameters that occur in all script files that compute Green's
functions are detailed here.

- `rhof`: mass density of fluid
- `rho`: mass density of solid
- `Gfr`: shear modulus of the framework of grains
- `eta`: fluid viscosity
- `k0`: medium permeability (static)
- `Kfr`: modulus of skeleton grains
- `Ks`: modulus of the framework of grains
- `Kf`: modulus of fluid
- `Concentr`: ionic concentration of fluid
- `bplus`: ionic mobility of anions
- `bmin`: ionic mobility of cations
- `porosity`: porosity
- `epsilonRF`: relative permittivity of fluid
- `epsilonRS`: relative permittivity of solid
- `alpha_inf`: tortuosity
- `pH`: acidity
- `similaritypar`: 8;
- `c0`: velocity of light in free-space
- `mu0`: free-space magnetic permeability
- `epsilon0`: free-space electric permittivity
- `e`: elementary charge
- `z_1`: ion valence
- `z_1c`: valence of the conjugate ion
- `NA`: Avogadro's constant
- `kb`: Boltzmann constant
- `T`: Temperature in Kelvin
- `epsilonR`: bulk relative permittivity
- `omegac`: critical frequency
- `zetap`: zeta potential
- `L0`: static coupling coefficient
- `L`: dynamic coupling coefficient
- `N`: bulk-ionic concentration of species i
- `sigmaF`: conductivity of the pore-fluid phase
- `sigmaE`: bulk electric conductivity (freq-(in)dependent)
- `rhoB`: effective density of the fluid
- `epsilon`: bulk absolute permittivity
- `Delta`: combination of the compression moduli
- `Kg`: Gassmann's bulk modulus
- `C`: bulk compression modulus
- `M`: fluid compression modulus
- `H`: bulk modulus for P-wave
- `S`: combination of moduli, not used
- `sigmaM`: magnetic conductivity
- `Kc`: combination of moduli; not used

Frequency dependent parameters

- `k`: dynamic permeability
- `L`: dynamic coupling coefficient
- `rhoE`: effective density
- `rhoc`: complex density
- `zeta`: generalized magnetic conductivity
- `etae`: generalized electric conductivity
- `varsigma`: generalized effective electric conductivity
- `chi`: parameter product
- `ypf`: spherical wavenumber of fast P-wave
- `yps`: spherical wavenumber of slow P-wave
- `ys`: spherical wavenumber of seismic S-wave
- `yem`: spherical wavenumber of diffusive EM-wave
- `dems`: scaling factor of transversal Green's function
- `dpspf`: scaling factor of longitudinal Green's function

[^1]: 2022-12-09: The code was uploaded to
  [github.com/emsig/3Dhom](https://github.com/emsig/3Dhom) by
  [@prisae](https://github.com/prisae),
  upon the request of the main author Evert Slob.
