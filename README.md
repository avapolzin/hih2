# `hih2`
Atomic-to-molecular hydrogen transition models for astrophysical simulations.

To install:
```bash
cd ~

git clone https://github.com/avapolzin/hih2.git

cd hih2

pip install .
```
or 

```bash
pip install hih2
```

Included models are:
- [Krumholz, McKee, & Tumlinson (2009b)](https://ui.adsabs.harvard.edu/abs/2009ApJ...699..850K/abstract) -- KMT09b
- [Gnedin & Kravtsov (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJ...728...88G/abstract) -- GK11
- [Krumholz (2013)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2747K/abstract) -- K13
- [Sternberg et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...790...10S/abstract) -- S14
- [Gnedin & Draine (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...795...37G/abstract) -- GD14
- [Polzin et al. (2024a)](https://ui.adsabs.harvard.edu/abs/2024ApJ...966..172P/abstract) -- P24

Volumetric models are in `hih2.vol` and projected models are in `hih2.proj`. All syntax is the same between those modules. It's very simple, so the docstrings are largely sufficient and limited documentation follows here:

All models take as arguments `nh` (hydrogen number density if `hih2.vol` -- cm<sup>-3</sup> by default -- or hydrogen column density if `hih2.proj` -- cm<sup>-2</sup> by default) and `met` (solar-scaled metallicity). `gk11()` and `p24()` require `uv` (UV field strength). All of the functions in `hih2.proj` also require an `scale` argument, is the size of the grid cells (by default in cm). `hih2.vol.gd14()`, `kmt09b`, `k13`, and `s14` also take a `scale` argument, which follows the convention of `hih2.proj`.

Both `k13()` and `s14()` functions also take a clumping factor argument, `fc`, and `k13()` additionally takes `rho_sd` (density of stars and dark matter, by default in M<sub>&#9737;</sub> pc<sup>-3</sup>), `iterate` (toggles whether values are computed iteratively), and `niter` (if `iterate = True`, sets number of interations). 

All functions take `dens_unit`, which takes `astropy.units` objects to set the units in which `nh` is reported. `k13()` has an equivalent for `rho_sd` called `sddens_unit`. `hih2.proj` and `hih2.vol.gd14()` take a `scale_unit` argument.

Each function returns the molecular hydrogen fraction (defined as H<sub>2</sub>/(HI + H<sub>2</sub>)) for the input gas density, metallicity, and UV field strength. The implicit assumption of these models is that all gas is neutral and that any ionized component is negligible.

**Why is this a package?** Good question. It's very easy to implement these HI-H<sub>2</sub> models yourself either from the literature or from code examples (see [avapolzin/hydrogen_models](https://github.com/avapolzin/hydrogen_models) for functions on their own), but some people may prefer a quick install/import. 

**Summary of required and optional parameters in models included in `hih2.vol`**
| Model | n<sub>H</sub> | metallicity | UV field | scale | f<sub>c</sub> | &rho;<sub>SD</sub> |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kmt09b` | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  |  $${\color{red}Required}$$| | | 
| `gk11` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  | | | 
| `k13` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{yellow}Default provided}$$ | $${\color{yellow}Default provided}$$| 
| `s14` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{yellow}Default provided}$$| | 
| `gd14` | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  | $${\color{red}Required}$$ | | | 
| `p24` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  | | | 

**Summary of required and optional parameters in models included in `hih2.proj`**
| Model | N<sub>H</sub> | metallicity | UV field | scale | f<sub>c</sub> | &rho;<sub>SD</sub> |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kmt09b` | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  |  $${\color{red}Required}$$| | | 
| `gk11` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  | $${\color{red}Required}$$ | | 
| `k13` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{yellow}Default provided}$$ | $${\color{yellow}Default provided}$$| 
| `s14` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{yellow}Default provided}$$| | 
| `gd14` | $${\color{red}Required}$$ | $${\color{red}Required}$$ |  | $${\color{red}Required}$$ | | | 
| `p24` | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | $${\color{red}Required}$$ | | | 

