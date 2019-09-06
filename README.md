## Radioactive Heating Rate and Macronova (Kilonova) Light Curve
### Authors: Kenta Hotokezaka and Ehud Nakar


- decay_chain_thermal_beta.ipynb      ~~~(beta-decay only)
- decay_chain_thermal_alpha.ipynb     ~~~(beta and alpha decay)
- decay_chain_theraml_fission.ipynb   (beta and spontaneous fission--currently only 254Cf--)

work with Jupyter Notebook. These codes calculate nuclear heating rates [erg/s/g]  of beta-decay, alpha-decay, and spontaneous fission of r-process nuclei, taking into account for thermalization of gamma-rays and charged decay products in r-process ejecta. The codes use the half-lives and injection energy spectra from an evaluated nuclear data library (ENDF/B-VII.1). The heating rate is computed for given abundances, ejecta mass, velocity, and density profile.

The codes also compute the bolometric light curve and the evolution of the effective temperature for given abundances, ejecta mass, velocity, and density profile assuming opacities independent of the wavelength.


These codes support both Python 2 and Python 3.


### Documentation
See [here](http://github.com/hotokezaka/HeatingRate) for details.
