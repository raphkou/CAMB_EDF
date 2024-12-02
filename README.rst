===================
CAMB_EDF
===================
This is a modified version of CAMB (https://github.com/cmbant/CAMB) implementing the Early Dark Fluid (EDF) model (https://arxiv.org/abs/2410.16185).

See main CAMB repository for installation.

Usage
=============================

The EDF model consists in the addition of 6 parameters to the standard LCDM model:

 - 4 parameters for the amplitudes of the density modes.
 - 2 sound speed parameters.

Those parameters can be used when calling the set_params function from CAMB by adding the following fields:

 - amp_delta : list of four amplitudes
 - cs2_1 : $c_1^2$
 - cs2_2 : $c_2^2$

Note this implementation assumes LCDM for the late-time expansion of the Universe.

Citation
=============================

If you use this code, please reference:
 - Flexible parameterization to test early physics solutions to the Hubble tension with future CMB data: https://arxiv.org/abs/2410.16185
 - Efficient Computation of Cosmic Microwave Background Anisotropies in Closed Friedmann-Robertson-Walker Models: https://arxiv.org/abs/astro-ph/9911177
 - CMB power spectrum parameter degeneracies in the era of precision cosmology: https://arxiv.org/abs/1201.3654
