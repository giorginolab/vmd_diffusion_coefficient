Important: please check the current version of software and 
documentation at its home-page:

http://www.multiscalelab.org/utilities/DiffusionCoefficientTool




Diffusion Coefficient Tool (alpha)

Toni Giorgino — © Universitat Pompeu Fabra 2011

The Diffusion Coefficient Tool is an analysis plugin for VMD that
computes one, two or three-dimensional MSD-based diffusion
coefficients of a chosen molecular species.

Usage

(Preliminary documentation)

The plugin is accessible from VMD in Extensions > Analysis > Diffusion
Coefficient Tool. The profile is computed for the currently loaded
trajectory in the top molecule. The atom selection must match one atom
per diffusing molecule; MSD values will be averaged over the matched
atoms.

Two buttons are provided to plot the MSD displacement and the
Diffusion coefficient. The latter is computed as 1/2D MSD(t)/t, where
D is the dimensionality of the system.

The Analysis interval boxes allow to specify a subsection of the
trajectory to be used for the MSD calculation. Only displacements in
the subspace spanned by the Diffusion along axes will be
considered. In other words, if z is deselected, motions along that
axes will be irrelevant for the calculation of MSD and D.

Upon completion, the location of a data file with the computed
profiles will be printed on the console.

Warnings:

    * the trajectory must not be wrapped;

    * a check is performed whether bonds exist between the selected
      atoms. This check may fail if wrong connectivity is inferred on
      molecule load.


Licensing

This is alpha software. By downloading the software you agree to
comply with the terms of GPL version 2. The Tool may be described in a
forthcoming publication. You are kindly requested check here if a
citation is available when publishing results obtained with the tool.

Acknowledgments

Work partially supported by the Generalitat de Catalunya. 

