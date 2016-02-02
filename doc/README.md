Diffusion Coefficient Tool (beta)
=================================

Toni Giorgino

-   2012- (ISIB-CNR) Consiglio Nazionale delle Ricerche

    2009-2011 (GRIB-IMIM-UPF) Computational Biophysics Laboratory

The [Diffusion Coefficient Tool](#) is an analysis plugin for VMD that
computes one, two or three-dimensional mean squared displacements
(MSD)-based diffusion coefficients of a chosen molecular species.

This is beta software. If you don't know how to compute the diffusion
coefficient yourself, probably you won't be able to use this software.
**NO support or warranty is available whatsoever.**

**Repeat.** The Diffusion Coefficient Tool is not a magic box. It only
computes mean squared displacements (MSD) at a variety of lag times (τ).
The MSD is predicted by Einstein's relation to grow linearly with τ.
Most often than not, this DOES NOT occur in practice, for a variety of
reasons, including poor sampling or non-diffusive behavior. It is YOUR
responsibility to understand why and, if you trust the linearity, to fit
the slope. Please consider the following paper to be **mandatory
reading**: [David Keffer, The Working Man's Guide to Obtaining Self
Diffusion Coefficients from Molecular Dynamics
Simulations](http://utkstair.org/clausius/docs/che548/pdf/selfD.pdf).

Installation
------------

This is beta software. By downloading the software you agree to comply
with the terms of GPL version 2.

**NOTE: This code comes without any warranty of fitness for any use. It
is UNSUPPORTED. After you download it, you are on your own.**

[Download](https://github.com/tonigi/vmd_diffusion_coefficient/releases/latest)
and follow the instructions therein.

You may or may not be able to get support by posting to the [External
tools](https://sourceforge.net/p/membplugin/discussion/external_tools/)
forum hosted at the [MEMBPLUGIN](http://membplugin.sourceforge.net)
site.

Citation and license
--------------------

Please cite the following publication:

Toni Giorgino, Computing diffusion coefficients in macromolecular
simulations: the Diffusion Coefficient Tool for VMD, Submitted (2014).
Available from
[](http://multiscalelab.org/utilities/DiffusionCoefficientTool)
(accessed xxx, xxx 2013).

Usage (GUI)
-----------

The plugin is accessible from VMD in *Extensions \> Analysis \>
Diffusion Coefficient Tool*. The profile is computed for the currently
loaded trajectory in the *top* molecule. The atom selection must match
one atom per diffusing molecule; MSD values will be averaged over the
matched atoms.

Two buttons are provided to plot the **MSD displacement** MSD(τ) and the
**Diffusion coefficient** D(τ). The latter is computed with the Einstein
relation as *D(τ)=MSD(τ)/2Eτ*, where *E* is the integer dimensionality
of the system (1, 2 or 3). (Note that it's not a linear fit.)

The **Analysis interval** boxes allow to specify a subsection of the
trajectory to be used for the MSD calculation. Only displacements in the
subspace spanned by the **Diffusion along** axes will be considered. In
other words, if *z* is deselected, motions along that axes will be
irrelevant for the calculation of MSD and D.

The results are given in Å^2^/ns. Conversion factors:

-   0.1 Å^2^/ns = 10^-12^ m^2^/s = 10^-8^ cm^2^/s = 1 μm^2^/s

-   1 cm²/s = 10^7^ Å²/ns

-   1 m²/s = 10^11^ Å²/ns

-   1 μm²/s = 0.1 Å²/ns

-   1 Å²/ns = 10^-7^ cm²/s = 10^-11^ m²/s = 10 μm²/s

Upon completion, the location of a data file with the computed profiles
will be printed on the console.

**Warnings**:

-   **Don't draw premature conclusions.** Computing converged diffusion
    coefficients is likely beyond your sampling capacity and patience,
    unless done for plenty of molecules (ie. solvent). Expect
    **microseconds** sampling.

-   The trajectory must not be wrapped.

-   A check is performed whether bonds exist between the selected atoms.
    This check may fail if wrong connectivity is inferred on molecule
    load;

Usage (command line)
--------------------

You need to load the package with
`package require diffusion_coefficient`. The plugin can perform the
following computations

-   MSD at a given lag time (in frames). The return value will be in Å,
    and its value is independent of the time units. This is the clearest
    (i.e. recommended) way to use the plugin.

-   MSD for an interval of lag times

-   D, computed as above, for an interval of lag times.

Invocation is self-explanatory, i.e.:

    VMD Diffusion Coefficient tool. Computes one, two or three-dimensional
    MSD-based diffusion coefficients of a chosen molecular species.

    Usage: diffusion_coefficient <options> <command>

    Command is one of:
    -msd <NN>     Compute mean squared displacement (MSD) at a tau of
                  NN frames; equivalent to msd_interval -from NN -to NN.
                  Returns a value as Angstrom^2 . This is the recommended 
                  way of using the plugin.
    -msd range    Compute MSD for taus between -from and -to (mandatory)
                  Returns two lists of {tau} {MSD(tau)}
    -d range      Compute D(tau)=MSD(tau)/(2*D*tau) between -from and
                  -to (mandatory). Returns two lists of {tau} {D(tau)}

Averaging
---------

![formulas.png](http://www1.multiscalelab.org//utilities/DiffusionCoefficientTool?action=AttachFile&do=get&target=formulas.png)

  ------------------ ------------------------ -------------------------
  **Symbol above**   **GUI parameter name**   **Command line option**
                                              

  w~f~               Analysis interval from   `-interval_from`
                                              

  w~t~               Analysis interval to     `-interval_to`
                                              

  w~s~               Analysis interval step   `-interval_stride`
                                              
  ------------------ ------------------------ -------------------------

Screenshot
----------

![gui.png](http://www1.multiscalelab.org//utilities/DiffusionCoefficientTool?action=AttachFile&do=get&target=gui.png)

Acknowledgments
---------------

Former support from the Agència de Gestió d'Ajuts Universitaris i de
Recerca - Generalitat de Catalunya is gratefully acknowledged.
