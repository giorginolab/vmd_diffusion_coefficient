# VMD Diffusion Coefficient Tool - demo and regression test

This directory contains a script intended both as 
a demonstration of use, and as a regression test.

It can be launched with

    vmd -dispdev none regtest.tcl

On first run, the script will download ~250MB of simulation
results part of the [case
study](https://sourceforge.net/projects/membplugin/files/Case%20Study/CaseStudy_v2.pdf/download)
for the [MEMBPLUGIN
tool](https://academic.oup.com/bioinformatics/article/30/10/1478/266901#supplementary-data).

The simulation data corresponds to 20 ns of simulation of a 1:1 molar
ratio POPC:Cholesterol bilayer. The bilayer was modeled at all-atom
resolution with the CHARMM36 forcefield.  The details of the system
preparation are reported in the [case
study](https://sourceforge.net/projects/membplugin/files/Case%20Study/CaseStudy_v2.pdf/download).

The script computes lateral MSD displacement for time intervals
ranging from 2 to 10 ns for the POPC and Cholesterol species
separately. Diffusion coefficients can be determined from *asymptotic
slopes* of the curves (usually fitting linear models).

