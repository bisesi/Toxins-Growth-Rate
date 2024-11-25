# Capabilities of the COMETS platform

COMETS can:

**Simulate a range of spatial conditions**

* Microbial growth in both 2D and 3D worlds
* Diffusive and convective propagation of biomass in space
* Substrate-dependent nutrient and biomass propagation
* Boundary conditions

**Simulate a range of biologically realistic details**

* Continuous (chemostat) and batch growth modes
* Lag-phases
* Parsimonious dFBA
* Cell death
* Neutral population drift
* Evolutionary processes like random mutation (reaction deletions)

The software is implemented in JAVA and can be run via the command-line or in a graphical user interface that includes visualization tools. In multi-CPU systems, dFBA is parallelized in a multi-threaded process for greater computational performance. Toolboxes in both `MATLAB` and `Python` are available to modify the input files for COMETS in a programmatic way. 