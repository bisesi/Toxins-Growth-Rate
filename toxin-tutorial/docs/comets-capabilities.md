# Capabilities of the COMETS platform

COMETS can:

**Simulate a range of spatial conditions:**

* Microbial growth in 2D and 3D environments
* Diffusive and convective spatial propagation of biomass
* Substrate-dependent nutrient and biomass propagation
* Boundary conditions

**Simulate a range of other biologically realistic details:**

* Chemostat and batch growth modes
* Parsimonious dFBA
* Cell death
* Neutral population drift
* Evolutionary processes like random mutation (reaction deletions)

The software is implemented in JAVA and can be run via the command-line or in a graphical user interface. In multi-CPU systems, dFBA is parallelized for greater computational performance. Toolboxes in both `MATLAB` and `Python` are available to modify the input files for COMETS in a programmatic way. 

**With the introduction of the signaling class in COMETS, we augment the capabilities of the platform to include:**

* Single signal-reaction relationships, where the presence of a metabolite changes the bounds of a metabolic reaction
* Multiple signal-single reaction relationships, where the presence of multiple metabolites additively change the bounds of a metabolic reaction
