# Getting started with COMETS

## Installation of COMETS prerequisites

There are two software prerequisites for COMETS. 

First, users must download and install [Java](https://www.java.com/en/download/manual.jsp). The minimum version of required 64-bit Java is 1.8.

Second, users should download and install [Gurobi](https://www.gurobi.com/downloads/gurobi-software/). The installation of Gurobi requires a license, which is [free for academics](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). **Follow the Gurobi instructions to activate your license *before* attempting to use COMETS.** 

If COMETS is struggling to connect to Gurobi, users can begin the troubleshooting process by making sure the `GUROBI_HOME` environment variable is set. How exactly this is done will depend on the operating system of the user's computer. More information can be found on the [COMETS website](https://www.runcomets.org/home).

## Installation of the COMETS software

COMETS can be installed in two ways: using the COMETS installer or by unpacking the .tar.gz file. The installer is available [here](https://www.runcomets.org/installation). Installation for Windows, MacOS, and Linux systems are all available. The installer will guide users through a standard GUI installation procedure.

Alternatively, users can use the `comets_x.x.x.tar.gz` file for custom installation. More information on this process is available [here](https://cometspy.readthedocs.io/en/latest/installation.html).

## Installation of the COMETS `Python` Toolbox

While both `MATLAB` and `Python` toolboxes exist for COMETS, this tutorial uses `Python` to implement secondary metabolites.

It is recommended that users install `Python` (version >=3.6) using the [Anaconda distribution](https://www.anaconda.com/download) or other method as preferred. Anaconda and related methods have the benefit of installing Jupyter Notebooks at the same time. 

The COMETS `Python` toolbox (`cometspy`) is available from the package manager `PyPI` using the `pip` command.

```
pip install cometspy
```

## Version issues

There are two important version issues to note when attempting to get started with COMETS and the COMETS Python Toolbox. 

* Gurobi versions >10 are no longer compatible with Java
* `pandas` versions >2 and <1.5 are not compatible with the COMETS `Python` toolbox

## Dependencies

All simulations in this tutorial will rely on 4 `Python` packages:

* `cobra`
* `random`
* `pandas`
* `cometspy`
