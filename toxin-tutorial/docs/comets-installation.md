# Getting started with COMETS

## Installation of COMETS prerequisites

There are two software prerequisites for the successful usage of COMETS. 

First, Java must be [downloaded and installed](https://www.java.com/en/download/manual.jsp). The minimum version of required 64-bit Java is 1.8.

Second, the Gurobi software should be [downloaded and installed](https://www.gurobi.com/downloads/gurobi-software/). The installation of Gurobi requires obtaining a license, which is [free for academics](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). **Follow the Gurobi instructions to activate your license *before* attempting to use COMETS.** 

If COMETS is struggling to connect to Gurobi, users should begin the troubleshooting process by making sure the `GUROBI_HOME` environment variable is set. In Windows, this can be done in the Control Panel. In Linux, the file `.bashrc` can be edited to contain the line `export GUROBI_HOME=/usr/gurobi/gurobi904/linux64` where `gurobi904` is replaced with the name and number of the appropriate installed version. For MacOS, the line `export GUROBI_HOME=/Library/gurobi904/mac64` should be added to the `.zshrc` file, where `gurobi904` is replaced with the identifying information for the installed version.

## Installation of the COMETS software

COMETS can be installed in two ways: using the COMETS installer or by unpacking the .tar.gz file. Use of the installer is recommended and is available [here](https://www.runcomets.org/installation). Users are asked to register, after which they can obtain the correct installer for their system. Installation for Windows, MacOS, and Linux systems are all available. The installer will guide users through a standard GUI installation procedure.

Alternatively, users can use the `comets_x.x.x.tar.gz` file for custome installation, best on a Linux system. The file should be unpacked in the directory where COMETS will be installed in order to create the COMETS installation directory:

```
$tar -xzvf comets_x.x.x.tar.gz   ./
```

If a user is completing custom installation on a Unix system (Linux or MacOS), they will also need to specify the `COMETS_HOME` environmental variable so that it points to the COMETS installation folder. This is done by adding the following line to the`.bashrc` file located in the home folder:

```
export COMETS_HOME = "/your/comets/installation/folder" 
```

## Installation of the COMETS `Python` Toolbox

While both `MATLAB` and `Python` toolboxes exist for COMETS, this tutorial uses `Python` to implement secondary metabolites.

It is recommended that users install `Python` (version >=3.6) using the [Anaconda distribution](https://www.anaconda.com/download) or other method as preferred. Anaconda and related methods have the benefit of installing Jupyter Notebooks at the same time. 

The COMETS `Python` toolbox (`cometspy`) is available from the package manager `PyPI` using the `pip` command. If `Python` has been installed through Anaconda, `pip` will have also been installed. Otherwise, in Linux systems, the `pip` command is installed through available repositories (e.g. `sudo apt-get install python3-pipin` Debian-based distributions). Once `pip` has been installed, `cometspy` can be downloaded by running the terminal command:

```
pip install cometspy
```

## Version issues

There are two important version issues to note when attempting to get started with COMETS and the COMETS Python Toolbox. 

* Gurobi versions >10 are no longer compatible with Java
* `pandas` versions >2 and <1.5 are not compatible with the COMETS `Python` toolbox

## Dependencies

All simulations in this tutorial will rely on 4 `Python` packages:

* `cobra`, especially `Metabolite`, `Reaction`, and `Model` 
* `random`
* `pandas`
* `cometspy`
