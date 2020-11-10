[![PyPI version](https://badge.fury.io/py/autopacmen-Paulocracy.svg)](https://badge.fury.io/py/autopacmen-Paulocracy)

# AutoPACMEN (Automatic integration of Protein Allocation Constraints for stoichiometric MEtabolic Networks)


## General Description

AutoPACMEN allows one to apply the sMOMENT method of automatically expanding a stoichiometric metabolic
model with protein allocation constraints (as described in further detail in [Bekiaris & Klamt, 2020](#autopacmens-publication)).

This repository of AutoPACMEN consists of three source code parts:

1) The AutoPACMEN Model Generator, implemented as the Python module "autopacmen-Paulocracy" module, which provides a mostly automated way to generate enzyme-constraint-enhanced stoichiometric metabolic models, including an automatic retrieval of k<sub>cat</sub> values. It is primarily dependent on [cobrapy](https://github.com/opencobra/cobrapy).
<br>→This module can be found in the "autopacmen" subfolder.

2) An optional mixed Python 3/MATLAB AutoPACMEN Model Calibrator whose MATLAB parts primarily use [CellNetAnalyzer](https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html) and which allows one to optimize given protein allocation constraints in order to get a better fit with <i>in vivo</i> data.
<br>→The Python parts can be found in the "autopacmen" subfolder as described in AutoPACMEN's manual,the MATLAB parts can be found in the subfolder "AutoPACMEN_Model_Calibrator_MATLAB".

3) The exemplary usage of autopacmen-Paulocracy and the Model Calibrator resulting in the enzyme-constraint-enhanced model iJO1366*, as described in Supplementary File 1 of AutoPACMEN's publication. The final iJO1366* model is stored in a ready-to-use SBML form as "./iJO1366star/ec_model_2019_06_25_output_optimization/iJO1366star.xml".
<br>→The relevant scripts and data can be found in the "iJO1366" subfolder. <b>Note:</b> As it has a huge file size, the obligatory downloaded complete BRENDA text file brenda_downloads.txt (as described in AutoPACMEN's manual) is not included here. Instead, it can be downloaded from BRENDA's web site.

## Documentation

The combined manual for autopacmen-Paulocracy and the Model Calibrator can be found as "manual.pdf" or LibreOffice-compatible "manual.odt" in the "docs" subfolder folder.
It explains the manual installation process (without pip, see next chapter for the installation with pip) and usage of AutoPACMEN in detail.

An additional HTML documentation of the source code of AutoPACMEN's Python modules can be found under the "autopacmen" folder's subfolder "html".
This HTML documentation was automatically generated using [pdoc3](https://pdoc3.github.io/pdoc/) (link accessed on Oct 30, 2019).
The HTML documentation's starting point is "index.html" in the "./autopacmen/html/autopacmen" subfolder.

If you are particularily interested in the generation of k<sub>cat</sub> databases from BRENDA and SABIO-RK, look up
the scripts "data_parse_brenda_textfile.py", "data_parse_sabio_rk_for_model.py" as well as the combining
script "data_create_combined_kcat_database.py". These scripts create JSON files with the k<sub>cat</sub> data from these
databases with EC number, organism and substrate information. An exemplary created database with BRENDA is
"kcat_database_brenda" in the main folder's  subfolder "./iJO1366star/ec_model_2019_06_25_output", a database with SABIO-RK is
"kcat_database_sabio_rk" in the same subfolder.


## Installation of autopacmen-Paulocracy using pip

You can install autopacmen-Paulocracy [from PyPI](https://pypi.org/project/autopacmen-Paulocracy/) using pip as follows:
<pre>
pip install autopacmen-Paulocracy
</pre>

autopacmen-Paulocracy requires Python >=3.7 for its Python parts, and MATLAB >=2017a for its optional Model Calibrator
MATLAB scripts.


## Structure of AutoPACMEN's source code

All relevant scripts of autopacmen-Paulocracy are in the "autopacmen" main folder.

In this main folder, the scripts which start with "analysis_", "data_" and "modeling_" are command-line interfaces (CLI)
for AutoPACMEN's Python modules. These Python modules can be found the main folder's "submodules" subfolder.

The subfolder "AutoPACMEN_Model_Calibrator_MATLAB" contains the Model Calibrator's MATLAB parts.

All scripts and folders within "iJO1366star" are part of the generation and analysis of iJO1366* (see AutoPACMEN's publication for more about it). From these scripts, these ones
starting with "./iJO1366star/ec_model_2019_06_25_figure" create either a full figure or data for a figure used in AutoPACMEN's publication.

The main script for the generation of the uncalibrated iJO1366* model is "./iJO1366star/ec_model_2019_06_25_sMOMENT_iJO_CREATION.py" in the "iJO1366" subfolder. This
main script uses AutoPACMEN's functionalities as Python modules. The commented steps in this script correspond to the steps described in supplementary
File 1 of (Bekiaris & Klamt, in submission).

The "iJO1366star" folder's subfolder "iJOstar_MCS_analysis_scripts" contains the scripts used for the computation and the analysis of the published Minimal Cut Set enumeration
with iJO1366 and iJO1366*.

AutoPACMEN creates a cache of SABIO-RK, NCBI TAXONOMY and UniProt data in the "_cache" main folder's subfolder.  In the current state, the cache is the one
from AutoPACMEN's run for iJO1366* around the 25th of June 2019.


## AutoPACMEN's Publication

* Bekiaris, P.S., Klamt, S. Automatic construction of metabolic models with enzyme constraints. <i>BMC Bioinformatics</i> <b>21</b>, 19 (2020). [https://doi.org/10.1186/s12859-019-3329-9](https://doi.org/10.1186/s12859-019-3329-9)

## License

This project is free and open-source, using the Apache License Version 2.0.


## External sources

External sources which are included in this package are given in the respective SOURCES.txt files.
