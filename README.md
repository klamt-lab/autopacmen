# AutoPACMEN (Automatic integration of Protein Allocation Constraints for stoichiometric MEtabolic Networks)


## General Description

AutoPACMEN allows one to apply the sMOMENT method of automatically expanding a stoichiometric metabolic
model with protein allocation constraints (Bekiaris & Klamt, in submission).
AutoPACMEN consists of 2 parts:
1) A Python 3 Model Generator which primarily uses [cobrapy](https://github.com/opencobra/cobrapy) and which applies the sMOMENT
   method on a stoichiometric metabolic model.
2) An optional mixed Python 3/MATLAB Model Calibrator whose MATLAB parts primarily use [CellNetAnalyzer](https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html)
   and which allows one to optimize given protein allocation constraints in order to get a better fit
   with <i>in vivo</i> data.


## Documentation

AutoPACMEN's manual can be found as "manual.pdf" or "manual.odt" in the "docs" subfolder of the "autopacmen" folder.
It explains the manual installation process (without pip, see next chapter for the installation with pip) and usage of AutoPACMEN in detail.

An HTML documentation of the source code of AutoPACMEN's Python modules can be found under the main folder's subfolder "html".
This HTML documentation was automatically generated using [pdoc3](https://pdoc3.github.io/pdoc/) (link accessed on Oct 30, 2019).
The HTML documentation's starting point is "index.html" in the "./html/autopacmen" subfolder.

In addition, the exemplary usage of AutoPACMEN with iJO1366* is explained in AutoPACMEN's publication Supplementary File 1.

<b>Note:</b> As it has a huge file size, the downloaded complete BRENDA text file brenda_downloads.txt (as described in
AutoPACMEN's manual) is not included here.

If you are particularily interested in the generation of k<sub>cat</sub> databases from BRENDA and SABIO-RK, look up
the scripts "data_parse_brenda_textfile.py", "data_parse_sabio_rk_for_model.py" as well as the combining
script "data_create_combined_kcat_database.py". These scripts create JSON files with the k<sub>cat</sub> data from these
databases with EC number, organism and substrate information. An exemplary created database with BRENDA is
"kcat_database_brenda" in the main folder's  subfolder "ec_model_2019_06_25_output", a database with SABIO-RK is
"kcat_database_sabio_rk" in the same subfolder.

AutoPACMEN requires Python >=3.7 for its Python parts, and MATLAB >=2017a for its optional Model Calibrator
MATLAB scripts.



## Installation of Python parts using pip

You can install the Python parts of AutoPACMEN's latest release [from PyPI]() using pip as follows:
<pre>
pip install autopacmen-Paulocracy
</pre>


## Structure of AutoPACMEN's source code

<b>All relevant scripts are in the "autopacmen" main folder.</b>

In this main folder, the scripts which start with "analysis_", "data_" and "modeling_" are command-line interfaces (CLI)
for AutoPACMEN's Python modules. These Python modules can be found the main folder's "submodules" subfolder.

The main folder's subfolder "AutoPACMEN_Model_Calibrator_MATLAB" contains the Model Calibrator's MATLAB parts.

All scripts and folders starting with "ec_model_2019_06_25" are part of the generation and analysis of iJO1366* (see AutoPACMEN's publication for more on it). From these scripts, these ones
starting with "ec_model_2019_06_25_figure" create either a full figure or data for a figure used in AutoPACMEN's publication.

The main script for the generation of the uncalibrated iJO1366* model is "ec_model_2019_06_25_sMOMENT_iJO_CREATION.py" in the main folder. This
main script uses AutoPACMEN's functionalities as Python modules. The commented steps in this script correspond to the steps described in supplementary
File 1 of (Bekiaris & Klamt, in submission).


The main folder's subfolder "iJOstar_MCS_analysis_scripts" contains the scripts used for the computation and the analysis of the published Minimal Cut Set enumeration
with iJO1366 and iJO1366*.

AutoPACMEN creates a cache of SABIO-RK, NCBI TAXONOMY and UniProt data in the "_cache" main folder's subfolder.  In the current state, the cache is the one
from AutoPACMEN's run for iJO1366* around the 25th of June 2019.


## AutoPACMEN's Publication

* Bekiaris PS & Klamt S; "Automatic Construction of Metabolic Models with Enzyme Constraints"; in submission

## License

This project is free and open-source, using the Apache License Version 2.0.


## External sources

External sources which are included in this package are given in the respective SOURCES.txt files.
