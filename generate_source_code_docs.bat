REM This Windows script generates AutoPACMEN's source code documentation
REM In order to run this script, you have to install pdoc3 first
REM and you have to remove all active Python scripts (such as all
REM beggining with "ec_model_")
pdoc3 --force --output-dir "./autopacmen/html/" "autopacmen"
