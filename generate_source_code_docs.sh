# This Unix script generates AutoPACMEN's source code documentation
# In order to run this script, you have to install pdoc3 first
# and you have to remove all active Python scripts (such as all
# beggining with "ec_model_")
pdoc3 --force --output-dir "./autopacmen/html/" "autopacmen"
