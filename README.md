# collaboration-bestcsp-experimental-data

Sharing of experimental data for validating computational methods as part of the BEST-CSP project, see <https://www.cost.eu/actions/CA22107> and <https://best-csp.eu/>

This repository contains individual folders for each of the systems considered in this Action, and an example "System_template" to illustrate the format

Within each folder are .csv files for each polymorph, for standardised recording of data. Keeping the experimental data
to a standard format will aid model building/ statistical assessment of results, and csv files are machine readable

There is also a python script for deriving consensus means from data generated in different labs, in the top directory, called stats.py.
It requires the modules statsmodules and matplotlib; it may be simpler to make a conda env and install those modules, but I have included
environment.yaml as a snapshot for the pictures used in the Paris meeting (20/06/2024)
