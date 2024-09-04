# collaboration-bestcsp-experimental-data

Sharing of experimental data for validating computational methods as part of the BEST-CSP project,
see <https://www.cost.eu/actions/CA22107> and <https://best-csp.eu/>

This repository contains individual folders for each of the compounds considered in this Action,
and an example "System_template" to illustrate the format.

Within each folder are .csv files for each polymorph, for standardised recording of data.
Keeping the experimental data to a standard format will aid in the statistical treatment of results,
and csv files are machine readable, facilitating plotting and tabulating the data.

There is also a Python script in the top directory for calculating consensus means from data generated in different labs called stats.py.
The script requires the modules statsmodels and matplotlib; it may be simpler to make a conda environment and install those modules,
but I have included environment.yaml as a snapshot for the pictures used in the Paris meeting (20/06/2024).

## Filling in CSV files

In writing the csv files, please consider the following guidance.
Each line can contain either a single data point, or the mean value and standard deviation of several data points from the same experiment/lab.
It is very important that the uncertainty or error bar written in the csv file is one sample standard deviation.
If you calculate it yourself, make sure to use the formula for sample standard deviation, with n-1 in the denominator.
In literature, other kinds of error bars are often used, and it is important that you try to figure out what kind of error bar it is.
If the paper doesn't specify, we assume that it is one std.
All data we have should be entered into the CSV. Outliers and questionable data should also be included for the sake of transparancy and completeness.
We then exclude such data by commenting out lines, and adding a statement about why that data point was rejected.
It is important that the entered data can be traced back to its origin. For literature, DOI or other reference that makes it possible
to find the study should be added.

What the Python script does:
If there are several data points from the same lab/name, the script will assume that they
are all from the same experiment and merge the data into a mean value and standard deviation.
If there are several literature studies that have only reported a single number without error bar,
the script will merge them and treat the data from all of them as one additonal laboratory.

DSC and other thermal data:
Melting temperature (temperature of fusion) should be the onset of the melting endotherm.
Values in J/g should be converted to kJ/mol. Take care in finding the correct molar mass.
We only use SI units. All temperatures should be in Kelvin.  

The first line in each CSV file contains a CSD reference specifying a crystal structure.
This is a prototypical structure that defines which compound and polymorph the data is for.
Make sure that the data you enter are for that polymorph.

Each line in the csv file then contaisn the following fields, separated by commas:

- Id: We number the data, so that it becomes easier to refer to a specific data point by specifying the filename and row number.
- Physical property: Some description of what was measured, the physical property.
- Value: This is the measured value for the physical property at that temperature, pressure. The units should be in the column heading
- Std: The sample standard deviation, if available.
- N: The number of data points, if a mean value was entered.
- Temperature: (Optional) if the Property is temperature dependent (i.e. Heat Capacity at 290K etc.,) then add a column here
- Name: The name of the scientist or the PI of the lab that did the measurements. The main importance is traceability, we need to keep track of where tha data came from.
For literature data, add the year it was published as Name_YY. For newly collected data, omit the year.
- Comment, Ref: Free text comments, DOI, experimental conditions or other info.
