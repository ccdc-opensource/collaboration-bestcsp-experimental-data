# Manual for filling in csv files

Alongside the Readme in this directory, this Manual is intended to show how the csv files should be filled in, in a standardised way. There are probably 3 cases for reporting a value:
- If data has been generated during the Cost action (i.e. unpublished), then its almost certainly easier to fill in each data point on each line, leaving standard deviation and Number columns blank.
The stats.py script will collect all data points labelled with the same "Name" in the name column, and automatically calculate the mean and standard deviation for that set of data.
- If its published data, there will likely be one value quoted, with standard deviation and the number of data points collected reported somewhere; these can be filled in on one line and the code will
collect the relevant data for that set in one.
- There may be a value reported without standard deviation or number of experiments reported. We don't want to lose this data, but we should treat them as single data points; therefore report them
with number and standard deviation left blank, and the code will collect all sets (labelled with different "names") with < 4 experiments conducted, into one set. This will be called "Literature", and
treated as one lab for deriving the consensus mean
