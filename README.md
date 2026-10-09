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
For repeated-measurement summaries, `Std` must be one sample standard deviation, not the standard error of the mean.
If you calculate it yourself, make sure to use the formula for sample standard deviation, with n-1 in the denominator.
In literature, other kinds of error bars are often used, and it is important that you try to figure out what kind of error bar it is.
Record the uncertainty type in the comments. For a single fitted or measured value, a reported one-sigma standard uncertainty
may also be stored in `Std` for weighting; this is not a sample SD estimated from one observation.
If a literature paper does not specify the uncertainty type, treating its error bar as one standard deviation is an explicit,
unverified assumption. Expanded uncertainties and confidence intervals require an appropriate conversion before entry.
All data we have should be entered into the CSV. Outliers and questionable data should also be included for the sake of transparancy and completeness.
We then exclude such data by commenting out lines, and adding a statement about why that data point was rejected.
It is important that the entered data can be traced back to its origin. For literature, DOI or other reference that makes it possible
to find the study should be added.

What the Python script does:
A single summarized row with `Std` enters as one laboratory/source effect.
For raw-only data from one name, at least two observations are merged into a mean and sample SD.
Single values without `Std` from at least two distinct names are pooled into one "Literature data" effect,
using the standard error of that independent pool's mean. An isolated raw value without a supplied uncertainty
cannot enter the metaanalysis. Multiple summarized rows, or a mixture of raw and summarized rows,
under the same name stop the analysis: their overlap and dependence must be resolved rather than
silently dropping observations or treating repeated summaries as independent effects.

### Statistical assumptions and uncertainty

The working model distinguishes dependence within a laboratory from
independence between laboratories:

- A supplied repeated-measurement `Std` gives the working variance `Std**2`,
  without division by `N`. Raw repeats from one name give
  `np.var(values, ddof=1)`, also without division by the group size.
  This deliberately uses repeatability as a laboratory-level uncertainty;
  repeated measurements do not automatically increase precision.
- A pool of `M` single values from distinct, independent laboratories gives
  `np.var(values, ddof=1) / M`, the estimated variance of the pool's mean.
  At least two independent values are needed. Repeated rows from one name
  are grouped as laboratory repeats first and never counted as independent
  contributions to this pool.

`ddof=1` specifies the conventional sample-variance definition with denominator
`N-1`; it does not correct for correlation between repeated measurements.

For independent, identically distributed within-lab measurements, the
sampling variance of a laboratory mean would normally be `Std**2 / N`.
That formula is deliberately not used for within-lab repeats here.
`statsmodels` normally expects the sampling variance of each effect; the
SD-based laboratory inputs are an intentional substitution of working
uncertainties. Avoiding the automatic reduction by `sqrt(N)` is intended to
be conservative about repeated measurements. Sample SD measures repeatability
and need not capture a shared calibration bias or other systematic errors,
so this choice does not guarantee coverage or an upper uncertainty bound.

`N` records a summary's measurement count where known; grouping counts
reported rows, not the sum of their `N` values. `N` is not a multiplier of
within-lab precision. A blank `N` does not assert that a summary had
one replicate. Rows with `Std` may therefore have unknown `N`, or `N=1` for
a supplied single-measurement uncertainty, without changing their weight.
A sample SD cannot itself be estimated from one observation. Provided `N`
must be a positive integer, and Name must identify the laboratory/source.
A single row without `Std` contributes one laboratory value to the
independent pool, including a reported mean with `N > 1`: its `N` does not
create additional independent contributions. Multiple rows without `Std`
from one name are treated as raw repeats only when their `N` is blank or `1`.
Multiple reported means without dispersions require their dependence and
dispersion to be resolved before combining. Zero, undefined, or non-finite working variances stop
the analysis rather than receiving infinite weight. These checks do not
modify any CSV.

### When a consensus is estimated

After exclusions, grouping and uncertainty checks, the script requires either:

- At least four distinct laboratory/source estimates with usable uncertainties;
  each must have a supplied positive `Std` or at least two raw repeats from
  which a positive sample SD can be calculated.
- At least three such laboratory estimates plus a pool of single values
  from at least two additional, distinct, independent laboratories/sources.

The pool counts as one effect in the metaanalysis, irrespective of its size.
Thus the second route uses at least five source laboratories but four fitted
effects. Two characterized laboratories plus a pool do not qualify, even if
many sources enter that pool. A dataset consisting only of uncharacterized
single values also does not qualify for this random-effects consensus.
Repeated measurements and different names for the same physical lab cannot
substitute for independent laboratory information.

If neither route is satisfied, the script prints a `Skipping` message with
the counts and continues to the next CSV without fitting or writing a new
consensus table or plot. Existing tables/plots are not overwritten or deleted
for skipped files; they are historical outputs, not new consensus estimates.
An isolated value without uncertainty is not silently assigned an error bar.

This is a pre-specified BEST-CSP eligibility policy, not a universal
statistical theorem or proof that four effects yield adequate precision.
Requiring at least four fitted effects supplies at least three residual
degrees of freedom for assessing excess between-effect scatter. With fewer
effects, a heterogeneity estimate and the resulting weights can be dominated
by one laboratory. Even four effects may give poorly determined heterogeneity;
the uncertainty qualifications below remain necessary.

### Scientific rationale

BEST-CSP's working premise is that a laboratory's few repeated measurements
often characterize repeatability better than total measurement error.
Repeats can share sample preparation, calibration, instrument and analysis
biases. A small within-lab SD therefore need not establish accuracy relative
to the true value. Using laboratory SD rather than SD/sqrt(N) deliberately
avoids making additional repeats alone justify increasingly strong weights.
This is an operational working-uncertainty choice, not a claim that sample SD
quantifies all systematic error.

Independent laboratories supply information about reproducibility that
within-laboratory repeats cannot provide. We model each laboratory/source
estimate as a common target value plus a laboratory-specific deviation and
measurement error. The random-effects variance represents excess scatter
beyond the supplied working variances. Paule–Mandel estimation and weights
`1 / (v_i + tau**2)` allow observed disagreement between laboratories to
moderate the influence of a laboratory reporting a very small uncertainty.
The original interlaboratory consensus framework of
[Paule and Mandel (1982)](https://doi.org/10.6028/jres.087.022) motivates
combining within- and between-group variability. It does not itself justify
BEST-CSP's substitution of sample SD squared for a mean's sampling variance.

If laboratory deviations are independent and centered on the true target,
averaging across laboratories can reduce their influence and improve
expected accuracy. Under equal variances, the variance contribution of
independent laboratory deviations to an equally weighted average decreases
as `1 / L`, where `L` is the number of independent laboratories. This is the
rationale for expecting consensus to be closer to the true value on average;
it does not establish that every consensus is more accurate than every lab.
A common calibration bias, shared methodological error, or unsuitable
experimental conditions can shift all results together and survives pooling.
Compatible properties, polymorphs and conditions, genuinely independent
sources, and external validation remain necessary for certified references.

With few effects, uncertainty in estimated heterogeneity is substantial.
[Röver, Knapp and Friede (2015)](https://doi.org/10.1186/s12874-015-0091-1)
show that small-study random-effects inference can require more cautious
interval procedures, including modified Knapp–Hartung inference. The
four-effect eligibility rule does not implement such a correction or validate
the coverage of the retained normal-based intervals.

### Estimator and reported uncertainty

The random-effects analysis retains the Paule–Mandel estimator (`method_re="pm"`).
For effect `i`, let `v_i` be the working variance described above.
Fixed-effect weights are proportional to `1 / v_i`; random-effects weights
are proportional to `1 / (v_i + tau**2)`. Here `tau**2` is additional
between-effect variance estimated conditional on the working variances;
it is not a separately measured laboratory bias. The output column
`model_sd` is `sqrt(v_i)` for an input effect and
`sqrt(1 / sum(1 / (v_i + tau**2)))` for the random-effects consensus.
For a literature input, `model_sd` is its pool's standard error; for a
repeated-measurement laboratory input, it is its supplied or estimated SD.

Distinct names are assumed to identify independent laboratories or sources.
Use the same name for dependent raw repeats; changing names or appending
years does not establish physical independence. `combine_effects` treats
effect rows as independent: multiple summaries from the same lab, overlapping
measurements, or shared calibrations can violate that assumption even with
SD-based working variances. Pool only compatible values for the same property,
polymorph, temperature, and conditions. The literature pool is a compressed
effect: its variance estimates uncertainty in an average of laboratory values,
but it does not retain individual study weights or allow separate heterogeneity
estimates for those studies in the subsequent metaanalysis.

A single effect cannot estimate between-effect heterogeneity. The retained
`statsmodels` defaults produce normal-based, nominal 95% model intervals;
they do not validate coverage for this working model or account for uncertainty
in estimated SDs. With few independent effects, assess small-sample behaviour
and systematic errors separately before certifying reference uncertainties.

Run the numerical regression tests with `python -m pytest -q tests` after
installing `pytest` in addition to the script's existing dependencies.

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
- Std: The sample SD for repeated measurements, or a documented one-sigma uncertainty for a single value, if available.
- N: The number of measurements contributing to a reported mean, if known; blank or 1 for a raw observation.
- Temperature: Put temperature in Property or comments (for example, Heat Capacity at 290 K). The current script expects Name as the sixth column; it does not support an inserted Temperature column or regression across temperatures.
- Name: The name of the scientist or the PI of the lab that did the measurements. The main importance is traceability, we need to keep track of where tha data came from.
For literature data, add the year it was published as Name_YY. For newly collected data, omit the year.
- Comment, Ref: Free text comments, DOI, experimental conditions or other info.

The dataset is provided under CC_BY_SA, please use for any purpose, but cite this repository or the relevant papers, and retain the same license if re-used.

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
