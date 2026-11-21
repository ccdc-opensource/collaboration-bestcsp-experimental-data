# imports
import os
import glob
import sys
import math
import numpy as np
from statsmodels.stats.meta_analysis import combine_effects
from statsmodels.graphics.dotplots import dot_plot
import matplotlib.pyplot as plt


def write_table_to_disk(table, name):
    lines = table.split('\n')
    newlines = []
    for line in lines:
        if line.split()[0] == 'eff':
            newlines.append('                         mean    model_sd      CI_low      CI_upp      weight')
        elif line.split()[0] == 'fixed':
            continue
        elif 'wls' in line[:19]:
            continue
        else:
            label = line[:19]
            ls = line[18:].split()
            newline = label + '{:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}'
            newline = newline.format(float(ls[0]), float(ls[1]), float(ls[2]), float(ls[3]), float(ls[5]))
            newlines.append(newline)
    f = open(filename.replace('.csv', '.dat'), 'w')
    f.write('\n'.join(newlines))
    f.close()

    return

folder = sys.argv[1]
files = glob.glob(os.path.join(folder, '*csv')) #'%s\\*.csv' % (sys.argv[1]))
print("Found", len(files), "csv files:")
for filename in files:
    print(filename)

for filename in files:
    mean_effect = np.array([])
    var_effect = np.array([])
    # nobs = np.array([]) # Number of observations for each lab.
    # And names/labels
    idx = []
    print("Reading", filename)
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    # Keep excluded rows out of every pass, not just the observation count.
    data = data[:2] + [line for line in data[2:]
                      if line.strip() and "#" not in line.split(',')[0]]
    ndict = dict()
    best_rep = data[0].split(',')[0]
    Units = data[1].split(',')[2].split()[-1].strip('()"')
    print(best_rep)
    for line in data[2:]:
        # print([s.strip() for s in line.split(',')])
        try:
            id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
        except ValueError as e:
            raise ValueError(f"{filename}: malformed row: {line.strip()}") from e
        # Find out how many data points there are for each lab...
        if "#" in id:
            continue
        try:
            if not author:
                raise ValueError("Name is required to identify dependent repeats and independent sources")
            if not np.isfinite(float(value)):
                raise ValueError("Value must be finite")
            if n and int(n) < 1:
                raise ValueError("N must be a positive integer when supplied")
            if s and (not np.isfinite(float(s)) or float(s) <= 0):
                raise ValueError("Std must be finite and positive")
        except ValueError as e:
            raise ValueError(f"{filename}: row {id}: {e}") from e
        if author not in ndict:
            ndict[author] = 1
        else:
            ndict[author] += 1
    combined_names = []
    combined_values = []
    for name in ndict:
        name_rows = [[s.strip() for s in line.split(',')][:6] for line in data[2:]
                     if line.split(',')[5].strip() == name]
        summaries = sum(bool(row[3]) for row in name_rows)
        if summaries > 1 or (summaries and len(name_rows) > 1):
            raise ValueError(f"{filename}: {name}: multiple summaries or mixed raw/summary rows; "
                             "resolve overlap and dependence before combining")
        if len(name_rows) > 1 and any(row[4] and int(row[4]) > 1 for row in name_rows):
            raise ValueError(f"{filename}: {name}: multiple reported means without dispersions; "
                             "resolve overlap and dependence before combining")
        values = []
        print(name, "has", ndict[name], "reported rows")
        if ndict[name] > 1:
            adding = True
            for line in data[2:]:
                id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
                if author == name:
                    if s.strip() != '':
                        print("has standard deviation")
                        mean_effect = np.append(mean_effect, np.array(float(value)))
                        var_effect = np.append(var_effect, np.array((float(s)**2)))
                        idx.append(name)
                        print(name, "OK, data added.")
                        adding = False
                    else:
                        values.append(float(value))
            if adding:
                mean_effect = np.append(mean_effect, np.mean(values))
                var_effect = np.append(var_effect, np.var(values, ddof=1))
                idx.append(name)
                print(name, "OK, data added")
        else:
            # Could be too few data or there is a standard deviation.
            for line in data[2:]:
                # print([s.strip() for s in line.split(',')][:6])
                id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
                if author == name:
                    if s.strip() != '':
                        print("has standard deviation")
                        mean_effect = np.append(mean_effect, np.array(float(value)))
                        var_effect = np.append(var_effect, np.array((float(s)**2)))
                        idx.append(name)
                        print(name, "OK, data added.")
                    else:
                        combined_names.append(name)
                        combined_values.append(float(value))
    n_lab_estimates = len(idx)
    n_pooled_labs = len(combined_values)
    if n_pooled_labs > 1:
        #        idx.append(', '.join(set(combined_names)))
        idx.append("Literature data")
        print("Combining data for", set(combined_names))
        mean_effect = np.append(mean_effect, np.mean(combined_values))
        # Each value here is from a distinct name/lab: variance of their mean.
        var_effect = np.append(var_effect, np.var(combined_values, ddof=1) / len(combined_values))

    if not np.all(np.isfinite(mean_effect)) or not np.all(np.isfinite(var_effect) & (var_effect > 0)):
        raise ValueError(f"{filename}: means must be finite and working variances finite and positive")

    # Eligibility policy: at least four effects with estimable uncertainties.
    # An independent literature pool counts as one effect, not M effects.
    consensus_eligible = (n_lab_estimates >= 4 or
                          (n_lab_estimates >= 3 and n_pooled_labs >= 2))
    skip_reason = ""
    if not consensus_eligible:
        skip_reason = (f"{n_lab_estimates} laboratory estimates and {n_pooled_labs} uncharacterized "
                       "independent laboratory values; require four laboratory estimates, or "
                       "three plus a pool of at least two additional independent laboratories")
        print(f"Skipping {filename}: {skip_reason}")
        continue

    print("mean_effect", mean_effect)
    print("var_effect", var_effect)
    print("Calculating consensus for", filename)
    # Working model: lab repeats use SD**2; independent literature means use s**2/M.
    results = combine_effects(mean_effect, var_effect, method_re="pm", row_names=idx)
    # results.conf_int_samples(nobs=nobs)
    print("Heterogeneity, tau-square: {:.3f}".format(results.tau2), "tau: {:.3f}".format(math.sqrt(results.tau2)))
    table = results.summary_frame()
    print(table)
    write_table_to_disk(table.to_string(), filename)

    #fig = results.plot_forest()
    # fig.tight_layout()
    #fig.savefig(filename.replace('.csv', '') + '.png')
    fig1, ax = plt.subplots()
    res_df = results.summary_frame()
    res_df.drop(res_df.tail(2).index, inplace=True)
    hw = np.abs(res_df[["ci_low", "ci_upp"]] - res_df[["eff"]].values)
    fig1 = dot_plot(points=res_df["eff"], intervals=hw,
                    lines=res_df.index, line_order=res_df.index, ax=ax)
    x_title = prop + " [" + Units + "]"
    ax.set_xlabel(x_title)
    fig1.tight_layout()
    fig1.savefig(filename.replace('.csv', '') + '_without_wls.png')

print("All done!")
