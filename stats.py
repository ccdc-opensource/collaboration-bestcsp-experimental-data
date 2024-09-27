# imports
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
            newlines.append('                         mean       stdev      CI_low      CI_upp      weight')
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


files = glob.glob('%s\\*.csv' % (sys.argv[1]))
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
    ndict = dict()
    best_rep = data[0].split(',')[0]
    Units = data[1].split(',')[2].split()[-1].strip(')').strip('(')
    print(best_rep)
    for line in data[2:]:
        # print([s.strip() for s in line.split(',')])
        try:
            id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
        except Exception as e:
            print("Some issue with line: ", line, " Error: ", e)
        # Find out how many data points there are for each lab...
        if "#" in id:
            continue
        if author not in ndict:
            if n == '':
                ndict[author] = 1
            else:
                ndict[author] = int(n)
        else:
            if n == '':
                ndict[author] += 1
            else:
                ndict[author] += int(n)
    combined_names = []
    combined_values = []
    for name in ndict:
        values = []
        print(name, "has", ndict[name], "data points")
        if ndict[name] > 3:
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
                var_effect = np.append(var_effect, np.var(np.array(values)))
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
    if len(combined_values) > 2:
        #        idx.append(', '.join(set(combined_names)))
        idx.append("Literature data")
        print("Combining data for", set(combined_names))
        mean_effect = np.append(mean_effect, np.mean(combined_values))
        var_effect = np.append(var_effect, np.var(combined_values))

    print("mean_effect", mean_effect)
    print("var_effect", var_effect)
    print("Calculating consensus for", filename)
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
    Title = sys.argv[1].strip('\\').strip('/') + " form " + best_rep
    x_title = prop + " (" + Units + ")"
    ax.set_xlabel(x_title)
    ax.set_title(Title)
    fig1.tight_layout()
    fig1.savefig(filename.replace('.csv', '') + '_without_wls.png')

print("All done!")
