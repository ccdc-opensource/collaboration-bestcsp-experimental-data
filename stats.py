# imports
import glob
import sys
import numpy as np
from statsmodels.stats.meta_analysis import combine_effects


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
    print ("Reading", filename)
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    ndict = dict()
    for line in data[2:]:
        #print([s.strip() for s in line.split(',')])
        id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
        # Find out how many data points there are for each lab... 
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
            for line in data[2:]:
                id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
                if author == name:
                    values.append(float(value))
            mean_effect = np.append(mean_effect, np.mean(values))        
            var_effect = np.append(var_effect, np.var(np.array(values)))        
            idx.append(name)
            print(name, "OK, data added")    
        else:
            # Could be too few data or there is a standard deviation. 
            for line in data[2:]:
                #print([s.strip() for s in line.split(',')][:6])
                id, prop, value, s, n, author = [s.strip() for s in line.split(',')][:6]
                if author == name:
                    if s != '':
                        print ("has standard deviation")
                        mean_effect = np.append(mean_effect, np.array(float(value)))
                        var_effect = np.append(var_effect, np.array((float(s)**2)))
                        idx.append(name)
                        print(name, "OK, data added.")
                    else:
                        combined_names.append(name)
                        combined_values.append(float(value))
    if len(combined_values) > 2:
        idx.append(', '.join(set(combined_names)))
        print("Combining data for", set(combined_names))
        mean_effect = np.append(mean_effect, np.mean(combined_values))
        var_effect = np.append(var_effect, np.var(combined_values))
 
    print ("mean_effect", mean_effect)
    print ("var_effect", var_effect)
    print("Calculating consensus for", filename)
    results = combine_effects(mean_effect, var_effect, method_re="pm", row_names=idx)
    #results.conf_int_samples(nobs=nobs)
    print("Heterogeneity, tau-square:", results.tau2)
    print(results.summary_frame())  
    print("Testing homogeneity:")
    homo = results.test_homogeneity()
    if homo.pvalue < 0.05:
        print("Warning, significant inhomogeneity detected!")
    print (homo)    
    fig = results.plot_forest()
    fig.tight_layout()
    fig.savefig(filename.replace('.csv', '') + '.png')
    # Plot results
