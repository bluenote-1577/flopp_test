import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from glob import glob
plt.style.use(['science','no-latex'])


covs = ['10x / hap','15x / hap','20x / hap']
higher_directories = sys.argv[1:]
for (i,high_dir) in enumerate(higher_directories):
    results_directories = glob(high_dir+"/*")
    print(high_dir)
    ploidy = sys.argv[1]
    all_results = []
    for results_directory in results_directories:
        results = pickle.load(open(results_directory+'/results_auto/results.pickle','rb'))
        all_results.append(results)

    num_results = len(all_results)

    plt.yscale('log')

    times = {'whp' : [], 'flopp' : [], 'hpop' : []}
    for alg in results:
        alg_results = results[alg]
        for params in alg_results:
            splitted = params.split('-')
            analysis_type = splitted[0]
            cov = splitted[1]
            #print(params,cov,alg)
            method = ""
            if 'whp' in alg:
                method = 'whp'
            if 'flopp' in alg:
                method = 'flopp'
            if 'hpop' in alg:
                method = 'hpop'

            if 'time' in params:
                times[method].append((int(cov),alg_results[params]))

    hpop_colour = '#e41a1c'
    whp_colour = '#377eb8'
    flopp_colour = '#4daf4a'

    print(times)
    for (alg,l) in times.items():
        if len(l) > 0:
            colour = ""
            if 'whp' in alg:
                colour = whp_colour
            if 'hpop' in alg:
                colour = hpop_colour
            if 'flopp' in alg:
                colour = flopp_colour

            x_vals = [x[0] for x in l]
            y_vals = [x[1] for x in l]
            if i == 0:
                if colour == flopp_colour:
                    plt.plot(x_vals,y_vals,'o--',label="flopp - 10 threads",c=colour)
                else:
                    plt.plot(x_vals,y_vals,'o-',label=alg,c=colour)
            else:
                if colour == flopp_colour:
                    plt.plot(x_vals,y_vals,'o--',c=colour)
                else:
                    plt.plot(x_vals,y_vals,'o-',c=colour)


    plt.title("Running time")

##Tested individually for the dataset specified in the paper. 

#flopp_1_3 = [[10,15,20],[9.4,16.7,24.3]]
#flopp_1_4 = [[10,15,20],[15.0,26.8,44.6]]
#flopp_1_5 = [[10,15,20],[22.25,41.2,67.5]]
#flopp_1_6 = [[10,15,20],[32.4,56.2,87.3]]
#flopp_1_res = [flopp_1_3,flopp_1_4,flopp_1_5,flopp_1_6]
#
#for (i,res) in enumerate(flopp_1_res):
#    if i == 0:
#        plt.plot(res[0],res[1],'o-',c=flopp_colour,label='flopp')
#    else:
#        plt.plot(res[0],res[1],'o-',c=flopp_colour)
#    plt.legend()
#
#plt.xlabel("Coverage per haplotype")
#plt.ylabel("Time (s)")
#plt.show()
