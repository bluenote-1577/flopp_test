import matplotlib.pyplot as plt
from copy import deepcopy
from glob import glob
import numpy as np
import pickle
import sys
plt.style.use(['science','no-latex'])

covs = ['10x / hap','15x / hap','20x / hap']
all_directories = sys.argv[1:]

fig, axs = plt.subplots(2,len(all_directories),figsize=(3, 9))
axs[0,0].set_yscale('log')

for (n,high_dir) in enumerate(all_directories):
    results_directories = glob(high_dir+"/*")
    #results_directories = sys.argv[2:]
    all_results = []
    for results_directory in results_directories:
        results = pickle.load(open(results_directory+'/results_auto/results.pickle','rb'))
        all_results.append(results)

    num_results = len(all_results)
    all_keys = all_results[0]['flopp'].keys()
    all_keys = sorted(all_keys)
    print(all_results)
    collapsed_ham_dict = dict()
    collapsed_swer_dict = dict()
    am_stdev_dicts = dict()
    swer_stdev_dicts = dict()

    for (i,result) in enumerate(all_results):
        ham_dicts = dict()
        swer_dicts = dict()
        for alg in result:
            if alg not in swer_dicts:
                swer_dicts[alg] = dict()
            if alg not in ham_dicts:
                ham_dicts[alg] = dict()
            for key in sorted(result[alg].keys()):
                type_read = key.split('-')[-1].split('.')[0]
                if 'swer' in key:
                    if type_read in swer_dicts[alg]:
                        swer_dicts[alg][type_read] = np.append(swer_dicts[alg][type_read],result[alg][key]/num_results)
                    else:
                        swer_dicts[alg][type_read] = np.array([result[alg][key]/ num_results])
                if 'ham' in key:
                    if type_read in ham_dicts[alg]:
                        ham_dicts[alg][type_read] = np.append(ham_dicts[alg][type_read], (result[alg][key] / num_results))
                    else:
                        ham_dicts[alg][type_read] = np.array([result[alg][key]/ num_results])


        print(ham_dicts)
        #print(swer_dicts)
        if i == 0:
            collapsed_ham_dict = ham_dicts
            collapsed_swer_dict = swer_dicts
            ham_stdev_dicts = deepcopy(ham_dicts)
            swer_stdev_dicts = deepcopy(swer_dicts)
        else:
            for alg in collapsed_ham_dict:
                for type_read in collapsed_ham_dict[alg]:
                    collapsed_swer_dict[alg][type_read] += swer_dicts[alg][type_read]
                    collapsed_ham_dict[alg][type_read] += ham_dicts[alg][type_read]
                    ham_stdev_dicts[alg][type_read] = np.append(ham_stdev_dicts[alg][type_read],ham_dicts[alg][type_read])
                    swer_stdev_dicts[alg][type_read] = np.append(swer_stdev_dicts[alg][type_read],swer_dicts[alg][type_read])


    for alg in collapsed_ham_dict:
        for type_read in collapsed_ham_dict[alg]:
            ham_stdev_dicts[alg][type_read] = ham_stdev_dicts[alg][type_read].reshape(len(all_results),len(covs))
            swer_stdev_dicts[alg][type_read] = swer_stdev_dicts[alg][type_read].reshape(len(all_results),len(covs))

    #        ham_stdev_dicts[alg][type_read] = num_results*np.std(ham_stdev_dicts[alg][type_read],axis=0) #originally scaled down above
    #        swer_stdev_dicts[alg][type_read] = num_results*np.std(swer_stdev_dicts[alg][type_read],axis=0)

            ham_stdev_dicts[alg][type_read] = [np.abs(num_results*np.amin(ham_stdev_dicts[alg][type_read],axis=0) - collapsed_ham_dict[alg][type_read])
                    ,num_results*np.amax(ham_stdev_dicts[alg][type_read],axis=0) - collapsed_ham_dict[alg][type_read]
    ] #originally scaled down above
            swer_stdev_dicts[alg][type_read] = [np.abs(num_results*np.amin(swer_stdev_dicts[alg][type_read],axis=0) - collapsed_swer_dict[alg][type_read])
                    ,num_results*np.amax(swer_stdev_dicts[alg][type_read],axis=0) - collapsed_swer_dict[alg][type_read]] #originally scaled down above


    print(ham_stdev_dicts)

    
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

    #print(times)
    #for (alg,l) in times.items():
    #    if len(l) > 0:
    #        colour = ""
    #        if 'whp' in alg:
    #            colour = whp_colour
    #        if 'hpop' in alg:
    #            colour = hpop_colour
    #        if 'flopp' in alg:
    #            colour = flopp_colour
    #
    #        x_vals = [x[0] for x in l]
    #        y_vals = [x[1] for x in l]
    #        axs[0].plot(x_vals,y_vals,'o-',label=alg,c=colour)
    #        axs[0].legend()
    #
    #axs[0].set_title("Time taken for run (note : flopp is run with 10x threads)")

    for j in range(0,2):
        d_to_use = dict()
        e_to_use = dict()
        axs[j,n].set_yscale('log')
        if j == 0:
            d_to_use = collapsed_ham_dict
            e_to_use = ham_stdev_dicts
            axs[j,n].set_ylabel("Hamming error rate")
            #axs[j,n].set_ylim(0.0005,0.15)
        else:
            d_to_use = collapsed_swer_dict
            e_to_use = swer_stdev_dicts
            axs[j,n].set_ylabel("Switch error rate")
            #axs[j,n].set_ylim(0.0020,0.05)

        x = np.arange(len(covs))  # the label locations
        width = 0.45  # the width of the bars
        
        rects1 = axs[j,n].bar(x - width/4, d_to_use['flopp']['pacbio'], width/4, label='flopp-pacbio',color=flopp_colour,edgecolor='black', yerr = e_to_use['flopp']['pacbio'])
        rects2 = axs[j,n].bar(x - width/2, d_to_use['flopp']['nanopore'], width/4, label='flopp-nanopore',color=flopp_colour,hatch='//',edgecolor='black',yerr = e_to_use['flopp']['nanopore'])
        rects4 = axs[j,n].bar(x + width/2, d_to_use['hpop']['pacbio'], width/4, label='hpop-pacbio',color=hpop_colour,edgecolor='black',yerr = e_to_use['hpop']['pacbio'])
        rects3 = axs[j,n].bar(x + width/4, d_to_use['hpop']['nanopore'], width/4, label='hpop-nanopore',color=hpop_colour,hatch='//',edgecolor='black',yerr = e_to_use['hpop']['nanopore'])

        if 'nanopore' and False in d_to_use['whp']:
            rects5 = axs[j,n].bar(x + width, d_to_use['whp']['nanopore'], width/4, label='whp-nanopore',color=whp_colour,hatch='//',edgecolor='black',yerr = e_to_use['whp']['nanopore'])
            rects6 = axs[j,n].bar(x + width*5/4, d_to_use['whp']['pacbio'], width/4, label='whp-pacbio',color=whp_colour,edgecolor='black',yerr = e_to_use['whp']['pacbio'])

        rects = [rects1,rects2,rects3,rects4]
        axs[j,n].set_xticks(x)
        axs[j,n].set_xticklabels(covs)
        axs[j,n].set_xlabel("Coverage")
        ploidy = n+3;
        if j == 1:
            title_str = "Switch error rate - "  + str(ploidy) + "x ploidy" 
            axs[j,n].set_title(title_str)
    #        axs[j].title(title_str)
        else:
            title_str = "Hamming error rate - "+str(ploidy)+"x ploidy"
            axs[j,n].set_title(title_str)
            if n == 0:
                axs[j,n].legend()
    #        axs[j].title(title_str)


for ax in axs.flat:
    ax.label_outer()

figure_name = 'long.svg' 
plt.savefig(figure_name, bbox_inches='tight')
plt.show()




##
#swer_rate-10-nanopore 10 hpop
#mean_ham_rate-10-nanopore 10 hpop



