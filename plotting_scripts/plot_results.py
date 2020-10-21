import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
plt.style.use(['science','no-latex'])

results_directories = sys.argv[1:]
all_results = []
for results_directory in results_directories:
    results = pickle.load(open(results_directory+'/results_auto/results.pickle','rb'))
    all_results.append(results)

num_results = len(all_results)
collapsed_results = {'whp':dict(),'flopp':dict(),'hpop':dict()}
kblock_keys = all_results[0]['flopp'].keys()

#print(kblock_keys)

for key in kblock_keys:
    for alg in collapsed_results.keys():
        for result in all_results:
            #print(key)
            if alg == 'whp':
                continue
            if key in collapsed_results[alg]:
                collapsed_results[alg][key] += np.array(result[alg][key])/num_results
            else:
                collapsed_results[alg][key] = np.array(result[alg][key])/num_results

results = collapsed_results

x_coords = []
covs = [10,15,20]
num_kblocks =9 
for i in range(num_kblocks):
    x_coords.append(50 + i*(100))

fig, axs = plt.subplots(3,sharey=True)
for (i,plt_cov) in enumerate(covs):
    title = "q-block error rate coverage %s" %(str(plt_cov))
    axs[i].set_title("title")
    for alg in results:
        alg_results = results[alg]
        print(alg_results)
        for params in alg_results:
            if 'kblock' not in params:
              #  print(alg,params,results[alg][params])
                continue
            print(params)
            splitted = params.split('-')
            analysis_type = splitted[0]
            cov = splitted[1]
            read_type = splitted[2]
            #label = alg+'-'+cov+'x-'+read_type
            label = alg+'-'+read_type
            percent_err_per_hap = np.divide(alg_results[params],x_coords)
            
            percent_err_per_hap = percent_err_per_hap[percent_err_per_hap > 0]
            str_xcoord = [str(x) for x in x_coords]
    #        plt.yscale("log")
            marker = '--'
            colour = 'r'

            print(cov,read_type)

            if 'flopp' in alg:
                colour = 'g'
            elif 'hpop' in alg:
                colour = 'r'
            elif 'whp' in alg:
                colour= 'b'

            if float(cov) != plt_cov:
                continue

            if read_type == 'nanopore':
                marker += "^"
            else:
                marker += 'o'

            axs[i].plot(str_xcoord[:len(percent_err_per_hap)],percent_err_per_hap,marker,label=label,c=colour)

            axs[i].legend()

#for ax in axs.flat:
#    ax.label_outer()
plt.xlabel("q")
plt.ylabel("q-block error rate")
plt.show()

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

#print(times)
#for (alg,l) in times.items():
#    if len(l) > 0:
#        x_vals = [x[0] for x in l]
#        y_vals = [x[1] for x in l]
#        plt.plot(x_vals,y_vals,'o-',label=alg)
#
#plt.title("Time taken for run (note : flopp is run with 10x threads)")
#plt.legend()
#plt.yscale('log')
#plt.show()
                





