import analysis_utils
import pickle
import sys
import glob

res_dir = sys.argv[1]
res_dir_auto = res_dir+'/results_auto/'
res_files =  glob.glob(res_dir_auto+"*.txt")
hap_file = glob.glob(res_dir_auto + 'haplo.p')[0]
ploidy = int(sys.argv[2])

print(res_dir_auto,res_files,hap_file)

results = dict()
results['whp'] = dict()
results['hpop'] = dict()
results['flopp'] = dict()


block_sizes = []
num_blocks = 9 
for i in range(0,num_blocks):
    block_sizes.append(100*(i)+50)

times = []
for file in res_files:
    print(file)
    if "times" in file:
        cov_of_file = file.split('_')[-1][:-4]
        f = open(file,'r')
        times = next(f).split(',')
        for time in times:
            key_name= 'time-'+cov_of_file
            software = time.split(':')[0]
            seconds = time.split(':')[1]
            results[software][key_name] = float(seconds)

    elif "txt" in file:
        for k in block_sizes:
            (err_per_block,blocks) = analysis_utils.kblock(file,k,ploidy,hap_file)
            type_of_read = ""
            cov = ""
            splitted = file[:-4].split('/')[-1].split('_')
            type_of_read = splitted[2]
            cov = splitted[1]
            kblock_key = 'kblock'+'-'+cov+'-'+type_of_read
            if blocks == 0:
                err_per_block = -1
            if 'whp' in file:
                if kblock_key in results['whp']:
                    results['whp'][kblock_key].append(err_per_block)
                else:
                    results['whp'][kblock_key] = [err_per_block]
            elif 'hpop' in file:
                if kblock_key in results['hpop']:
                    results['hpop'][kblock_key].append(err_per_block)
                else:
                    results['hpop'][kblock_key] = [err_per_block]
            elif 'flopp' in file:
                if kblock_key in results['flopp']:
                    results['flopp'][kblock_key].append(err_per_block)
                else:
                    results['flopp'][kblock_key] = [err_per_block]
        (N50,total_swer,total_wrongeno,mean_ham_rate,swer_rate) = analysis_utils.results(file,ploidy,hap_file)
        method = ""
        if 'whp' in file:
            method = 'whp'
        elif 'hpop' in file:
            method = 'hpop'
        elif 'flopp' in file:
            method = 'flopp'

        n50_key = 'n50' + '-'+cov+'-'+type_of_read
        swer_rate_key = 'swer_rate' + '-'+cov+'-'+type_of_read
        mean_hr_key = 'mean_ham_rate' + '-'+cov+'-'+type_of_read
        results[method][n50_key] = N50
        results[method][swer_rate_key] = swer_rate
        results[method][mean_hr_key] = mean_ham_rate
#        (N50,total_swer,total_wrongeno,avg_ham_rate_blocks) = analysis_utils.results(file,k,hap_file)

print(times)
print(results)

dump_str = res_dir_auto +"results.pickle"
print(dump_str)
pickle.dump(results,open(dump_str,"wb"))
