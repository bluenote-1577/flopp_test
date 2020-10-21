
# Testing pipeline for flopp paper

The following are required to reproduce the simulation pipeline described in our paper. 

1. [flopp](https://github.com/bluenote-1577/flopp) software built.
2. [Haplogenerator](https://github.com/EhsanMotazedi/Haplosim) script must be available. 
3. [NanoSim](https://github.com/bcgsc/NanoSim) must be useable and downloaded.

We include the following in this folder:

4. [H-PoPG](https://github.com/MinzhuXie/H-PoPG) jar file; we include a jar file in this folder. 
5. [PaSS](http://cgm.sjtu.edu.cn/PaSS/) PacBio simulator; we include the binary in this folder. 

We also require that ``python2`` and ``python3`` be in PATH as Haplogenerator requires python2 but NanoSim requires python3. 

### To simulate genomes and test flopp/H-PoP against simulated genomes:
1. In ``get_reads_bam_pipeline.py``, modify the strings at the top of the script indicating the locations of the various required scripts and binaries. 
2. ``python get_reads_bam_pipeline.py``
3. The results for for each ploidy are in the folder ``(ploidy)/*/results_auto``, where ploidy = 3,4,5,6 and each of the folders in the ploidy folder encapsulates one iteration of a run.

### To plot the results of the previous step
1. To plot switch errors and hamming errors, do ``python plotting_scripts/plot_swerham 3/ 4/ 5/ 6/``. 
2. To plot the q-block error rates for one ploidy, do ``python plotting_scripts/plot_results (ploidy)/*`` NOTE: there has to be a wildcard at the end here.
