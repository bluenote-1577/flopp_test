
# Testing pipeline for flopp paper

## Software dependencies
The following are required to reproduce the simulation pipeline described in our paper. 

1. [flopp](https://github.com/bluenote-1577/flopp) software built.
2. [NanoSim](https://github.com/bcgsc/NanoSim), the main script must be useable and the pre-trained models must be available as well. Note that you **must unzip the pre-trained models before using this pipeline**. The pre-trained model used is NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy.tar.gz.
3. [samtools](https://github.com/samtools/samtools) is expected to be in PATH.
4. [minimap2](https://github.com/lh3/minimap2) is expected to be in PATH.
5. [SciencePlots](https://github.com/garrettj403/SciencePlots) is used for plotting. 

## Python dependencies
We also require that ``python2`` and ``python3`` be in PATH as Haplogenerator requires python2 but NanoSim requires python3 (it seems that NanoSim only works with python3). The appropriate python packages must be installed for both versions. 

python3 packages for NanoSim:
* six  
* numpy (Tested with version 1.10.1 or above)
* HTSeq (Tested with version 0.9.1)  
* Pysam (Tested with version 0.13)  
* scipy (Tested with verson 1.0.0)  
* scikit-learn (Tested with version 0.20.0)

These must be available for the *python3* executable. 

python2 packages for Haplogenerator:
* biopython
* scipy

These must be available for the *python2* executable. 

## Included software
We include the following in this folder:

1. [H-PoPG](https://github.com/MinzhuXie/H-PoPG) jar file; we include a jar file in this folder. 
2. [PaSS](http://cgm.sjtu.edu.cn/PaSS/) PacBio simulator; we include the binary in this folder. 
3. [Haplogenerator](https://github.com/EhsanMotazedi/Haplosim) haplotype simulator; we include the script originally by Motazedi et al - Exploiting next-generation sequencing to solve the haplotyping puzzle in polyploids: a simulation study. We modified the haplogenerator.py script slightly so it would run faster on large genomes.


### To simulate genomes and test flopp/H-PoP against simulated genomes:
1. In ``get_reads_bam_pipeline.py``, modify the strings at the top of the script indicating the locations of the various required scripts and binaries. Note : the pre-trained models in NanoSim must be unzipped prior to usage.
2. ``python get_reads_bam_pipeline.py quick`` to run a small subset of the pipeline for testing purposes. Use ``python get_reads_bam_pipeline.py full`` to run the full pipeline (will take days to run). 
3. The results for for each ploidy are in the folder ``(ploidy)/*/results_auto``, where ploidy = 3,4,5,6 and each of the folders in the ploidy folder encapsulates one iteration of a run.

### To plot the results of the previous step
1. To plot switch errors and hamming errors, do ``python plotting_scripts/plot_swerham.py 3/ 4/`` after the output of quick pipeline or ``python plotting_scripts/plot_swerham.py 3/ 4/ 5/ 6/`` for the full pipeline. 
2. To plot the q-block error rates for one ploidy, do ``python plotting_scripts/plot_results (ploidy)/*`` NOTE: there has to be a wildcard at the end here.
