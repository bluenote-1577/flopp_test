# Appendix for flopp paper

The Appendix to the paper "Practical probabilistic and graphical formulations of long read haplotyping" by Jim Shaw and Yun William Yu is found in this repository and is titled appendix.pdf. The rest of the repository is devoted to recreating the tesing pipeline for the paper.  

# Testing pipeline for flopp paper

## Software dependencies
The following are required to reproduce the simulation pipeline described in our paper. 

1. [flopp](https://github.com/bluenote-1577/flopp) software built.
    1. On our Ubuntu 18.04.4 LTS test machine, the following instructions sufficed.
    1. Installing Rust:
        ```
        curl https://sh.rustup.rs -sSf | sh
        echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> $HOME/.profile
        source $HOME/.profile
        ```
    1. Installing Flopp
        ```
        mkdir -p $HOME/software
        cd $HOME/software
        git clone https://github.com/bluenote-1577/flopp
        cd flopp
        cargo build --release
        ```
2. [NanoSim](https://github.com/bcgsc/NanoSim), the main script must be useable and the pre-trained models must be available as well. Note that you **must unzip the pre-trained models before using this pipeline**. The pre-trained model used is NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy.tar.gz.
    1. Making sure Python3-pip is available
        ```
        sudo apt install python3-pip
        ```
    1. Installing NanoSim
        ```
        mkdir -p $HOME/software
        cd $HOME/software
        git clone https://github.com/bcgsc/NanoSim.git
        cd $HOME/software/NanoSim
        pip3 install -r requirements.txt
        ```
    1. Unzipping pre-trained models
        ```
        cd $HOME/software/NanoSim/pre-trained_models
        for file in *.tar.gz; do tar -xzf $file; done
        ```
3. [samtools](https://github.com/samtools/samtools) is expected to be in PATH.
    1. Installing samtools from repo
        ```
        sudo apt install samtools
        ```
4. [minimap2](https://github.com/lh3/minimap2) is expected to be in PATH.
    1. Installing minimap2 via their README
        ```
        mkdir -p $HOME/software
        cd $HOME/software
        git clone https://github.com/lh3/minimap2
        cd minimap2 && make
        echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> $HOME/.profile
        source $HOME/.profile
        ```
5. [SciencePlots](https://github.com/garrettj403/SciencePlots) is used for plotting. 
    1. Installing SciencePlots via Pip3
        ```
        pip3 install matplotlib
        pip3 install SciencePlots
        ```

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
* biopython<=1.76
* scipy
* matplotlib

These must be available for the *python2* executable. 

On our Ubuntu 18.04.4 LTS machine, the following lines achieved this:
```
sudo apt install python python3
pip2 install --upgrade pip
pip2 install biopython==1.76 scipy --upgrade

pip3 install --upgrade pip
pip3 install six numpy HTSeq Pysam scipy scikit-learn --upgrade
```
Note that unlike the previous `pip3 install -r requirements.txt` line when installing NanoSim, we try not to force particular versions. Make sure that the versions of scipy and numpy line up, however, or you may get odd `ModuleNotFoundError`s due to version number mismatch. Aside, biopython 1.76 is the last version that supports pip2, so we do need to force that version.

## Included software
We include the following in this folder:

1. [H-PoPG](https://github.com/MinzhuXie/H-PoPG) jar file; we include a jar file in this folder. 
2. [PaSS](http://cgm.sjtu.edu.cn/PaSS/) PacBio simulator; we include the binary in this folder. 
3. [Haplogenerator](https://github.com/EhsanMotazedi/Haplosim) haplotype simulator; we include the script originally by Motazedi et al - Exploiting next-generation sequencing to solve the haplotyping puzzle in polyploids: a simulation study. We modified the haplogenerator.py script slightly so it would run faster on large genomes.

To get these on your machine, simply clone this repository:
```
mkdir -p $HOME/software
cd $HOME/software
git clone https://github.com/bluenote-1577/flopp_test.git
cd flopp_test
```

### To simulate genomes and test flopp/H-PoP against simulated genomes:
1. In ``get_reads_bam_pipeline.py``, modify the strings at the top of the script indicating the locations of the various required scripts and binaries. Note : the pre-trained models in NanoSim must be unzipped prior to usage.
2. ``python get_reads_bam_pipeline.py quick`` to run a small subset of the pipeline for testing purposes. Use ``python get_reads_bam_pipeline.py full`` to run the full pipeline (will take days to run). 
3. The results for for each ploidy are in the folder ``(ploidy)/*/results_auto``, where ploidy = 3,4,5,6 and each of the folders in the ploidy folder encapsulates one iteration of a run.

### To plot the results of the previous step
1. To plot switch errors and hamming errors, do ``python plotting_scripts/plot_swerham.py 3/ 4/`` after the output of quick pipeline or ``python plotting_scripts/plot_swerham.py 3/ 4/ 5/ 6/`` for the full pipeline. 
2. To plot the q-block error rates for one ploidy, do ``python plotting_scripts/plot_results (ploidy)/*`` NOTE: there has to be a wildcard at the end here.
