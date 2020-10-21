#REQUIRES PYTHON 2, ALSO REQUIRES 'python3' to be a version of python3 in your path for NanoSim usage. 

##Testing script :
##Generates haplotype files, reads, tests haplotype phasers, and outputs results. 

###Set these strings to be locations where the binary/folder/script is located. 
#Samtools is assumed to be in PATH
hpop_bin = "~/summer_project_phasing/H-PoPG/H-PoPGv0.2.0.jar"
whp_bin = "whatshap polyphase"
flopp_bin = "~/flopp/target/debug/flopp"
haplo_script = "~/software/Haplosim/haplogenerator.py"
nanosim_bin = "/home/jshaw/software/NanoSim/src/simulator.py"
nanosim_model ="/home/jshaw/software/NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training"

##This is for short-read simulation, not needed.
art_folder = "~/software/art/"

import subprocess
import time
import os
from datetime import datetime
from os import path

num_iterations = 1

for ploidy in range(3,5):
    for iternum in range(0,num_iterations):
        
        out_name = 'pds'
        
        ##SIMULATION PARAMETERS
        # p_mut = 0.02
        # 0.02 corresponds to 67.2 avg between bases
        # p_mut = 0.03
        # 0.03 corresponds to 45 avg between bases
        p_mut = 0.03

        ##Allele dosage for potatoes. Change depending on what your ploidy is.

        #3x dosage for potatoes
        #dosage = "[0.69,0.31,0]"
        #4x dosage for potatoes
        dosage = "[0.57,0.26,0.17,0]"
        #5x dosage for potatoes
        #dosage = "[0.50,0.23,0.14,0.13,0]"
        #6x dosage for potatoes
        #dosage = "[0.50,0.23,0.14,0.13,0,0]"

        
        #Coverages of reads per haplotype to be tested.
        covs = [10,15,20]

        #Genome length in bp of the reference.
        #gn_length = 3.02*10**6
        gn_length = 100000

        #Mean length of the generated nanopore and pacbio reads. Empirically tested; DON"T CHANGE
        mean_len = 8400

        ##Ignore
        recomb_freq = 0.001

        execute = True

        ##Simulate haplotypes based on reference. This step is needed for downstream analysis.
        get_hap_files = True

        ##Generate pacbio reads from simulated haplotypes.
        gen_pacbio_reads = True

        ##Generate nanopore reads from simulated haplotypes.
        gen_nanopore_reads = True

        ##Test flopp and H-PoPG, test WHP is run_whp is true.
        get_results = True
        run_whp = False

        ##Get SWER/HER/K-block error rates 
        kblock_results = True

        ##Remove intermediate files.
        clean=True


        ##Ignore
        gen_corr_reads = False
        recomb_genome = False
        gen_short_reads = False
        if recomb_genome and ploidy != 4:
            print("RECOMB GENOME REQUIRES PLOIDY 4")
            exit()

        folder_name_str = "PDS"+"_"+str(p_mut) + "_" + str(ploidy)+"_"+"iter"+str(iternum)

        if recomb_genome:
            folder_name_str = "RECOMB_" +str(recomb_freq)+ folder_name_str

        out_folder_name = './%s/%s/' %(str(ploidy),folder_name_str)
        ref_file = out_folder_name + 'ref/potato_chr1_50000_noN.fa'

        ## Random standard mutation dict
        mut_dict = {'A':'C','G':'A','C':'T','T':'G'}
        ## Polyalleic case
        mut_dict_poly = {'A':'C', 'A':'G' ,'G':'A','C':'T','T':'G'}

       
        hap_files_folder = '%s/hap_files' %(out_folder_name)
        vcf_outfile = '%s/%s.vcf' %(out_folder_name,out_name)

        def call(s,execute = execute, check_code = True):
            if execute:
                code = subprocess.call(s,shell=True)
                if code > 0 and check_code:
                    error_string = "Return code from command %s is not zero. Exiting!" %(s)
                    print(error_string)
                    exit()
            else:
                print(s)

        ##

        s = 'mkdir -p %s' %(out_folder_name)
        call(s,check_code=False)

        if get_hap_files:
            if not path.exists(hap_files_folder):
                s = 'cp -r ref/ %s/ref' %(out_folder_name)
                call(s)
                s = 'samtools faidx %s' %(ref_file)
                call(s)
                s = 'mkdir %s' %(hap_files_folder)
                call(s,check_code=False)
                s = "python %s -f %s -o %s/%s.fa --model poisson -s \"[%s,0,0]\" -p %s -v -m \"%s\" --dosage \"%s\"" %(haplo_script,ref_file,hap_files_folder,out_name,p_mut,ploidy,mut_dict,dosage)
                call(s)

                if recomb_genome:
                    s = 'python scripts/simulate_recomb.py %s %s' %(hap_files_folder,recomb_freq)
                    call(s)

                haplo_file = '%s/%s.fa_varianthaplos.txt' %(hap_files_folder,out_name)
                s = 'cut %s  -f-1,6- > %s/%s_haplo_truth.txt' %(haplo_file,hap_files_folder,out_name)
                call(s)

                s = 'python scripts/get_vcf_from_haplotypes.py %s %s %s/*.fa' % (ref_file, vcf_outfile, hap_files_folder)
                call(s)
            else:
                print("not creating hap files - hap files already exist")


        if gen_short_reads:
            for cov in covs:

                sim_short_folder = '%s/simulated_reads_illumina' %(out_folder_name)
                aln_short_folder = '%s/aln_files_illumina' %(out_folder_name)
                s = 'mkdir %s' %(sim_short_folder)
                call(s,check_code=False)
                s = 'mkdir %s' %(aln_short_folder)
                call(s,check_code=False)

                for i in range(1,ploidy+1):
                   ##Get paired reads for each hap
                    s = "%s/art_illumina --paired --in %s/%s.fa_hap%s.fa --out %s/%s_hap%s_paired_%sx_hs2500_350.fa --len 150 --fcov %s --mflen 650 --sdev 3.5 -1 ~/software/art/Illumina_profiles/HiSeq2500L150R1.txt -2 ~/software/art/Illumina_profiles/HiSeq2500L150R2.txt" %(art_folder,hap_files_folder,out_name,str(i),sim_short_folder,out_name,str(i),cov,cov)
                    print(s)
                    process = call(s)

                for i in range(1,ploidy+1):
                   ##Make reads have the same name
                    r_file2 = "%s/%s_hap%s_paired_%sx_hs2500_350.fa2.fq" %(sim_short_folder,out_name,i,cov)
                    r_file1 = "%s/%s_hap%s_paired_%sx_hs2500_350.fa1.fq" %(sim_short_folder,out_name,i,cov)
                    s = "sed -z 's/-2\\n/-1\\n/g' %s > %s_tmp" %(r_file2,r_file2)
                    call(s)
                    s = "mv %s_tmp %s" %(r_file2,r_file2)
                    call(s)

                    ##Map pairs
                    s = "bwa index %s" %(ref_file)
                    call(s)
                    s = "bwa mem -t 10 %s %s %s > %s/%s_hap%s_%sx_hs2500_350.sam" %(ref_file,r_file1,r_file2,aln_short_folder,out_name,i,cov)
                    process = call(s)

                merged_bam_name = "%s/%s_merged_%sx_hs2500_350.bam" %(aln_short_folder,out_name,cov)
                sorted_merged_bam_name = "%s/sorted_%s_merged_%sx_hs2500_350.bam" %(aln_short_folder,out_name,cov)
                s = "samtools merge -f %s ./%s/%s_hap*_%sx*" %(merged_bam_name,aln_short_folder,out_name,cov)
                process = call(s)
                s = 'samtools sort %s > %s' %(merged_bam_name,sorted_merged_bam_name)
                call(s)

                s = "java -jar %s -p %s -v %s -b %s -o %s/frag_files/illumina_%sx_frags.txt -b2f" %(hpop_bin,ploidy,vcf_outfile,sorted_merged_bam_name,out_folder_name,cov)
                process = call(s)

        for cov in covs:
            if gen_nanopore_reads:

                type_of_read = 'nanopore'
                aln_folder = '%s/aln_files_%s' %(out_folder_name,type_of_read)
                read_folder = '%s/simulated_reads_%s'%(out_folder_name,type_of_read)
                mean_length = 8400

                s = 'mkdir %s' %(aln_folder)
                call(s,check_code = False)
                s = 'mkdir %s' %(read_folder)
                call(s,check_code = False)

                #python src/simulator.py genome  -dna_type linear -rg ~/read_data/potato/potato_dosage_genome/ref/potato_chr1_50000_noN.fa -c pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training

                for i in range(1,ploidy+1):
                    hap_file = '%s/%s.fa_hap%s.fa'%(hap_files_folder,out_name,i)
                    num_reads = int(gn_length/mean_length * cov)
                    out_file = '%s/%sx_%s_sim_%s_hap%s' %(read_folder,cov,type_of_read,out_name,i)
                    s = 'python3 %s genome -dna_type linear -rg %s -c %s -n %s -o %s -t 15 --fastq -b guppy' %(nanosim_bin,hap_file,nanosim_model,num_reads,out_file)
                    process = call(s)

                for i in range(1,ploidy+1):
                    in_file = '%s/%sx_%s_sim_%s_hap%s_aligned_reads.fastq' %(read_folder,cov,type_of_read,out_name,i)
                    sam_out = '%s/%s_hap%s_%sx_%s.sam' %(aln_folder,out_name,i,cov,type_of_read)
                    #s = "minimap2 -Q -ax map-pb %s %s > %s" %(ref_file,in_file,sam_out)
                    s = "minimap2 -ax map-ont %s %s > %s" %(ref_file,in_file,sam_out)
                    #s = "graphmap align -r %s -d %s -o %s" %(ref_file,in_file,sam_out)
                    #s = "blasr movie.subreads.bam ecoli_K12.fasta --sam --out alignments.bam"


                    process = call(s)

                out_merged_sam = '%s/%s_merged_%sx_%s.bam' %(aln_folder,out_name,cov,type_of_read)
                s = "samtools merge -f %s/%s_merged_%sx_%s.bam %s/%s_hap*_%sx*" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov)
                process = call(s)

                s = "samtools sort %s/%s_merged_%sx_%s.bam > %s/sorted_%s_merged_%sx_%s.bam" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov,type_of_read)
                process = call(s)

                s = 'mkdir %s/frag_files' %(out_folder_name)
                call(s,check_code = False)

                s = "java -jar %s -p %s -v %s -b %s/sorted_%s_merged_%sx_%s.bam -o %s/frag_files/%s_%sx_frags.txt -b2f"%(hpop_bin,ploidy,vcf_outfile,aln_folder,out_name,cov,type_of_read,out_folder_name,type_of_read,cov)
                process = call(s)

                s = "rm %s/*.sam" %(aln_folder)
                call(s)

            if gen_pacbio_reads:

                type_of_read = 'pacbio_RS'
                aln_folder = '%s/aln_files_%s' %(out_folder_name,type_of_read)
                read_folder = '%s/simulated_reads_%s'%(out_folder_name,type_of_read)
                mean_length = 8400

                s = 'mkdir %s' %(aln_folder)
                call(s,check_code = False)
                s = 'mkdir %s' %(read_folder)
                call(s,check_code = False)
                s = 'mkdir %s/intermediate_files' %(out_folder_name)
                call(s,check_code = False)

                for i in range(1,ploidy+1):
                    hap_file = '%s/%s.fa_hap%s.fa'%(hap_files_folder,out_name,i)
                    s = 'perl pacbio_mkindex.pl %s %s/intermediate_files/' %(hap_file,out_folder_name)
                    process = call(s)

                    num_reads = int(gn_length/mean_length * cov)
                    out_file = '%s/%sx_%s_sim_%s_hap%s' %(read_folder,cov,type_of_read,out_name,i)
                    s = 'PaSS/PaSS -list %s/intermediate_files/percentage.txt -index %s/intermediate_files/index -m %s -c PaSS/sim.config -r %s -t 10 -o %s' %(out_folder_name,out_folder_name,type_of_read,num_reads,out_file)
                    process = call(s)
                    #time.sleep(1.2)

                for i in range(1,ploidy+1):
                    in_file = '%s/%sx_%s_sim_%s_hap%s.fq' %(read_folder,cov,type_of_read,out_name,i)
                    sam_out = '%s/%s_hap%s_%sx_%s.sam' %(aln_folder,out_name,i,cov,type_of_read)
                    #s = "minimap2 -Q -ax map-pb %s %s > %s" %(ref_file,in_file,sam_out)
                    s = "minimap2 -ax map-pb %s %s > %s" %(ref_file,in_file,sam_out)
                    #s = "graphmap align -r %s -d %s -o %s" %(ref_file,in_file,sam_out)
                    #s = "blasr movie.subreads.bam ecoli_K12.fasta --sam --out alignments.bam"


                    process = call(s)

                s = 'cat %s/%sx_* > %s/%sx_all_reads.fq' %(read_folder,cov,read_folder,cov)
                process = call(s)

                out_merged_sam = '%s/%s_merged_%sx_%s.bam' %(aln_folder,out_name,cov,type_of_read)
                s = "samtools merge -f %s/%s_merged_%sx_%s.bam %s/%s_hap*_%sx*" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov)
                process = call(s)

                s = "samtools sort %s/%s_merged_%sx_%s.bam > %s/sorted_%s_merged_%sx_%s.bam" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov,type_of_read)
                process = call(s)

                s = 'mkdir %s/frag_files' %(out_folder_name)
                call(s,check_code = False)
                s = "java -jar %s -p %s -v %s -b %s/sorted_%s_merged_%sx_%s.bam -o %s/frag_files/%s_%sx_frags.txt -b2f"%(hpop_bin,ploidy,vcf_outfile,aln_folder,out_name,cov,type_of_read,out_folder_name,type_of_read,cov)
                process = call(s)

                s = "rm %s/*.sam" %(aln_folder)
                call(s)

            if gen_corr_reads:
                type_of_read = 'pacbio_RS_corrected'
                aln_folder = 'aln_files_%s' %(type_of_read)
                read_folder = 'simulated_reads_%s'%(type_of_read)
                mean_length = 8400

                s = 'mkdir %s' %(aln_folder)
                call(s,check_code = False)
                s = 'mkdir %s' %(read_folder)
                call(s,check_code = False)
                s = 'mkdir intermediate_files'
                call(s,check_code = False)

                for i in range(1,ploidy+1):
                    in_file = '%s/rat_10x_pacbio_RS_sim_pot_dos_hap%s.fq.fasta' %(read_folder,i)
                    sam_out = '%s/%s_hap%s_%sx_%s.sam' %(aln_folder,out_name,i,cov,type_of_read)
                    ref_file = 'ref/potato_chr1_50000_noN.fa'
                    s = "minimap2 -Q -ax map-pb %s %s > %s" %(ref_file,in_file,sam_out)
                    #s = "minimap2 -ax map-pb %s %s > %s" %(ref_file,in_file,sam_out)
                    #s = "graphmap align -r %s -d %s -o %s" %(ref_file,in_file,sam_out)
                    #s = "blasr movie.subreads.bam ecoli_K12.fasta --sam --out alignments.bam"


                    process = call(s)

                out_merged_sam = '%s/%s_merged_%sx_%s.bam' %(aln_folder,out_name,cov,type_of_read)
                s = "samtools merge -f %s/%s_merged_%sx_%s.bam %s/%s_hap*_%sx*" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov)
                process = call(s)

                s = "samtools sort %s/%s_merged_%sx_%s.bam > %s/sorted_%s_merged_%sx_%s.bam" %(aln_folder,out_name,cov,type_of_read,aln_folder,out_name,cov,type_of_read)
                process = call(s)

                s = 'mkdir frag_files'
                call(s,check_code=False)

                s = "java -jar %s -p %s -v potato_50000chr1.vcf -b %s/sorted_%s_merged_%sx_%s.bam -o frag_files/%s_%sx_frags.txt -b2f"%(hpop_bin,ploidy,aln_folder,out_name,cov,type_of_read,type_of_read,cov)
                process = call(s)

            if get_results:
                read_types = ['nanopore','pacbio_RS']
                for type_of_read in read_types:
                    results_folder = '%s/results'%(out_folder_name)
                    results_folder_auto = '%s/results_auto'%(out_folder_name)

                    s = 'mkdir %s' %(results_folder)
                    call(s,check_code=False)

                    s = 'mkdir %s' %(results_folder_auto)
                    call(s,check_code=False)

                    aln_folder = '%s/aln_files_%s' %(out_folder_name,type_of_read)
                    read_folder = '%s/simulated_reads_%s'%(out_folder_name,type_of_read)

                    ###RUNNING HPOP
                    startTime = datetime.now()
                    s = "java -jar %s -p %s -v %s -b %s/sorted_%s_merged_%sx_%s.bam -o %s/hpop_%s_%s.txt"%(hpop_bin,ploidy,vcf_outfile,aln_folder,out_name,cov,type_of_read,results_folder_auto,cov,type_of_read)
                    process = call(s)
                    hpop_time = datetime.now() - startTime

                    ###RUNNING WHP
                    startTime = datetime.now()
                    if run_whp:
                        s = 'samtools index %s/sorted_%s_merged_%sx_%s.bam' % (aln_folder,out_name,cov,type_of_read)
                        process=call(s)
                        s = "%s -p %s %s %s/sorted_%s_merged_%sx_%s.bam -o %s/whp_%s_%s.vcf --ignore-read-groups"%(whp_bin,ploidy,vcf_outfile,aln_folder,out_name,cov,type_of_read,results_folder_auto,cov,type_of_read)
                        process = call(s)
                        s = "python scripts/phased_vcf2blocks.py %s/whp_%s_%s.vcf %s" %(results_folder_auto,cov,type_of_read, ploidy)
                        process = call(s)
                    whp_time = datetime.now() - startTime

                    ###RUNNING FLOPP
                    startTime = datetime.now()
                    s = "%s -p %s -v %s -b %s/sorted_%s_merged_%sx_%s.bam -o %s/flopp_%s_%s.txt"%(flopp_bin,ploidy,vcf_outfile,aln_folder,out_name,cov,type_of_read,results_folder_auto,cov,type_of_read)
                    process = call(s)
                    flopp_time = datetime.now() - startTime

                    s = "python scripts/haplo_truth_2_pickle.py %s/%s.fa_varianthaplos.txt %s/haplo.p %s" % (hap_files_folder,out_name,results_folder_auto,ploidy)
                    call(s)

                    time_file = "%s/times_%s.txt" %(results_folder_auto,cov)
                    if execute:
                        tfile = open(time_file,'w')
                    out_str = ""
                    if run_whp:
                        out_str = "hpop:%s,flopp:%s,whp:%s\n" % (str(hpop_time.total_seconds()),str(flopp_time.total_seconds()),str(whp_time.total_seconds()))
                    else:
                        out_str = "hpop:%s,flopp:%s" % (str(hpop_time.total_seconds()),str(flopp_time.total_seconds()))
                    if execute:
                        tfile.write(out_str)

        if kblock_results:
            s = 'python scripts/get_results.py %s %s' %(out_folder_name,ploidy)
            call(s)

        if clean:
            s = 'rm -r %s/simulated_reads*' %(out_folder_name)
            call(s)


