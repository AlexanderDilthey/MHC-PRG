# MHC-PRG

##HLA type inference

If you are interested in HLA typing, please see the HLA*PRG readme file:

https://github.com/AlexanderDilthey/MHC-PRG/blob/master/HLA-PRG.md

## Overview

Welcome to MHC-PRG.

MHC-PRG is a program that utilizes a Population Reference Graph of the human MHC region to improve genome inference (http://www.nature.com/ng/journal/v47/n6/full/ng.3257.html or open-access preprint http://biorxiv.org/content/early/2014/07/08/006973).

MHC-PRG will not do HLA type inference. For HLA type inference, see HLA-PRG (https://github.com/AlexanderDilthey/MHC-PRG/blob/master/HLA-PRG.md).

We have tested MHC-PRG on samples with relatively high coverage (>30x). It may or may not work well on sample data with less coverage (we would certainly appreciate any feedback).

Also, although the program has been executed many times on our servers, many bugs, in particular on non-Oxford servers, surely remain. When submitting bug reports, please provide as many details as possible.

## Paper-related data

See https://github.com/AlexanderDilthey/MHC-PRG/blob/master/PaperData.md for paper-related data and access.

## Wiggle files for UCSC

We provide the results from our genome-wide VCF evaluation in bigWig format for use in the UCSC genome browser (computed in 200bp windows). Include the following URL as 'custom track':

http://oxfordhla.well.ox.ac.uk/VCF.bw

## Installation

To install the package, please follow this procedure:

1. If you want to reproduce results from our paper, download the MHC-PRG data package (http://www.well.ox.ac.uk/MHC-PRG.tar.gz - c. 213GB), which will also set up the right file structure. We assume that you have uncompressed the data into ~/MHC-PRG. 

   If you are only interested in HLA typing, ignore the MHC-PRG data package (there is a separate HLA*PRG data package to be downloaded later, see https://github.com/AlexanderDilthey/MHC-PRG/blob/master/HLA-PRG.md).
   
   Instead, manually create the directory structure for the build (you can of course install into any other directory than ~/MHC-PRG):

   mkdir ~/MHC-PRG  
   mkdir ~/MHC-PRG/src  
   mkdir ~/MHC-PRG/obj  
   mkdir ~/MHC-PRG/bin  
   mkdir ~/MHC-PRG/tmp  
   mkdir ~/MHC-PRG/tmp2  
   mkdir ~/MHC-PRG/tmp2/GS_nextGen  
   mkdir ~/MHC-PRG/tmp2/GS_nextGen/hla  

2. Now, cd into the directory ~/MHC-PRG/src. Use git to clone a copy of this repository into the src directory.

3. Install all other programs / libraries that MHC-PRG depends on (see list below). Note that this list is much reduced for HLA typing only.

4. Modify the makefile so that MHC-PRG compiles on your machine (see below).

5. Start using the software. Our data package comes with data for NA12878 in the right formats.

### Software dependencies

#### For HLA typing

If you are only installing MHC-PRG to run the HLA*PRG component:

Boost (1.52)

Bamtools (commit https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e)

#### For MHC-PRG

If you are installing MHC-PRG to reproduce the results from our Nature Genetics paper:

The pipeline requires Cortex (1.0.5.12), BWA (0.6.2-r126), Stampy (1.0.20), the PICARD_SAM2FASTQ from Picard (1.83), samtools (0.1.18), Platypus (0.2.0), Boost (1.52), Bamtools (commit https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e) and sometimes zlib. The numbers in brackets indicate the version numbers as currently used for testing.

The pipeline also requires access to the human genome reference, a Cortex graph of the human genome reference, the PGF haplotype, and a haplotype graph of: the eight xMHC reference haplotypes, the HLA alleles and known SNPs. These files are provided in our data package.

After installing the required programs, open nextGenInferenceVariGraph.pl and nextGenValidationVariGraph.pl and make sure that the absolute paths to the required programs match your local installation.

It is also a good idea to try to execute the required other programs in your user context, simply by calling the main program file (even if without valid input data). Any errors relating to shared libraries, permissions etc. should be fixed before proceeding.

### Modifying the makefile

For reasons I (Alexander) don't fully understand, different machines seem to require different makefile settings.

The bottom line is that some Boost and bamtools libraries need to be compiled into MHC-PRG. You might have to play around a bit with different options (provided in the makefile) until everything works out.

1. Open the makefile and modify the path to the Boost C++ library and to all the Boost library files.
2. Open the makefile and modify the path to Bamtools (incl. include, src and the Bamtools object files)
3. Sometimes required: insert a correct path to zlib.

A simple "make all" should then be sufficient to build the C++ components of the pipeline.

#### Compiler

If you run into problems, check that your compiler has support for openMP and C++11 (we use G++ 4.7.2).

We have received reports that there are compilation issues with G++ 5.1, but haven't been able to reproduce these yet - if you run into problems, please let us know!

## Running MHC-PRG

### Sample IDs and data files

Each sample to be processed needs a unique sample ID.

For each sample, we require a BAM file and a Cortex graph, using the same k (kMer length) as the xMHC graph (this is usually 31).

We use the BAM the re-map reads to the inferred xMHC "haplotypes", and we use the Cortex graph to extract kMer counts, which our haplotype graph model is based on.

For each sample, these files are expected in the directory ~/MHC-PRG/data/samples/CortexGraphs and ~/MHC-PRG/data/samples/BAMs. Filenames should be composed according to the rule SAMPLE_ID.bam and SAMPLE_ID_$k.ctx.

We provide sample input data for NA12878 under the sample identifier AA02O9Q_Z2 don't ask). The corresponding files are hence called AA02O9Q_Z2.bam and AA02O9Q_Z2_31.ctx.

If you want to change the standard lookup paths, you can manually edit nextGenInferenceVariGraph.pl and nextGenValidationVariGraph.pl.

### Running the pipeline

Our examples use the sample ID "AA02O9Q_Z2". Change as necessary.

1.

  Increase standard stack size (typically required): on Linux systems, the command 'ulimit -s 81920' should do.

2.

  Data analysis step I: Execute
  
  ./nextGenInferenceVariGraph.pl --graph ../tmp2/GS_nextGen/varigraph3 --sample AA02O9Q_Z2 --kmer 31
  
  to prepare data analysis.
  
  If data preparation was successful, the script will print another command you need to execute.
  
  This will require some RAM (maybe 50 - 100G) and make use of openMP multithreading.

3.

  To be explicit, if you haven't executed the command that the previous script printed, do so now.
  
  If you encounter an error, first check 'ulimit -a' and make sure that stack size was really set to 81920. If this was the case, send me an email.

4.
  Now execute three commands to generate data in different output formats:
  
  ./nextGenInferenceVariGraph.pl --graph ../tmp2/GS_nextGen/varigraph3 --sample AA02O9Q_Z2 --kmer 31 --collect 2 --vcfPos ../data/VCFsnpRefData.txt ;\
  ./nextGenInferenceVariGraph.pl --graph ../tmp2/GS_nextGen/varigraph3 --sample AA02O9Q_Z2 --kmer 31 --collect 2viterbi --vcfPos ../data/VCFsnpRefData.txt ;\
  ./nextGenInferenceVariGraph.pl --graph ../tmp2/GS_nextGen/varigraph3 --sample AA02O9Q_Z2 --kmer 31 --collect 3 --vcfPos ../data/VCFsnpRefData.txt --mapping_display_only_instructions 1
  
  The "collect" parameter decides what data exactly is collected and what output is produced.
  
  "--collect 2" collects data that was generated by sampling from the posterior distribution over haplotype paths through the haplotype graph. The output file ../tmp/kMerCount__GS_nextGen_varigraph3_AA02O9Q_Z2_31_required.binaryCount.haplotypes.VCF will have quality measures attached to alleles which are generated from taking the marginal allele probabilities at each position independently.
  
  "--collect 2viterbi" collects data that was generated by calculating the Viterbi path through the haplotype graph model. The name of the output file  is. ../tmp/kMerCount__GS_nextGen_varigraph3_AA02O9Q_Z2_31_required.binaryCount.viterbiHaplotypes.viterbiVCF.
  
  "--collect 3" (with --mapping_display_only_instructions 1) takes the Viterbi-generated haplotypes and carries out read re-mapping and de novo variant calling on top of the inferred haplotypes. Specifically, we create two "personalized reference genomes", by taking the original reference genome, removing the xMHC and inserting our two "haplotypes". We then (independently) re-map all sample reads to both "personalized reference genomes". We examine the "personalized xMHC" in both genomes. For each read in the union of reads that map to the "personalized xMHC" in either "personalized reference genome", we want to make a decision which of the genomes it "comes from". We use mapping quality to make this decision. When mapping quality is equal for both genomes, we make a random decision. We produce 2 SAM files, one for each "personalized xMHC", containing the reads we assigned to each "haplotype". We then carry out variant discovery (currently using Platypus), and modify the Viterbi haplotype VCF to take into account the variants discovered by Platypus. 
  
  Importantly, the read re-mapping step involved two rounds of whole-genome alignment. In many settings, it will be advantageous to use a local, cluster-based pipeline to achieve this. The switch --mapping_display_only_instructions 1 enables the printing of mapping instructions - it will tell you which read files need to be mapped against which references, and where the resulting BAMs need to be copied in order for nextGenInferenceVariGraph.pl to find them. We have tested this step with BWA, and we highly recommend that you use BWA as well. Mapping should be carried out in single-end mode, and the BAMs should be sorted and indexed (and the index files should be copied 

  If --mapping_display_only_instructions is not active, the program will try to control a local instance of BWA  - the default setting for the number of threads is 42 (you need to edit the Perl code if you want to modify this), and it might take a long time.

## Acknowledgements

This software uses code from Zamin Iqbal's and Mario Caccamo's Cortex (http://cortexassembler.sourceforge.net/index.html) and Erik K Garrison's intervaltree package (https://github.com/ekg/intervaltree).

