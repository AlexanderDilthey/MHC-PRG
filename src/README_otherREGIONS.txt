Application of the haplotype graph system is not limited to the xMHC. However, it has to be stressed that this initial version was tested only on the xMHC, and that the code has not been optimized to handle very large regions. The current implementation can probably deal with regions of up to maybe 10mb. The specifics depend on the computational resources you are willing to utilize.

Some parts of the pipeline can currently not be expected to work for regions outside the xMHC (namely nextGenInferenceVariGraph.pl). 

Note that it is a good idea to increase stack size before working with our system. 81920 works well for all xMHC-related applications, but as you increase the size of the region you are working on, you might want to experiment with higher values. The Linux command to set stack size to 81920 is 'ulimit -s 81920'.

All relative paths are relative to the src directory of PnPHaplograph2.

1. Overview of the pipeline

- Set up the required files (more on that in a bit) and store them in a separate directory. We usually use ../tmp2/GS_nextGen/GRAPH_NAME, and if you follow this convention, all helper scripts should automatically find the right files.

- Build a nucleotide graph for GRAPH_NAME:
  
  ../bin/PnPHaploGraph2 domode createConcatenatedVariationGraphs ../tmp2/GS_nextGen/GRAPH_NAME
  
- "kMerify" the nucleotide graph:

  ../bin/PnPHaploGraph2 domode kMerify ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt --kmer 31
  
  "kMerification" means that we transform the nucleotide graph, an object in which each edge is labeled with a single nucleotide, into an equivalent object in which each edge is labeled with a kMer. This is more or less a large forward search from each node in the nucleotide graph.
  
  This process becomes more computationally intensive as k becomes bigger (because the look-ahead space reachable from each node in the nucleotide graph becomes bigger). k = 31 has produces sensible results for the xMHC.
  
- Carry out some simulations:

  ../bin/PnPHaploGraph2 domode simulationSuite ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31 ../tmp/ SIMU_GRAPH --genotypingMode 8 > simulation_results.txt
  
  The number of simulation iterations is currently determined in code (NextGen/simulationSuite.cpp - sorry!). The output will contain some summary statistics as well as detailed results broken down by, for example, variant length. We do currently not have a detailed documentation of the output format.
  
- Apply the model to your sample data:

  See below. Unfortunately, we cannot apply nextGenInferenceVariGraph.pl.
  
2. Required files

The model obviously requires data that describe regional haplotypes and polymorphisms. This section contains information on what these data look like and how to create them. Unless otherwise specified, all files should be placed in the directory ../tmp2/GS_nextGen/GRAPH_NAME/.

The basic intuition behind the input files is as follows: you have simple text-based files, in which each line refers to a line in an alignment of the regional haplotypes. For example, you could produce an alignment of the eight xMHC reference haplotypes (from Horton et al. 2008). Each line in this alignment (including gaps!) would end up as one line in one of our input files.

Now we want to allow for SNPs, and it is entirely possible that there are SNPs which exist only relative to some underlying reference sequence, which is not necesarily present on all utilized reference haplotypes. Therefore, each line which contains a line from the haplotype alignment can be followed by a number of lines that specify SNPs. In these lines, you specify a place holder symbol if there is no polymorphism at this position, and the nucleotide if there is a SNP relative to the base haplotype sequence. Each line can only contain one allele per position, so that you will have two additional lines if you have a maximum of two alleles per SNP for a particular reference haplotype.

Some regions can be split into subregions with differing numbers of reference sequences. In the HLA, for example, we have eight general reference haplotypes, and in between there are the six classical HLA loci, for which we have many hundreds of "mini reference haplotypes" (i.e. the sequences of the known HLA alleles). To accommodate such situations, you simply create a new file whenever the number of references changes (listing of course only the relevant positions), and then specify a master file ("segments.txt") which lists these files in linear order.

This is not as complicated as it sounds. For our MHC reference graph, segments.txt looks like this:

segment_1.txt
segment_HLAA.txt
segment_2.txt
segment_HLAC.txt
segment_3.txt
segment_HLAB.txt
segment_4.txt
segment_HLADRB.txt
segment_5.txt
segment_HLADQA.txt
segment_6.txt
segment_HLADQB.txt
segment_7.txt

.. and this here is an excerpt from one of the individual files, segment_3.txt:

IndividualID S3_0_2537944 S3_1_2537945 S3_2_2537946 S3_3_2537947 S3_4_2537948 S3_5_2537949 S3_6_2537950 S3_7_2537951 S3_8_2537952 S3_9_2537953 S3_10_2537954 S3_11_2537955 S3_12_2537956 S3_13_2537957 S3_14_2537958 S3_15_2537959 S3_16_2537960 rs115573554_17_2537961 S3_18_2537962 S3_19_2537963 S3_20_2537964 S3_21_2537965 S3_22_2537966 S3_23_2537967 S3_24_2537968 S3_25_2537969 S3_26_2537970 S3_27_2537971 S3_28_2537972 S3_29_2537973 S3_30_2537974 S3_31_2537975 S3_32_2537976 S3_33_2537977 S3_34_2537978 S3_35_2537979 S3_36_2537980 S3_37_2537981 S3_38_2537982 S3_39_2537983 S3_40_2537984 S3_41_2537985 S3_42_2537986 S3_43_2537987 rs188509578_44_2537988 S3_45_2537989 S3_46_2537990 S3_47_2537991 S3_48_2537992 S3_49_2537993
pgf C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_pgf * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A * *   
mann C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_mann * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A *    
ssto C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_ssto * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A *    
apd * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
SNPs_1_apd * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
cox C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_cox * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A * *   
qbl C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_qbl * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A * *   
dbb C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_dbb * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A * *   
mcf C C T C T T C T C C T A C A C C A A G C A T C T T T G T C A C A C T G T G T G C C T G A G T C C T G
SNPs_1_mcf * * * * * * * * * * * * * * * * * T * * * * * * * * * * * * * * * * * * * * * * * * * * A * *   

In these files, gaps are represented by underscores (_), and unknown nucleotides (often 'N' in FASTA) are encoded as stars (*). Fields are separated by spaces.

The locus identifiers in the first line of each file ('IndividualID' always has to be there as first field, for legacy reasons) should be unique and contain three fields, separated by underscores: 1) An identifier, which can refer to segment number (e.g. "S1") or identify a genomic locus (e.g. "rs123"). 2) A position identifier which indexes positions in the current segment alignment (ie a column index into the alignment). 3) A position identifier which refers to that position's alignment position in the reference genome. This value only increases if the reference haplotype has a non-gap character at this column. Or, the other way around: it is possible for two (or more) subsequent positions to have an identical third number, if the reference haplotype has gaps at the correspondining positions in the alignment of all haplotypes.

The SNP-specifying lines are enumerated (i.e. the next SNP line for PGF would be named 'SNPs_2_pgf'). Non-polymorphic positions are filled by stars, and stars are also used if the reservoir of alleles at a particular position relative to a haplotype has already been exhausted by the previous SNP lines for one haplotype (this rule only applies within rows with the same main identifier. That is, each main haplotype, e.g. 'cox', has its own complete SNP definition).

Finally, we also need to create a file named 'positions.txt' (same folder as the other files). For all position identifers specified in all segment files, positions.txt needs to specify an absolute position. This can be relative to the beginning of the main alignment. Identifiers and positions are separated by spaces. The specified positions should be continuous. The first few lines of positions.txt for the xMHC graph look like that:

rs150673860_0_1000 1
S1_1_1001 2
S1_2_1002 3
S1_3_1003 4
S1_4_1004 5
rs189316702_5_1005 6
...

3. Building the graph

The steps required to create a kmerified haplotype graph are described in section 1 of this file.

4. Application of the model

The model operates on genome-wide kMer counts. Therefore, we need to find out
- which kMers are specified by the model
- which of these kMers appear in the rest of the reference genome (they will be ignored)
- how often each of the remaining kMers appears in the sample genome (from the Cortex file)

(Note that you need to replace GRAPH_NAME and kmer_XX with whatever graph name and kMer size you are using in the following examples).

Find out which kMers appear in the haplotype graph:

   ../bin/PnPHaploGraph2 domode determineRequiredKMers ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31  ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers

Find out which of these appear in the reference genome:
   
   ./readCortexCoverage.pl --cortex_bin PATH/TO/HUMANGENOME/GRAPH  --kMer_size 31 --interesting_kMers ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers
	
   (Not surprisingly, PATH/TO/HUMANGENOME/GRAPH should point to a Cortex graph of the whole human genome for the correct kMer size. You can download such a graph for k = 31 from our website - see README.txt).
   
   This produces a file ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers.binaryCount, which specifies how often each kMer appears in the human reference genome.

Find out how often each kMer appears in the sample genome:

	cp ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers ../tmp2/GS_nextGen/GRAPH_NAME/SAMPLENAME_graph.txt.kmers_31.requiredKMers
   ./readCortexCoverage.pl --cortex_bin PATH/TO/SAMPLENAME_CORTEX_GRAPH  --kMer_size 31 --interesting_kMers ../tmp2/GS_nextGen/GRAPH_NAME/SAMPLENAME_graph.txt.kmers_31.requiredKMers
   
   This produces a file ../tmp2/GS_nextGen/GRAPH_NAME/SAMPLENAME_graph.txt.kmers_31.requiredKMers.binaryCount.
   
   (The rationale for the copying is that readCortexCoverage.pl writes to a file whose name is identical to that of the file which specifies which kMers we are interested in, with the suffix ".binaryCount" appended. We need one such file for the reference genome and one such file for each sample.)
   
Find out how many of the reference kMers come from the region we are modelling:

	Take a FASTA file of the reference genome sequence of the region which your haplotype graph covers, and subtract the kMer counts (reverse-complemented and non-complemented kMers are counted as identical) observed in this file from the counts ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers.binaryCount. Store the corrected counts in ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers.binaryCount.corrected.
	
	(Obviously, we do not want to ignore kMers which come from our region of interest.)

Apply the model:

   ../bin/PnPHaploGraph2 domode nextGenInference ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31 ../tmp2/GS_nextGen/GRAPH_NAME/graph.txt.kmers_31.requiredKMers.binaryCount.corrected ../tmp2/GS_nextGen/GRAPH_NAME/SAMPLENAME_graph.txt.kmers_31.requiredKMers.binaryCount  --labelonly --genotypingMode 8
	
   This will produce, for example, a file ../tmp2/GS_nextGen/GRAPH_NAME/SAMPLENAME_graph.txt.kmers_31.requiredKMers.binaryCount.haplotypes, which you can analyze.
	
Automating the pipeline:

   At some point, we will probably create a fully automated pipeline that can deal with arbitrary genomic reasons. Until then, nextGenInferenceVariGraph.pl can serve as an example for how to post-process the output data to, for example, generate a VCF.
 

	
	

   

	













  
 



