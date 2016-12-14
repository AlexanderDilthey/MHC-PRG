# HLA-PRG

HLA\*PRG is an algorithm for HLA type inference from next-generation sequencing data.

## Introduction

HLA*PRG infers HLA types by

1. Aligning reads putatively originating from the HLA genes to a Population Reference Graph (PRG) of the HLA genes
2. Evaluating the graph-aligned reads in a simple likelihood framework - for each gene we find the pair of underlying alleles that maximizes the probability of observing the reads aligned to the locus

HLA types are inferred at 6-digit "G" resolution. 6-digit "G" resolution specifies the sequence of the exons encoding the antigen-binding site of the HLA protein - (i.e. exons 2 and 3 for HLA class I genes and exon 2 for HLA class II genes) - all alleles with identical sequences over these exons end up in the same "G" group. 6-digit "G" resolution is a useful metric to consider - most biological variability associated with the HLA genes is thought to come from differences in peptide binding, and this in turn is governed mostly by the amino acid sequence of the peptide binding site.

In our evaluations (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151), we find HLA\*PRG to be highly accurate if applied to high-quality (PCR-free, >30x coverage) whole-genome sequencing data. We don't recommend applying HLA*PRG to exome sequencing data.

The inference process always begins with a (GRCh37, see below) BAM file and consists of four steps:

1. Positive selection - extract read pairs that might come from the HLA genes (according to some k-Mer metrics)
2. Negative selection - remove read pairs from the set so-created that are likely to come from other genomic locations (again, using some k-Mer metrics)
3. Read-to-PRG alignment
4. HLA type inference

The rationale for splitting the read selection process into two separate components is that the two steps have different resource footprints - positive selection is typically applied to individual BAM files (e.g. on a cluster) and not very memory-intensive. Negative selection, by contrast, is memory-intensive and typically applied to batches of samples (if applicable; this is only for computational reasons - negative selection is independent-per-sample, statistically speaking).

### GRCh37

Currently HLA\*PRG should only be applied to B37-aligned BAM files. We plan to update the pipeline to deal with GRCh38-based BAMs.

Any B37-based reference can be used, as long as the MHC is represented on only one contig.

As a workaround, we recommend extracting MHC reads from your B38 file and temporarily re-mapping these to a B37 reference. See "Speeding up read extraction" below for read extraction example commands.

IMPORTANT: GRCh38 contains multiple MHC contigs, and sometimes (e.g. for bwakit) also some genomic HLA sequences. Make sure that you include reads from all of these regions as well!


### Computational considerations

Applying HLA\*PRG is currently computationally intensive - we are working on an optimized version of the algorithms.

Analyzing the NA12878 Illumina PLatinum data (2 x 100bp, 55x) takes on the order of 55 CPU hours and about 75G of RAM (peak usage).

Analyzing the NA12878 1000 Genomes longer-read data (2 x 250bp, 63) takes on the order of 215 CPU hours and about 75G of RAM (peak usage).

Of note, negative selection and read-to-PRG alignment are the RAM-intensive steps - positive selection typically doesn't take more than 15G.

## Installation

Prerequisites: Boost 1.5.2, Bamtools, a compiler with support for openMP and C++11 (we use g++ 4.7.2 and have successfully tested on the 4.8 branch). For more details, see the notes in the MHC*PRG main readme file.

1. Follow the make process described for MHC\*PRG - you do NOT need the MHC-PRG data package.
2. Download the HLA\*PRG data package (http://www.well.ox.ac.uk/HLA-PRG.tar.gz, c. 100G when extracted) and extract it into ../tmp2/GS_nextGen/hla (relative to the MHC-PRG/src directory, which contains the HLAtypeinference.pl script).
3. Verify that the MHC\*PRG binary is functional by executing ../bin/MHC-PRG from the MHC-PRG/src directory - it is functional if it then complains that you didn't specify a valid mode.
4. Verify that HLAtypeinference.pl is functional by typing perl -c HLAtypeinference.pl - on many systems all required modules will already be present, but this might not be the case for you.

Done!

## HLA type inference

HLA type inference is carried out on a per-sample basis.

For each (indexed and sorted) input BAM, you need to specify a unique sample ID. All data associated with this sample will be stored in ../tmp/hla/$SAMPLEID - each sample will consume between 0.5 and 2G of data in ../tmp/hla.

To generate HLA types, you need to execute the following 4 commands from the MHC-PRG/src directory:

./HLAtypeinference.pl --actions p --sampleIDs SAMPLEID --BAMs /path/to/indexed/bam.bam --referenceGenome   /path/to/referenceGenome/as/one/fasta/file  
./HLAtypeinference.pl --actions n --sampleIDs SAMPLEID  
./HLAtypeinference.pl --actions a --sampleIDs SAMPLEID  
./HLAtypeinference.pl --actions i --sampleIDs SAMPLEID  

You can combine the actions in one command:

./HLAtypeinference.pl --actions pnai --sampleIDs SAMPLEID --BAMs /path/to/indexed/bam.bam --referenceGenome /path/to/referenceGenome/as/one/fasta/file  

... or specify multiple BAM files - concatenate the paths with the ',' character (you then also need to specify multiple --sampleIDs).

Some remarks:
- p, n, a, i stand for positive selection, negative selection, alignment, inference.
- if you don't provide explicit sample IDs, the script will try to guess them from the name of the BAM file. It is better to explicitly specify sample IDs.
- the --referenceGenome parameter must point to a single fasta file that contains the reference genome used for creating the input BAM files.

If the inference process was successful, output best-guess HLA types are in ../tmp/hla/$SAMPLEID/R1_bestguess.txt. There is a separate quality score associated with each allele, but our evaluations show that these are not very well-calibrated. The algorithm that extracts two best-guess alleles from a list of allele pairs with associated probabilities is identical to the one used in HLA*IMP:02. Full allele pair probability distributions for locus $LOCUS can be found in the file R1_PP_$LOCUS_pairs.txt.

## HiSeq 2 x 250bp reads

If you are dealing with 2 x 250bp HiSeq reads, add a --HiSeq250bp 1 to all calls of HLAtypeinference.pl. E.g.:

./HLAtypeinference.pl --actions pnai --sampleIDs SAMPLEID --BAMs /path/to/indexed/bam.bam --referenceGenome /path/to/referenceGenome/as/one/fasta/file --HiSeq250bp 1

## Speeding up read extraction

It is possible to speed up the process of read extraction by limiting the process to the MHC and to unmapped reads. For an input file called input.bam, you can do the following (and then apply HLA*PRG to extract_for_HLAPRG.bam):

samtools view -F 4 -bo extract_xMHC.bam input.bam 6:2800000-34000000  
samtools view -f 4 -bo extract_unmapped.bam input.bam  
samtools merge extract_for_HLAPRG.bam extract_xMHC.bam extract_unmapped.bam  
samtools index extract_for_HLAPRG.bam

This approach should work for both 2 x 100bp and 2 x 250bp HiSeq reads (we have repeated the Platinum validation experiments using this approach and found no reduction in accuracy, and the --HiSeq250bp switch activates a similar heuristic).

