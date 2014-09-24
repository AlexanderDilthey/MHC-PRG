//============================================================================
// Name        : MHC-PRG.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <assert.h>
#include <omp.h>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "MHC-PRG.h"
#include "LocusCodeAllocation.h"
#include "Data/HaplotypePanel.h"
#include "Data/GenotypePanel.h"
#include "Graph/Graph.h"
#include "Graph/Node.h"
#include "Graph/HMM.h"
#include "Graph/MultiHMM.h"
#include "Graph/AlphaHMM.h"
#include "NextGen/NextGen.h"
#include "NextGen/simulationSuite.h"

#include "NextGen/Validation.h"
#include "NextGen/HLAtypes.h"

#include "GraphAligner/GraphAligner.h"
#include "GraphAligner/AlignerTests.h"
#include "GraphAlignerUnique/UniqueAlignerTests.h"

#include "readFilter/readFilter.h"

#include "Utilities.h"

using namespace std;

void testing();
Config CONFIG;
double epsilon = 1.0e-7;

struct NodePair {
	Node* n1;
	Node* n2;
};


int main(int argc, char *argv[])
{
	CONFIG.threads = 40;
	CONFIG.quiet = true;
	
	vector<string> arguments (argv + 1, argv + argc + !argc);
	if((arguments.size() > 0) && (arguments.at(0) != "domode"))
		errEx("Please go into domode!");
	
	if((arguments.size() > 0) && (arguments.at(1) == "loadGraph"))
	{

	}
	else if((arguments.size() > 0) && (arguments.at(1) == "simulateGraph"))
	{
		string graph_file = arguments.at(2);
		Graph* g = new Graph();
		g->readFromFile(graph_file);
		g->simulateHaplotypes(500000);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "gappifyGraph"))
	{
		string graph_file = arguments.at(2);
		Graph* g = new Graph();
		g->readFromFile(graph_file);
		g->makeEdgesGaps(0.1);
		g->writeToFile(graph_file+".gaps");
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "determineRequiredKMers"))
	{
		omp_set_num_threads(CONFIG.threads);

		string graph_file = arguments.at(2);
		string output_file = arguments.at(3);

		LargeGraph* kMerG = new LargeGraph();
		kMerG->readFromFile(graph_file);

		vector<string> requiredKMers = kMerG->requiredKMers();

		ofstream output;
		output.open (output_file.c_str(), ios::out | ios::trunc);
		if (output.is_open())
		{
			for(int i = 0; i < (int)requiredKMers.size(); i++)
			{
				output << requiredKMers.at(i);
				if(i != (requiredKMers.size()-1))
				{
					output << "\n";
				}
			}
		}
		else
		{
			errEx("Cannot open for writing: "+output_file);
		}
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "describekMerGraph"))
	{
		string graph_file = arguments.at(2);
		string temp_dir = arguments.at(3);
		string temp_label = arguments.at(4);

		bool output_kMer_levels = false;
		std::string referenceGenomeCortexGraph;

		int k = 31;
		
		for(unsigned int i = 5; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--output_kMer_levels")
			{
				output_kMer_levels = true;
			}
			if(arguments.at(i) == "--referenceGenomeCortexGraph")
			{
				referenceGenomeCortexGraph = arguments.at(i+1);
			}
			if(arguments.at(i) == "--k")
			{
				k = Utilities::StrtoI(arguments.at(i+1));
			}			
		}

		if(referenceGenomeCortexGraph.length())
		{
			if(! output_kMer_levels)
			{
				throw std::runtime_error("If you specify --referenceGenomeCortexGraph, please activate also --output_kMer_levels");
			}
		}

		if(k == 25)
		{
			describeGraph_25(graph_file, temp_dir, temp_label, output_kMer_levels, referenceGenomeCortexGraph);
		}
		else
		{
			assert(k == 31);
			describeGraph(graph_file, temp_dir, temp_label, output_kMer_levels, referenceGenomeCortexGraph);
		}
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "plotGraph"))
	{
		string graph_file = arguments.at(2);
		string start_level_str = arguments.at(3);
		string stop_level_str = arguments.at(4);
		string graphviz_output_file = arguments.at(5);

		int start_level = Utilities::StrtoI(start_level_str);
		int stop_level = Utilities::StrtoI(stop_level_str);
		assert(stop_level > start_level);

		Graph* g = new Graph();
		g->readFromFile(graph_file);
		g->graphViz(start_level, stop_level, graphviz_output_file);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "describeNucleotideGraph"))
	{
		string graph_file = arguments.at(2);
		string temp_dir = arguments.at(3);
		string temp_label = arguments.at(4);

		describeNucleotideGraph(graph_file, temp_dir, temp_label);
	}	
	else if((arguments.size() > 0) && (arguments.at(1) == "simulationSuite"))
	{
		int threads = 40;
		int genotypingMode = -1;

		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--threads")
			{
				threads = Utilities::StrtoI(arguments.at(i+1));
			}
			if(arguments.at(i) == "--genotypingMode")
			{
				genotypingMode = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		omp_set_num_threads(threads);

		CONFIG.threads=threads;

		string graph_file = arguments.at(2);
		string temp_dir = arguments.at(3);
		string temp_label = arguments.at(4);

		simulationSuite(graph_file, temp_dir, temp_label, genotypingMode);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "testGraphAligner"))
	{
		//std::cout << "domode testGraphAligner!\n" << std::flush;
		//testGraphAligner();
		//testDiagonalInverseExtendAlignment();
		//testChainFindingAndExtension();
		//testExtensionStop();
		//testSeedAndExtend_Algorithm();
		//test_Simple_longRange_SeedAndExtend();
		//test_Simple_longRange_SeedAndExtend_2();
		//test_Simple_longRange_SeedAndExtend_3();
		//test_Simple_longRange_SeedAndExtend_4();
		// testSeedAndExtend();
		//GraphAlignerUnique::tests::testChains();
		//GraphAlignerUnique::tests::testSeedAndExtend();
		//GraphAlignerUnique::tests::testSeedAndExtend_local();
		GraphAlignerUnique::tests::testSeedAndExtend_short();
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "testGraphAligner_realGraph"))
	{
		string graph_file = arguments.at(2);

		std::string filename_qualityMatrix = "GraphAlignerUnique/predefinedQualityMatrices/I100.txt";

		double insertSize_mean = 200;
		double insertSize_sd = 15;
		int read_length = 100;

		for(unsigned int i = 5; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--filename_qualityMatrix")
			{
				filename_qualityMatrix = arguments.at(i+1);
			}

			if(arguments.at(i) == "--insertSize_mean")
			{
				insertSize_mean = Utilities::StrtoD(arguments.at(i+1));
			}
			if(arguments.at(i) == "--insertSize_sd")
			{
				insertSize_sd = Utilities::StrtoD(arguments.at(i+1));
			}

			if(arguments.at(i) == "--readLength")
			{
				read_length = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		assert(insertSize_mean > 0);
		assert(insertSize_mean < 1000);
		assert(insertSize_sd >= 0);
		assert(insertSize_sd <= 20);

		GraphAlignerUnique::tests::testSeedAndExtend_local_realGraph(graph_file, read_length, insertSize_mean, insertSize_sd, filename_qualityMatrix);
	}

	else if((arguments.size() > 0) && (arguments.at(1) == "nextGenValidationTest"))
	{
		string viterbi_diploid_gS = "../tmp/kMerCount__GS_nextGen_varigraph_new_AA02O9Q_Z1_31_required.binaryCount.viterbiGenomeString";
		string deBruijnGraph = "Test/test_refgenome_31.ctx";
		string VCFfile = "Test/test_variants.vcf";
		string referenceGenomeFile = "Test/test_refgenome.fa";

		testValidation(viterbi_diploid_gS, deBruijnGraph, 31, VCFfile, referenceGenomeFile);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "nextGenValidation"))
	{
		bool labelOnly = false;
		omp_set_num_threads(40);
		CONFIG.threads=40;

		int cortex_height = 26;
		int cortex_width = 140;
		int kMer_size = 31;

		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--cortex_height")
			{
				cortex_height = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--cortex_width")
			{
				cortex_width = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		string graph_file = arguments.at(2);
		string sampleCount = arguments.at(3);
		int chromotypes_start = Utilities::StrtoI(arguments.at(4));
		int chromotypes_stop = Utilities::StrtoI(arguments.at(5));
		string vcfFile = arguments.at(6);
		int vcf_start = Utilities::StrtoI(arguments.at(7));
		int vcf_stop = Utilities::StrtoI(arguments.at(8));
		string referenceGenome = arguments.at(9);
		string deBruijnGraph = arguments.at(10);
		string outputDirectory = arguments.at(11);
		string contigsFile_fasta = arguments.at(12);
		std::string graphDir = arguments.at(13);

		std::cout << "Amended haplotypes versus VCF...\n==============================================================\n" << std::flush;
		
		string output_file_amendedHaplotypes = sampleCount+".amendedHaplotypes";
		string output_file_viterbi_genomeString = sampleCount+".viterbiGenomeString";

		validateAllChromotypesVsVCF(output_file_viterbi_genomeString, output_file_amendedHaplotypes, chromotypes_start, chromotypes_stop, vcfFile, vcf_start, vcf_stop, referenceGenome, deBruijnGraph, kMer_size, cortex_height, cortex_width, outputDirectory, graphDir);

		//		validateAmendedChromotypesVsVCF(output_file_amendedHaplotypes, chromotypes_start, chromotypes_stop, vcfFile, vcf_start, vcf_stop, referenceGenome, deBruijnGraph, kMer_size, cortex_height, cortex_width);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "alignShortReadsToHLAGraph"))
	{
		bool labelOnly = false;
		// omp_set_num_threads(40);
		CONFIG.threads=40;


		std::string input_FASTQ;
		std::string graph_dir;
		std::string referenceGenome;

		for(unsigned int i = 0; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--input_FASTQ")
			{
				input_FASTQ = arguments.at(i+1);
			}

			if(arguments.at(i) == "--graphDir")
			{
				graph_dir = arguments.at(i+1);
			}

			if(arguments.at(i) == "--referenceGenome")
			{
				referenceGenome = arguments.at(i+1);
			}
		}

		assert(input_FASTQ.length());
		assert(graph_dir.length());
		assert(referenceGenome.length());

		std::vector<std::pair<double, double>> inserSize_mean_sd_perFile;
		estimateInsertSizeFromGraph(input_FASTQ, graph_dir, inserSize_mean_sd_perFile);

		alignShortReadsToHLAGraph(input_FASTQ, graph_dir, referenceGenome, inserSize_mean_sd_perFile);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "simulateHLAreads"))
	{
		std::string graph_dir;
		std::string filename_qualityMatrix = "GraphAlignerUnique/predefinedQualityMatrices/I100.txt";
		std::string outputDirectory;

		int nIndividuals = 10;
		double insertSize_mean = 200;
		double insertSize_sd = 15;
		double haploidCoverage = 15;

		int read_length = 100;
		bool readError = true;
		bool perturbHaplotypes = false;
		bool exon23 = true;
		
		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--graphDir")
			{
				graph_dir = arguments.at(i+1);
			}

			if(arguments.at(i) == "--filename_qualityMatrix")
			{
				filename_qualityMatrix = arguments.at(i+1);
			}

			if(arguments.at(i) == "--outputDirectory")
			{
				outputDirectory = arguments.at(i+1);
			}

			if(arguments.at(i) == "--nIndividuals")
			{
				nIndividuals = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--insertSize_mean")
			{
				insertSize_mean = Utilities::StrtoD(arguments.at(i+1));
			}
			if(arguments.at(i) == "--insertSize_sd")
			{
				insertSize_sd = Utilities::StrtoD(arguments.at(i+1));
			}

			if(arguments.at(i) == "--haploidCoverage")
			{
				haploidCoverage = Utilities::StrtoD(arguments.at(i+1));
			}

			if(arguments.at(i) == "--readLength")
			{
				read_length = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--readError")
			{
				readError = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--perturbHaplotypes")
			{
				perturbHaplotypes = Utilities::StrtoI(arguments.at(i+1));
			}
			
			if(arguments.at(i) == "--exon23")
			{
				exon23 = Utilities::StrtoI(arguments.at(i+1));
			}			
		}

		assert(graph_dir.length());
		assert(outputDirectory.length());

		assert(insertSize_mean > 0);
		assert(insertSize_mean < 1000);
		assert(insertSize_sd >= 0);
		assert(insertSize_sd <= 20);
		assert(nIndividuals > 0);

		simulateHLAreads(graph_dir, nIndividuals, exon23, perturbHaplotypes, readError, outputDirectory, filename_qualityMatrix, read_length, insertSize_mean, insertSize_sd, haploidCoverage);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "HLATypeInference"))
	{
		bool labelOnly = false;
		CONFIG.threads=32;
		omp_set_num_threads(CONFIG.threads);
				
		std::string input_alignedReads;
		std::string graph_dir;
		std::string sampleID;

		for(unsigned int i = 0; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--input_alignedReads")
			{
				input_alignedReads = arguments.at(i+1);
			}

			if(arguments.at(i) == "--graphDir")
			{
				graph_dir = arguments.at(i+1);
			}

			if(arguments.at(i) == "--sampleID")
			{
				sampleID = arguments.at(i+1);
			}
		}

		assert(input_alignedReads.length());
		assert(graph_dir.length());

 		std::string loci_string;
		std::string starting_haplotypes_perLocus_1_str;
		std::string starting_haplotypes_perLocus_2_str;

		HLATypeInference(input_alignedReads, graph_dir, sampleID, true, loci_string, starting_haplotypes_perLocus_1_str, starting_haplotypes_perLocus_2_str);

 		// loci_string = "A";
		// starting_haplotypes_perLocus_1_str = "A*02:07:01;A*02:265";
		// starting_haplotypes_perLocus_2_str = "A*02:06:01";

 		// loci_string = "A";
		// starting_haplotypes_perLocus_1_str = "A*02:03:01;A*02:264";
		// starting_haplotypes_perLocus_2_str = "A*24:02:01:01;A*24:02:01:02L;A*24:02:03Q;A*24:02:10;A*24:09N;A*24:11N";
			
		HLAHaplotypeInference(input_alignedReads, graph_dir, sampleID, loci_string, starting_haplotypes_perLocus_1_str, starting_haplotypes_perLocus_2_str);
		
		// nodo activate
		
		HLATypeInference(input_alignedReads, graph_dir, sampleID, false, loci_string, starting_haplotypes_perLocus_1_str, starting_haplotypes_perLocus_2_str);
		
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "nextGenContigValidation"))
	{
		bool labelOnly = false;
		omp_set_num_threads(40);
		CONFIG.threads=40;

		int cortex_height = 26;
		int cortex_width = 95;
		int kMer_size = 31;

		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--cortex_height")
			{
				cortex_height = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--cortex_width")
			{
				cortex_width = Utilities::StrtoI(arguments.at(i+1));
			}

			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		string graph_file = arguments.at(2);
		string sampleCount = arguments.at(3);
		int chromotypes_start = Utilities::StrtoI(arguments.at(4));
		int chromotypes_stop = Utilities::StrtoI(arguments.at(5));
		string vcfFile = arguments.at(6);
		int vcf_start = Utilities::StrtoI(arguments.at(7));
		int vcf_stop = Utilities::StrtoI(arguments.at(8));
		string referenceGenome = arguments.at(9);
		string deBruijnGraph = arguments.at(10);
		string outputDir_contigs = arguments.at(11);
		string contigsFile_fasta = arguments.at(12);
		std::string graphDir = arguments.at(13);

		string output_file_amendedHaplotypes = sampleCount+".amendedHaplotypes";
		string output_file_viterbi_genomeString = sampleCount+".viterbiGenomeString";

		alignContigsToAllChromotypes(output_file_viterbi_genomeString, output_file_amendedHaplotypes, chromotypes_start, chromotypes_stop, vcfFile, vcf_start, vcf_stop, referenceGenome, deBruijnGraph, kMer_size, cortex_height, cortex_width, outputDir_contigs, contigsFile_fasta, graphDir);

		//		validateAmendedChromotypesVsVCF(output_file_amendedHaplotypes, chromotypes_start, chromotypes_stop, vcfFile, vcf_start, vcf_stop, referenceGenome, deBruijnGraph, kMer_size, cortex_height, cortex_width);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "nextGenInference"))
	{

		//omp_set_num_threads(CONFIG.threads);
		bool labelOnly = false;
		omp_set_num_threads(40);
		CONFIG.threads=40;

		/*
		map<int, int> edge1;
		map<int, int> edge2;
		edge1[1] = 3;

		edge2[1] = 1;
		edge2[2] = 1;
		edge2[3] = 1;

		map<int, int> thresholds;
		double P_above_all_edgePresent;
		double P_above_all_edgeNotPresent;

		computeThresholds(edge1, 0.05, 20, thresholds, P_above_all_edgePresent, P_above_all_edgeNotPresent);
		cout << "Edge1:\n";
		for(map<int, int>::iterator tIt = thresholds.begin(); tIt != thresholds.end(); tIt++)
		{
			cout << "\tThreshold for kMer " << tIt->first << ": " << tIt->second << "\n";
		}
		cout << "\tP_above_all_edgePresent: " << P_above_all_edgePresent << "\n";
		cout << "\tP_above_all_edgeNotPresent: " << P_above_all_edgeNotPresent << " => power: " << (1-P_above_all_edgeNotPresent) << "\n\n";

		computeThresholds(edge2, 0.05, 20, thresholds, P_above_all_edgePresent, P_above_all_edgeNotPresent);
		cout << "Edge2:\n";
		for(map<int, int>::iterator tIt = thresholds.begin(); tIt != thresholds.end(); tIt++)
		{
			cout << "\tThreshold for kMer " << tIt->first << ": " << tIt->second << "\n";
		}
		cout << "\tP_above_all_edgePresent: " << P_above_all_edgePresent << "\n";
		cout << "\tP_above_all_edgeNotPresent: " << P_above_all_edgeNotPresent << " => power: " << (1-P_above_all_edgeNotPresent) << "\n\n";

		errEx("OK");
		*/

		string graph_file = arguments.at(2);
		string wholeGenomeButGraphCount = arguments.at(3);
		string sampleCount = arguments.at(4);

		string test1 = "A||B";
		vector<string> test1_split = Utilities::split(test1, "||");
		assert(test1_split.size() == 2);
		assert(test1_split.at(0) == "A");
		assert(test1_split.at(1) == "B");

		vector<string> test1_split_2 = Utilities::split(test1, "*");
		assert(test1_split_2.size() == 1);


		string test2 = "A||B||C";
		vector<string> test2_split = Utilities::split(test2, "||");
		assert(test2_split.size() == 3);

		string test3 = "A||B||C||dtdtdtd|dmdmd||nd||a||dkdkd|ksks";
		vector<string> test3_split = Utilities::split(test3, "||");
		assert(test3_split.size() == 7);

		int genotypingMode = 5;
		for(unsigned int i = 5; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--labelonly")
			{
				labelOnly = (bool) Utilities::StrtoI(arguments.at(i));
			}

			if(arguments.at(i) == "--genotypingMode")
			{
				genotypingMode = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		LargeGraph* kMerG = new LargeGraph();
		kMerG->readFromFile(graph_file);
		// kMerG->stats();

		/* more complicated model

		MultiGraph* multiG = multiBeautify(kMerG, wholeGenomeButGraphCount);
		multiG->stats();

		MultiHMM* mHMM = new MultiHMM(multiG);
		map<string, long long> estimatedEmissions = mHMM->estimateEmissions(sampleCount, wholeGenomeButGraphCount);
		mHMM->fillForwardBackwardTable(estimatedEmissions);

		*/



		MultiGraph* multiG = multiBeautifyForAlpha2(kMerG, "");
		set<int> emptySet;

		map<string, long long> estimatedEmissions = AlphaHMM::estimateEmissions(sampleCount, wholeGenomeButGraphCount, multiG->kMerSize);

		MultiGraph* multiGsimple = simplifyAccordingToCoverage(multiG, estimatedEmissions, emptySet);

		// multiG->stats();

		AlphaHMM* aHMM = new AlphaHMM(multiGsimple, genotypingMode);
		aHMM->fillForwardBackwardTable(estimatedEmissions);
		// aHMM->debug(estimatedEmissions);

		string output_file = sampleCount+".inference";
		string output_file_haplotypes = sampleCount+".haplotypes";
		string output_file_viterbi_haplotypes = sampleCount+".viterbiHaplotypes";
		string output_file_viterbi_genomeString = sampleCount+".viterbiGenomeString";

		vector<multiHaploLabelPair> samples;
		vector<multiHaploLabelPair> viterbiSamples;
		
		vector<diploidNucleotidePath> samplesAlignedHaplotypes;
		vector<diploidNucleotidePath> viterbiSamplesAlignedHaplotypes;

		multiHaplotypePair_and_P viterbiSample = aHMM->retrieveViterbiSample(labelOnly);
		diploidEdgePointerPath viterbiEdgePointers = viterbiSample.diploidEdgePointers;
		viterbiSamples.push_back(multiGsimple->edgePointerPathToLabels(viterbiEdgePointers));
		diploidEdgePointerPath viterbi_sample_lG = multiGsimple->MultiGraphEdgesToLargeGraphEdges(viterbiEdgePointers);
		diploidNucleotidePath viterbiNucleotidePath = kMerG->diploidPathToAlignedNucleotides(viterbi_sample_lG);
		viterbiSamplesAlignedHaplotypes.push_back(viterbiNucleotidePath);	
		
		diploidGenomeString viterbiGenomeString = kMerG->diploidPathToGenomeString(viterbi_sample_lG);
		diploidGenomeString viterbiGenomeString_compressed = compressGenomeString(viterbiGenomeString);
		storeGenomeStringInFile(viterbiGenomeString_compressed, output_file_viterbi_genomeString);   

		cout << "Viterbi P: " << viterbiSample.P << "\n";

		for(int sI = 0; sI < 10; sI++)
		{
			diploidEdgePointerPath sample_mG = aHMM->sampleFromPosterior(estimatedEmissions);
			samples.push_back(multiGsimple->edgePointerPathToLabels(sample_mG));

			diploidEdgePointerPath sample_lG = multiGsimple->MultiGraphEdgesToLargeGraphEdges(sample_mG);
			diploidNucleotidePath nucleotidePath = kMerG->diploidPathToAlignedNucleotides(sample_lG);
			samplesAlignedHaplotypes.push_back(nucleotidePath);
		}

		set<string> loci;
		vector< map<string, string> > h1_found;
		vector< map<string, string> > h2_found;

		for(int sample = 0; sample < (int)samples.size(); sample++ )
		{
			map<string, string> h1_localfound;
			map<string, string> h2_localfound;

			for(int samplePart = 0; samplePart < (int)samples.at(sample).h1.size(); samplePart++)
			{
				string sampleString_out =  samples.at(sample).h1.at(samplePart);
				// cout << "h1: " << sampleString_out << "\n";

				vector<string> parts = Utilities::split(sampleString_out, "||");
				for(int sI = 0; sI < (int)parts.size(); sI++)
				{
					string sampleString = parts.at(sI);

					int found=sampleString.find("*");
					if (found==(int)string::npos)
						errEx("Cannot determine locus in "+sampleString);

					//string locus = "HLA"+sampleString.substr(0,found);
					//string allele_colons = sampleString.substr(found+1);
					//assert(allele_colons.length() >= 5);
					//string allele = allele_colons.substr(0,2)+allele_colons.substr(3,2);

					string locus = sampleString.substr(0,found);
					string allele = sampleString.substr(found+1);
					int found_sep_allele = allele.find("||");
					if(found_sep_allele != (int)string::npos)
					{
						cout << "Allele has || symbol -- bad!\n\nallele: " << allele << "\n" << "parts.size(): " << parts.size() << "\nOriginal samplestring_out: " << sampleString_out << "\n";
					}
					assert(found_sep_allele == (int)string::npos);

					loci.insert(locus);
					if(h1_localfound.count(locus) > 0)
					{
						errEx("Same locus "+locus+" appears more than once in output?\nsampleString_out: "+sampleString_out+"\nlocus: "+locus+"\nallele: "+allele);
					}
					h1_localfound[locus] = allele;
				}
			}
			for(int samplePart = 0; samplePart < (int)samples.at(sample).h2.size(); samplePart++)
			{
				string sampleString_out =  samples.at(sample).h2.at(samplePart);
				// cout << "h2: " << sampleString_out << "\n";

				vector<string> parts = Utilities::split(sampleString_out, "||");
				for(int sI = 0; sI < (int)parts.size(); sI++)
				{
					string sampleString = parts.at(sI);

					int found=sampleString.find("*");
					if (found==(int)string::npos)
						errEx("Cannot determine locus in "+sampleString);

					//string locus = "HLA"+sampleString.substr(0,found);
					//string allele_colons = sampleString.substr(found+1);
					//assert(allele_colons.length() >= 5);
					//string allele = allele_colons.substr(0,2)+allele_colons.substr(3,2);

					string locus = sampleString.substr(0,found);
					string allele = sampleString.substr(found+1);
					int found_sep_allele = allele.find("||");
					if(found_sep_allele != (int)string::npos)
					{
						cout << "Allele has || symbol -- bad!\n\nallele: " << allele << "\n" << "parts.size(): " << parts.size() << "\nOriginal samplestring_out: " << sampleString_out << "\n";
					}
					assert(found_sep_allele == (int)string::npos);

					loci.insert(locus);
					if(h2_localfound.count(locus) > 0)
					{
						errEx("Same locus "+locus+" appears more than once in output?\nsampleString_out: "+sampleString_out+"\nlocus: "+locus+"\nallele: "+allele);
					}
					h2_localfound[locus] = allele;
				}
			}
			h1_found.push_back(h1_localfound);
			h2_found.push_back(h2_localfound);
		}

		vector<string> loci_in_order(loci.begin(), loci.end());

		ofstream output;
		output.open (output_file.c_str(), ios::out | ios::trunc);
		if (output.is_open())
		{
			output << "IndividualID Chromosome " << Utilities::join(loci_in_order, " ") << "\n";
			for(int sample = 0; sample < (int)samples.size(); sample++ )
			{
				string individualID = sampleCount+"-x-"+Utilities::ItoStr(sample+1);
				vector<string> fields_f1;
				vector<string> fields_f2;

				fields_f1.push_back(individualID);
				fields_f2.push_back(individualID);
				fields_f1.push_back(Utilities::ItoStr(1));
				fields_f2.push_back(Utilities::ItoStr(2));

				for(int fieldI = 0; fieldI < (int)loci_in_order.size(); fieldI++)
				{
					string fieldName = loci_in_order.at(fieldI);
					if(h1_found.at(sample).count(fieldName) > 0)
					{
						fields_f1.push_back(h1_found.at(sample)[fieldName]);
					}
					else
					{
						fields_f1.push_back("?");
					}

					if(h2_found.at(sample).count(fieldName) > 0)
					{
						fields_f2.push_back(h2_found.at(sample)[fieldName]);
					}
					else
					{
						fields_f2.push_back("?");
					}
				}

				assert(fields_f1.size() == fields_f2.size());

				output << Utilities::join(fields_f1, " ") << "\n";
				output << Utilities::join(fields_f2, " ") << "\n";
			}
		}
		else
		{
			errEx("Cannot open for writing: "+output_file);
		}

		ofstream output_haplotypes;
		output_haplotypes.open (output_file_haplotypes.c_str(), ios::out | ios::trunc);
		if (output_haplotypes.is_open())
		{
			for(int sample = 0; sample < (int)samples.size(); sample++ )
			{
				string individualID = sampleCount+"-x-"+Utilities::ItoStr(sample+1);
				vector<string> fields_f1;
				vector<string> fields_f2;

				fields_f1.push_back(individualID);
				fields_f2.push_back(individualID);
				fields_f1.push_back(Utilities::ItoStr(1));
				fields_f2.push_back(Utilities::ItoStr(2));

				for(int fieldI = 0; fieldI < (int)samplesAlignedHaplotypes.at(sample).h1.size(); fieldI++)
				{
					fields_f1.push_back(samplesAlignedHaplotypes.at(sample).h1.at(fieldI));
					fields_f2.push_back(samplesAlignedHaplotypes.at(sample).h2.at(fieldI));
				}

				assert(fields_f1.size() == fields_f2.size());

				output_haplotypes << Utilities::join(fields_f1, " ") << "\n";
				output_haplotypes << Utilities::join(fields_f2, " ") << "\n";
			}
		}
		else
		{
			errEx("Cannot open for writing: "+output_file);
		}
		
		ofstream output_viterbi_haplotypes;
		output_viterbi_haplotypes.open (output_file_viterbi_haplotypes.c_str(), ios::out | ios::trunc);
		if (output_viterbi_haplotypes.is_open())
		{
			for(int sample = 0; sample < (int)viterbiSamples.size(); sample++ )
			{
				string individualID = sampleCount+"-x-"+Utilities::ItoStr(sample+1);
				vector<string> fields_f1;
				vector<string> fields_f2;

				fields_f1.push_back(individualID);
				fields_f2.push_back(individualID);
				fields_f1.push_back(Utilities::ItoStr(1));
				fields_f2.push_back(Utilities::ItoStr(2));

				for(int fieldI = 0; fieldI < (int)viterbiSamplesAlignedHaplotypes.at(sample).h1.size(); fieldI++)
				{
					fields_f1.push_back(viterbiSamplesAlignedHaplotypes.at(sample).h1.at(fieldI));
					fields_f2.push_back(viterbiSamplesAlignedHaplotypes.at(sample).h2.at(fieldI));
				}

				assert(fields_f1.size() == fields_f2.size());

				output_viterbi_haplotypes << Utilities::join(fields_f1, " ") << "\n";
				output_viterbi_haplotypes << Utilities::join(fields_f2, " ") << "\n";
			}
		}
		else
		{
			errEx("Cannot open for writing: "+output_file_viterbi_haplotypes);
		}		

	}
	else if((arguments.size() > 0) && (arguments.at(1) == "loadKMerGraphHMM"))
	{
		string graph_file = arguments.at(2);

		LargeGraph* kMerG = new LargeGraph();
		kMerG->readFromFile(graph_file);

		kMerG->stats();
		errEx("OK");

		//LargeHMM* HMM = new LargeHMM(kMerG, true);

	}
	else if((arguments.size() > 0) && (arguments.at(1) == "kMerifyAndSimulate"))
	{
		string graph_file = arguments.at(2);
		Graph* g = new Graph();
		g->readFromFile(graph_file);
		LargeGraph* kMerG = kMerify(g, true);
		kMerG->simulateHaplotypes(500000);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "kMerify"))
	{
		string graph_file = arguments.at(2);

		int kMer_size = 55;
		bool protectPGF = true;
		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
			
			if(arguments.at(i) == "--noPGFprotection")
			{
				protectPGF = false;
			}
		}

		cout << "Using graph " << graph_file << ".\n";
		Graph* g = new Graph();
		g->readFromFile(graph_file);
		g->printComplexity(graph_file+".complexity");

		string kMerified_graph_file = graph_file+".kmers";
		kMerified_graph_file += ("_" + Utilities::ItoStr(kMer_size));
		LargeGraph* kMerG = kMerify(g, false, kMer_size, protectPGF);
		kMerG->writeToFile(kMerified_graph_file);

		kMerG->stats();

	}
	else if((arguments.size() > 0) && (arguments.at(1) == "createConcatenatedVariationGraphs"))
	{
		omp_set_num_threads(CONFIG.threads);

		string haplotypes_files_dir = arguments.at(2);
		string positions_file = haplotypes_files_dir+"/positions.txt";
		string graph_output_file = haplotypes_files_dir+"/graph.txt";

		int kMer_size = 55;
		bool protectPGF = true;
		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
			if(arguments.at(i) == "--noPGFprotection")
			{
				protectPGF = false;
			}			
		}

		string haplotypes_files_list = haplotypes_files_dir+"/segments.txt";

		cout << "Create variation graph based on haplotype specified in list files " << haplotypes_files_list << ", with positions" << positions_file << "\n\n" << flush;

		Graph* combinedGraph = new Graph();
		int combinedLevel = -1;
		ifstream filesStream;
		filesStream.open (haplotypes_files_list.c_str(), ios::in);

		set<Node*> previousGraphLastNodes;
		if(filesStream.is_open())
		{
			string line;

			vector<string> files;

			while(filesStream.good())
			{
				getline (filesStream, line);
				Utilities::eraseNL(line);

				if(line.length() == 0)
					continue;

				files.push_back(line);
			}

			for(int fI = 0; fI < (int)files.size(); fI++)
			{
				string fN = files.at(fI);

				string haplotypes_file = haplotypes_files_dir+"/"+fN;

				std::cout << "Reading " << haplotypes_file << "\n" << std::flush;

				Graph* g = variationGraph(haplotypes_file, positions_file, protectPGF);
				string output_file = haplotypes_file+".graph";
				g->writeToFile(output_file);

				if(fI < ((int)files.size() - 1))
				{
					g->CODE.removeLocus("END_PUFFER");
				}
				
				combinedGraph->CODE.takeLocusData(g->CODE);

				int levels = g->NodesPerLevel.size();
				int lastLevel = (fI == ((int)files.size() - 1)) ? g->NodesPerLevel.size() - 1 : g->NodesPerLevel.size() - 2;
				for(int l = 0; l <= lastLevel; l++)
				{
					for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(l).begin(); nodeIt != g->NodesPerLevel.at(l).end(); nodeIt++)
					{
						Node* n = *nodeIt;
						if(l == lastLevel)
						{
							//assert(n->Outgoing_Edges.size() == 0);
							if (fI != (files.size()-1))
							{
								previousGraphLastNodes.insert(n);
							}
							else
							{
								if(nodeIt == g->NodesPerLevel.at(l).begin())
									combinedLevel++;
								n->level = combinedLevel;
								n->Outgoing_Edges.clear();
								combinedGraph->registerNode(n, combinedLevel);
							}
						}
						else
						{
							if(nodeIt == g->NodesPerLevel.at(l).begin())
								combinedLevel++;

							n->level = combinedLevel;
							combinedGraph->registerNode(n, combinedLevel);

							for(set<Edge*>::iterator outgoingEdgesIt = n->Outgoing_Edges.begin(); outgoingEdgesIt != n->Outgoing_Edges.end(); outgoingEdgesIt++)
							{
								Edge* outgoingEdge = *outgoingEdgesIt;
								combinedGraph->registerEdge(outgoingEdge);
							}
						}
					}

					if(l == 0)
					{
						assert(g->NodesPerLevel.at(0).size() == 1);
						Node* n0 = *(g->NodesPerLevel.at(0).begin());
						for(set<Node*>::iterator lastNodeIt = previousGraphLastNodes.begin(); lastNodeIt != previousGraphLastNodes.end(); lastNodeIt++)
						{
							Node* lastLevelNode = *lastNodeIt;
							for(set<Edge*>::iterator incomingEdgesIt = lastLevelNode->Incoming_Edges.begin(); incomingEdgesIt != lastLevelNode->Incoming_Edges.end(); incomingEdgesIt++)
							{
								Edge* incomingEdge = *incomingEdgesIt;
								incomingEdge->To = n0;
								n0->Incoming_Edges.insert(incomingEdge);
							}
						}
						previousGraphLastNodes.clear();
					}
				}
			}
			filesStream.close();
		}

		combinedGraph->checkConsistency(false);

		combinedGraph->writeToFile(graph_output_file);
		combinedGraph->printComplexity(graph_output_file+".complexity");

		exit(0);

		string kMerified_graph_file = graph_output_file+".kmers";
		LargeGraph* kMerG = kMerify(combinedGraph, false, kMer_size);
		kMerG->writeToFile(kMerified_graph_file);

		kMerG->stats();
	}

	else if((arguments.size() > 0) && (arguments.at(1) == "createMHCVariationGraph"))
	{
		omp_set_num_threads(CONFIG.threads);


		string haplotypes_file = arguments.at(2);
		string positions_file = arguments.at(3);
		string graph_output_file = arguments.at(4);

		int kMer_size = 55;
		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		cout << "Create variation graph based on haplotypes from " << haplotypes_file << ", with positions" << positions_file << "\n\n" << flush;

		Graph* g = variationGraph(haplotypes_file, positions_file);
		g->writeToFile(graph_output_file);
		g->printComplexity(graph_output_file+".complexity");

		exit(0);

		string kMerified_graph_file = graph_output_file+".kmers";
		LargeGraph* kMerG = kMerify(g, false, kMer_size);
		kMerG->writeToFile(kMerified_graph_file);

		kMerG->stats();
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "createHLAAlleleshortAlignmentGraph"))
	{
		omp_set_num_threads(CONFIG.threads);

		string hla_allele_dir = arguments.at(2);
		string locus = arguments.at(3);
		string graph_output_file = arguments.at(4);

		int kMer_size = 55;
		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--kmer")
			{
				kMer_size = Utilities::StrtoI(arguments.at(i+1));
			}
		}

		cout << "Create graph for short HLA allele alignments for locus " << locus << ", using data from directory " << hla_allele_dir << "\n\n" << flush;

		Graph* g = HLAalignmentGraph(hla_allele_dir, locus);
		g->writeToFile(graph_output_file);

		string kMerified_graph_file = graph_output_file+".kmers";
		LargeGraph* kMerG = kMerify(g, false, kMer_size);
		kMerG->writeToFile(kMerified_graph_file);

		kMerG->stats();

		//LargeHMM* HMM = new LargeHMM(kMerG, true);
	}
	else if((arguments.size() > 0) && (arguments.at(1) == "filterReads"))
	{
		std::string positiveFilter;
		std::string negativeFilter;

		double positiveThreshold = 0.3;
		double negativeThreshold = 0.3;

		vector<string> arguments (argv + 1, argv + argc + !argc);

		std::string input_BAM;
		std::string input_FASTQ;

		std::string output_FASTQ;

		int k = 25;

		bool positiveUnique = false;
		bool negativePreserveUnique = false;

		std::string uniqueness_base;
		std::string uniqueness_subtract;

		for(unsigned int i = 2; i < arguments.size(); i++)
		{
			if(arguments.at(i) == "--positiveFilter")
			{
				positiveFilter = arguments.at(i+1);
			}
			if(arguments.at(i) == "--negativeFilter")
			{
				negativeFilter = arguments.at(i+1);
			}
			if(arguments.at(i) == "--input_BAM")
			{
				input_BAM = arguments.at(i+1);
			}
			if(arguments.at(i) == "--input_FASTQ")
			{
				input_FASTQ = arguments.at(i+1);
			}
			if(arguments.at(i) == "--output_FASTQ")
			{
				output_FASTQ = arguments.at(i+1);
			}
			if(arguments.at(i) == "--positiveThreshold")
			{
				positiveThreshold = Utilities::StrtoD(arguments.at(i+1));
			}
			if(arguments.at(i) == "--negativeThreshold")
			{
				negativeThreshold = Utilities::StrtoD(arguments.at(i+1));
			}
			if(arguments.at(i) == "--k")
			{
				k = Utilities::StrtoI(arguments.at(i+1));
			}
			if(arguments.at(i) == "--positiveUnique")
			{
				positiveUnique = true;
			}
			if(arguments.at(i) == "--negativePreserveUnique")
			{
				negativePreserveUnique = true;
			}
			if(arguments.at(i) == "--uniqueness_base")
			{
				uniqueness_base = arguments.at(i+1);
			}
			if(arguments.at(i) == "--uniqueness_subtract")
			{
				uniqueness_subtract = arguments.at(i+1);
			}
		}

		readFilter F;
		F.positiveFilter = positiveFilter;
		F.negativeFilter = negativeFilter;
		F.input_BAM = input_BAM;
		F.input_FASTQ = input_FASTQ;
		F.output_FASTQ = output_FASTQ;
		F.positiveThreshold = positiveThreshold;
		F.negativeThreshold = negativeThreshold;
		F.k = k;
		F.negativePreserveUnique = negativePreserveUnique;
		F.positiveUnique = positiveUnique;
		F.uniqueness_base = uniqueness_base;
		F.uniqueness_subtract = uniqueness_subtract;

		F.doFilter();

	}
	else
	{
		errEx("Please specify valid mode.");
	}
	return 0;
}

void errEx(std::string message)
{
	cerr << message << "\n";
	exit(1);
}

void testing()
{
	LocusCodeAllocation CODE;
	unsigned char codedHLA = CODE.doCode("HLA", "0101");
	string originalAllele = CODE.deCode("HLA", codedHLA);
	// std::cout << "Coded HLA: " << codedHLA << " -- original allele: " << originalAllele << "\n";
	assert(originalAllele == "0101");

}
