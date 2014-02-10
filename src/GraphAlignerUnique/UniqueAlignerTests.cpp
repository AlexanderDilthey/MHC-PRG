/*
 * UniqueAlignerTests.cpp
 *
 *  Created on: 30.07.2013
 *      Author: AlexanderDilthey
 */

#include "UniqueAlignerTests.h"

#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <exception>
#include <stdexcept>
#include "../Graph/Graph.h"
#include "../Utilities.h"
#include "../NextGen/Validation.h"

#include "../GraphAligner/GraphAlignerAffine.h"
#include "../GraphAligner/AlignerTests.h"
#include "GraphAlignerUnique.h"
#include "../NextGen/readSimulator.h"

#include <omp.h>

namespace GraphAlignerUnique  {
namespace tests {

void testSeedAndExtend()
{
	// this thing uses "mutation" and extension

	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	size_t achievableMatches = 0;
	size_t achievedMatches = 0;
	double sum_scores = 0;

	size_t rightPositionCharacters_counted = 0;
	size_t rightPositionCharacters_correct = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(800);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAlignerUnique gA(gS_graph, aligner_kMerSize);
		GraphAlignerUnique gA2(gS_graph, aligner_kMerSize);

		gA.setThreads(4);
		gA2.setThreads(1);
		gA2.setIterationsMainRandomizationLoop(0);

		GraphAligner_affine gA_classical(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,

					10,
					15,

					2,
					4,

					2,
					3,

					3,
					7,

					30,
					40,

					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

//			int minimumAchievableScore = gA_classical.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);
			int minimumAchievableMatches = gA_classical.countMatchesInSequence(underlyingEdges, underlyingEdges_levels, randomString);

//			if(verbose || true)
//				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore <<  " and " << minimumAchievableMatches << " matches of sequence string: " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;
//

			achievableMatches += minimumAchievableMatches;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

//			if(verbose)
//				std::cout << "Start full-string alignment...\n" << std::flush;
//

			seedAndExtend_return wholeString_alignments = gA.seedAndExtend(randomString_noGaps);
			seedAndExtend_return wholeString_alignments_2 = gA2.seedAndExtend(randomString_noGaps);

//			 std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "\n" << std::flush;

			{
				seedAndExtend_return& thisAlignment = wholeString_alignments;

				int thisAlignmentScore = gA.score(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);
				int thisAlignmentScore_2 = gA2.score(wholeString_alignments_2.graph_aligned, wholeString_alignments_2.graph_aligned_levels, wholeString_alignments_2.sequence_aligned);

				int thisAlignmentMatches = gA_classical.countMatchesInSequence(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);



//				if(true || verbose)
//				{
//					std::cout << "\tAlignment [internal score " << thisAlignment.Score << "]>\n";
//					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
//					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
//					std::cout << "\t\tAffine NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;
//				}
				assert(thisAlignment.Score == thisAlignmentScore);
				assert(wholeString_alignments_2.Score == thisAlignmentScore_2);

				assert(wholeString_alignments.Score >= wholeString_alignments_2.Score);
				sum_scores += thisAlignment.Score;

//				assert(thisAlignmentMatches >= minimumAchievableMatches);

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

//				if(verbose || true)
//					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

				achievedMatches += thisAlignmentMatches;

				rightPositionCharacters_counted += validatable_noGap_characters;
				rightPositionCharacters_correct += validatable_noGap_characters_OK;

//				assert( 1 == 0);
			}

			individualTests++;
		}

		delete(gS_graph);
	}

	assert(achievableMatches > 0);
	assert(rightPositionCharacters_counted > 0);

	std::cout << "testSeedAndExtend(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
	std::cout << "\t Matches: " << achievedMatches  << " / " << achievableMatches << " => " << ((double)achievedMatches/(double)achievableMatches) << "\n";
	std::cout << "\t Positions: " << rightPositionCharacters_correct  << " / " << rightPositionCharacters_counted << " => " << (double(rightPositionCharacters_correct)/(double)rightPositionCharacters_counted) << "\n";
	std::cout << "\t Sum of scores: " << sum_scores << "\n";

	std::cout << std::flush;
}


void testSeedAndExtend_local()
{
	// this thing uses "mutation" and extension

	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	size_t achievableMatches = 0;
	size_t achievedMatches = 0;
	double sum_scores = 0;

	size_t rightPositionCharacters_counted = 0;
	size_t rightPositionCharacters_correct = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(800);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAlignerUnique gA(gS_graph, aligner_kMerSize);
		GraphAlignerUnique gA2(gS_graph, aligner_kMerSize);

		gA.setThreads(4);
		gA2.setThreads(1);
		gA2.setIterationsMainRandomizationLoop(0);

		GraphAligner_affine gA_classical(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "[LOCAL] Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,

					10,
					15,

					2,
					4,

					2,
					3,

					3,
					7,

					30,
					40,

					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

//			int minimumAchievableScore = gA_classical.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);
			int minimumAchievableMatches = gA_classical.countMatchesInSequence(underlyingEdges, underlyingEdges_levels, randomString);

//			if(verbose || true)
//				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore <<  " and " << minimumAchievableMatches << " matches of sequence string: " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;
//

			achievableMatches += minimumAchievableMatches;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

//			if(verbose)
//				std::cout << "Start full-string alignment...\n" << std::flush;
//

			std::vector<seedAndExtend_return_local> allBacktraces;
			seedAndExtend_return_local wholeString_alignments_local = gA.seedAndExtend_local(randomString_noGaps, allBacktraces);
//
//			seedAndExtend_return wholeString_alignments = gA.seedAndExtend(randomString_noGaps);
			seedAndExtend_return wholeString_alignments_2 = gA2.seedAndExtend(randomString_noGaps);

//			 std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "\n" << std::flush;

			{
				seedAndExtend_return_local& thisAlignment = wholeString_alignments_local;

				int thisAlignmentScore = gA.score(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);

//				int thisAlignmentScore_2 = gA2.score(wholeString_alignments_2.graph_aligned, wholeString_alignments_2.graph_aligned_levels, wholeString_alignments_2.sequence_aligned);

				int thisAlignmentMatches = gA_classical.countMatchesInSequence(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);


				if(true || verbose)
				{
					std::cout << "\tLOCAL Alignment [internal score " << thisAlignment.Score << "]>\n";
					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
					std::cout << "\t\tLocal NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;

					std::cout << "\t\tCompare with (one, max-only) global alignment [global score " << wholeString_alignments_2.Score << "]>:\n";
					std::cout << "\t\t\t" << wholeString_alignments_2.graph_aligned << "\n";
					std::cout << "\t\t\t" << wholeString_alignments_2.sequence_aligned << "\n" << std::flush;
				}

				assert(thisAlignment.Score == thisAlignmentScore);
//
//			assert(wholeString_alignments_2.Score == thisAlignmentScore_2);
//			assert(wholeString_alignments.Score >= wholeString_alignments_2.Score);
				sum_scores += thisAlignment.Score;

//				assert(thisAlignmentMatches >= minimumAchievableMatches);

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

//				if(verbose || true)
//					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

				achievedMatches += thisAlignmentMatches;

				rightPositionCharacters_counted += validatable_noGap_characters;
				rightPositionCharacters_correct += validatable_noGap_characters_OK;

//				assert( 1 == 0);
			}

			individualTests++;
		}

		delete(gS_graph);
	}

	assert(achievableMatches > 0);
	assert(rightPositionCharacters_counted > 0);

	std::cout << "testSeedAndExtend_local(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
	std::cout << "\t Matches: " << achievedMatches  << " / " << achievableMatches << " => " << ((double)achievedMatches/(double)achievableMatches) << "\n";
	std::cout << "\t Positions: " << rightPositionCharacters_correct  << " / " << rightPositionCharacters_counted << " => " << (double(rightPositionCharacters_correct)/(double)rightPositionCharacters_counted) << "\n";
	std::cout << "\t Sum of scores: " << sum_scores << "\n";

	std::cout << std::flush;
}

void testChains()
{

	bool verbose = false;
	int individualTests = 0;
	int individualTests_successful = 0;

	for(unsigned int graphIteration = 1; graphIteration <= 1000; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(400);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAlignerUnique gA(gS_graph, aligner_kMerSize);
		// gA.getGI().printIndex();

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::vector<Edge*> underlyingEdges;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleExactStringFromGraph(
					gS_graph,
					10,
					100,
					randomString,
					underlyingEdges
			);

			if(verbose)
				std::cout << "... sampling done:\n" << "\t" << randomString << "\n" << std::flush;

			std::cout << "... sampling done:\n" << "\t" << randomString << "\n" << std::flush;

			Node* firstNode = underlyingEdges.front()->From;
			Node* lastNode = underlyingEdges.back()->To;

			std::string randomString_noGaps;
			for(unsigned int sI = 0; sI < randomString.size(); sI++)
			{
				if(randomString.at(sI) != '_')
				{
					randomString_noGaps.push_back(randomString.at(sI));
				}
			}

			bool oneGoodChain = false;
			std::vector<kMerEdgeChain*> chains_for_sequence = gA.getGI().findChains(randomString_noGaps);
			std::cout << "Found chains: " << chains_for_sequence.size() << "\n" << std::flush;
			for(unsigned int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
			{
				kMerEdgeChain* chain = chains_for_sequence.at(chainI);
				Node* chain_firstNode = chain->traversedEdges.front()->From;
				Node* chain_lastNode = chain->traversedEdges.back()->To;

				/*
				std::cout << "chainI: " << chainI << "\n";
				std::cout << "\t" << "chain_firstNode: " << chain_firstNode << " level: " << chain_firstNode->level << "\n";
				std::cout << "\t" << "chain_lastNode: " << chain_lastNode << " level: " << chain_lastNode->level << "\n";
				std::cout << "\t" << "firstNode: " << firstNode << " level: " << firstNode->level << "\n";
				std::cout << "\t" << "lastNode: " << lastNode  << " level: " << lastNode->level <<  "\n\n" << std::flush;
				*/

				if((chain_firstNode == firstNode) && (chain_lastNode == lastNode))
				{
					oneGoodChain = true;
					std::string impliedSequence;
					for(unsigned int eI = 0; eI < chain->traversedEdges.size(); eI++)
					{
						Edge* e = chain->traversedEdges.at(eI);
						std::string edgeEmission = gS_graph->CODE.deCode(e->locus_id, e->emission);
						assert(edgeEmission.length() == 1);
						if(edgeEmission != "_")
						{
							impliedSequence.append(edgeEmission);
						}
					}

					std::string extractedSequence(randomString_noGaps.begin() + chain->sequence_begin, randomString_noGaps.begin() + chain->sequence_end + 1);

					assert(extractedSequence == impliedSequence);

					break;
				}
			}

			assert(oneGoodChain);

			for(unsigned int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
			{
				delete(chains_for_sequence.at(chainI));
			}

			individualTests++;
			individualTests_successful++;
		}

		gS_graph->freeMemory();
		delete(gS_graph);
	}

	std::cout << "GraphAlignerUnique::tests::testChains(): " << individualTests << " tests, of which " << individualTests_successful << " were succesful!\n" << std::flush;
}


void testSeedAndExtend_local_realGraph(std::string graph_filename, int read_length, double insertSize_mean, double insertSize_sd, std::string qualityMatrixFile)
{
	readSimulator rS(qualityMatrixFile);

	double haploidCoverage = 30;
	int aligner_kMerSize = 25;

	std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Loading graph.\n" << std::flush;

	std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Create GraphAlignerUnique.\n" << std::flush;
	
	int outerThreads = 32;
	omp_set_num_threads(outerThreads);
	
	std::vector<Graph*> graphs;  
	std::vector<GraphAlignerUnique*> graphAligners;
	graphs.resize(outerThreads);
	graphAligners.resize(outerThreads);
	
	#pragma omp parallel for
	for(int tI = 0; tI < outerThreads; tI++)
	{
		// std::cout << "Thread " << tI << "\n" << std::flush;
		
		Graph* g = new Graph();
		g->readFromFile(graph_filename);
			
		graphAligners.at(tI) = new GraphAlignerUnique(g, aligner_kMerSize);
		graphAligners.at(tI)->setIterationsMainRandomizationLoop(1);
		graphAligners.at(tI)->setThreads(1);	
		
		graphs.at(tI) = g;
	}
	

	std::cout << Utilities::timestamp() << "\tdone.\n" << std::flush;	

	std::string fileName_for_incorrectAlignments = "../tmp/incorrectAlignments.txt";
	std::ofstream incorrectAlignmentsStream;
	incorrectAlignmentsStream.open(fileName_for_incorrectAlignments.c_str());
	assert(incorrectAlignmentsStream.is_open());

	int simulateGenomePairs = 1;
	for(int genomePair = 1; genomePair <= simulateGenomePairs; genomePair++)
	{
		std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Iteration " << genomePair << " / " << simulateGenomePairs << "\n" << std::flush;

		diploidEdgePointerPath diploidPath = graphs.at(0)->simulateRandomDiploidPath();

		std::vector<oneReadPair> simulatedReadPairs_h1 = rS.simulate_paired_reads_from_edgePath(diploidPath.h1, haploidCoverage, insertSize_mean, insertSize_sd, true);
		std::vector<oneReadPair> simulatedReadPairs_h2 = rS.simulate_paired_reads_from_edgePath(diploidPath.h2, haploidCoverage, insertSize_mean, insertSize_sd, true);

		std::vector<oneReadPair> combinedPairs_for_alignment;
		combinedPairs_for_alignment.insert(combinedPairs_for_alignment.end(), simulatedReadPairs_h1.begin(), simulatedReadPairs_h1.end());
		combinedPairs_for_alignment.insert(combinedPairs_for_alignment.end(), simulatedReadPairs_h2.begin(), simulatedReadPairs_h2.end());

		std::cout << "\t" << "Simulated " << combinedPairs_for_alignment.size() << " read pairs." << "\n" << std::flush;

		std::vector<oneReadPair> combinedPairs_for_alignment_filtered;
		for(unsigned int pI = 0; pI < combinedPairs_for_alignment.size(); pI++)
		{
			if((pI % 5) != 0)
			{
				continue;
			}
					
			oneReadPair& p = combinedPairs_for_alignment.at(pI);
			if(	(p.reads.first.sequence.find("*") == std::string::npos) && (p.reads.first.sequence.find("N") == std::string::npos) &&
				(p.reads.second.sequence.find("*") == std::string::npos) && (p.reads.second.sequence.find("N") == std::string::npos))
			{
				combinedPairs_for_alignment_filtered.push_back(p);
			}
		}

		std::cout << "\t" << "After filtering out N / *, have " << combinedPairs_for_alignment_filtered.size() << " read pairs." << "\n" << std::flush;

		auto checkOneReadLevelCorrectness = [&](oneRead& r, seedAndExtend_return_local& r_aligned, int& levels_total, int& levels_OK) -> void {
		
			levels_total = 0;
			levels_OK = 0;

			assert(r.coordinates_edgePath.size() == r.sequence.length());
			int cI_in_unaligned_sequence = -1;
			std::vector<std::string> correctLevelsInAlignmentOrder;

			for(unsigned int cI = 0; cI < r_aligned.sequence_aligned.size(); cI++)
			{
				char sequenceCharacter_from_alignment = r_aligned.sequence_aligned.at(cI);
				if(sequenceCharacter_from_alignment == '_')
				{
					// ignore
					correctLevelsInAlignmentOrder.push_back("/");
				}
				else
				{
					cI_in_unaligned_sequence++;
					int specifiedEdgeLevel = r_aligned.graph_aligned_levels.at(cI);


					int cI_in_unaligned_sequence_correctlyAligned = cI_in_unaligned_sequence;
					if(r_aligned.reverse)
					{
						cI_in_unaligned_sequence_correctlyAligned = (r.sequence.length() - cI_in_unaligned_sequence - 1);
					}
					assert(cI_in_unaligned_sequence_correctlyAligned >= 0);
					assert(cI_in_unaligned_sequence_correctlyAligned < (int)r.sequence.length());

					int correctEdgeLevel = r.coordinates_edgePath.at(cI_in_unaligned_sequence_correctlyAligned);
					correctLevelsInAlignmentOrder.push_back(Utilities::ItoStr(correctEdgeLevel));

					levels_total++;
					if(specifiedEdgeLevel == correctEdgeLevel)
					{
						levels_OK++;
					}

					char sequenceCharacter_from_originalRead = r.sequence.at(cI_in_unaligned_sequence_correctlyAligned);
					if(r_aligned.reverse)
					{	
						sequenceCharacter_from_originalRead = Utilities::reverse_char_nucleotide(sequenceCharacter_from_originalRead);
					}
					
					if(!(sequenceCharacter_from_alignment == sequenceCharacter_from_originalRead))
					{
						std::cerr << "! (sequenceCharacter_from_alignment == sequenceCharacter_from_originalRead)" << "\n" << std::flush;
						std::vector<std::string> sequence_from_alignment;
						std::vector<std::string> sequence_from_originalRead;
						for(unsigned int cI2 = 0; cI2 < r_aligned.sequence_aligned.size(); cI2++)
						{
							int cI2_in_unaligned_sequence = -1;
							std::string cStr = r_aligned.sequence_aligned.substr(cI2, 1);
							if(cStr == "_")
							{
							}
							else
							{
								sequence_from_alignment.push_back(cStr);
								
								cI2_in_unaligned_sequence++;
								int cI2_in_unaligned_sequence_correctlyAligned = cI2_in_unaligned_sequence;
								if(r_aligned.reverse)
								{
									cI2_in_unaligned_sequence_correctlyAligned = (r.sequence.length() - cI2_in_unaligned_sequence - 1);									
								}								
								std::string cStr_from_originalRead = r.sequence.substr(cI2_in_unaligned_sequence_correctlyAligned, 1);
								if(r_aligned.reverse)
								{
									cStr_from_originalRead = Utilities::seq_reverse_complement(cStr_from_originalRead);
								}
								sequence_from_originalRead.push_back(cStr_from_originalRead);
							}								
						}
						
						std::cerr << "Sequence alignment: " << r_aligned.sequence_aligned << "\n";
						std::cerr << "Sequence original read: " << r.sequence << "\n";
						std::cerr << "Alignment reverse: " << r_aligned.reverse << "\n";
						std::cerr << "Extracted characters from alignment: " << Utilities::join(sequence_from_alignment, " ") << "\n";
						std::cerr << "Extracted characters from alignment: " << Utilities::join(sequence_from_originalRead, " ") << "\n";
						std::cerr << "\n\n" << std::flush;
						
					}
					assert(sequenceCharacter_from_alignment == sequenceCharacter_from_originalRead);
				}
			}

			if(levels_total != levels_OK)
			{
				incorrectAlignmentsStream << "Genome " << genomePair << " read " << r.name << ": " << levels_OK << " of " << levels_total << " OK\n";
				incorrectAlignmentsStream << "\t" << "  Raw read:" << "\t" << r.sequence << "\n";
				incorrectAlignmentsStream << "\t" << " Qualities:" << "\t" << r.quality << "\n\n";
				incorrectAlignmentsStream << "\t" << "Al. genome:" << "\t" << r_aligned.graph_aligned << "\n";
				incorrectAlignmentsStream << "\t" << "  Al. read:" << "\t" << r_aligned.sequence_aligned << "\n\n";
				incorrectAlignmentsStream << "\t" << "Al. levels:" << "\t" << Utilities::join(Utilities::ItoStr(r_aligned.graph_aligned_levels), " ") << "\n";
				incorrectAlignmentsStream << "\t" << "True levls:" << "\t" << Utilities::join(correctLevelsInAlignmentOrder, " ") << "\n\n";
				incorrectAlignmentsStream << std::flush;
			}
		};

		std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> > alignments_perThread;
		std::vector< std::vector<int> > alignments_readPairI_perThread;
		
		alignments_perThread.resize(outerThreads);
		alignments_readPairI_perThread.resize(outerThreads);
		
		
		auto alignReadPairs = [&](std::vector<oneReadPair>& readPairs, bool usePairing) -> void
		{
			unsigned int pairI = 0;
			unsigned int pairMax = readPairs.size();
			#pragma omp parallel for schedule(dynamic)
			for(pairI = 0; pairI < pairMax; pairI++)
			{
				int tI = omp_get_thread_num();
				assert(omp_get_num_threads() == outerThreads);
				assert((tI >= 0) && (tI < outerThreads));
				
				assert((pairI >= 0) && (pairI < readPairs.size()));
				oneReadPair rP = readPairs.at(pairI);				
				    
				assert((tI >= 0) && (tI < graphAligners.size()));

				std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_pair = graphAligners.at(tI)->seedAndExtend_local_paired(rP, usePairing, insertSize_mean, insertSize_sd);

				alignments_perThread.at(tI).push_back(alignment_pair);
				alignments_readPairI_perThread.at(tI).push_back(pairI);
				
				if(tI == 2)
				{
					std::cout  << Utilities::timestamp() << "\t\t" << "Thread " << tI << ": align pair " << pairI << "\n" << std::flush;
				}
			}
		};
		
		
		std::cout << Utilities::timestamp() << "\t\t" << "Start alignment.\n" << std::flush;
		
		alignReadPairs(combinedPairs_for_alignment_filtered, false);		
		
		std::cout  << Utilities::timestamp() << "\t\t" << "All pairs aligned - merge.\n" << std::flush;
		std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > alignments;
		std::vector< int > alignments_readPairI;
		for(unsigned int tI = 0; tI < outerThreads; tI++)
		{
			alignments.insert(alignments.end(), alignments_perThread.at(tI).begin(), alignments_perThread.at(tI).end());
			alignments_readPairI.insert(alignments_readPairI.end(), alignments_readPairI_perThread.at(tI).begin(), alignments_readPairI_perThread.at(tI).end());		
		}	
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Merging done.\n" << std::flush;
		
		auto evaluateAlignments = [&](std::vector<oneReadPair>& readPairs, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector< int >& alignments_readPairI) -> void {

			
			assert(readPairs.size() == alignments.size());
			assert(alignments.size() == alignments_readPairI.size());
			
			int summary_levels_evaluated = 0;
			int summary_levels_OK = 0;
			
			for(unsigned int alignmentI = 0; alignmentI < alignments.size(); alignmentI++)
			{
				int pairI = alignments_readPairI.at(alignmentI);
				oneReadPair& rP = readPairs.at(pairI);
				oneRead& r1 = rP.reads.first;
				oneRead& r2 = rP.reads.second;
			
				std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignments_thisPair = alignments.at(alignmentI);

				int r1_levels; int r1_levels_OK;
				int r2_levels; int r2_levels_OK;

				checkOneReadLevelCorrectness(r1, alignments_thisPair.first, r1_levels, r1_levels_OK);
				checkOneReadLevelCorrectness(r2, alignments_thisPair.second, r2_levels, r2_levels_OK);

				summary_levels_evaluated += (r1_levels + r2_levels);
				summary_levels_OK += (r1_levels_OK + r2_levels_OK);
			}

			std::cout << "\t\t" << "Summary: " << summary_levels_evaluated << " evaluated, of which " << summary_levels_OK << " were OK. ";
				std::cout << "(" <<  std::setw(5) << (double)summary_levels_OK/(double)summary_levels_evaluated << ")" << "\n\n" << std::flush;
		};

		evaluateAlignments(combinedPairs_for_alignment_filtered, alignments, alignments_readPairI);
	}

	incorrectAlignmentsStream.close();

	for(int tI = 0; tI < outerThreads; tI++)
	{
		delete(graphAligners.at(tI));
		delete(graphs.at(tI));
	}
}

void sampleExactStringFromGraph(Graph* g, int minLength_string, int maxLength_string, std::string& string_ret, std::vector<Edge*>& traversedEdges_ret)
{

	int levels = g->NodesPerLevel.size();
	assert(minLength_string < levels);

	assert(maxLength_string >= minLength_string);
	int lastStartLevel = levels - minLength_string - 1;
	assert(lastStartLevel >= 0);
	assert(lastStartLevel < levels);

	int startLevel = Utilities::randomNumber(lastStartLevel);
	std::set<Node*> startNodes_set = g->NodesPerLevel.at(startLevel);
	std::vector<Node*> startNodes(startNodes_set.begin(), startNodes_set.end());
	int selectedStartNodeIndex = Utilities::randomNumber(startNodes_set.size() - 1);

	Node* startNode = startNodes.at(selectedStartNodeIndex);

	bool error = false;
	auto sampleCharactersFromNode = [&](Node* n, int characters, int& have_nonGap_characters, std::vector<Edge*>& traversedEdges, std::string& edgeLabels) {
		traversedEdges.clear();
		edgeLabels.clear();

		assert(n != 0);
		Node* currentNode = n;
		have_nonGap_characters = 0;
		bool stop = false;
		while(((have_nonGap_characters != characters)) && (! error) && (! stop))
		{
			assert(currentNode != 0);

			if(!( currentNode->level <= (levels - 2)))
			{
				stop = true;
				break;
			}

			std::set<Edge*> availableEdges = currentNode->Outgoing_Edges;
			assert(availableEdges.size() > 0);
			std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

			int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
			Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);

			assert(selectedEdge != 0);
			traversedEdges.push_back(selectedEdge);
			string emission = g->CODE.deCode(selectedEdge->locus_id, selectedEdge->emission);
			assert(emission.length() == 1);

			edgeLabels.append(emission);
			currentNode = selectedEdge->To;

			if(emission != "_")
			{
				have_nonGap_characters++;
			}
		}
	};

	std::vector<Edge*> total_traversedEdges;
	std::string total_sequenceLabels;
	int have_nonGap_characters;

	int wantCharacters = minLength_string + Utilities::randomNumber(maxLength_string - minLength_string);

	sampleCharactersFromNode(startNode, wantCharacters, have_nonGap_characters, total_traversedEdges, total_sequenceLabels);

	if((have_nonGap_characters < minLength_string) || (total_sequenceLabels.at(0) == '_') || (total_sequenceLabels.at(total_sequenceLabels.size() - 1) == '_'))
	{
		sampleExactStringFromGraph(g, minLength_string, maxLength_string, string_ret, traversedEdges_ret);
	}
	else
	{
		string_ret = total_sequenceLabels;
		traversedEdges_ret = total_traversedEdges;
	}
}



};
};
