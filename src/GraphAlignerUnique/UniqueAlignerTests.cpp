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


void testSeedAndExtend_short()
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
			std::cout << "Create GraphAlignerUnique...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAlignerUnique gA(gS_graph, aligner_kMerSize);
		GraphAlignerUnique gA2(gS_graph, aligner_kMerSize);
//		GraphAlignerUnique gA3(gS_graph, aligner_kMerSize);

		gA.setThreads(1);
		gA2.setThreads(1);
		gA2.setIterationsMainRandomizationLoop(0);
//		gA3.setThreads(1);
//		gA3.setIterationsMainRandomizationLoop(2);


		GraphAligner_affine gA_classical(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "[SHORT] Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

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

					15,
					20,

					2,
					4,

					2,
					3,

					3,
					7,

					1,
					4,

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
			seedAndExtend_return_local wholeString_alignments_short = gA.seedAndExtend_short(randomString_noGaps, allBacktraces);
//
//			seedAndExtend_return_local wholeString_alignments_local = gA3.seedAndExtend_local(randomString_noGaps, allBacktraces);

//			seedAndExtend_return wholeString_alignments = gA.seedAndExtend(randomString_noGaps);
			seedAndExtend_return wholeString_alignments_2 = gA2.seedAndExtend(randomString_noGaps);

//			 std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "\n" << std::flush;

			{
				seedAndExtend_return_local& thisAlignment = wholeString_alignments_short;

				int thisAlignmentScore = gA.score(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);

//				int thisAlignmentScore_2 = gA2.score(wholeString_alignments_2.graph_aligned, wholeString_alignments_2.graph_aligned_levels, wholeString_alignments_2.sequence_aligned);

				int thisAlignmentMatches = gA_classical.countMatchesInSequence(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);


				if(true || verbose)
				{
					std::cout << "\tSHORT Alignment [internal score " << thisAlignment.Score << "]>\n";
					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
					std::cout << "\t\tLocal NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;

//					std::cout << "\t\tCompare with local alignment [score " << wholeString_alignments_local.Score << "]>:\n";
//					std::cout << "\t\t\t" << wholeString_alignments_local.graph_aligned << "\n";
//					std::cout << "\t\t\t" << wholeString_alignments_local.sequence_aligned << "\n\n" << std::flush;

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

	std::cout << "testSeedAndExtend_short(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
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
	int simulateGenomePairs = 1;
	int outerThreads = 4;
	int skipPairs_MOD = 50;
	bool evaluateWithoutPairing = false;
	bool useShort = true;
	
	// todo remove
	// boost::mt19937 rnd_gen;
		// auto seed = boost::random::random_device()();
		// rnd_gen.seed(seed);	
	// boost::random::uniform_int_distribution<> nucleotide_gen (0,3);
	// char nucleotides[4] = {'A', 'C', 'G', 'T'};
	// for(unsigned int i = 0; i < 100; i++)
	// {
		// std::cout << Utilities::randomNucleotide(rnd_gen) << " " << std::flush;
	// }
	// assert(1 == 0);
			
	
	// {
		// std::string read = "TTTGTTGACCTTTATTATGACATTCACCAGAAGTTGAAATTGTGTGTTTCTGGTTAATTTTTAATTTATATTTTTTATTTGTAATTCCTTTGAATTATTT";
		// read = "CGCCGTGGGCTACGTGGACGACACAGAGTTCGTGCGGTTCGACAGCGACTCCGTGAGTCCGAGGATGGAGCGGCGGGCGCCGTGGGTGGAGCAGGAGGGG";
		
		// Graph* g = new Graph();
		// g->readFromFile(graph_filename);
			
			
		// GraphAlignerUnique* gA = new GraphAlignerUnique(g, aligner_kMerSize);
		// gA->setIterationsMainRandomizationLoop(0);
		// gA->setThreads(1);	

		// std::vector<seedAndExtend_return_local> allBacktraces;
		// seedAndExtend_return_local wholeString_alignments_local = gA->seedAndExtend_local(read, allBacktraces);

		// int thisAlignmentScore = gA->score(wholeString_alignments_local.graph_aligned, wholeString_alignments_local.graph_aligned_levels, wholeString_alignments_local.sequence_aligned);
		// int thisAlignmentMatches = gA->countMatchesInSequence(wholeString_alignments_local.graph_aligned, wholeString_alignments_local.graph_aligned_levels, wholeString_alignments_local.sequence_aligned);

		// std::cout << "\tLOCAL Alignment [internal score " << wholeString_alignments_local.Score << "]>\n";
		// std::cout << "\t\t" << wholeString_alignments_local.graph_aligned << "\n";
		// std::cout << "\t\t" << wholeString_alignments_local.sequence_aligned << "\n";
		// std::cout << "\t\tLocal NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;
		
		// std::vector<kMerEdgeChain*> chains_for_sequence = gA->getGI().findChains(read);
		// std::cout << "Found chains: " << chains_for_sequence.size() << "\n" << std::flush;
		// for(unsigned int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
		// {
			// kMerEdgeChain* chain = chains_for_sequence.at(chainI);
			// Node* chain_firstNode = chain->traversedEdges.front()->From;
			// Node* chain_lastNode = chain->traversedEdges.back()->To;

			// std::cout << "chainI: " << chainI << "\n";
			// std::cout << "\t" << "chain_firstNode: " << chain_firstNode << " level: " << chain_firstNode->level << "\n";
			// std::cout << "\t" << "chain_lastNode: " << chain_lastNode << " level: " << chain_lastNode->level << "\n";
			// std::cout << "\t" << "sequence_begin: "  << chain->sequence_begin << "\n";
			// std::cout << "\t" << "sequence_end: "  << chain->sequence_end << "\n";
			
				// std::string impliedSequence;
				// for(unsigned int eI = 0; eI < chain->traversedEdges.size(); eI++)
				// {
					// Edge* e = chain->traversedEdges.at(eI);
					// std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					// assert(edgeEmission.length() == 1);
					// if(edgeEmission != "_")
					// {
						// impliedSequence.append(edgeEmission);
					// }
				// }
				
				// std::cout << "\t\t impliedSequence: " << impliedSequence << "\n" << std::flush;			
		// }
	// }	
	// assert(1 == 0);
	
	
	// {
		// std::string r1_seq = "TTCCCCACCCCAGGCGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCATGACAGATGCAAAATGCCTGAATTTTCTGACTCTTC";
		// std::string r2_seq = "CACCACCACAGCCGCCCACTTCTGGAAGGTTCCATCCCCTGCAGGCCTGGTCTCCACGAGCTCCGTGTGTGGGTCTGGTCCTCCCCATCCCGCTGCCAGG";
		// std::string r1_qualities = "`ccffdffhhhhgijjiijigj`ijihPjjjahjjjjf[jjgi\\jjjijejjjiihjbi^ijigdjjhfafhhfccfhg\\cddXjfdfSf^ddddddd^d";
		// std::string r2_qualities = "dedRddddfdbfdhdchddSecefhX^hbijjeji]jjdhjjigjjgjbX]djSjjfjiihijjhj^ijejiijjijajjjjieijiIhhgffffff``b";
		// int difference_starting_coordinates = 206;
		
		// oneRead r1("r1", r1_seq, r1_qualities);
		// oneRead r2("r2", r2_seq, r2_qualities);
		// oneReadPair rP(r1, r2, difference_starting_coordinates);
		
		// Graph* g = new Graph();
		// g->readFromFile(graph_filename);
			
		// GraphAlignerUnique* gA = new GraphAlignerUnique(g, aligner_kMerSize);
		// gA->setIterationsMainRandomizationLoop(4);
		// gA->setThreads(1);	

		// auto alignReadPair = [&](oneReadPair& rP, std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignment, bool usePairing) -> void
		// {		    
			// std::cout << "Align with usePairing = " << usePairing << "\n";
			// alignment = gA->seedAndExtend_local_paired(rP, usePairing, insertSize_mean, insertSize_sd);
		// };
		
		// std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_pairs;
		// std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_noPairs;
		// alignReadPair(rP, alignment_pairs, true);
		// alignReadPair(rP, alignment_noPairs, false);
		
		// auto printAlignment = [&](seedAndExtend_return_local& alignment) -> void {
			// std::cout <<  "\t\t" << "Al. genome:" << "\t" << alignment.graph_aligned << "\n";
			// std::cout << "\t\t" << "  Al. read:" << "\t" << alignment.sequence_aligned << "\n\n";
			// std::cout << "\t\t" << "Al. levels:" << "\t" << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), " ") << "\n";
		// };
		
		// std::cout << "Paired alignment:" << "\n\n";
		// std::cout << "\t" << "Read 1: " << "\n\n";
		// printAlignment(alignment_pairs.first);
		// std::cout << "\t" << "Read 2: " << "\n\n";
		// printAlignment(alignment_pairs.second);
		// std::cout << "\n";
		
		// std::cout << "Unpaired alignment:" << "\n\n";
		// std::cout << "\t" << "Read 1: " << "\n\n";
		// printAlignment(alignment_noPairs.first);
		// std::cout << "\t" << "Read 2: " << "\n\n";
		// printAlignment(alignment_noPairs.second);
		// std::cout << "\n" << std::flush;
				
	// }	
	// assert(1 == 0);	

	std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Loading graph.\n" << std::flush;

	std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Create GraphAlignerUnique.\n" << std::flush;

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
		graphAligners.at(tI)->setIterationsMainRandomizationLoop(4);
		graphAligners.at(tI)->setThreads(1);	
		
		graphs.at(tI) = g;
	}
	
	int assignedLevels_totalReads = 0;
	std::vector<int> levels_assigned_reads;
	std::vector<int> levels_assigned_reads_recovered;
	std::vector<std::string> levelIDs;

	levels_assigned_reads.resize(graphs.at(0)->NodesPerLevel.size(), 0);
	levels_assigned_reads_recovered.resize(graphs.at(0)->NodesPerLevel.size(), 0);
	for(unsigned int levelI = 0; levelI < (graphs.at(0)->NodesPerLevel.size() - 1); levelI++)
	{
		assert(graphs.at(0)->NodesPerLevel.at(levelI).size() > 0);
		Node* n = *(graphs.at(0)->NodesPerLevel.at(levelI).begin());
		
		assert(n->Outgoing_Edges.size() > 0);
		Edge* e = *(n->Outgoing_Edges.begin());
		
		std::string locus = e->locus_id;		
		levelIDs.push_back(locus);
	}
	
	std::cout << Utilities::timestamp() << "\tdone.\n" << std::flush;	

	std::string fileName_for_incorrectAlignments_unpaired = "../tmp/incorrectAlignments_unpaired.txt";
	std::string fileName_for_incorrectAlignments_paired = "../tmp/incorrectAlignments_paired.txt";
	std::ofstream incorrectAlignmentsStream_unpaired;
	incorrectAlignmentsStream_unpaired.open(fileName_for_incorrectAlignments_unpaired.c_str());
	assert(incorrectAlignmentsStream_unpaired.is_open());
	std::ofstream incorrectAlignmentsStream_paired;
	incorrectAlignmentsStream_paired.open(fileName_for_incorrectAlignments_paired.c_str());
	assert(incorrectAlignmentsStream_paired.is_open());
	
	std::string fileName_for_unpaired_better = "../tmp/alignments_unpaired_better.txt";
	std::ofstream unpairedAlignmentsBetterStream;
	unpairedAlignmentsBetterStream.open(fileName_for_unpaired_better.c_str());
	assert(unpairedAlignmentsBetterStream.is_open());

	std::string fileName_for_correctLevel_summary = "../tmp/summary_correct_perLevel.txt";

	for(int genomePair = 1; genomePair <= simulateGenomePairs; genomePair++)
	{
		std::cout << Utilities::timestamp() << "testSeedAndExtend_local_realGraph(..): Iteration " << genomePair << " / " << simulateGenomePairs << "\n" << std::flush;

		diploidEdgePointerPath diploidPath = graphs.at(0)->simulateRandomDiploidPath();

		std::vector<oneReadPair> simulatedReadPairs_h1 = rS.simulate_paired_reads_from_edgePath(diploidPath.h1, haploidCoverage, insertSize_mean, insertSize_sd, false);
		std::vector<oneReadPair> simulatedReadPairs_h2 = rS.simulate_paired_reads_from_edgePath(diploidPath.h2, haploidCoverage, insertSize_mean, insertSize_sd, false);

		std::vector<oneReadPair> combinedPairs_for_alignment;
		combinedPairs_for_alignment.insert(combinedPairs_for_alignment.end(), simulatedReadPairs_h1.begin(), simulatedReadPairs_h1.end());
		combinedPairs_for_alignment.insert(combinedPairs_for_alignment.end(), simulatedReadPairs_h2.begin(), simulatedReadPairs_h2.end());

		std::cout << "\t" << "Simulated " << combinedPairs_for_alignment.size() << " read pairs." << "\n" << std::flush;

		std::vector<oneReadPair> combinedPairs_for_alignment_filtered;
		for(unsigned int pI = 0; pI < combinedPairs_for_alignment.size(); pI++)
		{
			if((pI % skipPairs_MOD) != 0)
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

		auto checkOneReadLevelCorrectness = [&](oneRead& r, seedAndExtend_return_local& r_aligned, int& levels_total, int& levels_OK, int& positions_total, int& positions_matches, bool printToFile, std::ofstream& incorrectAlignmentsStream) -> void {
		
			assignedLevels_totalReads++;

			levels_total = 0;
			levels_OK = 0;

			positions_total = 0;
			positions_matches = 0;
			
			assert(r.coordinates_edgePath.size() == r.sequence.length());
			int cI_in_unaligned_sequence = -1;
			std::vector<std::string> correctLevelsInAlignmentOrder;

			for(unsigned int cI = 0; cI < r_aligned.sequence_aligned.size(); cI++)
			{
				char sequenceCharacter_from_alignment = r_aligned.sequence_aligned.at(cI);
				char graphCharacter_from_alignment = r_aligned.graph_aligned.at(cI);
				
				positions_total++;				
				if(sequenceCharacter_from_alignment == graphCharacter_from_alignment)
				{
					positions_matches++;
				}
				
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
					
					if(correctEdgeLevel != -1)
					{
						levels_assigned_reads.at(correctEdgeLevel)++;
						if(specifiedEdgeLevel == correctEdgeLevel)
						{
							levels_assigned_reads_recovered.at(correctEdgeLevel)++;
						}						
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

			if(printToFile && (levels_total != levels_OK) && (positions_matches != positions_total))
			{
				incorrectAlignmentsStream << "Genome " << genomePair << " read " << r.name << ": " << levels_OK << " of " << levels_total << " levels OK; " << positions_matches << " of " << positions_total << " positions OK.\n";
				incorrectAlignmentsStream << "\t" << "  Raw read:" << "\t" << r.sequence << "\n";
				incorrectAlignmentsStream << "\t" << " Qualities:" << "\t" << r.quality << "\n\n";
				incorrectAlignmentsStream << "\t" << "Al. genome:" << "\t" << r_aligned.graph_aligned << "\n";
				incorrectAlignmentsStream << "\t" << "  Al. read:" << "\t" << r_aligned.sequence_aligned << "\n\n";
				incorrectAlignmentsStream << "\t" << "Al. levels:" << "\t" << Utilities::join(Utilities::ItoStr(r_aligned.graph_aligned_levels), " ") << "\n";
				incorrectAlignmentsStream << "\t" << "True levls:" << "\t" << Utilities::join(correctLevelsInAlignmentOrder, " ") << "\n\n";
				incorrectAlignmentsStream << std::flush;
			}
		};
		
		auto alignReadPairs = [&](std::vector<oneReadPair>& readPairs, std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> >& alignments_perThread, std::vector< std::vector<int> >& alignments_readPairI_perThread, bool usePairing) -> void
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

				std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_pair = graphAligners.at(tI)->seedAndExtend_local_paired_or_short(rP, usePairing, useShort, insertSize_mean, insertSize_sd);

				alignments_perThread.at(tI).push_back(alignment_pair);
				alignments_readPairI_perThread.at(tI).push_back(pairI);
				
				if(tI == 0)
				{
					std::cout  << Utilities::timestamp() << "\t\t" << "Thread " << tI << ": align pair " << pairI << "\n" << std::flush;
				}
			}
		};
		
		std::cout << Utilities::timestamp() << "\t\t" << "Start alignment.\n" << std::flush;
		
		// no pairing
		
		std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> > noPairing_alignments_perThread;
		std::vector< std::vector<int> > noPairing_alignments_readPairI_perThread; 
		
		noPairing_alignments_perThread.resize(outerThreads);
		noPairing_alignments_readPairI_perThread.resize(outerThreads);
				
		if(evaluateWithoutPairing)
		{
			alignReadPairs(combinedPairs_for_alignment_filtered, noPairing_alignments_perThread, noPairing_alignments_readPairI_perThread, false);		
		}
		
		// with pairs
		
		std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> > withPairing_alignments_perThread;
		std::vector< std::vector<int> > withPairing_alignments_readPairI_perThread; 
		
		withPairing_alignments_perThread.resize(outerThreads);
		withPairing_alignments_readPairI_perThread.resize(outerThreads);
				
		alignReadPairs(combinedPairs_for_alignment_filtered, withPairing_alignments_perThread, withPairing_alignments_readPairI_perThread, true);		

		// merge
		
		std::cout  << Utilities::timestamp() << "\t\t" << "All pairs aligned - merge.\n" << std::flush;
		
		std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > noPairing_alignments;
		std::vector< int > noPairing_alignments_readPairI;
		for(unsigned int tI = 0; tI < outerThreads; tI++)
		{
			noPairing_alignments.insert(noPairing_alignments.end(), noPairing_alignments_perThread.at(tI).begin(), noPairing_alignments_perThread.at(tI).end());
			noPairing_alignments_readPairI.insert(noPairing_alignments_readPairI.end(), noPairing_alignments_readPairI_perThread.at(tI).begin(), noPairing_alignments_readPairI_perThread.at(tI).end());		
		}	
		

		std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > withPairing_alignments;
		std::vector< int > withPairing_alignments_readPairI;
		for(unsigned int tI = 0; tI < outerThreads; tI++)
		{
			withPairing_alignments.insert(withPairing_alignments.end(), withPairing_alignments_perThread.at(tI).begin(), withPairing_alignments_perThread.at(tI).end());
			withPairing_alignments_readPairI.insert(withPairing_alignments_readPairI.end(), withPairing_alignments_readPairI_perThread.at(tI).begin(), withPairing_alignments_readPairI_perThread.at(tI).end());		
		}	
		
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Merging done.\n" << std::flush;
		
		auto evaluateAlignments = [&](std::vector<oneReadPair>& readPairs, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector< int >& alignments_readPairI, bool printToFile, std::ofstream& incorrectAlignmentsStream, std::map<int, int>&  levelsOK_perPair) -> void {

			
			assert(readPairs.size() == alignments.size());
			assert(alignments.size() == alignments_readPairI.size());
			
			int summary_levels_evaluated = 0;
			int summary_levels_OK = 0;
			int positions_total = 0;
			int positions_matches = 0;			
			
			for(unsigned int alignmentI = 0; alignmentI < alignments.size(); alignmentI++)
			{
				int pairI = alignments_readPairI.at(alignmentI);
				oneReadPair& rP = readPairs.at(pairI);
				oneRead& r1 = rP.reads.first;
				oneRead& r2 = rP.reads.second;
			
				std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignments_thisPair = alignments.at(alignmentI);

				int r1_levels; int r1_levels_OK;
				int r2_levels; int r2_levels_OK;
				int r1_positions_total; int r1_positions_OK;
				int r2_positions_total; int r2_positions_OK;
				
				checkOneReadLevelCorrectness(r1, alignments_thisPair.first, r1_levels, r1_levels_OK, r1_positions_total, r1_positions_OK, printToFile, incorrectAlignmentsStream);
				checkOneReadLevelCorrectness(r2, alignments_thisPair.second, r2_levels, r2_levels_OK, r2_positions_total, r2_positions_OK, printToFile, incorrectAlignmentsStream);

				summary_levels_evaluated += (r1_levels + r2_levels);
				summary_levels_OK += (r1_levels_OK + r2_levels_OK);
				
				positions_total += (r1_positions_total + r2_positions_total);
				positions_matches += (r1_positions_OK + r2_positions_OK);
				
				assert(levelsOK_perPair.count(pairI) == 0);
				levelsOK_perPair[pairI] = (r1_levels_OK + r2_levels_OK);			
				
			}

			std::cout << "\t\t" << "Matches: " << positions_total << " evaluated, of which " << positions_matches << " were OK. ";
				std::cout << "(" <<  std::setw(5) << (double)positions_matches/(double)positions_total << ")" << "\n\n" << std::flush;				
			
			std::cout << "\t\t" << "Levels: " << summary_levels_evaluated << " evaluated, of which " << summary_levels_OK << " were OK. ";
				std::cout << "(" <<  std::setw(5) << (double)summary_levels_OK/(double)summary_levels_evaluated << ")" << "\n\n" << std::flush;
		};

		std::map<int, int> noPairing_levelsOK_perPair;
		std::map<int, int> withPairing_levelsOK_perPair;
		
		if(evaluateWithoutPairing)
		{
			std::cout << "Evaluation WITHOUT pairs:\n\n" << std::flush;
			evaluateAlignments(combinedPairs_for_alignment_filtered, noPairing_alignments, noPairing_alignments_readPairI, true, incorrectAlignmentsStream_unpaired, noPairing_levelsOK_perPair);
		}
		
		std::cout << "Evaluation WITH pairs:\n\n" << std::flush;
		evaluateAlignments(combinedPairs_for_alignment_filtered, withPairing_alignments, withPairing_alignments_readPairI, true, incorrectAlignmentsStream_paired, withPairing_levelsOK_perPair);
		
		// find alignments for which UNPAIRED is better (unexpectedly)
	
		if(evaluateWithoutPairing)
		{
			std::map<int, std::pair<seedAndExtend_return_local, seedAndExtend_return_local>*> noPairing_alignments_perPair;
			std::map<int, std::pair<seedAndExtend_return_local, seedAndExtend_return_local>*> withPairing_alignments_perPair;
			for(unsigned int i = 0; i < noPairing_alignments.size(); i++)
			{
				int pairI = noPairing_alignments_readPairI.at(i);
				assert(noPairing_alignments_perPair.count(pairI) == 0);
				noPairing_alignments_perPair[pairI] = &(noPairing_alignments.at(i));
			}
			for(unsigned int i = 0; i < withPairing_alignments.size(); i++)
			{
				int pairI = withPairing_alignments_readPairI.at(i);
				assert(withPairing_alignments_perPair.count(pairI) == 0);
				withPairing_alignments_perPair[pairI] = &(withPairing_alignments.at(i));
			}
			
			
			assert(noPairing_levelsOK_perPair.size() == withPairing_levelsOK_perPair.size());
			assert(noPairing_levelsOK_perPair.size() == combinedPairs_for_alignment_filtered.size());
			int paired_better = 0;
			int unpaired_better = 0;
			int paired_unpaired_equal = 0;
			
			for(std::map<int, int>::iterator pIt = noPairing_levelsOK_perPair.begin(); pIt != noPairing_levelsOK_perPair.end(); pIt++)
			{
				int pairI = pIt->first;
				int noPairing_OK = noPairing_levelsOK_perPair.at(pairI);
				int withPairing_OK = withPairing_levelsOK_perPair.at(pairI);
				
				oneReadPair& rP = combinedPairs_for_alignment_filtered.at(pairI);
				oneRead& r1 = rP.reads.first;
				oneRead& r2 = rP.reads.second;
			
				std::pair<seedAndExtend_return_local, seedAndExtend_return_local>* alignment_paired = withPairing_alignments_perPair.at(pairI);
				std::pair<seedAndExtend_return_local, seedAndExtend_return_local>* alignment_unpaired = noPairing_alignments_perPair.at(pairI);
							
				if(noPairing_OK > withPairing_OK)
				{
					unpaired_better++;
					
					auto printAlignmentVersusTruth = [&](oneRead& r, seedAndExtend_return_local& r_aligned, std::ofstream& outputStream, std::string indent) -> void {
					
						int cI_in_unaligned_sequence = -1;
						std::vector<std::string> correctLevelsInAlignmentOrder;

						int levels_OK = 0;
						int levels_total = 0;
						
						for(unsigned int cI = 0; cI < r_aligned.sequence_aligned.size(); cI++)
						{
							char sequenceCharacter_from_alignment = r_aligned.sequence_aligned.at(cI);
							char graphCharacter_from_alignment = r_aligned.graph_aligned.at(cI);
							
							if(sequenceCharacter_from_alignment == '_')
							{
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
								assert(sequenceCharacter_from_alignment == sequenceCharacter_from_originalRead);
							}
						}

						outputStream << indent << "Read " << r.name << ": " << levels_OK << " of " << levels_total << " levels OK\n";
						outputStream << indent << "\t" << "  Raw read:" << "\t" << r.sequence << "\n";
						outputStream << indent << "\t" << " Qualities:" << "\t" << r.quality << "\n\n";
						outputStream << indent << "\t" << "Al. genome:" << "\t" << r_aligned.graph_aligned << "\n";
						outputStream << indent << "\t" << "  Al. read:" << "\t" << r_aligned.sequence_aligned << "\n\n";
						outputStream << indent << "\t" << "Al. levels:" << "\t" << Utilities::join(Utilities::ItoStr(r_aligned.graph_aligned_levels), " ") << "\n";
						outputStream << indent << "\t" << "True levls:" << "\t" << Utilities::join(correctLevelsInAlignmentOrder, " ") << "\n\n";
						outputStream << std::flush;				
					};
					
					unpairedAlignmentsBetterStream << "Genome " << genomePair << ", read pair " << pairI << ": " << noPairing_OK << " no-pairs vs " << withPairing_OK << " with-pairs.\n";
						unpairedAlignmentsBetterStream << "\t" << "diff_starting_coordinates: " << rP.diff_starting_coordinates << "\n";
						unpairedAlignmentsBetterStream << "\t" << "Unpaired: " << "\n";
							unpairedAlignmentsBetterStream << "\t\t" << "Read 1: " << "\n";
								printAlignmentVersusTruth(r1, alignment_unpaired->first, unpairedAlignmentsBetterStream, "\t\t\t");
							unpairedAlignmentsBetterStream << "\t\t" << "Read 2: " << "\n";				
								printAlignmentVersusTruth(r2, alignment_unpaired->second, unpairedAlignmentsBetterStream, "\t\t\t");
								
						unpairedAlignmentsBetterStream << "\t" << "Paired: " << "\n";
							unpairedAlignmentsBetterStream << "\t\t" << "Read 1: " << "\n";
								printAlignmentVersusTruth(r1, alignment_paired->first, unpairedAlignmentsBetterStream, "\t\t\t");
							unpairedAlignmentsBetterStream << "\t\t" << "Read 2: " << "\n";				
								printAlignmentVersusTruth(r2, alignment_paired->second, unpairedAlignmentsBetterStream, "\t\t\t");
									
					unpairedAlignmentsBetterStream << "=====================================================\n\n\n" << std::flush;
				}
				else
				{
					if(noPairing_OK == withPairing_OK)
					{
						paired_unpaired_equal++;
					}
					else
					{
						assert(withPairing_OK > noPairing_OK);
						paired_better++;
					}
				}
			}
		
		
			std::cout << "\nPairing statistics:\n";
			std::cout << "\t" << "With pairing better : " << paired_better << "\n";
			std::cout << "\t" << "With/without equal  : " << paired_unpaired_equal << "\n";
			std::cout << "\t" << "W/out pairing better: " << unpaired_better << "\n";
			std::cout << "\t" << "Total               : " << (paired_better + paired_unpaired_equal + unpaired_better) << "\n" << std::flush;
		}
		
		std::ofstream correctLevel_summary_stream;
		correctLevel_summary_stream.open(fileName_for_correctLevel_summary.c_str());
		assert(correctLevel_summary_stream.is_open());
		correctLevel_summary_stream << "Genome iteration: " << genomePair << "\n" << std::flush;
		correctLevel_summary_stream << "Total reads evaluated: " << assignedLevels_totalReads << "\n" << std::flush;
		correctLevel_summary_stream << Utilities::join(levelIDs, " ") << "\n" << std::flush;
		correctLevel_summary_stream << Utilities::join(Utilities::ItoStr(levels_assigned_reads), " ") << "\n" << std::flush;
		correctLevel_summary_stream << Utilities::join(Utilities::ItoStr(levels_assigned_reads_recovered), " ") << "\n" << std::flush;
		
		correctLevel_summary_stream.close();
	}

	incorrectAlignmentsStream_unpaired.close();
	incorrectAlignmentsStream_paired.close();
	unpairedAlignmentsBetterStream.close();
	
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
