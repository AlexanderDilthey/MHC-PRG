/*
 * GraphAlignerUnique.cpp
 *
 *  Created on: 30.07.2013
 *      Author: AlexanderDilthey
 */

#include "GraphAlignerUnique.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <map>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <utility>
#include <omp.h>
#include <sstream>

#include <boost/math/distributions/normal.hpp>

#include "../Utilities.h"
#include "../hash/sequence/basic.h"
#include "VirtualNWUnique.h"

namespace GraphAlignerUnique {
GraphAlignerUnique::GraphAlignerUnique(Graph* graph, int k) : g(graph), kMerSize(k), gI(g, k){
	S_match = 2;
	S_mismatch = -5;
	S_gap = -2;
	S_graphGap = 0;

	S_openGap = -4;
	S_extendGap = -2;

	iterationsMainRandomizationLoop = 20;
	minimumChainUniqueness = 1;
	randomizationParameter = 3;

	verbose = false;

	threads = 4;

	// enumerate states at levels
	unsigned int levels = g->NodesPerLevel.size();
	nodesPerLevel_ordered.resize(levels);
	nodesPerLevel_ordered_rev.resize(levels);
	for(unsigned int levelI = 0; levelI < levels; levelI++)
	{
		nodesPerLevel_ordered.at(levelI) = std::vector<Node*>(g->NodesPerLevel.at(levelI).begin(), g->NodesPerLevel.at(levelI).end());
		for(unsigned int nodeI = 0; nodeI < nodesPerLevel_ordered.at(levelI).size(); nodeI++)
		{
			nodesPerLevel_ordered_rev.at(levelI)[nodesPerLevel_ordered.at(levelI).at(nodeI)] = nodeI;
		}
	}
}


double GraphAlignerUnique::scoreOneAlignment(oneRead& underlyingRead, seedAndExtend_return_local& alignment, int& totalMismatches)
{
	int indexIntoOriginalReadData = -1;

	double rate_deletions = log(0.01);
	double rate_insertions = log(0.01);
	double combined_log_likelihood = 0;

	totalMismatches = 0;
	
	bool verbose = false;
	
	if(verbose) std::cout << "SCORE\n";
	
	for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
	{
		std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
		std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
		int graphLevel = alignment.graph_aligned_levels.at(cI);

		if(sequenceCharacter != "_")
		{
			indexIntoOriginalReadData++;
			int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
			if(alignment.reverse)
			{
				indexIntoOriginalReadData_correctlyAligned = underlyingRead.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
			}
			assert(indexIntoOriginalReadData_correctlyAligned >= 0);
			assert(indexIntoOriginalReadData_correctlyAligned < underlyingRead.sequence.length());;

			std::string underlyingReadCharacter = underlyingRead.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
			if(alignment.reverse)
			{
				underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
			}
			assert(underlyingReadCharacter == sequenceCharacter);

			if(graphCharacter == "_")
			{
				// sequence non gap, graph gap -- insertion
				combined_log_likelihood += (rate_insertions + log(1.0/4.0));
				totalMismatches++;
			}
			else
			{
				// two well-defined characters
				char qualityCharacter = underlyingRead.quality.at(indexIntoOriginalReadData_correctlyAligned);
				double pCorrect = Utilities::PhredToPCorrect(qualityCharacter);
				assert((pCorrect > 0) && (pCorrect <= 1));
				if(sequenceCharacter == graphCharacter)
				{
					combined_log_likelihood += log(pCorrect);
				}
				else
				{
					double pIncorrect = 1 - pCorrect;
					assert((pIncorrect > 0) && (pIncorrect < 1));
					combined_log_likelihood += log(pIncorrect);
					totalMismatches++;
				}
			}

		}
		else
		{
			assert(sequenceCharacter == "_");
			if(graphCharacter == "_")
			{
				// sequence gap, graph gap - likelihood 1
				assert(graphLevel != -1);
			}
			else
			{
				// sequence gap, graph non gap - deletion in sequence
				combined_log_likelihood += rate_deletions;
				totalMismatches++;
			}
		}
		
		if(verbose) std::cout << "\t" << cI << "\t" << combined_log_likelihood << "\n";
	}

	return combined_log_likelihood;
}


void GraphAlignerUnique::printkMerProfile(std::string file)
{
	std::vector<std::string> kMersInGraph = gI.getIndexedkMers();
	unsigned int graphLevels = g->NodesPerLevel.size();

	std::vector<int> inducedkMers_perLevel;
	std::vector<int> inducedkMers_perLevel_unique;

	inducedkMers_perLevel.resize(graphLevels, 0);
	inducedkMers_perLevel_unique.resize(graphLevels, 0);

	for(unsigned int kI = 0; kI < kMersInGraph.size(); kI++)
	{
		std::string kMer = kMersInGraph.at(kI);
		std::vector<kMerInGraphSpec> kMerPositions = gI.queryIndex(kMer);
		for(unsigned int posI = 0; posI < kMerPositions.size(); posI++)
		{
			kMerInGraphSpec thisPosition = kMerPositions.at(posI);
			Edge* firstEdge = thisPosition.traversedEdges.at(0);
			unsigned int level = firstEdge->From->level;

			inducedkMers_perLevel.at(level)++;

			if(kMerPositions.size() == 1)
			{
				inducedkMers_perLevel_unique.at(level)++;
			}
		}
	}


	std::ofstream output;
	output.open (file.c_str(), ios::out);
	if(! output.is_open())
	{
		throw std::runtime_error("Cannot open file "+file+" for writing.\n");
	}
	assert(output.is_open());

	output 	<< "Level" << "\t"
			<< "IndexedkMers" << "\t"
			<< "IndexedkMers_Unique" << "\n";

	for(unsigned int level = 0; level < graphLevels; level++)
	{
		output 	<< level << "\t"
				<< inducedkMers_perLevel.at(level) << "\t"
				<< inducedkMers_perLevel_unique.at(level) << "\n";
	}

	output.close();
}

void GraphAlignerUnique::analyzeChainUniqueness(std::string& sequence, std::vector<kMerEdgeChain*>& uniquelyTrimmedChains, std::map<std::string, int>& kMer_sequence_occurrences, std::map<kMerEdgeChain*, int>& uniquelyTrimmedChains_doubleUniquekMers, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered)
{
	for(size_t chainI = 0; chainI < uniquelyTrimmedChains.size(); chainI++)
	{
		kMerEdgeChain* chain = uniquelyTrimmedChains.at(chainI);
		std::string impliedSequence(sequence.begin() + chain->sequence_begin, sequence.begin() + chain->sequence_end + 1);
		std::vector<std::string> kMers_impliedSequence = partitionStringIntokMers(impliedSequence, kMerSize);

		size_t kMers_doubleUnique = 0;
		for(size_t kMerI = 0; kMerI < kMers_impliedSequence.size(); kMerI++)
		{
			std::string kMer = kMers_impliedSequence.at(kMerI);
			if(iskMerDoubleUnique(kMer, kMer_sequence_occurrences))
			{
				kMers_doubleUnique++;
			}
		}

		assert(kMers_impliedSequence.size() > 0);

		if(verbose)
		{
			std::cout << "\t\t" << "analyzeChainUniqueness(..): " << chainI << "[" << chain << "]\t" << kMers_doubleUnique << "\t" << ((double)kMers_doubleUnique/(double)kMers_impliedSequence.size()) << "\n" << std::flush;
		}
		
		uniquelyTrimmedChains_doubleUniquekMers[chain] = kMers_doubleUnique;
		uniquelyTrimmedChains_ordered.insert(chain);
	}

	if(uniquelyTrimmedChains_ordered.size() > 0)
	{
		kMerEdgeChain* set_currentChain = *uniquelyTrimmedChains_ordered.begin();
		std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>::iterator set_nextChain_iterator = uniquelyTrimmedChains_ordered.begin();
		set_nextChain_iterator++;
		while(set_nextChain_iterator != uniquelyTrimmedChains_ordered.end())
		{
			kMerEdgeChain* nextChain = *set_nextChain_iterator;
			assert(uniquelyTrimmedChains_doubleUniquekMers.at(set_currentChain) >= uniquelyTrimmedChains_doubleUniquekMers.at(nextChain));
			set_currentChain = *set_nextChain_iterator;
			set_nextChain_iterator++;
		}
	}
}

void GraphAlignerUnique::buildkMerOccurrenceMap(std::string& S, std::map<std::string, int>& occurrencesInSequence)
{
	std::vector<std::string> kMers = partitionStringIntokMers(S, kMerSize);
	for(unsigned int kI = 0; kI < kMers.size(); kI++)
	{
		std::string kMer = kMers.at(kI);
		if(occurrencesInSequence.count(kMer) == 0)
		{
			occurrencesInSequence[kMer] = 0;
		}
		occurrencesInSequence[kMer]++;
	}
};

void  GraphAlignerUnique::seedAndExtend_init_occurrence_strand_etc(std::string& sequence_nonReverse, bool& useReverse, std::string& sequence, std::vector<std::string>& kMers_sequence, std::map<std::string, int>& kMer_sequence_occurrences)
{
	// build occurrence hashes

	std::string sequence_reverse = seq_reverse_complement(sequence_nonReverse);

	std::vector<std::string> kMers_sequence_nonReverse = partitionStringIntokMers(sequence_nonReverse, kMerSize);
	std::vector<std::string> kMers_sequence_reverse = partitionStringIntokMers(sequence_reverse, kMerSize);
	assert(kMers_sequence_nonReverse.size() == kMers_sequence_reverse.size());

	size_t kMers_sequence_inGraph = 0;
	size_t kMers_sequence_reverse_inGraph = 0;

	for(size_t kI = 0; kI < kMers_sequence_nonReverse.size(); kI++)
	{
		if(gI.queryIndex(kMers_sequence_nonReverse.at(kI)).size() > 0)
		{
			kMers_sequence_inGraph++;
		}
		if(gI.queryIndex(kMers_sequence_reverse.at(kI)).size() > 0)
		{
			kMers_sequence_reverse_inGraph++;
		}
	}

	useReverse = (kMers_sequence_reverse_inGraph > kMers_sequence_inGraph);

	if(verbose)
	{
		std::cout << Utilities::timestamp() << "Make decision on orientation.\n" << std::flush;
		std::cout << "\t\t\t" << "Total kMers in sequence: " << kMers_sequence_nonReverse.size() << "\n";
		std::cout << "\t\t\t" << "kMers_sequence_inGraph: " << kMers_sequence_inGraph << "\n";
		std::cout << "\t\t\t" << "kMers_sequence_reverse_inGraph: " << kMers_sequence_reverse_inGraph << "\n";
		std::cout << "\t\t\t" << "useReverse: " << useReverse << "\n" << std::flush;
	}
	
	sequence = (useReverse) ? sequence_reverse : sequence_nonReverse;
	kMers_sequence = (useReverse) ? kMers_sequence_reverse : kMers_sequence_nonReverse;
	buildkMerOccurrenceMap(sequence, kMer_sequence_occurrences);
}

bool GraphAlignerUnique::iskMerDoubleUnique(std::string& kMer, std::map<std::string, int>& occurrencesInSequence)
{
	if(occurrencesInSequence.count(kMer) && (occurrencesInSequence.at(kMer) == 1))
	{
		if(gI.queryIndex(kMer).size() == 1)
		{
			return true;
		}
	}
	return false;
}

void GraphAlignerUnique::trimChain(kMerEdgeChain* inputChain, std::string& sequence, int removeLeft, int removeRight)
{
	int existingLength = inputChain->sequence_end - inputChain->sequence_begin + 1;
	assert(existingLength > 0);
	assert(existingLength > (removeLeft + removeRight));

	auto getEdgeEmission = [&](Edge* e) -> std::string {
			return g->CODE.deCode(e->locus_id, e->emission);
		};

	int sequenceRemove_begin_characters = removeLeft;
	int begin_removed_characters = 0;
	while(begin_removed_characters < sequenceRemove_begin_characters)
	{
		Edge* firstEdge = inputChain->traversedEdges.at(0);
		while(getEdgeEmission(firstEdge) == "_")
		{
			assert(begin_removed_characters == 0);
			inputChain->traversedEdges.erase(inputChain->traversedEdges.begin());
			firstEdge = inputChain->traversedEdges.at(0);
		}

		std::string firstCharacter = getEdgeEmission(firstEdge);
		assert(firstCharacter == sequence.substr(inputChain->sequence_begin, 1));

		inputChain->traversedEdges.erase(inputChain->traversedEdges.begin());
		inputChain->sequence_begin++;

		Edge* afterEdge = inputChain->traversedEdges.at(0);
		while(getEdgeEmission(afterEdge) == "_")
		{
			inputChain->traversedEdges.erase(inputChain->traversedEdges.begin());
			afterEdge = inputChain->traversedEdges.at(0);
		}

		begin_removed_characters++;
	}

	int sequenceRemove_end_characters = removeRight;
	int end_removed_characters = 0;
	while(end_removed_characters < sequenceRemove_end_characters)
	{
		Edge* lastEdge = inputChain->traversedEdges.at(inputChain->traversedEdges.size() - 1);
		while(getEdgeEmission(lastEdge) == "_")
		{
			assert(end_removed_characters == 0);
			inputChain->traversedEdges.erase(inputChain->traversedEdges.end() - 1);
			lastEdge = inputChain->traversedEdges.at(inputChain->traversedEdges.size() - 1);
		}

		std::string lastCharacter = getEdgeEmission(lastEdge);
		assert(lastCharacter == sequence.substr(inputChain->sequence_end, 1));

		inputChain->traversedEdges.erase(inputChain->traversedEdges.end() - 1);
		inputChain->sequence_end--;

		Edge* beforeEdge = inputChain->traversedEdges.at(inputChain->traversedEdges.size() - 1);
		while(getEdgeEmission(beforeEdge) == "_")
		{
			inputChain->traversedEdges.erase(inputChain->traversedEdges.end() - 1);
			beforeEdge = inputChain->traversedEdges.at(inputChain->traversedEdges.size() - 1);
		}

		end_removed_characters++;
	}

	assert(inputChain->sequence_begin <= inputChain->sequence_end);

	std::string impliedSequence;
	for(unsigned int eI = 0; eI < inputChain->traversedEdges.size(); eI++)
	{
		Edge* e = inputChain->traversedEdges.at(eI);
		std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
		assert(edgeEmission.length() == 1);
		if(edgeEmission != "_")
		{
			impliedSequence.append(edgeEmission);
		}
	}
	std::string extractedSequence(sequence.begin() + inputChain->sequence_begin, sequence.begin() + inputChain->sequence_end + 1);
	if(!(extractedSequence == impliedSequence))
	{
		std::string extractedSequence_originalChain(sequence.begin() + inputChain->sequence_begin, sequence.begin() + inputChain->sequence_end + 1);

		std::cerr << "!(extractedSequence == impliedSequence)!\n" << std::flush;
		std::cerr << "extractedSequence: " << extractedSequence << "\n" << std::flush;
		std::cerr << "impliedSequence: " << impliedSequence << "\n" << std::flush;
		std::cerr << "extractedSequence_originalChain: " << extractedSequence_originalChain << "\n" << std::flush;

	}
	assert(extractedSequence == impliedSequence);

	for(unsigned int eI = 0; eI < inputChain->traversedEdges.size(); eI++)
	{
		Edge* e = inputChain->traversedEdges.at(eI);
		Node* e_To = e->To;
		Node* e_From = e->From;
		if(eI < (inputChain->traversedEdges.size() - 1))
		{
			std::set<Edge*> nextEdges = e_To->Outgoing_Edges;
			Edge* e_next = inputChain->traversedEdges.at(eI+1);
			assert(nextEdges.count(e_next));
		}
		if(eI > 0)
		{
			std::set<Edge*> previousEdges = e_From->Incoming_Edges;
			Edge* e_previous = inputChain->traversedEdges.at(eI-1);
			assert(previousEdges.count(e_previous));
		}
	}
}


std::vector<kMerEdgeChain*> GraphAlignerUnique::trimChainsForUniqueness(std::vector<kMerEdgeChain*>& inputChains, std::string& sequence, std::map<std::string, int>& occurrencesInSequence)
{
	std::vector<kMerEdgeChain*> uniquelyTrimmedChains;

	auto trimChainForUniquess = [&](kMerEdgeChain* inputChain) -> kMerEdgeChain* {
		std::string impliedSequence(sequence.begin() + inputChain->sequence_begin, sequence.begin() + inputChain->sequence_end + 1);
		std::vector<std::string> kMers_impliedSequence = partitionStringIntokMers(impliedSequence, kMerSize);

		int remove_kMers_beginning = 0;
		int remove_kMers_end = 0;

		for(size_t kMerI = 0; kMerI < kMers_impliedSequence.size(); kMerI++)
		{
			std::string kMer = kMers_impliedSequence.at(kMerI);
			if(iskMerDoubleUnique(kMer, occurrencesInSequence))
			{
				break;
			}
			else
			{
				remove_kMers_beginning++;
			}
		}


		for(long long kMerI = (kMers_impliedSequence.size() - 1); kMerI >= 0; kMerI--)
		{
			std::string kMer = kMers_impliedSequence.at(kMerI);
			if(iskMerDoubleUnique(kMer, occurrencesInSequence))
			{
				break;
			}
			else
			{
				remove_kMers_end++;
			}
		}

		if(0)
			std::cout << "Chain " << (&inputChain) << " remove " << remove_kMers_beginning << " beginning and " << remove_kMers_end <<  " ends.\n" << std::flush;

		if((remove_kMers_beginning + remove_kMers_end) < kMers_impliedSequence.size())
		{
			kMerEdgeChain* outputChain = new kMerEdgeChain(*inputChain);
			trimChain(outputChain, sequence, remove_kMers_beginning, remove_kMers_end);
			return outputChain;
		}
		else
		{
			return (new kMerEdgeChain());
		}

	};

	for(size_t chainI = 0; chainI < inputChains.size(); chainI++)
	{
		kMerEdgeChain* chain = inputChains.at(chainI);
		kMerEdgeChain* uniquelyTrimmedChain = trimChainForUniquess(chain);
		if(uniquelyTrimmedChain->traversedEdges.size() > 0)
		{
			uniquelyTrimmedChains.push_back(uniquelyTrimmedChain);
		}
	}

	return uniquelyTrimmedChains;

}

std::vector<std::pair<int, int> > GraphAlignerUnique::findGaps_chainCoverage(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered_input)
{
	std::vector<std::pair<int, int> > gaps;
	long long lastGapStart_sequence = -1;
	for(int seqI = 0; seqI < (int)sequence.length(); seqI++)
	{
		if(sequencePositions_covered_input.at(seqI))
		{
			if(lastGapStart_sequence != -1)
			{
				std::pair<int, int> gapCoordinates;
				gapCoordinates.first = lastGapStart_sequence;
				gapCoordinates.second = seqI-1;
				assert(gapCoordinates.second >= gapCoordinates.first);
				gaps.push_back(gapCoordinates);
			}

			lastGapStart_sequence = -1;
		}
		else
		{
			if(lastGapStart_sequence == -1)
			{
				lastGapStart_sequence = seqI;
			}
		}
	}
	if(lastGapStart_sequence != -1)
	{
		std::pair<int, int> gapCoordinates;
		gapCoordinates.first = lastGapStart_sequence;
		gapCoordinates.second = (sequence.length() - 1);
		gaps.push_back(gapCoordinates);
	}

	return gaps;
}

void GraphAlignerUnique::fixNonUniqueChains(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered, bool thisIterationRandomization, std::vector<kMerEdgeChain*>& allChains, std::vector<kMerEdgeChain*>& newChains)
{
	if(verbose)
		std::cout << Utilities::timestamp() << "GraphAlignerUnique::fixNonUniqueChain(..): Enter function with " << allChains.size() << " total chains.\n" << std::flush;

	int assignedChains;
	int iteration_gapsClosing = 0;
	do {
		iteration_gapsClosing++;
		std::vector<std::pair<int, int> > gaps = findGaps_chainCoverage(sequence, sequencePositions_covered);

		if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::fixNonUniqueChain(..): Process " << gaps.size() << " gaps.\n" << std::flush;

		assignedChains = 0;
		for(unsigned int gapI = 0; gapI < gaps.size(); gapI++)
		{
			closeOneGap_withNonUniqueChains(sequence, gaps.at(gapI), sequencePositions_covered, thisIterationRandomization, allChains, newChains, assignedChains);
		}

		if(verbose) std::cout << "GraphAlignerUnique::fixNonUniqueChain(..), iteration " << iteration_gapsClosing << ": " << assignedChains << " gaps closed.\n";
	} while(assignedChains > 0);

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::fixNonUniqueChain(..): Done.\n" << std::flush;
}

void GraphAlignerUnique::closeOneGap_withNonUniqueChains(std::string& sequence, std::pair<int, int> gapCoordinates, std::vector<kMerEdgeChain*>& sequencePositions_covered, bool thisIterationRandomization, std::vector<kMerEdgeChain*>& chainsToConsider, std::vector<kMerEdgeChain*>& newTrimmedChains, int& assignedChains)
{
	if(verbose) std::cout << Utilities::timestamp() << "closeOneGap_withNonUniqueChains(..): Examine gap from " << gapCoordinates.first << " to " << gapCoordinates.second << "\n" << std::flush;

	int gapLength = gapCoordinates.second - gapCoordinates.first + 1;
	assert(gapLength > 0);
	if(gapLength <= kMerSize)
	{
		return;
	}

	std::string gapSequence(sequence.begin() + gapCoordinates.first, sequence.begin() + gapCoordinates.second + 1);

	Node* leftBoundaryNode = (gapCoordinates.first > 0) ? sequencePositions_covered.at(gapCoordinates.first - 1)->traversedEdges.back()->To : 0;
	Node* rightBoundaryNode = (gapCoordinates.second < ((int)sequencePositions_covered.size() - 1)) ? sequencePositions_covered.at(gapCoordinates.second + 1)->traversedEdges.front()->From : 0;
	int gapBoundary_left_level = (leftBoundaryNode != 0) ? leftBoundaryNode->level : 0;
	int gapBoundary_right_level = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);
	assert(gapBoundary_right_level >= gapBoundary_left_level);
	if(gapBoundary_left_level == gapBoundary_right_level)
	{
		return;
	}
	// filter chains that are relevant to this bit of sequence
	std::vector<kMerEdgeChain*> relevantChains;
	for(unsigned int chainI = 0; chainI < chainsToConsider.size(); chainI++)
	{
		kMerEdgeChain* existingChain = chainsToConsider.at(chainI);
		bool chainIsInteresting = false;
		if((existingChain->sequence_begin <= gapCoordinates.first) && (existingChain->sequence_end >= gapCoordinates.second))
		{
			// overlap
			chainIsInteresting = true;
		}
		else if((existingChain->sequence_begin >= gapCoordinates.first) && (existingChain->sequence_begin <= gapCoordinates.second))
		{
			// start point contained
			chainIsInteresting = true;
		}
		else if((existingChain->sequence_end >= gapCoordinates.first) && (existingChain->sequence_end <= gapCoordinates.second))
		{
			// stop point contained
			chainIsInteresting = true;
		}

		int trimCharacters_left = 0;
		int trimCharacters_right = 0;
		if(chainIsInteresting)
		{
			// check that chain is long enough to contain kMers
			int chainOverLap_left = (existingChain->sequence_begin <= gapCoordinates.first) ? gapCoordinates.first : existingChain->sequence_begin;
			int chainOverLap_right = (existingChain->sequence_end >= gapCoordinates.second) ? gapCoordinates.second : existingChain->sequence_end;
			assert(chainOverLap_right >= chainOverLap_left);
			int chainOverlapLength = chainOverLap_right - chainOverLap_left + 1;
			assert(chainOverlapLength > 0);
			if(chainOverlapLength <= kMerSize)
			{
				chainIsInteresting = false;
			}

			trimCharacters_left = (existingChain->sequence_begin <= gapCoordinates.first) ? (gapCoordinates.first - existingChain->sequence_begin) : 0;
			trimCharacters_right = (existingChain->sequence_end >= gapCoordinates.second) ? (existingChain->sequence_end - gapCoordinates.second) : 0;
		}

		if(chainIsInteresting)
		{
			// now we need to trim the chain accordingly, and then find out whether it is still compatible
			kMerEdgeChain trimmedChain = *existingChain;
			// std::cout << "Chain " << chainI << ", going from " << existingChain->sequence_begin << " to " << existingChain->sequence_end << ": remove " << trimCharacters_left << " on the left and " << trimCharacters_right << " on the right.\n" << std::flush;

			trimChain(&trimmedChain, sequence, trimCharacters_left, trimCharacters_right);
			assert(trimmedChain.sequence_begin >= gapCoordinates.first);
			assert(trimmedChain.sequence_end <= gapCoordinates.second);

			int chain_graph_left = trimmedChain.traversedEdges.front()->From->level;
			int chain_graph_right = trimmedChain.traversedEdges.back()->To->level;

			if((chain_graph_left >= gapBoundary_left_level) && (chain_graph_right <= gapBoundary_right_level))
			{
				kMerEdgeChain* trimmedChain_pointer  = new kMerEdgeChain(trimmedChain);
				newTrimmedChains.push_back(trimmedChain_pointer);
				relevantChains.push_back(trimmedChain_pointer);
			}
		}
	}

	// kick out chains that cannot be reached from the two boundary edges
	cleanChainsAccordingToBoundaries(leftBoundaryNode, rightBoundaryNode, relevantChains);

	// rate chains according to their local optimality
	std::map<kMerEdgeChain*, int> chains_localUniqueness;
	std::function<bool(kMerEdgeChain*,kMerEdgeChain*)> cmpChainUniqueness = [&](kMerEdgeChain* lhs, kMerEdgeChain* rhs) -> bool {
		int rhs_uniqueness = chains_localUniqueness.at(rhs);
		int lhs_uniqueness = chains_localUniqueness.at(lhs);
		if(rhs_uniqueness == lhs_uniqueness)
		{
			return (rhs < lhs);
		}
		else
		{
			return (  rhs_uniqueness < lhs_uniqueness  );
		}
	};
	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> relevantChains_ordered(cmpChainUniqueness);

	std::pair<int, int> gapSpanInGraph = make_pair(gapBoundary_left_level, gapBoundary_right_level);
	analyzeLocalChainUniqueness(
			sequence,
			gapCoordinates,
			gapSpanInGraph,
			relevantChains,
			chains_localUniqueness,
			relevantChains_ordered
	);

	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> remainingChains = relevantChains_ordered;
	
	if(verbose) std::cout << Utilities::timestamp() << "closeOneGap_withNonUniqueChains(..): After analyzeLocalChainUniqueness(..), have  " << remainingChains.size() << " chains in remainingChains." << "\n" << std::flush;

	while((remainingChains.size() > 0) && (chains_localUniqueness.at(*remainingChains.begin()) > minimumChainUniqueness))
	{
		if(thisIterationRandomization)
		{
			std::vector<kMerEdgeChain*> candidateChains;
			std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>::iterator candidateIt = remainingChains.begin();
			assert(candidateIt != remainingChains.end());
			for(int chainI = 1; chainI <= randomizationParameter; chainI++)
			{
				kMerEdgeChain* candidateChain  = *candidateIt;
				if(chains_localUniqueness.at(*remainingChains.begin()) < minimumChainUniqueness)
				{
					break;
				}
				candidateChains.push_back(candidateChain);
				candidateIt++;
				if(candidateIt == remainingChains.end())
					break;
			}
			assert(candidateChains.size() > 0);

			int selectedCandidate = Utilities::randomNumber_nonCritical(candidateChains.size() - 1, &(rng_seeds.at(omp_get_thread_num())));
			kMerEdgeChain* selectedChain = candidateChains.at(selectedCandidate);

			remainingChains.erase(selectedChain);
			for(int seqI = selectedChain->sequence_begin; seqI <= selectedChain->sequence_end; seqI++)
			{
				assert(sequencePositions_covered.at(seqI) == 0);
				sequencePositions_covered.at(seqI) = selectedChain;
			}

			cleanChainsAccordingToSelectedChain(selectedChain, remainingChains);
		}
		else
		{
			kMerEdgeChain* selectedChain = *(remainingChains.begin());
			// std::cerr << "selected chain in sequence from " << selectedChain->sequence_begin << " to " << selectedChain->sequence_end << ", inclusive.\n" << std::flush;

			remainingChains.erase(selectedChain);

			for(int seqI = selectedChain->sequence_begin; seqI <= selectedChain->sequence_end; seqI++)
			{
				assert(sequencePositions_covered.at(seqI) == 0);
				sequencePositions_covered.at(seqI) = selectedChain;
			}

			cleanChainsAccordingToSelectedChain(selectedChain, remainingChains);
		}
		assignedChains++;
	}

	std::set<kMerEdgeChain*> selectedChains_afterGapClosed;
	for(int seqI = 0; seqI < (int)sequencePositions_covered.size(); seqI++)
	{
		if(sequencePositions_covered.at(seqI) != 0)
			selectedChains_afterGapClosed.insert(sequencePositions_covered.at(seqI));
	}
	selectedChainsConsistencyCheck(selectedChains_afterGapClosed);
}

bool GraphAlignerUnique::iskMerLocallyDoubleUnique(std::string& kMer, std::pair<int, int>& spannedRegionCoordinates, std::pair<int, int>& spannedRegionGraphBoundaries, std::map<std::string, int>& occurrencesInSequence)
{
	if(occurrencesInSequence.count(kMer) && (occurrencesInSequence.at(kMer) == 1))
	{
		std::vector<kMerInGraphSpec> kMerPositions = gI.queryIndex(kMer);
		if(gI.queryIndex(kMer).size() >= 1)
		{
			int eligibleCandidates = 0;
			for(unsigned int pI = 0; pI < kMerPositions.size(); pI++)
			{
				kMerInGraphSpec& pos = kMerPositions.at(pI);
				int firstLevel = pos.traversedEdges.front()->From->level;
				int lastLevel = pos.traversedEdges.back()->To->level;
				if((firstLevel >= spannedRegionGraphBoundaries.first) && (lastLevel <= spannedRegionGraphBoundaries.second))
				{
					eligibleCandidates++;
					if(eligibleCandidates > 1)
					{
						break;
					}
				}
			}
			if(eligibleCandidates == 1)
			{
				return true;
			}
		}
	}
	return false;
}

void GraphAlignerUnique::analyzeLocalChainUniqueness(std::string& sequence, std::pair<int, int>& spannedRegionCoordinates, std::pair<int, int>& spannedRegionGraphBoundaries, std::vector<kMerEdgeChain*>& availableChains, 	std::map<kMerEdgeChain*, int>& chains_localUniqueness, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered_output)
{
	assert(spannedRegionCoordinates.second >= spannedRegionCoordinates.first);
	assert(spannedRegionGraphBoundaries.second > spannedRegionGraphBoundaries.first);

	std::string coveredRegion(sequence.begin() + spannedRegionCoordinates.first, sequence.begin() + spannedRegionCoordinates.second + 1);
	std::map<std::string, int> local_kMer_occurrences;
	buildkMerOccurrenceMap(coveredRegion, local_kMer_occurrences);

	for(size_t chainI = 0; chainI < availableChains.size(); chainI++)
	{
		kMerEdgeChain* chain = availableChains.at(chainI);
		std::string impliedSequence(sequence.begin() + chain->sequence_begin, sequence.begin() + chain->sequence_end + 1);
		std::vector<std::string> kMers_impliedSequence = partitionStringIntokMers(impliedSequence, kMerSize);

		size_t kMers_doubleUnique = 0;
		for(size_t kMerI = 0; kMerI < kMers_impliedSequence.size(); kMerI++)
		{
			std::string kMer = kMers_impliedSequence.at(kMerI);
			if(iskMerLocallyDoubleUnique(kMer, spannedRegionCoordinates, spannedRegionGraphBoundaries, local_kMer_occurrences))
			{
				kMers_doubleUnique++;
			}
		}

		assert(kMers_impliedSequence.size() > 0);
		chains_localUniqueness[chain] = kMers_doubleUnique;
		uniquelyTrimmedChains_ordered_output.insert(chain);
	}
}

void GraphAlignerUnique::cleanChainsAccordingToBoundaries(Node* leftBoundaryNode, Node* rightBoundaryNode, std::vector<kMerEdgeChain*>& chains)
{
	int gapBoundary_left_level = (leftBoundaryNode != 0) ? leftBoundaryNode->level : 0;
	int gapBoundary_right_level = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::cleanChainsAccordingToBoundaries(..): Enter function with " << chains.size() << " chains, level from " << gapBoundary_left_level << " to " << gapBoundary_right_level << "\n" << std::flush;

	std::map<Node*, std::set<kMerEdgeChain*> > chainStartNodes;
	std::map<Node*, std::set<kMerEdgeChain*> > chainStopNodes;

	for(unsigned int chainI = 0; chainI < chains.size(); chainI++)
	{
		kMerEdgeChain* chain = chains.at(chainI);
		assert(chain != 0);
		assert(chain->traversedEdges.size() > 0);
		Edge* firstEdge = chain->traversedEdges.front();
		Edge* lastEdge = chain->traversedEdges.back();
		chainStartNodes[firstEdge->From].insert(chain);
		chainStopNodes[lastEdge->To].insert(chain);
	}

	std::set<kMerEdgeChain*> firstEdgeValidated;
	std::set<kMerEdgeChain*> lastEdgeValidated;

	std::map<Node*, std::set<Node*> > runningReachableNodes;
	std::set<Node*> leftBoundaryNodes;
	std::set<Node*> rightBoundaryNodes;
	for(int lI = gapBoundary_left_level; lI <= gapBoundary_right_level; lI++)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::set<Node*> > runningReachableNodes_thisLevel;
		std::set<Edge*> edgesToThisLevel;

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(lI == gapBoundary_left_level)
			{
				if((n == leftBoundaryNode) || ((leftBoundaryNode == 0) && (lI == 0)))
				{
					runningReachableNodes_thisLevel[n].insert(n);
					leftBoundaryNodes.insert(n);
				}
				else
				{
					runningReachableNodes_thisLevel[n] = std::set<Node*>();
				}
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					Node* nodeFrom = e->From;
					assert(runningReachableNodes.count(nodeFrom) > 0);
					runningReachableNodes_thisLevel[n].insert(runningReachableNodes.at(nodeFrom).begin(), runningReachableNodes.at(nodeFrom).end());
				}
			}

			if(chainStopNodes.count(n))
			{
				runningReachableNodes_thisLevel.at(n).insert(n);
			}

			if(chainStartNodes.count(n))
			{
				for(std::set<Node*>::iterator leftBoundaryNodeIt = leftBoundaryNodes.begin(); leftBoundaryNodeIt != leftBoundaryNodes.end(); leftBoundaryNodeIt++)
				{
					Node* thisLeftBoundaryNode = *leftBoundaryNodeIt;
					if(runningReachableNodes_thisLevel.at(n).count(thisLeftBoundaryNode))
					{
						for(std::set<kMerEdgeChain*>::iterator chainStartHereIt = chainStartNodes.at(n).begin(); chainStartHereIt != chainStartNodes.at(n).end(); chainStartHereIt++)
						{
							kMerEdgeChain* chainStartingHere = *chainStartHereIt;
							firstEdgeValidated.insert(chainStartingHere);
						}
						break;
					}
				}
			}

			if(lI == gapBoundary_right_level)
			{
				if((n == rightBoundaryNode) || ((rightBoundaryNode == 0) && (lI == (g->NodesPerLevel.size() - 1))))
				{
					for(unsigned int chainI = 0; chainI < chains.size(); chainI++)
					{
						kMerEdgeChain* chain = chains.at(chainI);
						Node* lastChainNode = chain->traversedEdges.back()->To;
						if(runningReachableNodes_thisLevel.at(n).count(lastChainNode))
						{
							lastEdgeValidated.insert(chain);
						}
					}
				}
				rightBoundaryNodes.insert(n);
			}
		}

		runningReachableNodes = runningReachableNodes_thisLevel;
	}

	assert(leftBoundaryNodes.size() > 0);
	assert(rightBoundaryNodes.size() > 0);

	std::vector<kMerEdgeChain*> remainingChains;
	int removedChains = 0;
	for(unsigned int chainI = 0; chainI < chains.size(); chainI++)
	{
		kMerEdgeChain* chain = chains.at(chainI);
		if(firstEdgeValidated.count(chain) && lastEdgeValidated.count(chain))
		{
			remainingChains.push_back(chain);
		}
		else
		{
			removedChains++;
		}
	}
	chains = remainingChains;

	if(verbose) std::cout << "\t" << "GraphAlignerUnique::cleanChainsAccordingToBoundaries(..): Leave function with " << chains.size() << " chains, removed " <<  removedChains << " chains." << Utilities::timestamp() << "\n" << std::flush;

}

seedAndExtend_return GraphAlignerUnique::seedAndExtend(std::string sequence_nonReverse)
{
	seedAndExtend_return forReturn;

	std::cout << Utilities::timestamp() << " Enter GraphAlignerUnique::seedAndExtend(..)!\n" << std::flush;

	bool useReverse;
	std::string sequence;
	std::vector<std::string> kMers_sequence;
	std::map<std::string, int> kMer_sequence_occurrences;

	seedAndExtend_init_occurrence_strand_etc(
			sequence_nonReverse,
			useReverse,
			sequence,
			kMers_sequence,
			kMer_sequence_occurrences
	);

	if(verbose) std::cout << Utilities::timestamp() << "Find chains.\n" << std::flush;

	std::vector<kMerEdgeChain*> chains_for_sequence = gI.findChains(sequence);

	if(verbose) std::cout << Utilities::timestamp() << "Make " << chains_for_sequence.size() << " chains uniquely ending.\n" << std::flush;

	std::vector<kMerEdgeChain*> uniquelyTrimmedChains = trimChainsForUniqueness(chains_for_sequence, sequence, kMer_sequence_occurrences);

	if(verbose) std::cout << Utilities::timestamp() << "Analyze " << uniquelyTrimmedChains.size() << " trimmed chains.\n" << std::flush;

	std::map<kMerEdgeChain*, int> uniquelyTrimmedChains_doubleUniquekMers;

	std::function<bool(kMerEdgeChain*,kMerEdgeChain*)> cmpChainUniqueness = [&](kMerEdgeChain* lhs, kMerEdgeChain* rhs) -> bool {
		int rhs_uniqueness = uniquelyTrimmedChains_doubleUniquekMers.at(rhs);
		int lhs_uniqueness = uniquelyTrimmedChains_doubleUniquekMers.at(lhs);
		if(lhs_uniqueness == rhs_uniqueness)
		{
			return (rhs < lhs);
		}
		else
		{
			return (  rhs_uniqueness < lhs_uniqueness  );
		}
	};

	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> uniquelyTrimmedChains_ordered(cmpChainUniqueness);

	analyzeChainUniqueness(
			sequence,
			uniquelyTrimmedChains,
			kMer_sequence_occurrences,
			uniquelyTrimmedChains_doubleUniquekMers,
			uniquelyTrimmedChains_ordered
	);
	
	std::vector<seedAndExtend_return> possibleBacktraces;
	std::vector<double> possibleBacktraces_scores;

	std::vector<std::map<std::string, double> > certainty_alignment_sequence;
	std::vector<std::map<std::string, double> > certainty_alignment_graph;
	double certainty_totalIterations = 0;
	certainty_alignment_sequence.resize(sequence.length());
	certainty_alignment_graph.resize(g->NodesPerLevel.size() - 1);

	rng_seeds.resize(threads);
	srand(time(NULL));
	for(unsigned int tI = 0; tI < threads; tI++)
	{
		rng_seeds.at(tI) = rand();
	}

	// very hacky, code duplication
	if(threads > 1)
	{
		omp_set_num_threads(threads);
	
		#pragma omp parallel for
		for(int iI = 0; iI <= iterationsMainRandomizationLoop; iI++)
		{
			assert(omp_get_num_threads() == threads);

			bool thisIterationRandomization = (iI != 0);

			std::set<kMerEdgeChain*> selectedChains;
			std::vector<kMerEdgeChain*> sequencePositions_covered;
			std::vector<kMerEdgeChain*> moreChains;

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)			
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Fix properly unique chains.\n" << std::flush;
				}
			}
			fixUniqueChains(sequence, thisIterationRandomization, uniquelyTrimmedChains_ordered, selectedChains, sequencePositions_covered, uniquelyTrimmedChains_doubleUniquekMers, false);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)			
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Fix non-unique chains.\n" << std::flush;
				}
			}

			fixNonUniqueChains(sequence, sequencePositions_covered, thisIterationRandomization, chains_for_sequence, moreChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)			
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Add chains to vNW.\n" << std::flush;
				}
			}

			VirtualNWTable_Unique vNW(this, &sequence);
			std::map<kMerEdgeChain*, int> currentChains_start;
			std::map<kMerEdgeChain*, NWPath*> chains2Paths;
			kMerEdgeChains2vNW(vNW, sequencePositions_covered, currentChains_start, chains2Paths);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Complete gaps in vNW.\n" << std::flush;
				}
			}

			double finalScore;
			int finalScore_z;
			NWEdge* finalScore_backtrack;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startNormal;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

			vNW_completeRemainingGaps_and_score(
					sequence,
					vNW,
					sequencePositions_covered,
					chains2Paths,
					currentChains_start,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap,
					finalScore,
					finalScore_z,
					finalScore_backtrack
			);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)			
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Backtrack!\n" << std::flush;
				}
			}

			std::string reconstructedSequence;
			std::string reconstructedGraph;
			std::vector<int> reconstructedGraph_levels;

			seedAndExtend_backtrack(
					vNW,
					sequence,
					finalScore,
					finalScore_z,
					finalScore_backtrack,
					reconstructedSequence,
					reconstructedGraph,
					reconstructedGraph_levels,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap
			);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Backtracking done.\n" << std::flush;
				}
			}
			
			seedAndExtend_return thisBacktrace;
			thisBacktrace.Score = finalScore;
			thisBacktrace.graph_aligned = reconstructedGraph;
			thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
			thisBacktrace.sequence_aligned = reconstructedSequence;

			if(omp_get_thread_num() == 2)
			{
				// std::cout << Utilities::timestamp() << "Free some memory.\n" << std::flush;
			}
			
			vNW.freeMemory();
			for(unsigned int moreChainI = 0; moreChainI < moreChains.size(); moreChainI++)
			{
				delete(moreChains.at(moreChainI));
			}
			if(omp_get_thread_num() == 2)
			{
				// std::cout << Utilities::timestamp() << "Freeing done.\n" << std::flush;
			}

			// continue;
			
			#pragma omp critical
			{
				possibleBacktraces.push_back(thisBacktrace);
				possibleBacktraces_scores.push_back(finalScore);

				// std::cerr << Utilities::timestamp() << "Thread " <<  omp_get_thread_num()<< ", iteration " << iI << ", score: " << finalScore << "\n" << std::flush;

				unsigned int alignment_seqI = 0;
				unsigned int alignment_graphI = 0;
				for(unsigned int alignmentI = 0; alignmentI < reconstructedSequence.size(); alignmentI++)
				{
					std::string S = reconstructedSequence.substr(alignmentI, 1);
					std::string G = reconstructedGraph.substr(alignmentI, 1);
					int GL = reconstructedGraph_levels.at(alignmentI);

					std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
					std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
					if(S != "_")
					{
						if(certainty_alignment_sequence.at(alignment_seqI).count(G_Identifier) == 0)
						{
							certainty_alignment_sequence.at(alignment_seqI)[G_Identifier] = 0;
						}
						certainty_alignment_sequence.at(alignment_seqI)[G_Identifier]++;

						alignment_seqI++;
					}
					if(GL != -1)
					{
						if(certainty_alignment_graph.at(alignment_graphI).count(S_Identifier) == 0)
						{
							certainty_alignment_graph.at(alignment_graphI)[S_Identifier] = 0;
						}
						certainty_alignment_graph.at(alignment_graphI)[S_Identifier]++;

						alignment_graphI++;
					}
				}
				certainty_totalIterations++;
			}
		}
	}
	else
	{
		assert(threads == 1);
		
		for(int iI = 0; iI <= iterationsMainRandomizationLoop; iI++)
		{
			bool thisIterationRandomization = (iI != 0);

			std::set<kMerEdgeChain*> selectedChains;
			std::vector<kMerEdgeChain*> sequencePositions_covered;
			std::vector<kMerEdgeChain*> moreChains;

			fixUniqueChains(sequence, thisIterationRandomization, uniquelyTrimmedChains_ordered, selectedChains, sequencePositions_covered, uniquelyTrimmedChains_doubleUniquekMers, false);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			fixNonUniqueChains(sequence, sequencePositions_covered, thisIterationRandomization, chains_for_sequence, moreChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			VirtualNWTable_Unique vNW(this, &sequence);
			std::map<kMerEdgeChain*, int> currentChains_start;
			std::map<kMerEdgeChain*, NWPath*> chains2Paths;
			kMerEdgeChains2vNW(vNW, sequencePositions_covered, currentChains_start, chains2Paths);

			double finalScore;
			int finalScore_z;
			NWEdge* finalScore_backtrack;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startNormal;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

			vNW_completeRemainingGaps_and_score(
					sequence,
					vNW,
					sequencePositions_covered,
					chains2Paths,
					currentChains_start,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap,
					finalScore,
					finalScore_z,
					finalScore_backtrack
			);

			std::string reconstructedSequence;
			std::string reconstructedGraph;
			std::vector<int> reconstructedGraph_levels;

			seedAndExtend_backtrack(
					vNW,
					sequence,
					finalScore,
					finalScore_z,
					finalScore_backtrack,
					reconstructedSequence,
					reconstructedGraph,
					reconstructedGraph_levels,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap
			);

			seedAndExtend_return thisBacktrace;
			thisBacktrace.Score = finalScore;
			thisBacktrace.graph_aligned = reconstructedGraph;
			thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
			thisBacktrace.sequence_aligned = reconstructedSequence;

			vNW.freeMemory();
			for(unsigned int moreChainI = 0; moreChainI < moreChains.size(); moreChainI++)
			{
				delete(moreChains.at(moreChainI));
			}
			
			{
				possibleBacktraces.push_back(thisBacktrace);
				possibleBacktraces_scores.push_back(finalScore);

				// std::cerr << Utilities::timestamp() << "Thread " <<  omp_get_thread_num()<< ", iteration " << iI << ", score: " << finalScore << "\n" << std::flush;

				unsigned int alignment_seqI = 0;
				unsigned int alignment_graphI = 0;
				for(unsigned int alignmentI = 0; alignmentI < reconstructedSequence.size(); alignmentI++)
				{
					std::string S = reconstructedSequence.substr(alignmentI, 1);
					std::string G = reconstructedGraph.substr(alignmentI, 1);
					int GL = reconstructedGraph_levels.at(alignmentI);

					std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
					std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
					if(S != "_")
					{
						if(certainty_alignment_sequence.at(alignment_seqI).count(G_Identifier) == 0)
						{
							certainty_alignment_sequence.at(alignment_seqI)[G_Identifier] = 0;
						}
						certainty_alignment_sequence.at(alignment_seqI)[G_Identifier]++;

						alignment_seqI++;
					}
					if(GL != -1)
					{
						if(certainty_alignment_graph.at(alignment_graphI).count(S_Identifier) == 0)
						{
							certainty_alignment_graph.at(alignment_graphI)[S_Identifier] = 0;
						}
						certainty_alignment_graph.at(alignment_graphI)[S_Identifier]++;

						alignment_graphI++;
					}
				}
				certainty_totalIterations++;
			}
		}		
	}

	std::cerr << std::flush;
	if(verbose)
	{
		std::cout << std::flush << Utilities::timestamp() << "All threads have returned.\n" << std::flush;
	}
	
	// assert( 1 == 0);

	std::pair<double, unsigned int> selectedMaximum = Utilities::findVectorMaxP_nonCritical(possibleBacktraces_scores, &(rng_seeds.at(omp_get_thread_num())));
	seedAndExtend_return selectedBacktrace = possibleBacktraces.at(selectedMaximum.second);

	std::vector<double> selectedBacktrace_certainty_alignment_sequence;
	std::vector<double> selectedBacktrace_certainty_alignment_graph;
	selectedBacktrace_certainty_alignment_sequence.resize(sequence.length(), 0);
	selectedBacktrace_certainty_alignment_graph.resize(g->NodesPerLevel.size() - 1, 0);

	unsigned int alignment_seqI = 0;
	unsigned int alignment_graphI = 0;
	double certainty_sum_sequenceCertainty = 0;
	double certainty_sum_graphCertainty = 0;

	for(unsigned int alignmentI = 0; alignmentI < selectedBacktrace.sequence_aligned.size(); alignmentI++)
	{
		std::string S = selectedBacktrace.sequence_aligned.substr(alignmentI, 1);
		std::string G = selectedBacktrace.graph_aligned.substr(alignmentI, 1);
		int GL = selectedBacktrace.graph_aligned_levels.at(alignmentI);

		std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
		std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
		if(S != "_")
		{
			selectedBacktrace_certainty_alignment_sequence.at(alignment_seqI) = ( certainty_alignment_sequence.at(alignment_seqI).at(G_Identifier) / certainty_totalIterations );
			certainty_sum_sequenceCertainty += selectedBacktrace_certainty_alignment_sequence.at(alignment_seqI);
			alignment_seqI++;
		}
		if(GL != -1)
		{
			selectedBacktrace_certainty_alignment_graph.at(alignment_graphI) = ( certainty_alignment_graph.at(alignment_graphI).at(S_Identifier) / certainty_totalIterations );
			certainty_sum_graphCertainty +=selectedBacktrace_certainty_alignment_graph.at(alignment_graphI);
			alignment_graphI++;
		}
	}


	for(int i = 0; i < chains_for_sequence.size(); i++)
	{
		kMerEdgeChain* c = chains_for_sequence.at(i);
		delete(c);
	}

	for(int i = 0; i < uniquelyTrimmedChains.size(); i++)
	{
		kMerEdgeChain* c = uniquelyTrimmedChains.at(i);
		delete(c);
	}

	double certainty_average_sequenceCertainty = certainty_sum_sequenceCertainty / (double)sequence.length();
	double certainty_average_graphCertainty = certainty_sum_graphCertainty / (double)(g->NodesPerLevel.size() - 1);

	int matches_alignment = countMatchesInSequence(selectedBacktrace.graph_aligned, selectedBacktrace.graph_aligned_levels, selectedBacktrace.sequence_aligned);

	// std::cout << Utilities::timestamp() << " Finished GraphAlignerUnique::seedAndExtend(..)! Average certainty " << certainty_average_sequenceCertainty << " (sequence) / " << certainty_average_graphCertainty << " (graph).; matching positions in sequence: " << matches_alignment << " / " << sequence.length() << ".\n" << std::flush;

	return selectedBacktrace;
}   


std::pair<seedAndExtend_return_local, seedAndExtend_return_local> GraphAlignerUnique::seedAndExtend_local_paired_or_short(oneReadPair readPair, bool usePairing, bool use_short, double insertSize_mean, double insertSize_sd)
{

	// bool verbose = true;
	   
	assert(readPair.reads.first.sequence.find("_") == std::string::npos);
	assert(readPair.reads.first.sequence.find("*") == std::string::npos);
	assert(readPair.reads.first.sequence.find("N") == std::string::npos);
	assert(readPair.reads.second.sequence.find("_") == std::string::npos);
	assert(readPair.reads.second.sequence.find("*") == std::string::npos);
	assert(readPair.reads.second.sequence.find("N") == std::string::npos);

	std::vector<seedAndExtend_return_local> read1_backtraces;
	std::vector<seedAndExtend_return_local> read2_backtraces;

	double minusInfinity = -1 * numeric_limits<double>::max();
	boost::math::normal rnd_InsertSize(insertSize_mean, insertSize_sd);

	//boost::random::normal_distribution<> rnd_InsertSize (insertSize_mean, insertSize_sd);

	if(usePairing)
	{
		if(use_short)
		{
			seedAndExtend_short(readPair.reads.first.sequence, read1_backtraces);
			seedAndExtend_short(readPair.reads.second.sequence, read2_backtraces);
		}
		else
		{
			seedAndExtend_local(readPair.reads.first.sequence, read1_backtraces);
			seedAndExtend_local(readPair.reads.second.sequence, read2_backtraces);
		}

		if(verbose) std::cout << "Read 1 alternatives:\n" << std::flush;
		std::vector<double> likelihoods_read1_alternatives;
		for(unsigned int i = 0; i < read1_backtraces.size(); i++)
		{
			seedAndExtend_return_local& thisBacktrace = read1_backtraces.at(i);
			int ignore;
			double LL = scoreOneAlignment(readPair.reads.first, thisBacktrace, ignore);
			likelihoods_read1_alternatives.push_back(LL);
			
			if(verbose) std::cout << "\t" << i << " " << LL << " [" << thisBacktrace.graph_aligned_levels.front() << " - " << thisBacktrace.graph_aligned_levels.back() << "]\n" << std::flush;
		}
		assert(likelihoods_read1_alternatives.size() == read1_backtraces.size());

		if(verbose) std::cout << "Read 2 alternatives:\n" << std::flush;		
		std::vector<double> likelihoods_read2_alternatives;
		for(unsigned int i = 0; i < read2_backtraces.size(); i++)
		{
			seedAndExtend_return_local& thisBacktrace = read2_backtraces.at(i);
			int ignore;
			double LL = scoreOneAlignment(readPair.reads.second, thisBacktrace, ignore);
			likelihoods_read2_alternatives.push_back(LL);
			
			if(verbose) std::cout << "\t" << i << " " << LL << " [" << thisBacktrace.graph_aligned_levels.front() << " - " << thisBacktrace.graph_aligned_levels.back() << "]\n" << std::flush;
		}
		assert(likelihoods_read2_alternatives.size() == read2_backtraces.size());

		

		std::vector<double> combinedScores;
		std::vector<std::pair<unsigned int, unsigned int> > combinedScores_indices;
		if(verbose) std::cout << "Read combinations alternatives:\n" << std::flush;
		for(unsigned int aI1 = 0; aI1 < read1_backtraces.size(); aI1++)
		{
			for(unsigned int aI2 = 0; aI2 < read2_backtraces.size(); aI2++)
			{
				double combinedScore = likelihoods_read1_alternatives.at(aI1) + likelihoods_read2_alternatives.at(aI2);

				int distance_graph_levels = read2_backtraces.at(aI2).graph_aligned_levels.front() - read1_backtraces.at(aI1).graph_aligned_levels.back();
				double distance_graph_levels_P = boost::math::pdf(rnd_InsertSize, distance_graph_levels);

				assert((distance_graph_levels_P >= 0) && (distance_graph_levels_P <= 1));

				combinedScore += log(distance_graph_levels_P);

				if(read1_backtraces.at(aI1).reverse == read2_backtraces.at(aI2).reverse)
				{
					combinedScore = minusInfinity;
				}

				combinedScores.push_back(combinedScore);
				combinedScores_indices.push_back(std::make_pair(aI1, aI2));
				
			
				if(verbose)
				{
					std::cout << "\t" << aI1 << "/" << aI2 << ": ";
					std::cout << "LL1: " << likelihoods_read1_alternatives.at(aI1) << " LL2: " << likelihoods_read2_alternatives.at(aI2) << " ";
					std::cout << "REV1: " << read1_backtraces.at(aI1).reverse << " REV2:" << read2_backtraces.at(aI2).reverse << " ";
					std::cout << "Distance: " <<  read1_backtraces.at(aI1).graph_aligned_levels.back() << " to " << read2_backtraces.at(aI2).graph_aligned_levels.front() << ", i.e. " << distance_graph_levels << " (" << log(distance_graph_levels_P) << ")" << "  == > " << combinedScore <<  "\n" << std::flush;
				}
			}
		}

		std::pair<double, unsigned int> bestCombination = Utilities::findVectorMax(combinedScores);

		if(verbose)
		{
			std::cout << "CHOOSE alternative " << combinedScores_indices.at(bestCombination.second).first << " / " << combinedScores_indices.at(bestCombination.second).second << "\n" << std::flush;
		}
		
		std::pair<seedAndExtend_return_local, seedAndExtend_return_local> forReturn;
		forReturn.first = read1_backtraces.at(combinedScores_indices.at(bestCombination.second).first);
		forReturn.second = read2_backtraces.at(combinedScores_indices.at(bestCombination.second).second);
		return forReturn;
	}
	else
	{
		std::pair<seedAndExtend_return_local, seedAndExtend_return_local> forReturn;
		forReturn.first = (use_short)? seedAndExtend_short(readPair.reads.first.sequence, read1_backtraces) : seedAndExtend_local(readPair.reads.first.sequence, read1_backtraces);
		forReturn.second = (use_short) ? seedAndExtend_short(readPair.reads.second.sequence, read2_backtraces) : seedAndExtend_local(readPair.reads.second.sequence, read2_backtraces);
		return forReturn;
	}
}

seedAndExtend_return_local GraphAlignerUnique::seedAndExtend_local(std::string sequence_nonReverse, std::vector<seedAndExtend_return_local>& allBacktraces)
{
	seedAndExtend_return_local forReturn;
	
	verbose = false;
	
	// std::cout << Utilities::timestamp() << " Enter GraphAlignerUnique::seedAndExtend(..)!\n" << std::flush;

	bool useReverse;
	std::string sequence;
	std::vector<std::string> kMers_sequence;
	std::map<std::string, int> kMer_sequence_occurrences;

	seedAndExtend_init_occurrence_strand_etc(
			sequence_nonReverse,
			useReverse,
			sequence,
			kMers_sequence,
			kMer_sequence_occurrences
	);

	if(verbose) std::cout << Utilities::timestamp() << "Find chains.\n" << std::flush;

	std::vector<kMerEdgeChain*> chains_for_sequence = gI.findChains(sequence);

	if(verbose) std::cout << Utilities::timestamp() << "Make chain " << chains_for_sequence.size() << " uniquely ending.\n" << std::flush;

	std::vector<kMerEdgeChain*> uniquelyTrimmedChains = trimChainsForUniqueness(chains_for_sequence, sequence, kMer_sequence_occurrences);

	if(verbose) std::cout << Utilities::timestamp() << "Analyze " << uniquelyTrimmedChains.size() << " trimmed chains.\n" << std::flush;

	bool rescureNonUniqueChains = false;
	if(uniquelyTrimmedChains.size() == 0)
	{
		if(verbose) std::cout << Utilities::timestamp() << "Uniqueness-trimming has removed all chains, rescure by accepting all.\n" << std::flush;
		for(unsigned int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
		{
			kMerEdgeChain* chainCopy = new kMerEdgeChain(*(chains_for_sequence.at(chainI)));
			uniquelyTrimmedChains.push_back(chainCopy);
		}
		if(verbose) std::cout << Utilities::timestamp() << "Rescured, now have " << uniquelyTrimmedChains.size() << " chains.\n" << std::flush;
		rescureNonUniqueChains = true;
	}
	
	std::map<kMerEdgeChain*, int> uniquelyTrimmedChains_doubleUniquekMers;

	std::function<bool(kMerEdgeChain*,kMerEdgeChain*)> cmpChainUniqueness = [&](kMerEdgeChain* lhs, kMerEdgeChain* rhs) -> bool {
		int rhs_uniqueness = uniquelyTrimmedChains_doubleUniquekMers.at(rhs);
		int lhs_uniqueness = uniquelyTrimmedChains_doubleUniquekMers.at(lhs);
		if(lhs_uniqueness == rhs_uniqueness)
		{
			int sequence_length_rhs = rhs->sequence_end - rhs->sequence_begin;
			assert(sequence_length_rhs > 0);
			int sequence_length_lhs = lhs->sequence_end - lhs->sequence_begin;
			assert(sequence_length_lhs > 0);
			if(sequence_length_rhs == sequence_length_lhs)
			{
				return (rhs < lhs);
			}
			else
			{
				return (sequence_length_rhs < sequence_length_lhs);
			}
		}
		else
		{
			return (  rhs_uniqueness < lhs_uniqueness  );
		}
	};
	

	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> uniquelyTrimmedChains_ordered(cmpChainUniqueness);

	analyzeChainUniqueness(
			sequence,
			uniquelyTrimmedChains,
			kMer_sequence_occurrences,
			uniquelyTrimmedChains_doubleUniquekMers,
			uniquelyTrimmedChains_ordered
	);

	if(verbose) std::cout << Utilities::timestamp() << "After analyzeChainUniqueness(..), have " << uniquelyTrimmedChains_ordered.size() << " elements in uniquelyTrimmedChains_ordered trimmed chains.\n" << std::flush;

	
	std::vector<seedAndExtend_return> possibleBacktraces;
	std::vector<double> possibleBacktraces_scores;

	std::vector<std::map<std::string, double> > certainty_alignment_sequence;
	std::vector<std::map<std::string, double> > certainty_alignment_graph;
	double certainty_totalIterations = 0;
	certainty_alignment_sequence.resize(sequence.length());
	certainty_alignment_graph.resize(g->NodesPerLevel.size() - 1);

	int required_rng_seeds;
	if(threads == 1)
	{
		required_rng_seeds = omp_get_num_threads();
		if(required_rng_seeds < 1)
		{
			required_rng_seeds = 1;
		}
	}
	else
	{
		required_rng_seeds = threads;
	}
	
	rng_seeds.resize(required_rng_seeds);
	srand(time(NULL));
	for(unsigned int tI = 0; tI < required_rng_seeds; tI++)
	{
		rng_seeds.at(tI) = rand();
	}

	// very hacky, code duplication
	if(threads > 1)
	{
		// todo remove later
		assert(1 == 0);
		omp_set_num_threads(threads);
	
		#pragma omp parallel for
		for(int iI = 0; iI <= iterationsMainRandomizationLoop; iI++)
		{
			assert(omp_get_num_threads() == threads);

			bool thisIterationRandomization = (iI != 0);

			std::set<kMerEdgeChain*> selectedChains;
			std::vector<kMerEdgeChain*> sequencePositions_covered;
			std::vector<kMerEdgeChain*> moreChains;

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Fix properly unique chains.\n" << std::flush;
				}
			}
			
			fixUniqueChains(sequence, thisIterationRandomization, uniquelyTrimmedChains_ordered, selectedChains, sequencePositions_covered, uniquelyTrimmedChains_doubleUniquekMers, rescureNonUniqueChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Fix non-unique chains.\n" << std::flush;
				}
			}

			fixNonUniqueChains(sequence, sequencePositions_covered, thisIterationRandomization, chains_for_sequence, moreChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Add chains to vNW.\n" << std::flush;
				}
			}

			VirtualNWTable_Unique vNW(this, &sequence);
			std::map<kMerEdgeChain*, int> currentChains_start;
			std::map<kMerEdgeChain*, NWPath*> chains2Paths;
			kMerEdgeChains2vNW(vNW, sequencePositions_covered, currentChains_start, chains2Paths);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Complete gaps in vNW.\n" << std::flush;
				}
			}

			double finalScore;
			int finalScore_z;
			NWEdge* finalScore_backtrack;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startNormal;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

			vNW_completeRemainingGaps_and_score_local(
					sequence,
					vNW,
					sequencePositions_covered,
					chains2Paths,
					currentChains_start,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap,
					finalScore,
					finalScore_z,
					finalScore_backtrack
			);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Backtrack!\n" << std::flush;
				}
			}

			std::string reconstructedSequence;
			std::string reconstructedGraph;
			std::vector<int> reconstructedGraph_levels;

			seedAndExtend_backtrack_local(
					vNW,
					sequence,
					finalScore,
					finalScore_z,
					finalScore_backtrack,
					reconstructedSequence,
					reconstructedGraph,
					reconstructedGraph_levels,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap
			);

			if(verbose)
			{
				//pragma omp critical
				if(omp_get_thread_num() == 2)
				{
					std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << "; Backtracking done.\n" << std::flush;
				}
			}

			seedAndExtend_return thisBacktrace;
			thisBacktrace.Score = finalScore;
			thisBacktrace.graph_aligned = reconstructedGraph;
			thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
			thisBacktrace.sequence_aligned = reconstructedSequence;

			if(omp_get_thread_num() == 2)
			{
				// std::cout << Utilities::timestamp() << "Free some memory.\n" << std::flush;
			}

			vNW.freeMemory();
			for(unsigned int moreChainI = 0; moreChainI < moreChains.size(); moreChainI++)
			{
				delete(moreChains.at(moreChainI));
			}
			if(omp_get_thread_num() == 2)
			{
				// std::cout << Utilities::timestamp() << "Freeing done.\n" << std::flush;
			}

			// continue;

			#pragma omp critical
			{
				possibleBacktraces.push_back(thisBacktrace);
				possibleBacktraces_scores.push_back(finalScore);

				// std::cerr << Utilities::timestamp() << "Thread " <<  omp_get_thread_num()<< ", iteration " << iI << ", score: " << finalScore << "\n" << std::flush;

				unsigned int alignment_seqI = 0;
				unsigned int alignment_graphI = 0;
				for(unsigned int alignmentI = 0; alignmentI < reconstructedSequence.size(); alignmentI++)
				{
					std::string S = reconstructedSequence.substr(alignmentI, 1);
					std::string G = reconstructedGraph.substr(alignmentI, 1);
					int GL = reconstructedGraph_levels.at(alignmentI);

					std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
					std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
					if(S != "_")
					{
						if(certainty_alignment_sequence.at(alignment_seqI).count(G_Identifier) == 0)
						{
							certainty_alignment_sequence.at(alignment_seqI)[G_Identifier] = 0;
						}
						certainty_alignment_sequence.at(alignment_seqI)[G_Identifier]++;

						alignment_seqI++;
					}
					if(GL != -1)
					{
						if(certainty_alignment_graph.at(alignment_graphI).count(S_Identifier) == 0)
						{
							certainty_alignment_graph.at(alignment_graphI)[S_Identifier] = 0;
						}
						certainty_alignment_graph.at(alignment_graphI)[S_Identifier]++;

						alignment_graphI++;
					}
				}
				certainty_totalIterations++;
			}
		}
	}
	else
	{
		assert(threads == 1);

		for(int iI = 0; iI <= iterationsMainRandomizationLoop; iI++)
		{
			if(verbose)
			{
				std::cout << "Main iterations loop: " << iI << "\n" << std::flush;
			}
			bool thisIterationRandomization = (iI != 0);

			std::set<kMerEdgeChain*> selectedChains;
			std::vector<kMerEdgeChain*> sequencePositions_covered;
			std::vector<kMerEdgeChain*> moreChains;

			fixUniqueChains(sequence, thisIterationRandomization, uniquelyTrimmedChains_ordered, selectedChains, sequencePositions_covered, uniquelyTrimmedChains_doubleUniquekMers, rescureNonUniqueChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);


			fixNonUniqueChains(sequence, sequencePositions_covered, thisIterationRandomization, chains_for_sequence, moreChains);
			if(iI == 0) printSequenceChainCoverageStats(sequence, sequencePositions_covered);

			VirtualNWTable_Unique vNW(this, &sequence);
			std::map<kMerEdgeChain*, int> currentChains_start;
			std::map<kMerEdgeChain*, NWPath*> chains2Paths;
			kMerEdgeChains2vNW(vNW, sequencePositions_covered, currentChains_start, chains2Paths);

			double finalScore;
			int finalScore_z;
			NWEdge* finalScore_backtrack;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startNormal;
			std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
			std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

			vNW_completeRemainingGaps_and_score_local(
					sequence,
					vNW,
					sequencePositions_covered,
					chains2Paths,
					currentChains_start,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap,
					finalScore,
					finalScore_z,
					finalScore_backtrack
			);

			std::string reconstructedSequence;
			std::string reconstructedGraph;
			std::vector<int> reconstructedGraph_levels;

			seedAndExtend_backtrack_local(
					vNW,
					sequence,
					finalScore,
					finalScore_z,
					finalScore_backtrack,
					reconstructedSequence,
					reconstructedGraph,
					reconstructedGraph_levels,
					lastPositionDistances_perZ_startNormal,
					lastPositionDistances_perZ_startAffineGap,
					NWedges_graphDist_startNormal,
					NWedges_graphDist_startAffineGap
			);

			seedAndExtend_return thisBacktrace;
			thisBacktrace.Score = finalScore;
			thisBacktrace.graph_aligned = reconstructedGraph;
			thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
			thisBacktrace.sequence_aligned = reconstructedSequence;

			vNW.freeMemory();
			for(unsigned int moreChainI = 0; moreChainI < moreChains.size(); moreChainI++)
			{
				delete(moreChains.at(moreChainI));
			}

			{
				possibleBacktraces.push_back(thisBacktrace);
				possibleBacktraces_scores.push_back(finalScore);

				// std::cerr << Utilities::timestamp() << "Thread " <<  omp_get_thread_num()<< ", iteration " << iI << ", score: " << finalScore << "\n" << std::flush;

				unsigned int alignment_seqI = 0;
				unsigned int alignment_graphI = 0;
				for(unsigned int alignmentI = 0; alignmentI < reconstructedSequence.size(); alignmentI++)
				{
					std::string S = reconstructedSequence.substr(alignmentI, 1);
					std::string G = reconstructedGraph.substr(alignmentI, 1);
					int GL = reconstructedGraph_levels.at(alignmentI);

					std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
					std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
					if(S != "_")
					{
						if(certainty_alignment_sequence.at(alignment_seqI).count(G_Identifier) == 0)
						{
							certainty_alignment_sequence.at(alignment_seqI)[G_Identifier] = 0;
						}
						certainty_alignment_sequence.at(alignment_seqI)[G_Identifier]++;

						alignment_seqI++;
					}
					if(GL != -1)
					{
						if(certainty_alignment_graph.at(alignment_graphI).count(S_Identifier) == 0)
						{
							certainty_alignment_graph.at(alignment_graphI)[S_Identifier] = 0;
						}
						certainty_alignment_graph.at(alignment_graphI)[S_Identifier]++;

						alignment_graphI++;
					}
				}
				certainty_totalIterations++;
			}
		}		
	}
	std::cerr << std::flush;
	if(verbose)
	{
		std::cout << std::flush << Utilities::timestamp() << "All threads have returned.\n" << std::flush;
	}

	// assert( 1 == 0);

	std::pair<double, unsigned int> selectedMaximum = Utilities::findVectorMaxP_nonCritical(possibleBacktraces_scores, &(rng_seeds.at(omp_get_thread_num())));
	seedAndExtend_return selectedBacktrace = possibleBacktraces.at(selectedMaximum.second);

	std::vector<double> selectedBacktrace_certainty_alignment_sequence;
	std::vector<double> selectedBacktrace_certainty_alignment_graph;
	selectedBacktrace_certainty_alignment_sequence.resize(sequence.length(), 0);
	selectedBacktrace_certainty_alignment_graph.resize(g->NodesPerLevel.size() - 1, 0);

	unsigned int alignment_seqI = 0;
	unsigned int alignment_graphI = 0;
	double certainty_sum_sequenceCertainty = 0;
	double certainty_sum_graphCertainty = 0;

	for(unsigned int alignmentI = 0; alignmentI < selectedBacktrace.sequence_aligned.size(); alignmentI++)
	{
		std::string S = selectedBacktrace.sequence_aligned.substr(alignmentI, 1);
		std::string G = selectedBacktrace.graph_aligned.substr(alignmentI, 1);
		int GL = selectedBacktrace.graph_aligned_levels.at(alignmentI);

		std::string S_Identifier = Utilities::ItoStr(alignment_seqI)+"-"+S;
		std::string G_Identifier = Utilities::ItoStr(GL)+"-"+G;
		if(S != "_")
		{
			selectedBacktrace_certainty_alignment_sequence.at(alignment_seqI) = ( certainty_alignment_sequence.at(alignment_seqI).at(G_Identifier) / certainty_totalIterations );
			certainty_sum_sequenceCertainty += selectedBacktrace_certainty_alignment_sequence.at(alignment_seqI);
			alignment_seqI++;
		}
		if(GL != -1)
		{
			selectedBacktrace_certainty_alignment_graph.at(alignment_graphI) = ( certainty_alignment_graph.at(alignment_graphI).at(S_Identifier) / certainty_totalIterations );
			certainty_sum_graphCertainty +=selectedBacktrace_certainty_alignment_graph.at(alignment_graphI);
			alignment_graphI++;
		}
	}

	// compute uniqueness
	std::vector<std::string> kMers_impliedSequence = partitionStringIntokMers(sequence, kMerSize);

	size_t kMers_unique = 0;
	for(size_t kMerI = 0; kMerI < kMers_impliedSequence.size(); kMerI++)
	{
		std::string kMer = kMers_impliedSequence.at(kMerI);
		if(iskMerDoubleUnique(kMer, kMer_sequence_occurrences))
		{
			kMers_unique++;
		}
	}

	size_t kMers_unique_utilized = 0;
	if((int)selectedBacktrace.graph_aligned.length() >= kMerSize)
	{
		for(int windowStop = (int)(kMerSize - 1); windowStop < (int)selectedBacktrace.graph_aligned.length(); windowStop++)
		{
			int windowStart = windowStop - kMerSize + 1;
			assert(windowStart >= 0);

			std::string kMer_sequence = selectedBacktrace.sequence_aligned.substr(windowStart, windowStop - windowStart + 1);
			std::string kMer_graph = selectedBacktrace.graph_aligned.substr(windowStart, windowStop - windowStart + 1);
			kMer_sequence.erase(std::remove_if(kMer_sequence.begin(),kMer_sequence.end(), [&](char c){return ((c == '_') ? true : false);}), kMer_sequence.end());
			kMer_graph.erase(std::remove_if(kMer_graph.begin(),kMer_graph.end(), [&](char c){return ((c == '_') ? true : false);}), kMer_graph.end());

			if((kMer_sequence == kMer_graph) && ((int)kMer_graph.length() == kMerSize))
			{
				if(iskMerDoubleUnique(kMer_sequence, kMer_sequence_occurrences))
				{
					kMers_unique_utilized++;
				}
			}
		}
	}

	allBacktraces.clear();
	for(unsigned int bI = 0; bI < possibleBacktraces.size(); bI++)
	{
		seedAndExtend_return& possibleBacktrace = possibleBacktraces.at(bI);
		seedAndExtend_return_local possibleBacktrace_forReturn;

		possibleBacktrace_forReturn.Score = possibleBacktrace.Score;
		possibleBacktrace_forReturn.graph_aligned = possibleBacktrace.graph_aligned;
		possibleBacktrace_forReturn.sequence_aligned = possibleBacktrace.sequence_aligned;
		possibleBacktrace_forReturn.graph_aligned_levels = possibleBacktrace.graph_aligned_levels;
		possibleBacktrace_forReturn.reverse = useReverse;
		allBacktraces.push_back(possibleBacktrace_forReturn);
	}

	seedAndExtend_return_local selectedBacktrace_forReturn;
	selectedBacktrace_forReturn.Score = selectedBacktrace.Score;
	selectedBacktrace_forReturn.graph_aligned = selectedBacktrace.graph_aligned;
	selectedBacktrace_forReturn.sequence_aligned = selectedBacktrace.sequence_aligned;
	selectedBacktrace_forReturn.graph_aligned_levels = selectedBacktrace.graph_aligned_levels;
	selectedBacktrace_forReturn.certainty_sequence2Graph = selectedBacktrace_certainty_alignment_sequence;
	selectedBacktrace_forReturn.kMers_total = kMers_impliedSequence.size();
	selectedBacktrace_forReturn.kMers_unique_total = kMers_unique;
	selectedBacktrace_forReturn.kMers_unique_utilized = kMers_unique_utilized;
	selectedBacktrace_forReturn.reverse = useReverse;

	double certainty_average_sequenceCertainty = certainty_sum_sequenceCertainty / (double)sequence.length();
	double certainty_average_graphCertainty = certainty_sum_graphCertainty / (double)(g->NodesPerLevel.size() - 1);

	int matches_alignment = countMatchesInSequence(selectedBacktrace.graph_aligned, selectedBacktrace.graph_aligned_levels, selectedBacktrace.sequence_aligned);

	// std::cout << Utilities::timestamp() << " Finished GraphAlignerUnique::seedAndExtend_local(..)! Average certainty " << certainty_average_sequenceCertainty << " (sequence) / " << certainty_average_graphCertainty << " (graph).; matching positions in sequence: " << matches_alignment << " / " << sequence.length() << ".\n" << std::flush;

	for(unsigned int i = 0; i < chains_for_sequence.size(); i++)
	{
		kMerEdgeChain* c = chains_for_sequence.at(i);
		delete(c);
	}

	for(unsigned int i = 0; i < uniquelyTrimmedChains.size(); i++)
	{
		kMerEdgeChain* c = uniquelyTrimmedChains.at(i);
		delete(c);
	}

	return selectedBacktrace_forReturn;
}


seedAndExtend_return_local GraphAlignerUnique::seedAndExtend_short(std::string sequence_nonReverse, std::vector<seedAndExtend_return_local>& allBacktraces)
{
	seedAndExtend_return_local forReturn;

	verbose = true;
	int seedAndExtend_short_maximumBacktraces = 5;

	// just to remind everyone that this is for short reads
	assert(sequence_nonReverse.length() < 200);

	if(verbose)
		std::cout << Utilities::timestamp() << " Enter GraphAlignerUnique::seedAndExtend_short(..)!\n" << std::flush;

	bool useReverse;
	std::string sequence;
	std::vector<std::string> kMers_sequence;
	std::map<std::string, int> kMer_sequence_occurrences;

	seedAndExtend_init_occurrence_strand_etc(
			sequence_nonReverse,
			useReverse,
			sequence,
			kMers_sequence,
			kMer_sequence_occurrences
	);

	if(verbose) std::cout << Utilities::timestamp() << "Find chains.\n" << std::flush;

	std::vector<kMerEdgeChain*> chains_for_sequence = gI.findChains(sequence);

	std::function<bool(kMerEdgeChain*,kMerEdgeChain*)> cmpChainLength = [&](kMerEdgeChain* lhs, kMerEdgeChain* rhs) -> bool {
		int sequence_length_rhs = rhs->sequence_end - rhs->sequence_begin;
		assert(sequence_length_rhs > 0);
		int sequence_length_lhs = lhs->sequence_end - lhs->sequence_begin;
		assert(sequence_length_lhs > 0);
		if(sequence_length_rhs == sequence_length_lhs)
		{
			return (rhs < lhs);
		}
		else
		{
			return (sequence_length_rhs < sequence_length_lhs);
		}
	};

	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> chains_orderedByLength(cmpChainLength);
	for(unsigned int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
	{
		chains_orderedByLength.insert(chains_for_sequence.at(chainI));
	}
	assert(chains_orderedByLength.size() == chains_for_sequence.size());

	int required_rng_seeds;
	if(threads == 1)
	{
		required_rng_seeds = omp_get_num_threads();
		if(required_rng_seeds < 1)
		{
			required_rng_seeds = 1;
		}
	}
	else
	{
		required_rng_seeds = threads;
	}

	rng_seeds.resize(required_rng_seeds);
	for(unsigned int tI = 0; tI < required_rng_seeds; tI++)
	{
		rng_seeds.at(tI) = 0;
	}


	std::vector<seedAndExtend_return> possibleBacktraces;
	std::vector<double> possibleBacktraces_scores;

	assert(threads == 1);
	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>::iterator currentChain = chains_orderedByLength.begin();

	double best_chain_totalScore;
	double best_chain_subOptimality;
	int best_chain_initialMatches;

	for(int chainI = 0; chainI < chains_for_sequence.size(); chainI++)
	{
		assert(currentChain != chains_orderedByLength.end());
		kMerEdgeChain* thisChain = *currentChain;

		int matches_this_chain = thisChain->sequence_end - thisChain->sequence_begin + 1;
		assert(matches_this_chain > 0);

		if(verbose)
					std::cout << "\tChain " << chainI << "/" << chains_for_sequence.size() << ", going from " << thisChain->sequence_begin << " to " << thisChain->sequence_end << "\n" << std::flush;

		if(chainI > 0)
		{
			int matches_less_than_best_chain = best_chain_initialMatches - matches_this_chain;
			assert(matches_less_than_best_chain >= 0);
			double forgone_score = matches_less_than_best_chain * S_match;
			if((1*forgone_score) > best_chain_subOptimality)
			{
				if(verbose)
				{
					std::cout << "\t\tABORT. This chain has " << matches_this_chain << ", which is " << matches_less_than_best_chain << " less than the best chain. We assume that " << 0.5*forgone_score << " of this are lost, and the best chain is " << best_chain_subOptimality << " away from optimality.\n" << std::flush;
				}
				break;
			}
		}
		std::set<kMerEdgeChain*> selectedChains;
		std::vector<kMerEdgeChain*> sequencePositions_covered;
		sequencePositions_covered.resize(sequence.length());

		selectedChains.insert(thisChain);

		for(int seqI = thisChain->sequence_begin; seqI <= thisChain->sequence_end; seqI++)
		{
			assert(sequencePositions_covered.at(seqI) == 0);
			sequencePositions_covered.at(seqI) = thisChain;
		}

		VirtualNWTable_Unique vNW(this, &sequence);
		std::map<kMerEdgeChain*, int> currentChains_start;
		std::map<kMerEdgeChain*, NWPath*> chains2Paths;
		kMerEdgeChains2vNW(vNW, sequencePositions_covered, currentChains_start, chains2Paths);

		double finalScore;
		int finalScore_z;
		NWEdge* finalScore_backtrack;
		std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startNormal;
		std::vector<std::map<NWEdge*, graphPointDistance> > lastPositionDistances_perZ_startAffineGap;
		std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
		std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

		vNW_completeRemainingGaps_and_score_local(
				sequence,
				vNW,
				sequencePositions_covered,
				chains2Paths,
				currentChains_start,
				lastPositionDistances_perZ_startNormal,
				lastPositionDistances_perZ_startAffineGap,
				NWedges_graphDist_startNormal,
				NWedges_graphDist_startAffineGap,
				finalScore,
				finalScore_z,
				finalScore_backtrack
		);

		std::string reconstructedSequence;
		std::string reconstructedGraph;
		std::vector<int> reconstructedGraph_levels;

		seedAndExtend_backtrack_local(
				vNW,
				sequence,
				finalScore,
				finalScore_z,
				finalScore_backtrack,
				reconstructedSequence,
				reconstructedGraph,
				reconstructedGraph_levels,
				lastPositionDistances_perZ_startNormal,
				lastPositionDistances_perZ_startAffineGap,
				NWedges_graphDist_startNormal,
				NWedges_graphDist_startAffineGap
		);

		seedAndExtend_return thisBacktrace;
		thisBacktrace.Score = finalScore;
		thisBacktrace.graph_aligned = reconstructedGraph;
		thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
		thisBacktrace.sequence_aligned = reconstructedSequence;

		possibleBacktraces.push_back(thisBacktrace);
		possibleBacktraces_scores.push_back(finalScore);

		vNW.freeMemory();

		if(verbose)
			std::cout << "\t\tScore " << finalScore << "\n" << std::flush;

		int score_for_matches_from_chain = matches_this_chain * S_match;
		int optimal_chain_score = S_match * sequence.length();
		int alignment_suboptimality = optimal_chain_score - finalScore;
		assert(alignment_suboptimality >= 0);


		if((chainI == 0) || (best_chain_totalScore < finalScore))
		{
			best_chain_totalScore = finalScore;
			best_chain_subOptimality = alignment_suboptimality;
			best_chain_initialMatches = matches_this_chain;
		}

		currentChain++;
	}
	// assert(currentChain == chains_orderedByLength.end());

	std::cerr << std::flush;
	if(verbose)
	{
		std::cout << std::flush << Utilities::timestamp() << "All chains have been examined have returned.\n" << std::flush;
	}

	std::pair<double, unsigned int> selectedMaximum = Utilities::findVectorMax(possibleBacktraces_scores);
	seedAndExtend_return selectedBacktrace = possibleBacktraces.at(selectedMaximum.second);

	std::vector<std::pair<double, unsigned int>> backtraces_scores_and_positions;
	for(unsigned int i = 0; i < possibleBacktraces.size(); i++)
	{
		std::pair<double, unsigned int> thisBacktrace_score_and_position = std::make_pair(possibleBacktraces_scores.at(i), i);
		backtraces_scores_and_positions.push_back(thisBacktrace_score_and_position);
	}
	std::sort(backtraces_scores_and_positions.begin(), backtraces_scores_and_positions.end(), [](std::pair<double, unsigned int> a, std::pair<double, unsigned int> b) {
		return (a.first < b.first);
	});
	std::reverse(backtraces_scores_and_positions.begin(), backtraces_scores_and_positions.end());
	if(backtraces_scores_and_positions.size() > 1)
	{
		assert(possibleBacktraces_scores.at(backtraces_scores_and_positions.at(0).second) >= possibleBacktraces_scores.at(backtraces_scores_and_positions.at(1).second));
	}
	assert((possibleBacktraces_scores.size() == 0) || (selectedMaximum.first == possibleBacktraces_scores.at(backtraces_scores_and_positions.at(0).second)));

	allBacktraces.clear();
	unsigned int maxNumberReturnedBacktraces = (backtraces_scores_and_positions.size() > seedAndExtend_short_maximumBacktraces) ? seedAndExtend_short_maximumBacktraces : backtraces_scores_and_positions.size();
	for(unsigned int bI = 0; bI < maxNumberReturnedBacktraces; bI++)
	{
		seedAndExtend_return& possibleBacktrace = possibleBacktraces.at(backtraces_scores_and_positions.at(bI).second);
		seedAndExtend_return_local possibleBacktrace_forReturn;

		possibleBacktrace_forReturn.Score = possibleBacktrace.Score;
		possibleBacktrace_forReturn.graph_aligned = possibleBacktrace.graph_aligned;
		possibleBacktrace_forReturn.sequence_aligned = possibleBacktrace.sequence_aligned;
		possibleBacktrace_forReturn.graph_aligned_levels = possibleBacktrace.graph_aligned_levels;
		possibleBacktrace_forReturn.reverse = useReverse;

		allBacktraces.push_back(possibleBacktrace_forReturn);
	}

	seedAndExtend_return_local selectedBacktrace_forReturn;
	selectedBacktrace_forReturn.Score = selectedBacktrace.Score;
	selectedBacktrace_forReturn.graph_aligned = selectedBacktrace.graph_aligned;
	selectedBacktrace_forReturn.sequence_aligned = selectedBacktrace.sequence_aligned;
	selectedBacktrace_forReturn.graph_aligned_levels = selectedBacktrace.graph_aligned_levels;
	selectedBacktrace_forReturn.reverse = useReverse;

	for(unsigned int i = 0; i < chains_for_sequence.size(); i++)
	{
		kMerEdgeChain* c = chains_for_sequence.at(i);
		delete(c);
	}

	return selectedBacktrace_forReturn;
}

void GraphAlignerUnique::vNW_completeRemainingGaps_and_score(std::string& sequence, VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, NWPath*>& chains2Paths, std::map<kMerEdgeChain*, int>& currentChains_start, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, double& finalScore, int& finalScore_z, NWEdge*& finalScore_backtrack)
{
	bool verbose = false;
	bool superquiet = false;

	double minusInfinity = -1 * numeric_limits<double>::max();

	auto scoreEdgeJump = [&](NWEdge* entryEdge, NWEdge* exitEdge, graphPointDistance& graphDistance_startAffineGap, graphPointDistance& graphDistance_startNormal, bool& jumpComingFromSequenceGap) -> double {

		if(exitEdge == 0)
		{
			if(verbose)
				std::cout << "\t" << "THE COORDINATE ORIGIN" << "\n\n" << std::flush;
		}
		else
		{
			if(verbose)
			{
				std::cout << "\t" << "exitEdge->to_x: " << exitEdge->to_x << "\n";
				std::cout << "\t" << "exitEdge->to_y: " << exitEdge->to_y << "\n";
				std::cout << "\t" << "exitEdge->to_z: " << exitEdge->to_z << "\n\n" << std::flush;
			}
		}

		bool debug = false;
		if((void*)entryEdge == (void*)0x817fa5d8)
		{
			if(exitEdge != 0)
			{
				if((exitEdge->to_x == 54) && (exitEdge->to_y == 31))
				{
					debug = true;
				}

			}
		}


		double previousScore = (exitEdge == 0) ? 0 : exitEdge->getScoreAfterEdge();

		bool exitEdge_isAffineSequenceGap = (exitEdge == 0) ? false : exitEdge->isSequenceGap_affine();
		bool exitEdge_isGraphGap = (exitEdge == 0) ? false : exitEdge->isGraphGap();
		// bool entryEdge_isAffineSequenceGap = (entryEdge == 0) ? false : entryEdge->isSequenceGap_affine();
		// bool entryEdge_isGraphGap = (entryEdge == 0) ? false : entryEdge->isGraphGap();

		int distance_along_y;
		int raw_distance_along_x;

		if((entryEdge == 0) && (exitEdge == 0))
		{
			distance_along_y = sequence.length();
			raw_distance_along_x = g->NodesPerLevel.size() - 1;
		}
		else if(entryEdge == 0)
		{
			assert(exitEdge != 0);
			distance_along_y = sequence.length() - exitEdge->to_y;
			raw_distance_along_x = g->NodesPerLevel.size() - 1 - exitEdge->to_x;
		}
		else if(exitEdge == 0)
		{
			assert(entryEdge != 0);
			distance_along_y = entryEdge->from_y;
			raw_distance_along_x = entryEdge->from_x;
		}
		else
		{
			assert(entryEdge != 0);
			assert(exitEdge != 0);
			distance_along_y = entryEdge->from_y - exitEdge->to_y;
			raw_distance_along_x = entryEdge->from_x - exitEdge->to_x;
		}

		assert(distance_along_y >= 0);
		assert(graphDistance_startAffineGap.Score_endAffinely <= 0);
		assert(graphDistance_startAffineGap.Score_endInAnything <= 0);
		assert(graphDistance_startAffineGap.Score_endInAnything >= graphDistance_startAffineGap.Score_endAffinely);
		assert(graphDistance_startAffineGap.Score_endInAnything > minusInfinity);

		assert(graphDistance_startNormal.Score_endAffinely <= 0);
		assert(graphDistance_startNormal.Score_endInAnything <= 0);
		assert(graphDistance_startNormal.Score_endInAnything >= graphDistance_startNormal.Score_endAffinely);
		assert(graphDistance_startNormal.Score_endInAnything > minusInfinity);

		if(raw_distance_along_x == 0)
		{
			assert(graphDistance_startAffineGap.Score_endInAnything == 0);
			assert(graphDistance_startNormal.Score_endInAnything == 0);
		}

		double combinedScore;
		if((raw_distance_along_x != 0) && (distance_along_y != 0))
		{
			// we have two alternatives: either the sequence gap comes first and the graph gap follows, or the other way around
			std::vector<double> scores;
			std::vector<bool> scores_fromAffineSequenceGap;

			{
				// case 1: sequence gap comes first
				{
					// case 1.1: the sequence gap ends in anything

					graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);
					double sequenceGapScore = utilizedDistance.Score_endInAnything;
					double graphGapScore = S_openGap + distance_along_y * S_extendGap;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(false);

				}
				{
					// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
				}
			}
			{
				// case 2: graph gap comes first
				{
					// case 2.1: sequence gap ends in anything
					double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
					double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
					double sequenceGapScore = graphDistance_startNormal.Score_endInAnything;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(false);

				}
				{
					// case 2.2: sequence gap ends affinely
					double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
					double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
					double sequenceGapScore = graphDistance_startNormal.Score_endAffinely;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(true);
				}
			}

			std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
			combinedScore = scoreAlternativesMax.first;
			jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);

		}
		else if((raw_distance_along_x == 0) && (distance_along_y != 0))
		{
			// graph gap
			double gapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
			double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
			combinedScore = gapOpenPenalty + distance_along_y * S_extendGap + scoreFromEntryEdge;
			jumpComingFromSequenceGap = false;

		}
		else if((raw_distance_along_x != 0) && (distance_along_y == 0))
		{
			// sequence gap
			graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);

			std::vector<double> scores;
			std::vector<bool> scores_fromAffineSequenceGap;

			// case 1: we want to end in an affine gap
			{
				double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
				double scoreGraphDistance = utilizedDistance.Score_endAffinely;
				double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
				scores.push_back(local_combinedScore);
				scores_fromAffineSequenceGap.push_back(true);

			}

			// case 2: end in anything
			{
				double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
				double scoreGraphDistance = utilizedDistance.Score_endInAnything;
				double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
				scores.push_back(local_combinedScore);
				scores_fromAffineSequenceGap.push_back(false);

			}

			std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
			combinedScore = scoreAlternativesMax.first;
			jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);

		}
		else
		{
			assert((raw_distance_along_x == 0) && (distance_along_y == 0));
			double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(exitEdge_isAffineSequenceGap, exitEdge_isGraphGap, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
			combinedScore = scoreFromEntryEdge;
			jumpComingFromSequenceGap = (exitEdge != 0) ? exitEdge->isSequenceGap_endsFree() : false;
		}

		if(verbose)
		{
			std::cout << "\t" << "graph jump scores:\n";
			std::cout << "\t\t" << " affine start, affine end: " << graphDistance_startAffineGap.Score_endAffinely << "\n";
			std::cout << "\t\t" << " non-affine start, affine end: " << graphDistance_startNormal.Score_endAffinely << "\n";
			std::cout << "\t\t" << " affine start, arbitrary end: " << graphDistance_startAffineGap.Score_endInAnything << "\n";
			std::cout << "\t\t" << " non-affine start, arbitrary end: " << graphDistance_startNormal.Score_endInAnything << "\n";
			std::cout << "\t" << "distance_along_y: " << distance_along_y << "\n" << std::flush;
		}

		if((void*)entryEdge == (void*)0x817fa5d8)
		{
			if(exitEdge != 0)
			{
				if((exitEdge->to_x == 54) && (exitEdge->to_y == 31))
				{
//					std::cerr << "Entry edge: " << entryEdge << "\n";
//					std::cerr << "Exit edge: " << exitEdge << "\n";
//					std::cerr << "\tto_x: " << exitEdge->to_x << "\n";
//					std::cerr << "\tto_y: " << exitEdge->to_y << "\n";
//					std::cerr << "\tto_z: " << exitEdge->to_z << "\n" << std::flush;
//					std::cerr << "\texitEdge->getScoreAfterEdge(): " << exitEdge->getScoreAfterEdge() << "\n" << std::flush;
//					std::cerr << "This jump score: " << combinedScore << "\n" << std::flush;
				}

			}
//			std::cerr << "\n\n==================\n\n" << (previousScore + combinedScore) << "\n\n" << previousScore << "\n\n" << combinedScore << "\n\n=========================\n\n" << std::flush;
		}

		return (previousScore + combinedScore);
	};

	std::vector<std::pair<int, int> > gaps = findGaps_chainCoverage(sequence, sequencePositions_covered);
	bool haveStartGap = ((gaps.size() > 0) && (gaps.at(0).first == 0));
	bool haveStopGap = ((gaps.size() > 0) && (gaps.at(gaps.size() - 1).second == (sequence.length() - 1)));
	if(! haveStartGap)
	{
		gaps.insert(gaps.begin(), make_pair(0, -1));
	}
	if(! haveStopGap)
	{
		gaps.insert(gaps.end(), make_pair(sequence.length(), sequence.length() - 1));
	}

	for(int seqI = 1; seqI < (int)sequencePositions_covered.size(); seqI++)
	{
		kMerEdgeChain* chain2 = sequencePositions_covered.at(seqI);
		kMerEdgeChain* chain1 = sequencePositions_covered.at(seqI - 1);
		if((chain1 != 0) && (chain2 != 0) && (chain2 != chain1))
		{
			std::pair<int, int> pseudoGap;
			pseudoGap.first = seqI;
			pseudoGap.second = seqI - 1;
			gaps.push_back(pseudoGap);
		}
	}

//	std::cerr << "Process " << gaps.size() << " gaps (or pseudo-gaps).\n" << std::flush;

//	std::cerr << "sequencePositions_covered.at(0): " << sequencePositions_covered.at(0) << "\n" << std::flush;

	for(int gapI = 0; gapI < (int)gaps.size(); gapI++)
	{
		std::pair<int, int>& gapCoordinates = gaps.at(gapI);

		kMerEdgeChain* leftChain = (gapCoordinates.first > 0) ? sequencePositions_covered.at(gapCoordinates.first - 1) : 0;
		kMerEdgeChain* rightChain = (gapCoordinates.second < ((int)sequencePositions_covered.size() - 1)) ? sequencePositions_covered.at(gapCoordinates.second + 1) : 0;

		NWPath* leftPath = (leftChain != 0) ? chains2Paths.at(leftChain) : 0;
		NWPath* rightPath = (rightChain != 0) ? chains2Paths.at(rightChain) : 0;

		Node* leftBoundaryNode = (leftChain != 0) ? leftChain->traversedEdges.back()->To : 0;
		Node* rightBoundaryNode = (rightChain != 0) ? rightChain->traversedEdges.front()->From : 0;


//		std::cerr << " Gap " << gapI << ", from " << gapCoordinates.first << " to " << gapCoordinates.second << "\n" << std::flush;
		// std::cerr << "\trichtChain: " << rightChain << "\n";
//		if((gapCoordinates.first >= 0) and (gapCoordinates.second >= 0))
//		{
//			std::string strCovered(sequence.begin() + gapCoordinates.first, sequence.begin() + gapCoordinates.second + 1);
//			std::cerr << "\t Sequence: " << strCovered << "\n" << std::flush;
//		}

		int gapBoundary_left_level = (leftBoundaryNode != 0) ? leftBoundaryNode->level : 0;
		int gapBoundary_right_level = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);

		int graphDistance_gap = gapBoundary_right_level - gapBoundary_left_level;
		assert(graphDistance_gap >= 0);

		int sequencePosition_left = (leftChain) ? leftChain->sequence_end : 0;
		int sequencePosition_right = (rightChain) ? rightChain->sequence_begin : sequence.length();
		int sequenceDistance_gap = sequencePosition_right - sequencePosition_left - 1;
		assert(sequenceDistance_gap >= -1);

//		std::cout << std::flush;
//		std::cerr << "\t" << "leftChain: " << leftChain << "\n";
//		std::cerr << "\t" << "rightChain: " << rightChain << "\n";
//		std::cerr << "\t" << "leftPath: " << leftPath << "\n";
//		std::cerr << "\t" << "rightPath: " << rightPath << "\n";
//		std::cerr << "\t" << "leftBoundaryNode: " << leftBoundaryNode << "\n";
//		std::cerr << "\t" << "rightBoundaryNode: " << rightBoundaryNode << "\n";
//		std::cerr << "\t" << "gapBoundary_left_level: " << gapBoundary_left_level << "\n";
//		std::cerr << "\t" << "gapBoundary_right_level: " << gapBoundary_right_level << "\n" << std::flush;
//		std::cerr << "\t" << "sequenceDistance_gap: " << sequenceDistance_gap << "\n" << std::flush;
//		std::cerr << "\t" << "graphDistance_gap: " << graphDistance_gap << "\n" << std::flush;

		if((sequenceDistance_gap > 0) && (graphDistance_gap > 0))
		{
			if(leftPath != 0)
			{
				std::set<NWEdge*> existingExitEdges = leftPath->exit_edges;
				assert(existingExitEdges.size() == 1);

				NWEdge* exitEdgeToExtend = *(existingExitEdges.begin());

				vNW.removePath(leftPath);

				assert((sequencePosition_left+1) == exitEdgeToExtend->to_y);

//				std::cerr << "vNW_completeRemainingGaps_and_score:\n\tsequence.length: " << sequence.length() << "\n";
//				std::cerr << "\tsequencePosition_right: " << sequencePosition_right << "\n" << std::flush;

				std::vector<localExtension_pathDescription> forwardExtensions = fullNeedleman_diagonal_extension(
						sequence,
						exitEdgeToExtend->to_y,
						exitEdgeToExtend->to_x,
						exitEdgeToExtend->to_z,
						gapBoundary_right_level,
						sequencePosition_right,
						-11,
						&vNW,
						true,
						false
				);

//				if(rightPath == 0)
//				{
//					std::cout << "Now printing all paths for right extension towards the end!\n" << std::flush;
//					for(unsigned int fI = 0; fI < forwardExtensions.size(); fI++)
//					{
//						std::cout << "Extension #" << fI << "\n" << std::flush;
//						forwardExtensions.at(fI)._printExtension();
//					}
//
//					assert( 1 == 0 );
//				}

				for(unsigned int fwI = 0; fwI < forwardExtensions.size(); fwI++)
				{
					localExtension_pathDescription& this_fwE = forwardExtensions.at(fwI);
					leftPath->takeInExtensionPath(&this_fwE, 1, 0);
				}

				leftPath->recalculateFirstLast();

				vNW.addPath(leftPath);
			}

			if(rightPath != 0)
			{

				std::set<NWEdge*> existingEntryEdges = rightPath->entry_edges;
				assert(existingEntryEdges.size() == 1);

				NWEdge* entryEdgeToExtend = *(existingEntryEdges.begin());

				vNW.removePath(rightPath);

				std::vector<localExtension_pathDescription> backwardExtensions = fullNeedleman_diagonal_extension(
						sequence,
						entryEdgeToExtend->from_y,
						entryEdgeToExtend->from_x,
						entryEdgeToExtend->from_z,
						gapBoundary_left_level,
						((leftChain) ? (sequencePosition_left+1) : 0),
						-11,
						&vNW,
						false,
						false
				);

				for(unsigned int bwI = 0; bwI < backwardExtensions.size(); bwI++)
				{
					localExtension_pathDescription& this_bwE = backwardExtensions.at(bwI);
					rightPath->takeInExtensionPath(&this_bwE, -1, 0);
				}

				rightPath->recalculateFirstLast();

				vNW.addPath(rightPath);
			}

			if(leftPath && rightPath)
			{
				for(std::set<NWEdge*>::iterator entryIt = rightPath->entry_edges.begin(); entryIt != rightPath->entry_edges.end(); entryIt++)
				{
					NWEdge* entryEdge = *entryIt;

					int entry_x = entryEdge->from_x;
					int entry_y = entryEdge->from_y;

					std::set<NWEdge*> exitEdges = leftPath->exit_edges;
					std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths_X;
					std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths_Y;

					for(std::set<NWEdge*>::iterator potentialExitEdgeIt = exitEdges.begin(); potentialExitEdgeIt != exitEdges.end(); potentialExitEdgeIt++)
					{
						NWEdge* potentialExitEdge = *potentialExitEdgeIt;
						if((entryEdge->from_x == potentialExitEdge->to_x) && (potentialExitEdge->to_y <= entry_y))
						{
							if(potentialExitEdge->path != entryEdge->path)
							{
								if((hitPaths_X.count(potentialExitEdge->path) == 0) || (hitPaths_X.at(potentialExitEdge->path).first < potentialExitEdge->to_y))
								{
									hitPaths_X[potentialExitEdge->path] = std::pair<int, NWEdge*>(potentialExitEdge->to_y, potentialExitEdge);
								}
							}
						}

						if((entryEdge->from_y == potentialExitEdge->to_y) && (potentialExitEdge->to_x <= entry_x))
						{
							if(potentialExitEdge->path != entryEdge->path)
							{
								if((hitPaths_Y.count(potentialExitEdge->path) == 0) || (hitPaths_Y.at(potentialExitEdge->path).first < potentialExitEdge->to_x))
								{
									hitPaths_Y[potentialExitEdge->path] = std::pair<int, NWEdge*>(potentialExitEdge->to_x, potentialExitEdge);
								}
							}
						}
					}

					for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths_X.begin(); makeEdgeExitIt != hitPaths_X.end(); makeEdgeExitIt++)
					{
						assert(makeEdgeExitIt->first != entryEdge->path);
						assert(makeEdgeExitIt->second.second != entryEdge);
						makeEdgeExitIt->second.second->makeExitEdge();
					}

					for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths_Y.begin(); makeEdgeExitIt != hitPaths_Y.end(); makeEdgeExitIt++)
					{
						assert(makeEdgeExitIt->first != entryEdge->path);
						assert(makeEdgeExitIt->second.second != entryEdge);
						makeEdgeExitIt->second.second->makeExitEdge();
					}
				}
			}
		}

		std::vector<NWEdge*> entryEdges_vector;
		std::vector<NWEdge*> exitEdges_vector;

		if(leftPath != 0)
		{
			exitEdges_vector.insert(exitEdges_vector.end(), leftPath->exit_edges.begin(), leftPath->exit_edges.end());
		}
		if(rightPath != 0)
		{
			entryEdges_vector.insert(entryEdges_vector.end(), rightPath->entry_edges.begin(), rightPath->entry_edges.end());
		}


		int gapBoundary_left_level_forDistances = (leftBoundaryNode != 0) ? (leftBoundaryNode->level - 1): 0;
		int gapBoundary_right_level_forDistances = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);

		// std::cerr << "Calculate distances from " << gapBoundary_left_level_forDistances << " to " << gapBoundary_right_level_forDistances << "\n" << std::flush;

		affinelyCalculateGraphDistancesForVirtualNW(
				gapBoundary_left_level_forDistances,
				gapBoundary_right_level_forDistances,
				entryEdges_vector,
				exitEdges_vector,
				lastPositionDistances_perZ_startAffineGap,
				lastPositionDistances_perZ_startNormal,
				NWedges_graphDist_startAffineGap,
				NWedges_graphDist_startNormal
		);
	}
	
	
	if(vNW.getEntryEdges().size() == 0)
	{
		assert(vNW.getNumPaths() == 0);
		int maxY = sequence.length();
		int fixX = (g->NodesPerLevel.size() - 1);
		
		NWEdge* entryEdge = 0;
		NWEdge* exitEdge = 0;
		assert(maxY > 1);
		NWPath* p = new NWPath();
		for(unsigned int y = 0; y < maxY; y++)
		{
			int entryExitStatus = 0;
			if(y == 0)
			{
				entryExitStatus = -1;
			}
			if(y == (maxY-1))
			{
				entryExitStatus = 1;
			}

			NWEdge* createdEdge = p->createAndAddEdge(fixX, y, 0, fixX, y+1, 0, 0, 0, entryExitStatus);
			if(entryExitStatus == -1)
			{
				entryEdge = createdEdge;
			}
			if(entryExitStatus == 1)
			{
				exitEdge = createdEdge;
			}
		}
		
		p->recalculateFirstLast();	
		vNW.addPath(p);

		assert(entryEdge && exitEdge);
		
		std::vector<NWEdge*> entryEdges_vector;
		std::vector<NWEdge*> exitEdges_vector;

		entryEdges_vector.push_back(entryEdge);
		exitEdges_vector.push_back(exitEdge);
		
		affinelyCalculateGraphDistancesForVirtualNW(
				0,
				(g->NodesPerLevel.size() - 1),
				entryEdges_vector,
				exitEdges_vector,
				lastPositionDistances_perZ_startAffineGap,
				lastPositionDistances_perZ_startNormal,
				NWedges_graphDist_startAffineGap,
				NWedges_graphDist_startNormal
		);		
	}	


//	std::set<NWPath*> paths = vNW.paths;
//	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
//	{
//		std::cout << "Path " << *pathIt << "\n==============\n";
//		(*pathIt)->_printPath();
//		std::cout << std::flush;
//		std::cerr << std::flush;
//	}

	std::vector<NWEdge*> allEntryEdges = vNW.getEntryEdges();
	std::sort(allEntryEdges.begin(), allEntryEdges.end(), [](NWEdge* i, NWEdge* j){
		if(i->from_y != j->from_y)
		{
			return (i->from_y < j->from_y);
		}
		else
		{
			return (i->from_x < j->from_x);
		}
	});

	for(unsigned int entryPointI = 0; entryPointI < allEntryEdges.size(); entryPointI++)
	{
		NWEdge* entryEdge = allEntryEdges.at(entryPointI);
		if(! NWedges_graphDist_startNormal.count(entryEdge))
		{
			std::cerr << "No distance information for entryEdge " << entryEdge << "\n" << std::flush;
		}
		assert(NWedges_graphDist_startNormal.count(entryEdge) > 0);

		if(verbose)
		{
			assert(entryEdge != 0);
			std::cout << "entryPointI " << entryPointI << " " << entryEdge << "\n" << std::flush;
			std::cout << "\t" << "from_x: " << entryEdge->from_x << "\n";
			std::cout << "\t" << "from_y: " << entryEdge->from_y << "\n";
			std::cout << "\t" << "from_z: " << entryEdge->from_z << "\n\n" << std::flush;
		}
//
//		if((! superquiet) && ((entryPointI % 1000) == 0))
//			std::cout  << "\r" << "entryPointI: " << entryPointI << " / " << allEntryEdges.size() << std::flush;

		std::vector<NWEdge*> scores_backtrack;
		std::vector<double> scores_achieved;
		std::vector<bool> scores_jumpFromAffineSequenceGap;

		for(std::map<NWEdge*, graphPointDistance>::iterator comeFromIt = NWedges_graphDist_startNormal.at(entryEdge).begin(); comeFromIt != NWedges_graphDist_startNormal.at(entryEdge).end(); comeFromIt++)
		{
			NWEdge* exitEdge = comeFromIt->first;

			if(verbose)
				std::cout << "\n\t" << "-------\n\tConsider source " << exitEdge << "\n" << std::flush;

			bool jumpComingFromSequenceGap;
			double jumpScore = scoreEdgeJump(entryEdge, exitEdge, NWedges_graphDist_startAffineGap.at(entryEdge).at(exitEdge), comeFromIt->second, jumpComingFromSequenceGap);
			scores_backtrack.push_back(exitEdge);
			scores_achieved.push_back(jumpScore);
			scores_jumpFromAffineSequenceGap.push_back(jumpComingFromSequenceGap);

		}

		if(scores_achieved.size() > 0)
		{
			std::pair<double, unsigned int> maxOrigin = Utilities::findVectorMaxP_nonCritical(scores_achieved, &(rng_seeds.at(omp_get_thread_num())));
			assert(entryEdge != 0);
	//		 std::cout << "entryPoint" << entryPointI << ": selected score " << maxOrigin.first << "\n" << std::flush;

			entryEdge->takeScore_endsFree_nonCritical(maxOrigin.first, scores_jumpFromAffineSequenceGap.at(maxOrigin.second), scores_backtrack.at(maxOrigin.second), S_match, S_mismatch, S_openGap, S_extendGap, &(rng_seeds.at(omp_get_thread_num())));
		}
		else
		{
			entryEdge->takeScore_endsFree_nonCritical(minusInfinity, true, 0, S_match, S_mismatch, S_openGap, S_extendGap, &(rng_seeds.at(omp_get_thread_num())));
		}
	}


	std::vector<double> finalScore_alternatives;
	std::vector<NWEdge*> finalBacktrack_alternatives;
	std::vector<int> finalBacktrack_Zs;

	std::set<Node*> finalNodes = g->NodesPerLevel.at(g->NodesPerLevel.size() - 1);
	for(std::set<Node*>::iterator nodeIt = finalNodes.begin(); nodeIt != finalNodes.end(); nodeIt++)
	{
		Node* finalNode = *nodeIt;
		int z = nodesPerLevel_ordered_rev.at(g->NodesPerLevel.size() - 1).at(finalNode);
		for(std::map<NWEdge*, graphPointDistance>::iterator comeFromIt = lastPositionDistances_perZ_startNormal.at(z).begin(); comeFromIt != lastPositionDistances_perZ_startNormal.at(z).end(); comeFromIt++)
		{
			NWEdge* exitEdge = comeFromIt->first;

			if(verbose)
				std::cout << "\n\t" << "-------\n\tConsider source " << exitEdge << "\n" << std::flush;

			bool jumpComingFromSequenceGap;
			double jumpScore = scoreEdgeJump(0, exitEdge, lastPositionDistances_perZ_startAffineGap.at(z).at(exitEdge), comeFromIt->second, jumpComingFromSequenceGap);

			finalScore_alternatives.push_back(jumpScore);
			finalBacktrack_alternatives.push_back(exitEdge);
			finalBacktrack_Zs.push_back(z);

			if(exitEdge != 0)
			{
//				std::cout << "Final z alternatives " << z << " / " << exitEdge  << " coming from " << exitEdge->to_x << ", " << exitEdge->to_y << ", " << exitEdge->to_z << ", distance " << comeFromIt->second.Score_endInAnything << ": " << jumpScore << "\n";
			}
			else
			{
//				std::cout << "Final z alternatives " << z << " / " << "ORIGIN"  << " coming from " << 0 << ", " << 0 << ", " << "?" << ", distance " << comeFromIt->second.Score_endInAnything << ": " << jumpScore << "\n";
			}
		}
	}

	std::pair<double, unsigned int> maxFinalScore = Utilities::findVectorMaxP_nonCritical(finalScore_alternatives, &(rng_seeds.at(omp_get_thread_num())));
	finalScore = maxFinalScore.first;
	finalScore_z = finalBacktrack_Zs.at(maxFinalScore.second);
	finalScore_backtrack = finalBacktrack_alternatives.at(maxFinalScore.second);
}


void GraphAlignerUnique::vNW_completeRemainingGaps_and_score_local(std::string& sequence, VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, NWPath*>& chains2Paths, std::map<kMerEdgeChain*, int>& currentChains_start, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, double& finalScore, int& finalScore_z, NWEdge*& finalScore_backtrack)
{
	bool verbose = false;
	bool superquiet = true;

	double minusInfinity = -1 * numeric_limits<double>::max();
	double Infinity = numeric_limits<double>::max();

	auto scoreEdgeJump = [&](NWEdge* entryEdge, NWEdge* exitEdge, graphPointDistance& graphDistance_startAffineGap, graphPointDistance& graphDistance_startNormal, bool& jumpComingFromSequenceGap) -> double {

		if(exitEdge == 0)
		{
			if(verbose)
				std::cout << "\t" << "THE COORDINATE ORIGIN" << "\n\n" << std::flush;
		}
		else
		{
			if(verbose)
			{
				std::cout << "\t" << "exitEdge->to_x: " << exitEdge->to_x << "\n";
				std::cout << "\t" << "exitEdge->to_y: " << exitEdge->to_y << "\n";
				std::cout << "\t" << "exitEdge->to_z: " << exitEdge->to_z << "\n\n" << std::flush;
			}
		}

		bool debug = false;
		if((void*)entryEdge == (void*)0x817fa5d8)
		{
			if(exitEdge != 0)
			{
				if((exitEdge->to_x == 54) && (exitEdge->to_y == 31))
				{
					debug = true;
				}

			}
		}


		double previousScore = (exitEdge == 0) ? 0 : exitEdge->getScoreAfterEdge();

		bool exitEdge_isAffineSequenceGap = (exitEdge == 0) ? false : exitEdge->isSequenceGap_affine();
		bool exitEdge_isGraphGap = (exitEdge == 0) ? false : exitEdge->isGraphGap();
		// bool entryEdge_isAffineSequenceGap = (entryEdge == 0) ? false : entryEdge->isSequenceGap_affine();
		// bool entryEdge_isGraphGap = (entryEdge == 0) ? false : entryEdge->isGraphGap();

		int distance_along_y;
		int raw_distance_along_x;

		if((entryEdge == 0) && (exitEdge == 0))
		{
			distance_along_y = sequence.length();
			raw_distance_along_x = g->NodesPerLevel.size() - 1;
		}
		else if(entryEdge == 0)
		{
			assert(exitEdge != 0);
			distance_along_y = sequence.length() - exitEdge->to_y;
			raw_distance_along_x = g->NodesPerLevel.size() - 1 - exitEdge->to_x;
		}
		else if(exitEdge == 0)
		{
			assert(entryEdge != 0);
			distance_along_y = entryEdge->from_y;
			raw_distance_along_x = entryEdge->from_x;
		}
		else
		{
			assert(entryEdge != 0);
			assert(exitEdge != 0);
			distance_along_y = entryEdge->from_y - exitEdge->to_y;
			raw_distance_along_x = entryEdge->from_x - exitEdge->to_x;
		}

		assert(distance_along_y >= 0);
		assert(graphDistance_startAffineGap.Score_endAffinely <= 0);
		if(!(graphDistance_startAffineGap.Score_endInAnything <= 0))
		{
			std::cerr << "!(graphDistance_startAffineGap.Score_endInAnything <= 0)" << "\n";
			std::cerr << graphDistance_startAffineGap.Score_endInAnything << "\n" << std::flush;
		}
		assert(graphDistance_startAffineGap.Score_endInAnything <= 0);
		assert(graphDistance_startAffineGap.Score_endInAnything >= graphDistance_startAffineGap.Score_endAffinely);
		assert(graphDistance_startAffineGap.Score_endInAnything >= minusInfinity);

		assert(graphDistance_startNormal.Score_endAffinely <= 0);
		assert(graphDistance_startNormal.Score_endInAnything <= 0);
		assert(graphDistance_startNormal.Score_endInAnything >= graphDistance_startNormal.Score_endAffinely);
		assert(graphDistance_startNormal.Score_endInAnything >= minusInfinity);

		if(raw_distance_along_x == 0)
		{
			assert((graphDistance_startAffineGap.Score_endInAnything == 0) || (graphDistance_startAffineGap.Score_endInAnything == minusInfinity));
			assert(graphDistance_startNormal.Score_endInAnything == 0);
		}

		double combinedScore;
		if((raw_distance_along_x != 0) && (distance_along_y != 0))
		{
			// we have two alternatives: either the sequence gap comes first and the graph gap follows, or the other way around
			std::vector<double> scores;
			std::vector<bool> scores_fromAffineSequenceGap;

			{
				// case 1: sequence gap comes first
				{
					// case 1.1: the sequence gap ends in anything

					graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);
					double sequenceGapScore = utilizedDistance.Score_endInAnything;
					double graphGapScore = S_openGap + distance_along_y * S_extendGap;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(false);

				}
				{
					// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
				}
			}
			{
				// case 2: graph gap comes first
				{
					// case 2.1: sequence gap ends in anything
					double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
					double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
					double sequenceGapScore = graphDistance_startNormal.Score_endInAnything;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(false);

				}
				{
					// case 2.2: sequence gap ends affinely
					double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
					double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
					double sequenceGapScore = graphDistance_startNormal.Score_endAffinely;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
					scores.push_back(local_combinedScore);
					scores_fromAffineSequenceGap.push_back(true);
				}
			}

			std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
			combinedScore = scoreAlternativesMax.first;
			jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);

		}
		else if((raw_distance_along_x == 0) && (distance_along_y != 0))
		{
			// graph gap
			double gapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
			double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
			combinedScore = gapOpenPenalty + distance_along_y * S_extendGap + scoreFromEntryEdge;
			jumpComingFromSequenceGap = false;

		}
		else if((raw_distance_along_x != 0) && (distance_along_y == 0))
		{
			// sequence gap
			graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);

			std::vector<double> scores;
			std::vector<bool> scores_fromAffineSequenceGap;

			// case 1: we want to end in an affine gap
			{
				double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
				double scoreGraphDistance = utilizedDistance.Score_endAffinely;
				double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
				scores.push_back(local_combinedScore);
				scores_fromAffineSequenceGap.push_back(true);

			}

			// case 2: end in anything
			{
				double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
				double scoreGraphDistance = utilizedDistance.Score_endInAnything;
				double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
				scores.push_back(local_combinedScore);
				scores_fromAffineSequenceGap.push_back(false);

			}

			std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
			combinedScore = scoreAlternativesMax.first;
			jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);

		}
		else
		{
			assert((raw_distance_along_x == 0) && (distance_along_y == 0));
			double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(exitEdge_isAffineSequenceGap, exitEdge_isGraphGap, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
			combinedScore = scoreFromEntryEdge;
			jumpComingFromSequenceGap = (exitEdge != 0) ? exitEdge->isSequenceGap_endsFree() : false;
		}

		if(verbose)
		{
			std::cout << "\t" << "graph jump scores:\n";
			std::cout << "\t\t" << " affine start, affine end: " << graphDistance_startAffineGap.Score_endAffinely << "\n";
			std::cout << "\t\t" << " non-affine start, affine end: " << graphDistance_startNormal.Score_endAffinely << "\n";
			std::cout << "\t\t" << " affine start, arbitrary end: " << graphDistance_startAffineGap.Score_endInAnything << "\n";
			std::cout << "\t\t" << " non-affine start, arbitrary end: " << graphDistance_startNormal.Score_endInAnything << "\n";
			std::cout << "\t" << "distance_along_y: " << distance_along_y << "\n" << std::flush;
		}

		if((void*)entryEdge == (void*)0x817fa5d8)
		{
			if(exitEdge != 0)
			{
				if((exitEdge->to_x == 54) && (exitEdge->to_y == 31))
				{
//					std::cerr << "Entry edge: " << entryEdge << "\n";
//					std::cerr << "Exit edge: " << exitEdge << "\n";
//					std::cerr << "\tto_x: " << exitEdge->to_x << "\n";
//					std::cerr << "\tto_y: " << exitEdge->to_y << "\n";
//					std::cerr << "\tto_z: " << exitEdge->to_z << "\n" << std::flush;
//					std::cerr << "\texitEdge->getScoreAfterEdge(): " << exitEdge->getScoreAfterEdge() << "\n" << std::flush;
//					std::cerr << "This jump score: " << combinedScore << "\n" << std::flush;
				}

			}
//			std::cerr << "\n\n==================\n\n" << (previousScore + combinedScore) << "\n\n" << previousScore << "\n\n" << combinedScore << "\n\n=========================\n\n" << std::flush;
		}

		return (previousScore + combinedScore);
	};

	std::vector<std::pair<int, int> > gaps = findGaps_chainCoverage(sequence, sequencePositions_covered);
	bool haveStartGap = ((gaps.size() > 0) && (gaps.at(0).first == 0));
	bool haveStopGap = ((gaps.size() > 0) && (gaps.at(gaps.size() - 1).second == (sequence.length() - 1)));

	// todo deal with this later

//	assert(haveStartGap);
//	assert(haveStopGap);

//	if(haveStopGap)
//	{
//		gaps.erase(gaps.end() - 1);
//	}
//	if(haveStartGap)
//	{
//		gaps.erase(gaps.begin());
//	}

	if(! haveStartGap)
	{
		gaps.insert(gaps.begin(), make_pair(0, -1));
	}
	if(! haveStopGap)
	{
		gaps.insert(gaps.end(), make_pair(sequence.length(), sequence.length() - 1));
	}

	for(int seqI = 1; seqI < (int)sequencePositions_covered.size(); seqI++)
	{
		kMerEdgeChain* chain2 = sequencePositions_covered.at(seqI);
		kMerEdgeChain* chain1 = sequencePositions_covered.at(seqI - 1);
		if((chain1 != 0) && (chain2 != 0) && (chain2 != chain1))
		{
			std::pair<int, int> pseudoGap;
			pseudoGap.first = seqI;
			pseudoGap.second = seqI - 1;
			gaps.push_back(pseudoGap);
		}
	}


	if(verbose)
	{
		std::cout << "vNW_completeRemainingGaps_and_score_local(..): Process " << gaps.size() << " gaps (or pseudo-gaps).\n" << std::flush;
	}
	
//	std::cerr << "sequencePositions_covered.at(0): " << sequencePositions_covered.at(0) << "\n" << std::flush;

	for(int gapI = 0; gapI < (int)gaps.size(); gapI++)
	{
		std::pair<int, int>& gapCoordinates = gaps.at(gapI);

		kMerEdgeChain* leftChain = (gapCoordinates.first > 0) ? sequencePositions_covered.at(gapCoordinates.first - 1) : 0;
		kMerEdgeChain* rightChain = (gapCoordinates.second < ((int)sequencePositions_covered.size() - 1)) ? sequencePositions_covered.at(gapCoordinates.second + 1) : 0;

		NWPath* leftPath = (leftChain != 0) ? chains2Paths.at(leftChain) : 0;
		NWPath* rightPath = (rightChain != 0) ? chains2Paths.at(rightChain) : 0;

		Node* leftBoundaryNode = (leftChain != 0) ? leftChain->traversedEdges.back()->To : 0;
		Node* rightBoundaryNode = (rightChain != 0) ? rightChain->traversedEdges.front()->From : 0;

		int gapBoundary_left_level = (leftBoundaryNode != 0) ? leftBoundaryNode->level : 0;
		int gapBoundary_right_level = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);

		int graphDistance_gap = gapBoundary_right_level - gapBoundary_left_level;
		assert(graphDistance_gap >= 0);

		int sequencePosition_left = (leftChain) ? leftChain->sequence_end : 0;
		int sequencePosition_right = (rightChain) ? rightChain->sequence_begin : sequence.length();   
		int sequenceDistance_gap = sequencePosition_right - sequencePosition_left - 1;
		assert(sequenceDistance_gap >= -1);

		bool isLeftBeginning = (leftChain == 0);
		bool isRightEnd = (rightChain == 0);

		bool extendedLeftPathToZero = false;
		bool extendedRightPathToZero = false;

		if(verbose)
		{
			std::cout << std::flush;
			std::cout << " Gap " << gapI << ", from " << gapCoordinates.first << " to " << gapCoordinates.second << "\n" << std::flush;
			std::cout << "\t" << "leftChain: " << leftChain << "\n";
			std::cout << "\t" << "rightChain: " << rightChain << "\n";
			std::cout << "\t" << "leftPath: " << leftPath << "\n";
			std::cout << "\t" << "rightPath: " << rightPath << "\n";
			std::cout << "\t" << "leftBoundaryNode: " << leftBoundaryNode << "\n";
			std::cout << "\t" << "rightBoundaryNode: " << rightBoundaryNode << "\n";
			std::cout << "\t" << "gapBoundary_left_level: " << gapBoundary_left_level << "\n";
			std::cout << "\t" << "gapBoundary_right_level: " << gapBoundary_right_level << "\n" << std::flush;
			std::cout << "\t" << "isLeftBeginning: " << isLeftBeginning << "\n";
			std::cout << "\t" << "isRightEnd: " << isRightEnd << "\n\n" << std::flush;
			std::cout << "\t" << "sequencePosition_left: " << sequencePosition_left << "\n";
			std::cout << "\t" << "sequencePosition_right: " << sequencePosition_right << "\n";
			std::cout << "\t" << "sequenceDistance_gap: " << sequenceDistance_gap << "\n";
			std::cout << "\t" << "graphDistance_gap: " << graphDistance_gap << "\n"<< std::flush;
		}

		// I don't know exactly what this -1 in sequenceDistance_gap is supposed to achieve. The following
		// condition has a special case for !leftChain and start position of chain == 1.
		// We might have to do the same for the graph.
		if(((sequenceDistance_gap > 0) || (isLeftBeginning && (sequencePosition_right == 1))) && (graphDistance_gap > 0))
		{
			if(leftPath != 0)
			{
				std::set<NWEdge*> existingExitEdges = leftPath->exit_edges;
				assert(existingExitEdges.size() == 1);

				NWEdge* exitEdgeToExtend = *(existingExitEdges.begin());

				vNW.removePath(leftPath);

				assert((sequencePosition_left+1) == exitEdgeToExtend->to_y);

				if(verbose)
				{
					std::cout << "\t" << "fullNeedleman_diagonal_extension from left end (inwards)\n"; 
				}
				
				std::vector<localExtension_pathDescription> forwardExtensions = fullNeedleman_diagonal_extension(
						sequence,
						exitEdgeToExtend->to_y,
						exitEdgeToExtend->to_x,
						exitEdgeToExtend->to_z,
						gapBoundary_right_level,
						sequencePosition_right,
						-11,
						&vNW,
						true,
						false
				);

				if(forwardExtensions.size() > 0)
				{
					for(unsigned int fwI = 0; fwI < forwardExtensions.size(); fwI++)
					{
						localExtension_pathDescription& this_fwE = forwardExtensions.at(fwI);
						int extensionsUntilBorder = 0; // (isRightEnd) ? 1 : 0;
						leftPath->takeInExtensionPath(&this_fwE, 1, extensionsUntilBorder);
						// extendedLeftPathToZero = true;
					}
				}
				else
				{
//					if(isRightEnd)
//					{
//						leftPath->extendToZero(1);
//						extendedLeftPathToZero = true;
//					}
				}

				leftPath->recalculateFirstLast();

				vNW.addPath(leftPath);
			}

			if(rightPath != 0)
			{
				std::set<NWEdge*> existingEntryEdges = rightPath->entry_edges;
				assert(existingEntryEdges.size() == 1);

				NWEdge* entryEdgeToExtend = *(existingEntryEdges.begin());

				vNW.removePath(rightPath);

				if(verbose)
				{
					std::cout << "\t" << "fullNeedleman_diagonal_extension from right end (inwards)\n"; 
				}
				
				
				std::vector<localExtension_pathDescription> backwardExtensions = fullNeedleman_diagonal_extension(
						sequence,
						entryEdgeToExtend->from_y,
						entryEdgeToExtend->from_x,
						entryEdgeToExtend->from_z,
						gapBoundary_left_level,
						((leftChain) ? (sequencePosition_left+1) : 0),
						-11,
						&vNW,
						false,
						false
				);

				if(backwardExtensions.size() > 0)
				{
					for(unsigned int bwI = 0; bwI < backwardExtensions.size(); bwI++)
					{
						localExtension_pathDescription& this_bwE = backwardExtensions.at(bwI);
						int extensionsUntilBorder = 0; //(isLeftBeginning) ? -1 : 0;
						rightPath->takeInExtensionPath(&this_bwE, -1, extensionsUntilBorder);
						// extendedRightPathToZero = true;
					}
				}
				else
				{
//					if(isLeftBeginning)
//					{
//						rightPath->extendToZero(-1);
//						extendedRightPathToZero = true;
//					}
				}

				rightPath->recalculateFirstLast();

				vNW.addPath(rightPath);
			}

			if(leftPath && rightPath)
			{
				assert(! isLeftBeginning);
				assert(! isRightEnd);

				for(std::set<NWEdge*>::iterator entryIt = rightPath->entry_edges.begin(); entryIt != rightPath->entry_edges.end(); entryIt++)
				{
					NWEdge* entryEdge = *entryIt;

					int entry_x = entryEdge->from_x;
					int entry_y = entryEdge->from_y;

					std::set<NWEdge*> exitEdges = leftPath->exit_edges;
					std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths_X;
					std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths_Y;

					for(std::set<NWEdge*>::iterator potentialExitEdgeIt = exitEdges.begin(); potentialExitEdgeIt != exitEdges.end(); potentialExitEdgeIt++)
					{
						NWEdge* potentialExitEdge = *potentialExitEdgeIt;
						if((entryEdge->from_x == potentialExitEdge->to_x) && (potentialExitEdge->to_y <= entry_y))
						{
							if(potentialExitEdge->path != entryEdge->path)
							{
								if((hitPaths_X.count(potentialExitEdge->path) == 0) || (hitPaths_X.at(potentialExitEdge->path).first < potentialExitEdge->to_y))
								{
									hitPaths_X[potentialExitEdge->path] = std::pair<int, NWEdge*>(potentialExitEdge->to_y, potentialExitEdge);
								}
							}
						}

						if((entryEdge->from_y == potentialExitEdge->to_y) && (potentialExitEdge->to_x <= entry_x))
						{
							if(potentialExitEdge->path != entryEdge->path)
							{
								if((hitPaths_Y.count(potentialExitEdge->path) == 0) || (hitPaths_Y.at(potentialExitEdge->path).first < potentialExitEdge->to_x))
								{
									hitPaths_Y[potentialExitEdge->path] = std::pair<int, NWEdge*>(potentialExitEdge->to_x, potentialExitEdge);
								}
							}
						}
					}

					for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths_X.begin(); makeEdgeExitIt != hitPaths_X.end(); makeEdgeExitIt++)
					{
						assert(makeEdgeExitIt->first != entryEdge->path);
						assert(makeEdgeExitIt->second.second != entryEdge);
						makeEdgeExitIt->second.second->makeExitEdge();
					}

					for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths_Y.begin(); makeEdgeExitIt != hitPaths_Y.end(); makeEdgeExitIt++)
					{
						assert(makeEdgeExitIt->first != entryEdge->path);
						assert(makeEdgeExitIt->second.second != entryEdge);
						makeEdgeExitIt->second.second->makeExitEdge();
					}
				}
			}
		}

		if(isLeftBeginning && (! extendedRightPathToZero) && (rightPath))
		{
			vNW.removePath(rightPath);
			rightPath->extendToZero(-1);
			extendedRightPathToZero = true;
			rightPath->recalculateFirstLast();
			vNW.addPath(rightPath);
		}

		if(isRightEnd && (! extendedLeftPathToZero) && (leftPath))
		{
//			std::cerr << "A\n" << std::flush;
			vNW.removePath(leftPath);
//			std::cerr << "B\n" << std::flush;
			assert(leftPath != 0);
			leftPath->extendToZero(1, sequence.length());
//			std::cerr << "C\n" << std::flush;
			extendedLeftPathToZero = true;
//			std::cerr << "D\n" << std::flush;
			leftPath->recalculateFirstLast();
			vNW.addPath(leftPath);
//			std::cerr << "E\n" << std::flush;
		}

		std::vector<NWEdge*> entryEdges_vector;
		std::vector<NWEdge*> exitEdges_vector;

		if(leftPath != 0)
		{
			exitEdges_vector.insert(exitEdges_vector.end(), leftPath->exit_edges.begin(), leftPath->exit_edges.end());
		}
		if(rightPath != 0)
		{
			entryEdges_vector.insert(entryEdges_vector.end(), rightPath->entry_edges.begin(), rightPath->entry_edges.end());
		}


		int gapBoundary_left_level_forDistances = (leftBoundaryNode != 0) ? (leftBoundaryNode->level - 1): 0;
		int gapBoundary_right_level_forDistances = (rightBoundaryNode != 0) ? rightBoundaryNode->level : (g->NodesPerLevel.size() - 1);

		// std::cerr << "Calculate distances from " << gapBoundary_left_level_forDistances << " to " << gapBoundary_right_level_forDistances << "\n" << std::flush;

		// assert(!(isLeftBeginning && isRightEnd));

		if(isLeftBeginning && isRightEnd)
		{
			lastPositionDistances_perZ_startNormal.resize(nodesPerLevel_ordered_rev.back().size());
			lastPositionDistances_perZ_startAffineGap.resize(nodesPerLevel_ordered_rev.back().size());

			for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
			{
				NWEdge* exitEdge = exitEdges_vector.at(eI);
				assert(exitEdge->to_y == sequence.length());

				graphPointDistance startNormal;
				startNormal.Score_endAffinely = 0;
				startNormal.Score_endInAnything = 0;
				startNormal.start_in_affine_sequenceGap = false;

				graphPointDistance startAffinely;
				startAffinely.Score_endAffinely = minusInfinity;
				startAffinely.Score_endInAnything = minusInfinity;
				startAffinely.start_in_affine_sequenceGap = true;

				// NWedges_graphDist_startNormal[exitEdge][0] = startNormal;
				// NWedges_graphDist_startAffineGap[exitEdge][0] = startAffinely;

				for(unsigned int zI = 0; zI < nodesPerLevel_ordered_rev.back().size(); zI++)
				{
					lastPositionDistances_perZ_startNormal.at(zI)[exitEdge] = startNormal;
					lastPositionDistances_perZ_startAffineGap.at(zI)[exitEdge] = startAffinely;
				}

			}		
		}
		else if(isLeftBeginning)
		{
			for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
			{
				NWEdge* entryEdge = entryEdges_vector.at(eI);
				if(!(entryEdge->from_y == 0))
				{
					std::cerr << "Entry edge has from_y != 0, but should have! Value: " << entryEdge->from_y << "\n" << std::flush;
 				}
				assert(entryEdge->from_y == 0);

				graphPointDistance startNormal;
				startNormal.Score_endAffinely = minusInfinity;
				startNormal.Score_endInAnything = 0;
				startNormal.start_in_affine_sequenceGap = false;

				graphPointDistance startAffinely;
				startAffinely.Score_endAffinely = minusInfinity;
				startAffinely.Score_endInAnything = minusInfinity;
				startAffinely.start_in_affine_sequenceGap = true;

				NWedges_graphDist_startNormal[entryEdge][0] = startNormal;
				NWedges_graphDist_startAffineGap[entryEdge][0] = startAffinely;
			}
//			std::cerr << "Passed left beginning tests!\n" << std::flush;
		}
		else if(isRightEnd)
		{
			lastPositionDistances_perZ_startNormal.resize(nodesPerLevel_ordered_rev.back().size());
			lastPositionDistances_perZ_startAffineGap.resize(nodesPerLevel_ordered_rev.back().size());

			for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
			{
				NWEdge* exitEdge = exitEdges_vector.at(eI);
				assert(exitEdge->to_y == sequence.length());

				graphPointDistance startNormal;
				startNormal.Score_endAffinely = 0;
				startNormal.Score_endInAnything = 0;
				startNormal.start_in_affine_sequenceGap = false;

				graphPointDistance startAffinely;
				startAffinely.Score_endAffinely = minusInfinity;
				startAffinely.Score_endInAnything = minusInfinity;
				startAffinely.start_in_affine_sequenceGap = true;

				// NWedges_graphDist_startNormal[exitEdge][0] = startNormal;
				// NWedges_graphDist_startAffineGap[exitEdge][0] = startAffinely;

				for(unsigned int zI = 0; zI < nodesPerLevel_ordered_rev.back().size(); zI++)
				{
					lastPositionDistances_perZ_startNormal.at(zI)[exitEdge] = startNormal;
					lastPositionDistances_perZ_startAffineGap.at(zI)[exitEdge] = startAffinely;
				}

			}
//			std::cerr << "Passed right end tests!\n" << std::flush;
		}
		else
		{
			affinelyCalculateGraphDistancesForVirtualNW_local(
					gapBoundary_left_level_forDistances,
					gapBoundary_right_level_forDistances,
					entryEdges_vector,
					exitEdges_vector,
					lastPositionDistances_perZ_startAffineGap,
					lastPositionDistances_perZ_startNormal,
					NWedges_graphDist_startAffineGap,
					NWedges_graphDist_startNormal
			);
		}
	}


//	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
//	{
//		std::cerr << "Path " << *pathIt << "\n==============\n";
//		(*pathIt)->_printPath();
//		std::cerr << std::flush;
//		std::cerr << std::flush;
//	}
	
	if(vNW.getEntryEdges().size() == 0)
	{
		assert(vNW.getNumPaths() == 0);
		int maxY = sequence.length();
		int fixX = (g->NodesPerLevel.size() - 1);
		
		NWEdge* entryEdge = 0;
		NWEdge* exitEdge = 0;
		assert(maxY > 1);
		NWPath* p = new NWPath();
		for(unsigned int y = 0; y < maxY; y++)
		{
			int entryExitStatus = 0;
			if(y == 0)
			{
				entryExitStatus = -1;
			}
			if(y == (maxY-1))
			{
				entryExitStatus = 1;
			}

			NWEdge* createdEdge = p->createAndAddEdge(fixX, y, 0, fixX, y+1, 0, 0, 0, entryExitStatus);
			if(entryExitStatus == -1)
			{
				entryEdge = createdEdge;
			}
			if(entryExitStatus == 1)
			{
				exitEdge = createdEdge;
			}
		}
		
		p->recalculateFirstLast();	
		vNW.addPath(p);


		assert(entryEdge && exitEdge);
			
		// add entry edge distances
		
		if(!(entryEdge->from_y == 0))
		{
			std::cerr << "Entry edge has from_y != 0, but should have! Value: " << entryEdge->from_y << "\n" << std::flush;
		}
		assert(entryEdge->from_y == 0);

		graphPointDistance startNormal_1;
		startNormal_1.Score_endAffinely = minusInfinity;
		startNormal_1.Score_endInAnything = 0;
		startNormal_1.start_in_affine_sequenceGap = false;

		graphPointDistance startAffinely_1;
		startAffinely_1.Score_endAffinely = minusInfinity;
		startAffinely_1.Score_endInAnything = minusInfinity;
		startAffinely_1.start_in_affine_sequenceGap = true;

		NWedges_graphDist_startNormal[entryEdge][0] = startNormal_1;
		NWedges_graphDist_startAffineGap[entryEdge][0] = startAffinely_1;

		// add exit edge distances
		assert(exitEdge->to_y == sequence.length());

		graphPointDistance startNormal_2;
		startNormal_2.Score_endAffinely = 0;
		startNormal_2.Score_endInAnything = 0;
		startNormal_2.start_in_affine_sequenceGap = false;

		graphPointDistance startAffinely_2;
		startAffinely_2.Score_endAffinely = minusInfinity;
		startAffinely_2.Score_endInAnything = minusInfinity;
		startAffinely_2.start_in_affine_sequenceGap = true;

		for(unsigned int zI = 0; zI < nodesPerLevel_ordered_rev.back().size(); zI++)
		{
			lastPositionDistances_perZ_startNormal.at(zI)[exitEdge] = startNormal_2;
			lastPositionDistances_perZ_startAffineGap.at(zI)[exitEdge] = startAffinely_2;
		}		
	}
	
	std::set<NWPath*> paths = vNW.paths;
	
	
	std::vector<NWEdge*> allEntryEdges = vNW.getEntryEdges();
	std::sort(allEntryEdges.begin(), allEntryEdges.end(), [](NWEdge* i, NWEdge* j){
		if(i->from_y != j->from_y)
		{
			return (i->from_y < j->from_y);
		}
		else
		{
			return (i->from_x < j->from_x);
		}
	});

	assert(allEntryEdges.size() > 0); // I only think that this is correct!
	for(unsigned int entryPointI = 0; entryPointI < allEntryEdges.size(); entryPointI++)
	{
		NWEdge* entryEdge = allEntryEdges.at(entryPointI);
		if(! NWedges_graphDist_startNormal.count(entryEdge))
		{
			std::cerr << "No distance information for entryEdge " << entryEdge << "\n" << std::flush;
		}
		assert(NWedges_graphDist_startNormal.count(entryEdge) > 0);

		if(verbose)
		{
			assert(entryEdge != 0);
			std::cout << "entryPointI " << entryPointI << " " << entryEdge << "\n" << std::flush;
			std::cout << "\t" << "from_x: " << entryEdge->from_x << "\n";
			std::cout << "\t" << "from_y: " << entryEdge->from_y << "\n";
			std::cout << "\t" << "from_z: " << entryEdge->from_z << "\n\n" << std::flush;
		}
//
//		if((! superquiet) && ((entryPointI % 1000) == 0))
//			std::cout  << "\r" << "entryPointI: " << entryPointI << " / " << allEntryEdges.size() << std::flush;

		std::vector<NWEdge*> scores_backtrack;
		std::vector<double> scores_achieved;
		std::vector<bool> scores_jumpFromAffineSequenceGap;

		for(std::map<NWEdge*, graphPointDistance>::iterator comeFromIt = NWedges_graphDist_startNormal.at(entryEdge).begin(); comeFromIt != NWedges_graphDist_startNormal.at(entryEdge).end(); comeFromIt++)
		{
			NWEdge* exitEdge = comeFromIt->first;

			if(verbose)
				std::cout << "\n\t" << "-------\n\tConsider source " << exitEdge << "\n" << std::flush;

			bool jumpComingFromSequenceGap;
			double jumpScore = scoreEdgeJump(entryEdge, exitEdge, NWedges_graphDist_startAffineGap.at(entryEdge).at(exitEdge), comeFromIt->second, jumpComingFromSequenceGap);
			scores_backtrack.push_back(exitEdge);
			scores_achieved.push_back(jumpScore);
			scores_jumpFromAffineSequenceGap.push_back(jumpComingFromSequenceGap);

		}

		if(scores_achieved.size() > 0)
		{
			std::pair<double, unsigned int> maxOrigin = Utilities::findVectorMaxP_nonCritical(scores_achieved, &(rng_seeds.at(omp_get_thread_num())));
			assert(entryEdge != 0);
	//		 std::cout << "entryPoint" << entryPointI << ": selected score " << maxOrigin.first << "\n" << std::flush;

			entryEdge->takeScore_endsFree_nonCritical(maxOrigin.first, scores_jumpFromAffineSequenceGap.at(maxOrigin.second), scores_backtrack.at(maxOrigin.second), S_match, S_mismatch, S_openGap, S_extendGap, &(rng_seeds.at(omp_get_thread_num())));
		}
		else
		{
			entryEdge->takeScore_endsFree_nonCritical(minusInfinity, true, 0, S_match, S_mismatch, S_openGap, S_extendGap, &(rng_seeds.at(omp_get_thread_num())));
		}
	}


	std::vector<double> finalScore_alternatives;
	std::vector<NWEdge*> finalBacktrack_alternatives;
	std::vector<int> finalBacktrack_Zs;

	std::set<Node*> finalNodes = g->NodesPerLevel.at(g->NodesPerLevel.size() - 1);
	for(std::set<Node*>::iterator nodeIt = finalNodes.begin(); nodeIt != finalNodes.end(); nodeIt++)
	{
		Node* finalNode = *nodeIt;
		int z = nodesPerLevel_ordered_rev.at(g->NodesPerLevel.size() - 1).at(finalNode);
		for(std::map<NWEdge*, graphPointDistance>::iterator comeFromIt = lastPositionDistances_perZ_startNormal.at(z).begin(); comeFromIt != lastPositionDistances_perZ_startNormal.at(z).end(); comeFromIt++)
		{
			NWEdge* exitEdge = comeFromIt->first;

			if(verbose)
				std::cout << "\n\t" << "-------\n\tConsider source " << exitEdge << "\n" << std::flush;

			bool jumpComingFromSequenceGap;
			double jumpScore = scoreEdgeJump(0, exitEdge, lastPositionDistances_perZ_startAffineGap.at(z).at(exitEdge), comeFromIt->second, jumpComingFromSequenceGap);

			finalScore_alternatives.push_back(jumpScore);
			finalBacktrack_alternatives.push_back(exitEdge);
			finalBacktrack_Zs.push_back(z);

			if(exitEdge != 0)
			{
//				std::cerr << "Final z alternatives " << z << " / " << exitEdge  << " coming from " << exitEdge->to_x << ", " << exitEdge->to_y << ", " << exitEdge->to_z << ", distance " << comeFromIt->second.Score_endInAnything << ": " << jumpScore << "\n" << std::flush;
			}
			else
			{
//				std::cerr << "Final z alternatives " << z << " / " << "ORIGIN"  << " coming from " << 0 << ", " << 0 << ", " << "?" << ", distance " << comeFromIt->second.Score_endInAnything << ": " << jumpScore << "\n" << std::flush;
			}
		}
	}

	std::pair<double, unsigned int> maxFinalScore = Utilities::findVectorMaxP_nonCritical(finalScore_alternatives, &(rng_seeds.at(omp_get_thread_num())));
	finalScore = maxFinalScore.first;
	finalScore_z = finalBacktrack_Zs.at(maxFinalScore.second);
	finalScore_backtrack = finalBacktrack_alternatives.at(maxFinalScore.second);
}


void GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW(int fromLevel, int toLevel, std::vector<NWEdge*>& entryEdges_vector, std::vector<NWEdge*>& exitEdges_vector, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal)
{
	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();

	bool verbose = false;
	bool debug = false;

	if(fromLevel != 0)
	{
		assert(exitEdges_vector.size() > 0);
	}
	if(toLevel != (g->NodesPerLevel.size() - 1))
	{
		assert(entryEdges_vector.size() > 0);
	}

	std::map<int, std::map<int, std::set<NWEdge*> > >  entryEdges_from_coordinates;
	std::map<int, std::map<int, std::set<NWEdge*> > >  exitEdges_to_coordinates;
	for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
	{
		NWEdge* edge = entryEdges_vector.at(eI);
		entryEdges_from_coordinates[edge->from_x][edge->from_z].insert(edge);
	}

	for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
	{
		NWEdge* edge = exitEdges_vector.at(eI);
		exitEdges_to_coordinates[edge->to_x][edge->to_z].insert(edge);
	}

	if(verbose)
		std::cout << "GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW(..) for " << entryEdges_vector.size() <<  " entry NW edges and " << exitEdges_vector.size() << " exit NW edges, going from level " << fromLevel << " to " << toLevel << "\n" << std::flush;

	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAnything;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAffineGap;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAnything;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAffineGap;

	// std::set<Node*> nodes_l0 = g->NodesPerLevel.at(fromLevel);

	for(int lI = fromLevel; lI <= toLevel; lI++)
	{
		if(debug) std::cout << "lI: " << lI << "\n" << std::flush;

		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

		if(lI == (g->NodesPerLevel.size() - 1))
		{
			lastPositionDistances_perZ_startInAffineGap.resize(nodes_thisLevel.size());
			lastPositionDistances_perZ_startNormal.resize(nodes_thisLevel.size());
		}

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(debug) std::cout << "Node: " << n << "\n" << std::flush;

			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == fromLevel)
			{
				if(debug) std::cout << "a1" << "\n" << std::flush;

				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n][0] = 0;

				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n][0] = minusInfinity;

				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n][0] = 0;

				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n][0] = minusInfinity;

				if(debug) std::cout << "a2" << "\n" << std::flush;
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				if(debug) std::cout << "b1" << "\n" << std::flush;


				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();

				if(debug) std::cout << "b2" << "\n" << std::flush;


				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					if(debug) std::cout << "c1" << "\n" << std::flush;

					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;

					if(debug) std::cout << "c2" << "\n" << std::flush;

					for(std::map<NWEdge*, double>::iterator edgeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).begin(); edgeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).end(); edgeDistIt++)
					{
						NWEdge* interestingEdge = edgeDistIt->first;

						if(debug) std::cout << "d1" << "\n" << std::flush;

						double previousDistance_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).at(interestingEdge);
						double previousDistance_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap.at(nodeFrom).at(interestingEdge);
						double previousDistance_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything.at(nodeFrom).at(interestingEdge);
						double previousDistance_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap.at(nodeFrom).at(interestingEdge);

						if(debug) std::cout << "d2" << "\n" << std::flush;


						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								double newDistance_notStartInAffineGap_affine_extend;
								double newDistance_notStartInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start = minusInfinity;
									if(previousDistance_notStartInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_notStartInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_graphGap;
									}
								}
								else
								{
									newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_extendGap;
									newDistance_notStartInAffineGap_affine_start = previousDistance_notStartInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								double newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start > newDistance_notStartInAffineGap_affine_extend) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if((runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge) < newDistance_notStartInAffineGap_affine_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[interestingEdge] = newDistance_notStartInAffineGap_affine_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								double newDistance_startInAffineGap_affine_extend;
								double newDistance_startInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start = minusInfinity;
									if(previousDistance_startInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_startInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_graphGap;
									}
								}
								else
								{
									newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_extendGap;
									newDistance_startInAffineGap_affine_start = previousDistance_startInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								double newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start > newDistance_startInAffineGap_affine_extend) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if((runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge) < newDistance_startInAffineGap_affine_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[interestingEdge] = newDistance_startInAffineGap_affine_maximum;
								}
							}
						}

						if(debug) std::cout << "d3" << "\n" << std::flush;

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								double alternativeScore_endInAffine = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_notStartInAffineGap_endInAnything;
								}
								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[interestingEdge] = endInAnything_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								double alternativeScore_endInAffine = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_startInAffineGap_endInAnything;
								}
								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[interestingEdge] = endInAnything_maximum;
								}
							}
						}
					}

					if(debug) std::cout << "c3" << "\n" << std::flush;



					if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
					{
						std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
						for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
						{
							NWEdge* e = *nwEdgeIt;

							runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[e] = 0;
							runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[e] = minusInfinity;

							runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[e] = 0;
							runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[e] = 0;
						}
					}

					if(debug) std::cout << "c4" << "\n" << std::flush;

				}

				if(debug) std::cout << "b3" << "\n" << std::flush;

			}

			if(lI != (g->NodesPerLevel.size() -1 ))
			{
				if(entryEdges_from_coordinates.count(lI) && entryEdges_from_coordinates.at(lI).count(z))
				{
					std::set<NWEdge*> nwEdges = entryEdges_from_coordinates.at(lI).at(z);

					for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
					{
						NWEdge* entryEdge = *nwEdgeIt;

						if(NWedges_graphDist_startAffineGap.count(entryEdge) == 0)
						{
							NWedges_graphDist_startAffineGap[entryEdge] = std::map<NWEdge*, graphPointDistance>();
							NWedges_graphDist_startNormal[entryEdge] = std::map<NWEdge*, graphPointDistance>();
						}

						for(std::map<NWEdge*, double>::iterator exitEdgeIt = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).begin(); exitEdgeIt != runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).end(); exitEdgeIt++)
						{
							NWEdge* exitEdge = exitEdgeIt->first;
							assert(entryEdge != 0);
							if((exitEdge == 0) || (entryEdge->from_y >= exitEdge->to_y))
							{
								{
									double score = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									assert(score <= 0);

									graphPointDistance distanceSpecifier;
									distanceSpecifier.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.start_in_affine_sequenceGap = true;
									NWedges_graphDist_startAffineGap[entryEdge][exitEdge] = distanceSpecifier;
								}
								{
									double score = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									assert(score <= 0);

									graphPointDistance distanceSpecifier;
									distanceSpecifier.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.start_in_affine_sequenceGap = false;
									NWedges_graphDist_startNormal[entryEdge][exitEdge] = distanceSpecifier;
								}
							}
						}
					}
				}
			}
			else
			{
				for(std::map<NWEdge*, double>::iterator edgeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).begin(); edgeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).end(); edgeDistIt++)
				{
					// entryEdge is now the virtual edge beyond the final field
					NWEdge* exitEdge = edgeDistIt->first;
					{
						graphPointDistance distanceSpecifier;
						distanceSpecifier.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.start_in_affine_sequenceGap = true;
						assert(lastPositionDistances_perZ_startInAffineGap.at(z).count(exitEdge) == 0);
						lastPositionDistances_perZ_startInAffineGap.at(z)[exitEdge] = distanceSpecifier;
					}

					{
						graphPointDistance distanceSpecifier;
						distanceSpecifier.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.start_in_affine_sequenceGap = false;
						assert(lastPositionDistances_perZ_startNormal.at(z).count(exitEdge) == 0);
						lastPositionDistances_perZ_startNormal.at(z)[exitEdge] = distanceSpecifier;
					}
				}
			}
		}

		runningNodeDistances_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		runningNodeDistances_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

	}

	size_t calculated_distances = 0;
	for(std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >::iterator edgeIt = NWedges_graphDist_startNormal.begin(); edgeIt != NWedges_graphDist_startNormal.end(); edgeIt++)
	{
		for(std::map<NWEdge*, graphPointDistance>::iterator edgeIt2 = edgeIt->second.begin(); edgeIt2 != edgeIt->second.end(); edgeIt2++)
		{
			calculated_distances += 1;
		}
	}

	if(verbose) std::cout << "GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW(..): Have " << calculated_distances << " distances.\n" << std::flush;
}


void GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW_local(int fromLevel, int toLevel, std::vector<NWEdge*>& entryEdges_vector, std::vector<NWEdge*>& exitEdges_vector, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal)
{
	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();

	bool verbose = false;
	bool debug = false;

	if(fromLevel != 0)
	{
		assert(exitEdges_vector.size() > 0);
	}
	if(toLevel != (g->NodesPerLevel.size() - 1))
	{
		assert(entryEdges_vector.size() > 0);
	}

	std::map<int, std::map<int, std::set<NWEdge*> > >  entryEdges_from_coordinates;
	std::map<int, std::map<int, std::set<NWEdge*> > >  exitEdges_to_coordinates;
	for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
	{
		NWEdge* edge = entryEdges_vector.at(eI);
		entryEdges_from_coordinates[edge->from_x][edge->from_z].insert(edge);
	}

	for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
	{
		NWEdge* edge = exitEdges_vector.at(eI);
		exitEdges_to_coordinates[edge->to_x][edge->to_z].insert(edge);
	}

	if(verbose)
		std::cout << "GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW(..) for " << entryEdges_vector.size() <<  " entry NW edges and " << exitEdges_vector.size() << " exit NW edges, going from level " << fromLevel << " to " << toLevel << "\n" << std::flush;

	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAnything;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAffineGap;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAnything;
	std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAffineGap;

	// std::set<Node*> nodes_l0 = g->NodesPerLevel.at(fromLevel);

	for(int lI = fromLevel; lI <= toLevel; lI++)
	{
		if(debug) std::cout << "lI: " << lI << "\n" << std::flush;

		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<NWEdge*, double> > runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

		if(lI == (g->NodesPerLevel.size() - 1))
		{
			lastPositionDistances_perZ_startInAffineGap.resize(nodes_thisLevel.size());
			lastPositionDistances_perZ_startNormal.resize(nodes_thisLevel.size());
		}

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(debug) std::cout << "Node: " << n << "\n" << std::flush;

			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == fromLevel)
			{
				if(debug) std::cout << "a1" << "\n" << std::flush;

				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n][0] = 0;

				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n][0] = minusInfinity;

				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n][0] = 0;

				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				if(lI == 0)
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n][0] = minusInfinity;

				if(debug) std::cout << "a2" << "\n" << std::flush;
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				if(debug) std::cout << "b1" << "\n" << std::flush;


				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<NWEdge*, double>();
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<NWEdge*, double>();

				if(debug) std::cout << "b2" << "\n" << std::flush;


				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					if(debug) std::cout << "c1" << "\n" << std::flush;

					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;

					if(debug) std::cout << "c2" << "\n" << std::flush;

					for(std::map<NWEdge*, double>::iterator edgeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).begin(); edgeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).end(); edgeDistIt++)
					{
						NWEdge* interestingEdge = edgeDistIt->first;

						if(debug) std::cout << "d1" << "\n" << std::flush;

						double previousDistance_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).at(interestingEdge);
						double previousDistance_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap.at(nodeFrom).at(interestingEdge);
						double previousDistance_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything.at(nodeFrom).at(interestingEdge);
						double previousDistance_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap.at(nodeFrom).at(interestingEdge);

						if(debug) std::cout << "d2" << "\n" << std::flush;


						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								double newDistance_notStartInAffineGap_affine_extend;
								double newDistance_notStartInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start = minusInfinity;
									if(previousDistance_notStartInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_notStartInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_graphGap;
									}
								}
								else
								{
									newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_extendGap;
									newDistance_notStartInAffineGap_affine_start = previousDistance_notStartInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								double newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start > newDistance_notStartInAffineGap_affine_extend) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if((runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge) < newDistance_notStartInAffineGap_affine_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[interestingEdge] = newDistance_notStartInAffineGap_affine_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								double newDistance_startInAffineGap_affine_extend;
								double newDistance_startInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start = minusInfinity;
									if(previousDistance_startInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_startInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_graphGap;
									}
								}
								else
								{
									newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_extendGap;
									newDistance_startInAffineGap_affine_start = previousDistance_startInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								double newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start > newDistance_startInAffineGap_affine_extend) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if((runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge) < newDistance_startInAffineGap_affine_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[interestingEdge] = newDistance_startInAffineGap_affine_maximum;
								}
							}
						}

						if(debug) std::cout << "d3" << "\n" << std::flush;

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								double alternativeScore_endInAffine = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_notStartInAffineGap_endInAnything;
								}
								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[interestingEdge] = endInAnything_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								double alternativeScore_endInAffine = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_startInAffineGap_endInAnything;
								}
								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[interestingEdge] = endInAnything_maximum;
								}
							}
						}
					}

					if(debug) std::cout << "c3" << "\n" << std::flush;



					if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
					{
						std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
						for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
						{
							NWEdge* e = *nwEdgeIt;

							runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[e] = 0;
							runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[e] = minusInfinity;

							runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[e] = 0;
							runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[e] = 0;
						}
					}

					if(debug) std::cout << "c4" << "\n" << std::flush;

				}

				if(debug) std::cout << "b3" << "\n" << std::flush;

			}

			if(lI != (g->NodesPerLevel.size() -1 ))
			{
				if(entryEdges_from_coordinates.count(lI) && entryEdges_from_coordinates.at(lI).count(z))
				{
					std::set<NWEdge*> nwEdges = entryEdges_from_coordinates.at(lI).at(z);

					for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
					{
						NWEdge* entryEdge = *nwEdgeIt;

						if(NWedges_graphDist_startAffineGap.count(entryEdge) == 0)
						{
							NWedges_graphDist_startAffineGap[entryEdge] = std::map<NWEdge*, graphPointDistance>();
							NWedges_graphDist_startNormal[entryEdge] = std::map<NWEdge*, graphPointDistance>();
						}

						for(std::map<NWEdge*, double>::iterator exitEdgeIt = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).begin(); exitEdgeIt != runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).end(); exitEdgeIt++)
						{
							NWEdge* exitEdge = exitEdgeIt->first;
							assert(entryEdge != 0);
							if((exitEdge == 0) || (entryEdge->from_y >= exitEdge->to_y))
							{
								{
									double score = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									assert(score <= 0);

									graphPointDistance distanceSpecifier;
									distanceSpecifier.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.start_in_affine_sequenceGap = true;
									NWedges_graphDist_startAffineGap[entryEdge][exitEdge] = distanceSpecifier;
								}
								{
									double score = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									assert(score <= 0);

									graphPointDistance distanceSpecifier;
									distanceSpecifier.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
									distanceSpecifier.start_in_affine_sequenceGap = false;
									NWedges_graphDist_startNormal[entryEdge][exitEdge] = distanceSpecifier;
								}
							}
						}
					}
				}
			}
			else
			{
				for(std::map<NWEdge*, double>::iterator edgeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).begin(); edgeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).end(); edgeDistIt++)
				{
					// entryEdge is now the virtual edge beyond the final field
					NWEdge* exitEdge = edgeDistIt->first;
					{
						graphPointDistance distanceSpecifier;
						distanceSpecifier.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.start_in_affine_sequenceGap = true;
						assert(lastPositionDistances_perZ_startInAffineGap.at(z).count(exitEdge) == 0);
						lastPositionDistances_perZ_startInAffineGap.at(z)[exitEdge] = distanceSpecifier;
					}

					{
						graphPointDistance distanceSpecifier;
						distanceSpecifier.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(exitEdge);
						distanceSpecifier.start_in_affine_sequenceGap = false;
						assert(lastPositionDistances_perZ_startNormal.at(z).count(exitEdge) == 0);
						lastPositionDistances_perZ_startNormal.at(z)[exitEdge] = distanceSpecifier;
					}
				}
			}
		}

		runningNodeDistances_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		runningNodeDistances_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

	}

	size_t calculated_distances = 0;
	for(std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >::iterator edgeIt = NWedges_graphDist_startNormal.begin(); edgeIt != NWedges_graphDist_startNormal.end(); edgeIt++)
	{
		for(std::map<NWEdge*, graphPointDistance>::iterator edgeIt2 = edgeIt->second.begin(); edgeIt2 != edgeIt->second.end(); edgeIt2++)
		{
			calculated_distances += 1;
		}
	}

	if(verbose) std::cout << "GraphAlignerUnique::affinelyCalculateGraphDistancesForVirtualNW(..): Have " << calculated_distances << " distances.\n" << std::flush;
}




void GraphAlignerUnique::seedAndExtend_backtrack(VirtualNWTable_Unique& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap)
{
	int completed_x = g->NodesPerLevel.size() - 1;
	int completed_y = sequence.length();
	int completed_z = finalScore_z;
	NWEdge* currentEdge = 0; // = finalScore_backtrack;

	assert(finalScore_backtrack != 0);

	vNW2.checkConsistency();

	std::set<NWEdge*> knownEdges = vNW2.getAllEdges();

	bool verboseBacktrack = (omp_get_thread_num() == 2);
	verboseBacktrack = false;

//	if(verbose || verboseBacktrack)
//		std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Start backtrace...\n" << std::flush;

	bool firstStep = true;
	do {
		assert(completed_x >= 0);
		assert(completed_y >= 0);

		if(verboseBacktrack)
		{
			#pragma omp critical
			{
				std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Backtrace step...\n";
				std::cout << "\t" << "completed_x: " << completed_x << "\n";
				std::cout << "\t" << "completed_y: " << completed_y << "\n";
				std::cout << "\t" << "completed_z: " << completed_z << "\n";
				std::cout << "\t" << "firstStep: " << firstStep << "\n" << std::flush;
			}
		}

		if(! firstStep)
		{
			// std::cout << "\t" << "currentEdge: " << currentEdge << " with score after edge: " << currentEdge->scoreAfterEdge << "\n" << std::flush;
			assert(currentEdge != 0);
			currentEdge->checkConsistency();
			assert(knownEdges.count(currentEdge));
		}
		else
		{
			if(verboseBacktrack)
			{
				#pragma omp critical
				{
					std::cout << "\t" << "finalScore: " << finalScore << "\n" << std::flush;
				}
			}
		}

		NWEdge* nextEdge_goingBack;
		if(firstStep)
		{
			nextEdge_goingBack = finalScore_backtrack;
		}
		else
		{
			nextEdge_goingBack = currentEdge->scoreBacktrack;
		}


		if(verboseBacktrack)
		{
			std::cout << "\t" << "nextEdge_goingBack: " << nextEdge_goingBack << "\n" << std::flush;
		}


		if((nextEdge_goingBack != 0) && (completed_x == nextEdge_goingBack->to_x) && (completed_y == nextEdge_goingBack->to_y) && (completed_z == nextEdge_goingBack->to_z))
		{
			// switch from the next edge to the current edge
			if(verboseBacktrack) std::cout << "\t" << "switch from " << currentEdge << " to " <<  nextEdge_goingBack << " as currentEdge.\n" << std::flush;
			currentEdge = nextEdge_goingBack;

		}
		else if((! firstStep) && (completed_x == currentEdge->to_x) && (completed_y == currentEdge->to_y) && (completed_z == currentEdge->to_z))
		{
			// output the current edge
			int diff_current_x = currentEdge->to_x - currentEdge->from_x;
			int diff_current_y = currentEdge->to_y - currentEdge->from_y;
			assert((diff_current_x >= 0) && (diff_current_x <= 1));
			assert((diff_current_y >= 0) && (diff_current_y <= 1));
			assert(diff_current_x || diff_current_y);

			std::string sequence_character;
			std::string graph_character;
			int graph_level;

			if((diff_current_x == 1) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);

				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 1) && (diff_current_y == 0))
			{
				// gap in sequence
				sequence_character = "_";
				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 0) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);
				graph_character = "_";
				graph_level = -1;
			}

			reconstructedSequence.append(sequence_character);
			reconstructedGraph.append(graph_character);
			reconstructedGraph_levels.push_back(graph_level);

			if(verboseBacktrack) std::cout << "\t" << "output " << graph_character << " / " << sequence_character << " from current edge.\n" << std::flush;
			if(verboseBacktrack) std::cout << "\t\t" << "score *after* this edge " << currentEdge << " (in normal coordinates): " << currentEdge->scoreAfterEdge <<  " (i.e. after consuming the " << currentEdge->to_x << "-th graph level and the " << currentEdge->to_y << "-th sequence character)\n" << std::flush;

			completed_x = currentEdge->from_x;
			completed_y = currentEdge->from_y;
			completed_z = currentEdge->from_z;
		}
		else if(firstStep || ((completed_x == currentEdge->from_x) && (completed_y == currentEdge->from_y) && (completed_z == currentEdge->from_z)))
		{
			// complete the gap behind the current edge

			int need_gaps_graph;
			int need_gaps_sequence;

			if(nextEdge_goingBack == 0)
			{
				// going back to origin
				need_gaps_sequence = completed_x;
				need_gaps_graph = completed_y;
			}
			else
			{
				// going back to other normal edge
				need_gaps_sequence = completed_x - nextEdge_goingBack->to_x;
				need_gaps_graph = completed_y - nextEdge_goingBack->to_y;
			}

			if(nextEdge_goingBack == 0)
			{
				if(verboseBacktrack) std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << "\t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << 0 << " / " << 0  << "\n" << std::flush;
			}
			else
			{
				if(verboseBacktrack) std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " \t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << nextEdge_goingBack->to_x << " / " << nextEdge_goingBack->to_y  << "\n" << std::flush;
			}

			assert(need_gaps_sequence >= 0);
			assert(need_gaps_graph >= 0);

			// gap in graph - ie we take sequence, and graph has gaps
			std::string gaps_in_graph_graphCharacters;
			std::string gaps_in_graph_sequenceCharacters;
			gaps_in_graph_graphCharacters.resize(need_gaps_graph, '_');
			std::vector<int> gaps_in_graph_levels;

			if(need_gaps_graph > 0)
			{
				gaps_in_graph_levels.resize(need_gaps_graph, -1);
				int stopYForGapInGraph = (nextEdge_goingBack != 0) ? nextEdge_goingBack->to_y : 0;
				for(int seqI = completed_y - 1; seqI >= stopYForGapInGraph; seqI--)
				{
					std::string sequence_character = sequence.substr(seqI, 1);
					gaps_in_graph_sequenceCharacters.append(sequence_character);
				}
			}

			std::string gaps_in_sequence_graphCharacters;
			std::string gaps_in_sequence_sequenceCharacters;
			std::vector<int> gaps_in_sequence_levels;
			std::vector<Edge*> gaps_in_sequence_usedEdges;

			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startAffinely;
			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startNormally;

			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startAffinely_MTM;
			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startNormally_MTM;


			if(need_gaps_sequence > 0)
			{
				int start_x_graph = completed_x;
				int stop_x_graph = (nextEdge_goingBack == 0) ? 0 : nextEdge_goingBack->to_x;

				int start_z_graph = completed_z;
				int stop_z_graph = (nextEdge_goingBack == 0) ? -1 : nextEdge_goingBack->to_z;

				gaps_in_sequence_sequenceCharacters.resize(need_gaps_sequence, '_');

				bool testMTM = false;

				if(verboseBacktrack)
				{
					//pragma omp critical
//					{
//						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) with:\n";
//						std::cout << "\t testMTM: " << testMTM << "\n";
//						std::cout << "\t start_x_graph: " << start_x_graph << "\n";
//						std::cout << "\t start_z_graph: " << start_z_graph << "\n";
//						std::cout << "\t stop_x_graph: " << stop_x_graph << "\n";
//						std::cout << "\t stop_z_graph: " << stop_z_graph << "\n" << std::flush;
//					}
				}

				if(testMTM)
				{
					findShortGappedGraphConnection_affine(
							stop_x_graph,
							stop_z_graph,
							start_x_graph,
							start_z_graph,
							distances_startAffinely,
							distances_startNormally
					);

					findShortGappedGraphConnection_affine_MTM(
							stop_x_graph,
							stop_z_graph,
							start_x_graph,
							start_z_graph,
							distances_startAffinely_MTM,
							distances_startNormally_MTM
					);

					assert(distances_startAffinely.size() == distances_startAffinely_MTM.size());
					assert(distances_startNormally.size() == distances_startNormally_MTM.size());

					int MTM_comparisons = 0;
					for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator outerZIt = distances_startAffinely.begin(); outerZIt != distances_startAffinely.end(); outerZIt++)
					{
						std::map<int, graphPointDistance_withBacktrack>& innerZMap = outerZIt->second;

						assert(distances_startAffinely_MTM.count(outerZIt->first));
						std::map<int, graphPointDistance_withBacktrack>& innerZMap_MTM = distances_startAffinely_MTM.at(outerZIt->first);

						assert(innerZMap.size() == innerZMap_MTM.size());
						for(std::map<int, graphPointDistance_withBacktrack>::iterator innerZIt = innerZMap.begin(); innerZIt != innerZMap.end(); innerZIt++)
						{
							assert(distances_startAffinely_MTM.at(outerZIt->first).count(innerZIt->first));
							assert(distances_startAffinely.at(outerZIt->first).at(innerZIt->first).endAffinely.S == distances_startAffinely_MTM.at(outerZIt->first).at(innerZIt->first).endAffinely.S);
							assert(distances_startAffinely.at(outerZIt->first).at(innerZIt->first).endInAnything.S == distances_startAffinely_MTM.at(outerZIt->first).at(innerZIt->first).endInAnything.S);
							MTM_comparisons++;
						}
					}

					for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator outerZIt = distances_startNormally.begin(); outerZIt != distances_startNormally.end(); outerZIt++)
					{
						std::map<int, graphPointDistance_withBacktrack>& innerZMap = outerZIt->second;

						assert(distances_startNormally_MTM.count(outerZIt->first));
						std::map<int, graphPointDistance_withBacktrack>& innerZMap_MTM = distances_startNormally_MTM.at(outerZIt->first);
						assert(innerZMap.size() == innerZMap_MTM.size());
						for(std::map<int, graphPointDistance_withBacktrack>::iterator innerZIt = innerZMap.begin(); innerZIt != innerZMap.end(); innerZIt++)
						{
							assert(distances_startNormally_MTM.at(outerZIt->first).count(innerZIt->first));
							assert(distances_startNormally.at(outerZIt->first).at(innerZIt->first).endAffinely.S == distances_startNormally_MTM.at(outerZIt->first).at(innerZIt->first).endAffinely.S);
							assert(distances_startNormally.at(outerZIt->first).at(innerZIt->first).endInAnything.S == distances_startNormally_MTM.at(outerZIt->first).at(innerZIt->first).endInAnything.S);
							MTM_comparisons++;
						}
					}
				}
				else
				{
					findShortGappedGraphConnection_affine_MTM(
							stop_x_graph,
							stop_z_graph,
							start_x_graph,
							start_z_graph,
							distances_startAffinely,
							distances_startNormally
					);
				}

				if(verboseBacktrack) 
				{				
					//pragma omp critical
//					{
//						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) done.\n" << std::flush;
//					}
				}
				
				if((start_z_graph != -1) && (stop_z_graph != -1) && (currentEdge != 0) && (nextEdge_goingBack != 0))
				{
					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endAffinely.S == NWedges_graphDist_startAffineGap.at(currentEdge).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endInAnything.S == NWedges_graphDist_startAffineGap.at(currentEdge).at(nextEdge_goingBack).Score_endInAnything);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endAffinely.S == NWedges_graphDist_startNormal.at(currentEdge).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endInAnything.S == NWedges_graphDist_startNormal.at(currentEdge).at(nextEdge_goingBack).Score_endInAnything);
				}

				if(firstStep)
				{
					assert(currentEdge == 0);
					// assert(stop_z_graph != 0);
					assert(start_z_graph != -1);
					assert(nextEdge_goingBack != 0);

					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endAffinely.S == lastPositionDistances_perZ_startAffineGap.at(start_z_graph).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endInAnything.S == lastPositionDistances_perZ_startAffineGap.at(start_z_graph).at(nextEdge_goingBack).Score_endInAnything);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endAffinely.S == lastPositionDistances_perZ_startNormal.at(start_z_graph).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endInAnything.S == lastPositionDistances_perZ_startNormal.at(start_z_graph).at(nextEdge_goingBack).Score_endInAnything);

				}
			}

			if(need_gaps_graph && need_gaps_sequence)
			{
				// we need both sequence gap and graph gap. In order to find the optimal sequence, we need
				// to find out which one to place first.
				// this code strongly resembles some code from above, where we compute the actual optimal score.

				bool exitEdge_isAffineSequenceGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isSequenceGap_affine());
				bool exitEdge_isGraphGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isGraphGap());

				std::vector<double> scores;
				std::vector<backtrackBookkeeper> backtracks;
				std::vector<bool> sequenceGapFirst;

				if(verboseBacktrack) 
				{				
					#pragma omp critical
					{
						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) - (need_gaps_graph && need_gaps_sequence)! Now scan " << distances_startNormally.size() << " alternatives.\n" << std::flush;
					}
				}
				
				assert((nextEdge_goingBack == 0) || (distances_startNormally.size() == 1));
				for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
				{
					int zStart = zStarterIt->first;
					std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
					assert(zStops.size() == 1);

					for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
					{				
						if(verboseBacktrack) 
						{				
							std::cout << "." << std::flush;
						}					
						
						int zStop = zStopIt->first;

						{
							// case 1: sequence gap comes first
							{
								// case 1.1: the sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = (exitEdge_isAffineSequenceGap) ? distances_startAffinely.at(zStart).at(zStop) : distances_startNormally.at(zStart).at(zStop);
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double graphGapScore = S_openGap + need_gaps_graph * S_extendGap;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(true);
							}

							{
								// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
							}
						}
						{
							// case 2: graph gap comes first
							{
								// case 2.1: sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(false);
							}
							{
								// case 2.2: sequence gap ends affinely
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endAffinely.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
								sequenceGapFirst.push_back(false);
							}
						}
					}
				}

				assert(scores.size() > 0);
				if(verboseBacktrack) 
				{				
					std::cout << " " << scores.size() << " " << std::flush;
				}					
				std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
				if(verboseBacktrack) 
				{				
					std::cout << " " << "a" << " " << std::flush;
				}			
				
				backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);
				bool selectedSequenceGapFirst = sequenceGapFirst.at(maxScore_sequenceGap.second);

				gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
				gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

				if(verboseBacktrack) 
				{				
					std::cout << " " << "b" << " " << std::flush;
				}			
				
				std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
				std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				if(verboseBacktrack) std::cout << "\t\t\t" << "selectedSequenceGapFirst: " << selectedSequenceGapFirst << "\n" << std::flush;


				if(selectedSequenceGapFirst)
				{
					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());

					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());
				}
				else
				{
					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());
				}
				
				if(verboseBacktrack) 
				{				
					std::cout << " " << "---" << " \n" << std::flush;
				}							

			}
			else
			{
				// we either need the graph gap or the sequence gap. The other fields are going to be empty. The orde
				// of concatentation is thus irrelevant.

				if(verboseBacktrack) 
				{				
					std::cout << " " << "	!(need_gaps_graph && need_gaps_sequence) " << " \n" << std::flush;
				}			
				
				if(need_gaps_sequence > 0)
				{
					// if we go to the origin, we may end up in multiple z values. We traverse all z values to
					// find the best alternative.

					std::vector<double> scores;
					std::vector<backtrackBookkeeper> backtracks;

					assert((nextEdge_goingBack == 0) || (distances_startNormally.size() == 1));

					if(verboseBacktrack) 
					{				
						std::cout << " " << " distances_startNormally.size(): " << distances_startNormally.size() << " \n" << std::flush;
					}	
				
					for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
					{
						int zStart = zStarterIt->first;
						std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
						assert(zStops.size() == 1);

						for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
						{
							if(verboseBacktrack) 
							{				
								std::cout << " . " << std::flush;
							}	
							int zStop = zStopIt->first;

							graphPointDistance_withBacktrack utilizedBacktrack = zStopIt->second;
							if((nextEdge_goingBack != 0) && (nextEdge_goingBack->isSequenceGap_affine()))
							{
								utilizedBacktrack = distances_startAffinely.at(zStart).at(zStop);
							}

							// case 1: we want to end in an affine gap
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endAffinely.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
							}

							// case 2: end in anything
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endInAnything.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
							}

						}
					}

					if(verboseBacktrack) 
					{				
						std::cout << " " << " scores.size(): " << scores.size() << " \n" << std::flush;
					}	
				
				
					assert(scores.size() > 0);
					std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
					backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);

					
					if(verboseBacktrack) 
					{				
						std::cout << " . " << std::flush;
					}	
				
				
					gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
					gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

					std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
					std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

					if(verboseBacktrack) 
					{				
						std::cout << "--\n" << std::flush;
					}						
				}

				if(verboseBacktrack) 
				{				
					std::cout << " append1 " << "\n" << "\n";
				}				
				reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
				reconstructedGraph.append(gaps_in_sequence_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				if(verboseBacktrack) 
				{				
					std::cout << " append2 " << "\n" << "\n";
				}	
				
				reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
				reconstructedGraph.append(gaps_in_graph_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());

				if(verboseBacktrack) 
				{				
					std::cout << " append " << "\n" << "\n";
				}			
			}

			// gap in
			if(nextEdge_goingBack == 0)
			{
				completed_x = 0;
				completed_y = 0;
			}
			else
			{
				completed_x = nextEdge_goingBack->to_x;
				completed_y = nextEdge_goingBack->to_y;
				completed_z = nextEdge_goingBack->to_z;
			}

			if(firstStep)
			{
				currentEdge = nextEdge_goingBack;
			}
		}
		else
		{
			// this case is unforeseen!
			std::cout << "\tUNEXPECTED!\n";
			std::cout << "\t\t" << "currentEdge: " << currentEdge << "\n";
			std::cout << "\t\t" << "currentEdge->from_x: " << currentEdge->from_x << "\n";
			std::cout << "\t\t" << "currentEdge->from_y: " << currentEdge->from_y << "\n";
			std::cout << "\t\t" << "currentEdge->from_z: " << currentEdge->from_z << "\n";
			std::cout << "\t\t" << "currentEdge->to_x: " << currentEdge->to_x << "\n";
			std::cout << "\t\t" << "currentEdge->to_y: " << currentEdge->to_y << "\n";
			std::cout << "\t\t" << "currentEdge->to_z: " << currentEdge->to_z << "\n\n" << std::flush;
			if(nextEdge_goingBack != 0)
			{
				std::cout << "\t\t\t" << "nextEdge: " << nextEdge_goingBack << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_x: " << nextEdge_goingBack->from_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_y: " << nextEdge_goingBack->from_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_z: " << nextEdge_goingBack->from_z << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_x: " << nextEdge_goingBack->to_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_y: " << nextEdge_goingBack->to_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_z: " << nextEdge_goingBack->to_z << "\n\n" << std::flush;
			}

			assert( 1 == 0 );
		}

		firstStep = false;

	} while((completed_x != 0) || (completed_y != 0));

	if(verboseBacktrack)	std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Backtrace done!\n" << std::flush;

	assert(reconstructedGraph.length() == reconstructedSequence.length());
	assert(reconstructedGraph_levels.size() == reconstructedGraph.length());

	std::reverse(reconstructedGraph.begin(), reconstructedGraph.end());
	std::reverse(reconstructedSequence.begin(), reconstructedSequence.end());
	std::reverse(reconstructedGraph_levels.begin(), reconstructedGraph_levels.end());
}



void GraphAlignerUnique::seedAndExtend_backtrack_local(VirtualNWTable_Unique& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap)
{
	double minusInfinity = -1 * numeric_limits<double>::max();

	int completed_x = g->NodesPerLevel.size() - 1;
	int completed_y = sequence.length();
	int completed_z = finalScore_z;
	NWEdge* currentEdge = 0; // = finalScore_backtrack;

	assert(finalScore_backtrack != 0);

	vNW2.checkConsistency();

	std::set<NWEdge*> knownEdges = vNW2.getAllEdges();

	bool verboseBacktrack = (omp_get_thread_num() == 2);
	verboseBacktrack = false;

//	if(verbose || verboseBacktrack)
//		std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Start backtrace...\n" << std::flush;

	bool firstStep = true;
	do {
		assert(completed_x >= 0);
		assert(completed_y >= 0);

		if(verboseBacktrack)
		{
			#pragma omp critical
			{
				std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Backtrace step...\n";
				std::cout << "\t" << "completed_x: " << completed_x << "\n";
				std::cout << "\t" << "completed_y: " << completed_y << "\n";
				std::cout << "\t" << "completed_z: " << completed_z << "\n";
				std::cout << "\t" << "firstStep: " << firstStep << "\n" << std::flush;
			}
		}

		if(! firstStep)
		{
			// std::cout << "\t" << "currentEdge: " << currentEdge << " with score after edge: " << currentEdge->scoreAfterEdge << "\n" << std::flush;
			assert(currentEdge != 0);
			currentEdge->checkConsistency();
			assert(knownEdges.count(currentEdge));
		}
		else
		{
			if(verboseBacktrack)
			{
				#pragma omp critical
				{
					std::cout << "\t" << "finalScore: " << finalScore << "\n" << std::flush;
				}
			}
		}

		NWEdge* nextEdge_goingBack;
		if(firstStep)
		{
			nextEdge_goingBack = finalScore_backtrack;
		}
		else
		{
			nextEdge_goingBack = currentEdge->scoreBacktrack;
		}


		if(verboseBacktrack)
		{
			std::cout << "\t" << "nextEdge_goingBack: " << nextEdge_goingBack << "\n" << std::flush;
		}


		if((nextEdge_goingBack != 0) && (completed_x == nextEdge_goingBack->to_x) && (completed_y == nextEdge_goingBack->to_y) && (completed_z == nextEdge_goingBack->to_z))
		{
			// switch from the next edge to the current edge
			if(verboseBacktrack) std::cout << "\t" << "switch from " << currentEdge << " to " <<  nextEdge_goingBack << " as currentEdge.\n" << std::flush;
			currentEdge = nextEdge_goingBack;

		}
		else if((! firstStep) && (completed_x == currentEdge->to_x) && (completed_y == currentEdge->to_y) && (completed_z == currentEdge->to_z))
		{
			// output the current edge
			int diff_current_x = currentEdge->to_x - currentEdge->from_x;
			int diff_current_y = currentEdge->to_y - currentEdge->from_y;
			assert((diff_current_x >= 0) && (diff_current_x <= 1));
			assert((diff_current_y >= 0) && (diff_current_y <= 1));
			assert(diff_current_x || diff_current_y);

			std::string sequence_character;
			std::string graph_character;
			int graph_level;

			if((diff_current_x == 1) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);

				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 1) && (diff_current_y == 0))
			{
				// gap in sequence
				sequence_character = "_";
				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 0) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);
				graph_character = "_";
				graph_level = -1;
			}

			reconstructedSequence.append(sequence_character);
			reconstructedGraph.append(graph_character);
			reconstructedGraph_levels.push_back(graph_level);

			if(verboseBacktrack) std::cout << "\t" << "output " << graph_character << " / " << sequence_character << " from current edge.\n" << std::flush;
			if(verboseBacktrack) std::cout << "\t\t" << "score *after* this edge " << currentEdge << " (in normal coordinates): " << currentEdge->scoreAfterEdge <<  " (i.e. after consuming the " << currentEdge->to_x << "-th graph level and the " << currentEdge->to_y << "-th sequence character)\n" << std::flush;

			completed_x = currentEdge->from_x;
			completed_y = currentEdge->from_y;
			completed_z = currentEdge->from_z;
		}
		else if(firstStep || ((completed_x == currentEdge->from_x) && (completed_y == currentEdge->from_y) && (completed_z == currentEdge->from_z)))
		{
			// complete the gap behind the current edge

			int need_gaps_graph;
			int need_gaps_sequence;

			if(nextEdge_goingBack == 0)
			{
				// going back to origin
				need_gaps_sequence = completed_x;
				need_gaps_graph = completed_y;
			}
			else
			{
				// going back to other normal edge
				need_gaps_sequence = completed_x - nextEdge_goingBack->to_x;
				need_gaps_graph = completed_y - nextEdge_goingBack->to_y;
			}

			if(nextEdge_goingBack == 0)
			{
				if(verboseBacktrack) std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << "\t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << 0 << " / " << 0  << "\n" << std::flush;
			}
			else
			{
				if(verboseBacktrack) std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " \t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << nextEdge_goingBack->to_x << " / " << nextEdge_goingBack->to_y  << "\n" << std::flush;
			}

			assert(need_gaps_sequence >= 0);
			assert(need_gaps_graph >= 0);

			// gap in graph - ie we take sequence, and graph has gaps
			std::string gaps_in_graph_graphCharacters;
			std::string gaps_in_graph_sequenceCharacters;
			gaps_in_graph_graphCharacters.resize(need_gaps_graph, '_');
			std::vector<int> gaps_in_graph_levels;

			if(need_gaps_graph > 0)
			{
				gaps_in_graph_levels.resize(need_gaps_graph, -1);
				int stopYForGapInGraph = (nextEdge_goingBack != 0) ? nextEdge_goingBack->to_y : 0;
				for(int seqI = completed_y - 1; seqI >= stopYForGapInGraph; seqI--)
				{
					std::string sequence_character = sequence.substr(seqI, 1);
					gaps_in_graph_sequenceCharacters.append(sequence_character);
				}
			}

			std::string gaps_in_sequence_graphCharacters;
			std::string gaps_in_sequence_sequenceCharacters;
			std::vector<int> gaps_in_sequence_levels;
			std::vector<Edge*> gaps_in_sequence_usedEdges;

			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startAffinely;
			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startNormally;

			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startAffinely_MTM;
			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startNormally_MTM;


			if(need_gaps_sequence > 0)
			{
				int start_x_graph = completed_x;
				int stop_x_graph = (nextEdge_goingBack == 0) ? 0 : nextEdge_goingBack->to_x;

				int start_z_graph = completed_z;
				int stop_z_graph = (nextEdge_goingBack == 0) ? -1 : nextEdge_goingBack->to_z;

				gaps_in_sequence_sequenceCharacters.resize(need_gaps_sequence, '_');

				if(verboseBacktrack)
				{
					//pragma omp critical
//					{
//						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) with:\n";
//						std::cout << "\t testMTM: " << testMTM << "\n";
//						std::cout << "\t start_x_graph: " << start_x_graph << "\n";
//						std::cout << "\t start_z_graph: " << start_z_graph << "\n";
//						std::cout << "\t stop_x_graph: " << stop_x_graph << "\n";
//						std::cout << "\t stop_z_graph: " << stop_z_graph << "\n" << std::flush;
//					}
				}

				if(firstStep)
				{
					// std::cerr << "\nA firstStep\n" << std::flush;
					
					for(unsigned int zI = 0; zI < g->NodesPerLevel.back().size(); zI++)
					{
						graphPointDistance_withBacktrack distance_startAffinely;
						distance_startAffinely.start_in_affine_sequenceGap = true;
						distance_startAffinely.endAffinely.S = minusInfinity;
						distance_startAffinely.endInAnything.S = minusInfinity;

						graphPointDistance_withBacktrack distance_startNormally;
						distance_startNormally.start_in_affine_sequenceGap = false;
						distance_startNormally.endAffinely.S = minusInfinity;
						distance_startNormally.endInAnything.S = 0;

						distances_startAffinely[zI][start_z_graph] = distance_startAffinely;
						distances_startNormally[zI][start_z_graph] = distance_startNormally;
					}

					gaps_in_sequence_sequenceCharacters = "";
				}
				else if(nextEdge_goingBack == 0)
				{
					// std::cerr << "\nA nextEdge_goingBack == 0\n" << std::flush;
				
					// going back to origin
					graphPointDistance_withBacktrack distance_startAffinely;
					distance_startAffinely.start_in_affine_sequenceGap = true;
					distance_startAffinely.endAffinely.S = minusInfinity;
					distance_startAffinely.endInAnything.S = minusInfinity;

					graphPointDistance_withBacktrack distance_startNormally;
					distance_startNormally.start_in_affine_sequenceGap = false;
					distance_startNormally.endAffinely.S = minusInfinity;
					distance_startNormally.endInAnything.S = 0;

					distances_startAffinely[stop_z_graph][start_z_graph] = distance_startAffinely;
					distances_startNormally[stop_z_graph][start_z_graph] = distance_startNormally;

					gaps_in_sequence_sequenceCharacters = "";
				}
				else
				{
				
					// std::cerr << "\nA Call MTM\n" << std::flush;
					// std::cerr << "\t" << start_x_graph << " - " << stop_x_graph << "\n";
					// std::cerr << "\t" << start_z_graph << " - " << stop_z_graph << "\n";
		
					findShortGappedGraphConnection_affine_MTM(
							stop_x_graph,
							stop_z_graph,
							start_x_graph,
							start_z_graph,
							distances_startAffinely,
							distances_startNormally
					);
					
					// std::cerr << "\t" << distances_startAffinely.size() << "\n";
					// std::cerr << "\t" << distances_startNormally.size() << "\n" << std::flush;
				}

				if(verboseBacktrack)
				{
					//pragma omp critical
//					{
//						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) done.\n" << std::flush;
//					}
				}

				if((start_z_graph != -1) && (stop_z_graph != -1) && (currentEdge != 0) && (nextEdge_goingBack != 0))
				{
					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endAffinely.S == NWedges_graphDist_startAffineGap.at(currentEdge).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endInAnything.S == NWedges_graphDist_startAffineGap.at(currentEdge).at(nextEdge_goingBack).Score_endInAnything);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endAffinely.S == NWedges_graphDist_startNormal.at(currentEdge).at(nextEdge_goingBack).Score_endAffinely);
					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endInAnything.S == NWedges_graphDist_startNormal.at(currentEdge).at(nextEdge_goingBack).Score_endInAnything);
				}

				if(firstStep)
				{
					assert(currentEdge == 0);
					// assert(stop_z_graph != 0);
					assert(start_z_graph != -1);
					assert(nextEdge_goingBack != 0);

//					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endAffinely.S == lastPositionDistances_perZ_startAffineGap.at(start_z_graph).at(nextEdge_goingBack).Score_endAffinely);
//					assert(distances_startAffinely.at(stop_z_graph).at(start_z_graph).endInAnything.S == lastPositionDistances_perZ_startAffineGap.at(start_z_graph).at(nextEdge_goingBack).Score_endInAnything);
//					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endAffinely.S == lastPositionDistances_perZ_startNormal.at(start_z_graph).at(nextEdge_goingBack).Score_endAffinely);
//					assert(distances_startNormally.at(stop_z_graph).at(start_z_graph).endInAnything.S == lastPositionDistances_perZ_startNormal.at(start_z_graph).at(nextEdge_goingBack).Score_endInAnything);

				}
			}

			if(need_gaps_graph && need_gaps_sequence)
			{
				// we need both sequence gap and graph gap. In order to find the optimal sequence, we need
				// to find out which one to place first.
				// this code strongly resembles some code from above, where we compute the actual optimal score.

				bool exitEdge_isAffineSequenceGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isSequenceGap_affine());
				bool exitEdge_isGraphGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isGraphGap());

				std::vector<double> scores;
				std::vector<backtrackBookkeeper> backtracks;
				std::vector<bool> sequenceGapFirst;

				if(verboseBacktrack)
				{
					#pragma omp critical
					{
						std::cout << Utilities::timestamp() << " Thread " << omp_get_thread_num() << ", findShortGappedGraphConnection_nonAffine(..) - (need_gaps_graph && need_gaps_sequence)! Now scan " << distances_startNormally.size() << " alternatives.\n" << std::flush;
					}
				}

				assert((nextEdge_goingBack == 0) || (distances_startNormally.size() == 1));
				for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
				{
					int zStart = zStarterIt->first;
					std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
					assert(zStops.size() == 1);

					for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
					{
						if(verboseBacktrack)
						{
							std::cout << "." << std::flush;
						}

						int zStop = zStopIt->first;

						{
							// case 1: sequence gap comes first
							{
								// case 1.1: the sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = (exitEdge_isAffineSequenceGap) ? distances_startAffinely.at(zStart).at(zStop) : distances_startNormally.at(zStart).at(zStop);
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double graphGapScore = S_openGap + need_gaps_graph * S_extendGap;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(true);
							}

							{
								// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
							}
						}
						{
							// case 2: graph gap comes first
							{
								// case 2.1: sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(false);
							}
							{
								// case 2.2: sequence gap ends affinely
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endAffinely.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
								sequenceGapFirst.push_back(false);
							}
						}
					}
				}

				assert(scores.size() > 0);
				if(verboseBacktrack)
				{
					std::cout << " " << scores.size() << " " << std::flush;
				}
				std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
				if(verboseBacktrack)
				{
					std::cout << " " << "a" << " " << std::flush;
				}

				backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);
				bool selectedSequenceGapFirst = sequenceGapFirst.at(maxScore_sequenceGap.second);

				gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
				gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

				if(verboseBacktrack)
				{
					std::cout << " " << "b" << " " << std::flush;
				}

				std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
				std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				if(verboseBacktrack) std::cout << "\t\t\t" << "selectedSequenceGapFirst: " << selectedSequenceGapFirst << "\n" << std::flush;


				if(selectedSequenceGapFirst)
				{
					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());

					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());
				}
				else
				{
					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());
				}

				if(verboseBacktrack)
				{
					std::cout << " " << "---" << " \n" << std::flush;
				}

			}
			else
			{
				// we either need the graph gap or the sequence gap. The other fields are going to be empty. The orde
				// of concatentation is thus irrelevant.

				if(verboseBacktrack)
				{
					std::cout << " " << "	!(need_gaps_graph && need_gaps_sequence) " << " \n" << std::flush;
				}

				if(need_gaps_sequence > 0)
				{
					// if we go to the origin, we may end up in multiple z values. We traverse all z values to
					// find the best alternative.

					std::vector<double> scores;
					std::vector<backtrackBookkeeper> backtracks;

					if(!((nextEdge_goingBack == 0) || firstStep || (distances_startNormally.size() == 1)))
					{
						std::cout << "\n\n";
						std::cout << "firstStep: " << firstStep << "\n";						
						std::cout << "nextEdge_goingBack: " << nextEdge_goingBack << "\n";
						std::cout << "distances_startNormally.size(): " << distances_startNormally.size() << "\n\n" << std::flush;
					}
					assert((nextEdge_goingBack == 0) || firstStep || (distances_startNormally.size() == 1));

					if(verboseBacktrack)
					{
						std::cout << " " << " distances_startNormally.size(): " << distances_startNormally.size() << " \n" << std::flush;
					}

					for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
					{
						int zStart = zStarterIt->first;
						std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
						assert(zStops.size() == 1);

						for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
						{
							if(verboseBacktrack)
							{
								std::cout << " . " << std::flush;
							}
							int zStop = zStopIt->first;

							graphPointDistance_withBacktrack utilizedBacktrack = zStopIt->second;
							if((nextEdge_goingBack != 0) && (nextEdge_goingBack->isSequenceGap_affine()))
							{
								utilizedBacktrack = distances_startAffinely.at(zStart).at(zStop);
							}

							// case 1: we want to end in an affine gap
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endAffinely.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
							}

							// case 2: end in anything
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endInAnything.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
							}

						}
					}

					if(verboseBacktrack)
					{
						std::cout << " " << " scores.size(): " << scores.size() << " \n" << std::flush;
					}


					assert(scores.size() > 0);
					std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMaxP_nonCritical(scores, &(rng_seeds.at(omp_get_thread_num())));
					backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);


					if(verboseBacktrack)
					{
						std::cout << " . " << std::flush;
					}


					gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
					gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

					std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
					std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

					if(verboseBacktrack)
					{
						std::cout << "--\n" << std::flush;
					}
				}

				if(verboseBacktrack)
				{
					std::cout << " append1 " << "\n" << "\n";
				}
				reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
				reconstructedGraph.append(gaps_in_sequence_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				if(verboseBacktrack)
				{
					std::cout << " append2 " << "\n" << "\n";
				}

				reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
				reconstructedGraph.append(gaps_in_graph_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());

				if(verboseBacktrack)
				{
					std::cout << " append " << "\n" << "\n";
				}
			}

			// gap in
			if(nextEdge_goingBack == 0)
			{
				completed_x = 0;
				completed_y = 0;
			}
			else
			{
				completed_x = nextEdge_goingBack->to_x;
				completed_y = nextEdge_goingBack->to_y;
				completed_z = nextEdge_goingBack->to_z;
			}

			if(firstStep)
			{
				currentEdge = nextEdge_goingBack;
			}
		}
		else
		{
			// this case is unforeseen!
			std::cout << "\tUNEXPECTED!\n";
			std::cout << "\t\t" << "currentEdge: " << currentEdge << "\n";
			std::cout << "\t\t" << "currentEdge->from_x: " << currentEdge->from_x << "\n";
			std::cout << "\t\t" << "currentEdge->from_y: " << currentEdge->from_y << "\n";
			std::cout << "\t\t" << "currentEdge->from_z: " << currentEdge->from_z << "\n";
			std::cout << "\t\t" << "currentEdge->to_x: " << currentEdge->to_x << "\n";
			std::cout << "\t\t" << "currentEdge->to_y: " << currentEdge->to_y << "\n";
			std::cout << "\t\t" << "currentEdge->to_z: " << currentEdge->to_z << "\n\n" << std::flush;
			if(nextEdge_goingBack != 0)
			{
				std::cout << "\t\t\t" << "nextEdge: " << nextEdge_goingBack << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_x: " << nextEdge_goingBack->from_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_y: " << nextEdge_goingBack->from_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_z: " << nextEdge_goingBack->from_z << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_x: " << nextEdge_goingBack->to_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_y: " << nextEdge_goingBack->to_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_z: " << nextEdge_goingBack->to_z << "\n\n" << std::flush;
			}

			assert( 1 == 0 );
		}

		firstStep = false;

	} while((completed_x != 0) || (completed_y != 0));

	if(verboseBacktrack)	std::cout << Utilities::timestamp()  << " Thread " << omp_get_thread_num() << " Backtrace done!\n" << std::flush;

	assert(reconstructedGraph.length() == reconstructedSequence.length());
	assert(reconstructedGraph_levels.size() == reconstructedGraph.length());

	std::reverse(reconstructedGraph.begin(), reconstructedGraph.end());
	std::reverse(reconstructedSequence.begin(), reconstructedSequence.end());
	std::reverse(reconstructedGraph_levels.begin(), reconstructedGraph_levels.end());
}




void GraphAlignerUnique::findShortGappedGraphConnection_affine(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally)
{
	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();
	int thisThread = omp_get_thread_num();
	assert(stop_x > start_x);
	assert(start_x  >= 0);
	assert(stop_x < (int)g->NodesPerLevel.size());

	assert((start_z == -1) || ((start_z >= 0) && (start_z < (int)g->NodesPerLevel.at(start_x).size())));
	assert((stop_z == -1) || ((stop_z >= 0) && (stop_z < (int)g->NodesPerLevel.at(stop_x).size())));

	// bool verbose = true;
	distances_startAffinely.clear();
	distances_startNormally.clear();

	std::map<Node*, std::map<Node*, backtrackBookkeeper> >* runningNodeDistances_notStartInAffineGap_endInAnything = 0;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> >* runningNodeDistances_notStartInAffineGap_endInAffineGap = 0;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> >* runningNodeDistances_startInAffineGap_endInAnything = 0;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> >* runningNodeDistances_startInAffineGap_endInAffineGap = 0;

	bool printMe = (((stop_x - start_x) > 100000) && (omp_get_thread_num() == 2));
	// printMe = false;
	
	std::ostringstream print_buffer;
	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine(..): From " << start_x << " to " << stop_x << "\n" << std::flush;	
	}
	
	std::vector<std::map<Node*, std::map<Node*, backtrackBookkeeper> >> runningNodeDistances_notStartInAffineGap_endInAnything_allLevels;
	std::vector<std::map<Node*, std::map<Node*, backtrackBookkeeper> >> runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels;
	std::vector<std::map<Node*, std::map<Node*, backtrackBookkeeper> >> runningNodeDistances_startInAffineGap_endInAnything_allLevels;
	std::vector<std::map<Node*, std::map<Node*, backtrackBookkeeper> >> runningNodeDistances_startInAffineGap_endInAffineGap_allLevels;
	runningNodeDistances_notStartInAffineGap_endInAnything_allLevels.resize(stop_x - start_x + 1);
	runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels.resize(stop_x - start_x + 1);
	runningNodeDistances_startInAffineGap_endInAnything_allLevels.resize(stop_x - start_x + 1);
	runningNodeDistances_startInAffineGap_endInAffineGap_allLevels.resize(stop_x - start_x + 1);

	int lastPrintedDistance;
	for(int lI = start_x; lI <= stop_x; lI++)
	{
		if(printMe)
		{
			if((lI == start_x) || ((lI - lastPrintedDistance) >= 50000))
			{
				std::cout << "\t" << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", level " << lI << "\n" << std::flush;
				lastPrintedDistance = lI;
			}
		}
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<Node*, backtrackBookkeeper> >& runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel = runningNodeDistances_notStartInAffineGap_endInAnything_allLevels.at(lI - start_x);
		std::map<Node*, std::map<Node*, backtrackBookkeeper> >& runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel = runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels.at(lI - start_x);
		std::map<Node*, std::map<Node*, backtrackBookkeeper> >& runningNodeDistances_startInAffineGap_endInAnything_thisLevel = runningNodeDistances_startInAffineGap_endInAnything_allLevels.at(lI - start_x);
		std::map<Node*, std::map<Node*, backtrackBookkeeper> >& runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel = runningNodeDistances_startInAffineGap_endInAffineGap_allLevels.at(lI - start_x);

//		if(lI == start_x)
//		{
//			lastPositionDistances_perZ_startInAffineGap.resize(nodes_thisLevel.size());
//			lastPositionDistances_perZ_startNormal.resize(nodes_thisLevel.size());
//		}

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{

			Node* n = *nIt;
			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == start_x)
			{
				// jump over node if start_z specified
				bool takeNode = true;
				if(start_z != -1)
				{
					if(z != start_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n][n].S = 0;
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n][n].S = minusInfinity;

					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n][n].S = 0;
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();;
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n][n].S = (lI == 0) ? minusInfinity : 0;
				}
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;

					assert(runningNodeDistances_notStartInAffineGap_endInAnything != 0);
					if(! runningNodeDistances_notStartInAffineGap_endInAnything->count(nodeFrom))
					{
						continue;
					}

					for(std::map<Node*, backtrackBookkeeper>::iterator nodeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything->at(nodeFrom).begin(); nodeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything->at(nodeFrom).end(); nodeDistIt++)
					{
						Node* interestingNode = nodeDistIt->first;

						backtrackBookkeeper& previousDistance_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything->at(nodeFrom).at(interestingNode);
						backtrackBookkeeper& previousDistance_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap->at(nodeFrom).at(interestingNode);
						backtrackBookkeeper& previousDistance_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything->at(nodeFrom).at(interestingNode);
						backtrackBookkeeper& previousDistance_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap->at(nodeFrom).at(interestingNode);

						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								backtrackBookkeeper newDistance_notStartInAffineGap_affine_start; // = previousDistance_notStartInAffineGap_endInAnything;
								backtrackBookkeeper newDistance_notStartInAffineGap_affine_extend; // = previousDistance_notStartInAffineGap_endInAffineGap;

								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start.S =           minusInfinity;

									newDistance_notStartInAffineGap_affine_extend.S =          previousDistance_notStartInAffineGap_endInAffineGap.S;
									newDistance_notStartInAffineGap_affine_extend.backtrack = &previousDistance_notStartInAffineGap_endInAffineGap;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_notStartInAffineGap_affine_start.S =           previousDistance_notStartInAffineGap_endInAnything.S + (S_openGap + S_extendGap);
									newDistance_notStartInAffineGap_affine_start.backtrack =  &previousDistance_notStartInAffineGap_endInAnything;
									newDistance_notStartInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);

									newDistance_notStartInAffineGap_affine_extend.S =          previousDistance_notStartInAffineGap_endInAffineGap.S + S_extendGap;
									newDistance_notStartInAffineGap_affine_extend.backtrack = &previousDistance_notStartInAffineGap_endInAffineGap;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}

								backtrackBookkeeper newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start.S > newDistance_notStartInAffineGap_affine_extend.S) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if((runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode).S < newDistance_notStartInAffineGap_affine_maximum.S))
								{
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[interestingNode] = newDistance_notStartInAffineGap_affine_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								backtrackBookkeeper newDistance_startInAffineGap_affine_start; // = previousDistance_startInAffineGap_endInAnything;
								backtrackBookkeeper newDistance_startInAffineGap_affine_extend; // = previousDistance_startInAffineGap_endInAffineGap;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start.S = minusInfinity;

									newDistance_startInAffineGap_affine_extend.S = previousDistance_startInAffineGap_endInAffineGap.S;
									newDistance_startInAffineGap_affine_extend.backtrack = &previousDistance_startInAffineGap_endInAffineGap;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_startInAffineGap_affine_start.S = previousDistance_startInAffineGap_endInAnything.S + (S_openGap + S_extendGap);
									newDistance_startInAffineGap_affine_start.backtrack = &previousDistance_startInAffineGap_endInAnything;
									newDistance_startInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);

									newDistance_startInAffineGap_affine_extend.S = previousDistance_startInAffineGap_endInAffineGap.S + S_extendGap;
									newDistance_startInAffineGap_affine_extend.backtrack = &previousDistance_startInAffineGap_endInAffineGap;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}

								backtrackBookkeeper newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start.S > newDistance_startInAffineGap_affine_extend.S) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if((runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode).S < newDistance_startInAffineGap_affine_maximum.S))
								{
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[interestingNode] = newDistance_startInAffineGap_affine_maximum;
								}
							}
						}

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								backtrackBookkeeper alternativeScore_endInAffine = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode);

								backtrackBookkeeper followGraphGap; // = previousDistance_notStartInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap.S = previousDistance_notStartInAffineGap_endInAnything;
									followGraphGap.S = previousDistance_notStartInAffineGap_endInAnything.S;
									followGraphGap.backtrack = &previousDistance_notStartInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine.S) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(interestingNode).S < endInAnything_maximum.S))
								{
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[interestingNode] = endInAnything_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								backtrackBookkeeper alternativeScore_endInAffine = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode);
								backtrackBookkeeper followGraphGap;// = previousDistance_startInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap = previousDistance_startInAffineGap_endInAnything;
									followGraphGap.S = previousDistance_startInAffineGap_endInAnything.S;
									followGraphGap.backtrack = &previousDistance_startInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine.S) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(interestingNode).S < endInAnything_maximum.S))
								{
									runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[interestingNode] = endInAnything_maximum;
								}
							}
						}
					}
				}
			}

			if(lI == stop_x)
			{
				bool takeNode = true;
				if(stop_z != -1)
				{
					if(z != stop_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{

					auto backtrackThroughGraph = [&](backtrackBookkeeper start) -> backtrackBookkeeper {

						std::vector<Edge*> usedEdges;
						std::string graphSequence;
						std::vector<int> graphSequence_levels;

						backtrackBookkeeper* current = &start;
						assert(current != 0);

						while(current != 0)
						{
							if(current->usedEdges.size() > 0)
							{
								if(!(current->usedEdges.size() == 1))
								{
									std::cerr << "start_x: " << start_x << "\n";
									std::cerr << "stop_x: " << stop_x << "\n";
									std::cerr << "current: " << current << "\n";
									std::cerr << "current->S: " << current->S << "\n";

									std::cerr << "current->usedEdges.size(): " << current->usedEdges.size() << "\n" << std::flush;
								}
								assert(current->usedEdges.size() == 1);
								assert(current->graphSequence.size() == 1);
								assert(current->graphSequence_levels.size() == 1);

								usedEdges.push_back(current->usedEdges.front());
								graphSequence.append(current->graphSequence.substr(0, 1));
								graphSequence_levels.push_back(current->graphSequence_levels.front());
							}
							current = current->backtrack;
						}

						std::reverse(usedEdges.begin(), usedEdges.end());
						std::reverse(graphSequence.begin(), graphSequence.end());
						std::reverse(graphSequence_levels.begin(), graphSequence_levels.end());

						assert(usedEdges.size() == graphSequence.size());
						assert(graphSequence.size() == graphSequence_levels.size());

						int firstEdge_level = usedEdges.front()->From->level;
						int lastEdge_level = usedEdges.back()->To->level;

						assert(firstEdge_level == start_x);
						assert(lastEdge_level == stop_x);

						backtrackBookkeeper forReturn;
						forReturn.S = start.S;
						forReturn.usedEdges = usedEdges;
						forReturn.graphSequence = graphSequence;
						forReturn.graphSequence_levels = graphSequence_levels;

						return forReturn;
					};

					for(std::map<Node*, backtrackBookkeeper>::iterator startNodeIt = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).begin(); startNodeIt != runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).end(); startNodeIt++)
					{
						Node* startNode = startNodeIt->first;
						Node* stopNode = n;

						int z_startNode = nodesPerLevel_ordered_rev.at(start_x).at(startNode);
						int z_stopNode = nodesPerLevel_ordered_rev.at(stop_x).at(stopNode);

						{
							// start affinely

							graphPointDistance_withBacktrack distanceSpecifier;

							distanceSpecifier.endAffinely = backtrackThroughGraph(runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(stopNode).at(startNode));
							distanceSpecifier.endInAnything = backtrackThroughGraph(runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(stopNode).at(startNode));
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startAffinely.count(z_startNode) == 0) || (distances_startAffinely.at(z_startNode).count(z_stopNode) == 0));
							distances_startAffinely[z_startNode][z_stopNode] = distanceSpecifier;
							if(printMe)
							{
								std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", distances_startAffinely.size() " << distances_startAffinely.size() << " / " << distances_startAffinely.at(z_startNode).size() << "\n" << std::flush;	
							}							
						}

						{
							// start arbitrarily

							graphPointDistance_withBacktrack distanceSpecifier;

							distanceSpecifier.endAffinely = backtrackThroughGraph(runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(stopNode).at(startNode));
							distanceSpecifier.endInAnything = backtrackThroughGraph(runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(stopNode).at(startNode));
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startNormally.count(z_startNode) == 0) || (distances_startNormally.at(z_startNode).count(z_stopNode) == 0));
							distances_startNormally[z_startNode][z_stopNode] = distanceSpecifier;
						}
					}
				}
			}
		}
	
		runningNodeDistances_notStartInAffineGap_endInAnything = &runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_notStartInAffineGap_endInAffineGap = &runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		runningNodeDistances_startInAffineGap_endInAnything	= &runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_startInAffineGap_endInAffineGap = &runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;
	}

	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine(..): Now delete!\n" << std::flush;
	}
	
		runningNodeDistances_notStartInAffineGap_endInAnything_allLevels.clear();
		runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels.clear();
		runningNodeDistances_startInAffineGap_endInAnything_allLevels.clear();
		runningNodeDistances_startInAffineGap_endInAffineGap_allLevels.clear();
	
	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine(..): Done\n" << std::flush;	
	}
		
}



void GraphAlignerUnique::findShortGappedGraphConnection_affine_MTM(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally)
{
	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();
	int thisThread = omp_get_thread_num();
	assert(stop_x > start_x);
	assert(start_x  >= 0);
	assert(stop_x < (int)g->NodesPerLevel.size());

	assert((start_z == -1) || ((start_z >= 0) && (start_z < (int)g->NodesPerLevel.at(start_x).size())));
	assert((stop_z == -1) || ((stop_z >= 0) && (stop_z < (int)g->NodesPerLevel.at(stop_x).size())));

	// bool verbose = true;
	distances_startAffinely.clear();
	distances_startNormally.clear();

	bool printMe = (((stop_x - start_x) > 100000) && (omp_get_thread_num() == 2));
	printMe = false;

	std::ostringstream print_buffer;
	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine_MTM(..): From " << start_x << " to " << stop_x << "\n" << std::flush;
	}


	auto node2Int = [&](Node* n, int level) -> unsigned int {
		assert(level < (int)nodesPerLevel_ordered_rev.size());
		assert(nodesPerLevel_ordered_rev.at(level).count(n));
		return nodesPerLevel_ordered_rev.at(level).at(n);
	};

	// memory pre-allocation
	std::vector<backtrackBookkeeper_UGA> v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels;
	std::vector<backtrackBookkeeper_UGA> v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels;
	std::vector<backtrackBookkeeper_UGA> v_runningNodeDistances_startInAffineGap_endInAnything_allLevels;
	std::vector<backtrackBookkeeper_UGA> v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels;

	std::vector<bool> v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels_assigned;
	std::vector<bool> v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels_assigned;
	std::vector<bool> v_runningNodeDistances_startInAffineGap_endInAnything_allLevels_assigned;
	std::vector<bool> v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels_assigned;


	int memory_levels = ( stop_x - start_x + 1);

	unsigned int memory_keepTrackOfNodes = 0;
	std::set<Node*> keepTrackOfNodes;
	std::set<Node*> nodes_firstLevel = g->NodesPerLevel.at(start_x);
	std::map<Node*, unsigned int> keepTrackNode_2_int;
	std::map<unsigned int, Node*> keepTrackNode_2_int_rev;
	for(std::set<Node*>::iterator nIt = nodes_firstLevel.begin(); nIt != nodes_firstLevel.end(); nIt++)
	{
		Node* n = *nIt;
		int z = nodesPerLevel_ordered_rev.at(start_x).at(n);
		bool takeNode = true;
		if(start_z != -1)
		{
			if(z != start_z)
			{
				takeNode = false;
			}
		}
		if(takeNode)
		{
			keepTrackOfNodes.insert(n);
			keepTrackNode_2_int_rev[memory_keepTrackOfNodes] = n;
			keepTrackNode_2_int[n] = memory_keepTrackOfNodes;
			assert((memory_keepTrackOfNodes >= 0) && (memory_keepTrackOfNodes < nodes_firstLevel.size()));
			memory_keepTrackOfNodes++;
		}
	}

	size_t totalNodes_covered_levels = 0;
	std::vector<size_t> nodesBeforeLevels;
	nodesBeforeLevels.resize(memory_levels, 0);
	for(int lI = start_x; lI <= stop_x; lI++)
	{
		nodesBeforeLevels.at(lI - start_x) = totalNodes_covered_levels;
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);
		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;
			unsigned int nInt = node2Int(n, lI);
			assert((nInt >= 0) && (nInt < nodes_thisLevel.size()));
		}
		totalNodes_covered_levels += nodes_thisLevel.size();
	}

	size_t requiredCells = totalNodes_covered_levels * memory_keepTrackOfNodes;

	// std::cerr << "findShortGappedGraphConnection_affine_MTM(..): " << requiredCells << " required cells.\n" << std::flush;

	v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels.resize(requiredCells);
	v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels.resize(requiredCells);
	v_runningNodeDistances_startInAffineGap_endInAnything_allLevels.resize(requiredCells);
	v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels.resize(requiredCells);

	v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels_assigned.resize(requiredCells, false);
	v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels_assigned.resize(requiredCells, false);
	v_runningNodeDistances_startInAffineGap_endInAnything_allLevels_assigned.resize(requiredCells, false);
	v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels_assigned.resize(requiredCells, false);

	auto getAssignedVector = [&](const std::vector<backtrackBookkeeper_UGA>* runningDistances) -> std::vector<bool>* {
		if(runningDistances == &v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels)
		{
			return &v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels_assigned;
		}
		else if(runningDistances == &v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels)
		{
			return &v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels_assigned;
		}
		else if(runningDistances == &v_runningNodeDistances_startInAffineGap_endInAnything_allLevels)
		{
			return &v_runningNodeDistances_startInAffineGap_endInAnything_allLevels_assigned;
		}
		else
		{
			assert(runningDistances == &v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels);
			return &v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels_assigned;
		}
	};

	auto getRunningValue = [&](const std::vector<backtrackBookkeeper_UGA>& runningDistances, Node* thisLevelNode, Node* loggedNode, int level) -> const backtrackBookkeeper_UGA*  {
		assert((int)thisLevelNode->level == level);
		assert(keepTrackOfNodes.count(loggedNode));

		size_t thisLevelNode_index = nodesBeforeLevels.at(level-start_x) + node2Int(thisLevelNode, level);
		size_t loggedNode_index = keepTrackNode_2_int.at(loggedNode);
		size_t combinedIndex = (thisLevelNode_index * keepTrackNode_2_int.size()) + loggedNode_index;

		assert(combinedIndex >= 0);
		assert(combinedIndex < runningDistances.size());

		std::vector<bool>* assignedNess = getAssignedVector(&runningDistances);
		assert(assignedNess->size() == runningDistances.size());

		assert(assignedNess->at(combinedIndex));

		return &(runningDistances.at(combinedIndex));
	};

	auto haveRunningValue = [&](const std::vector<backtrackBookkeeper_UGA>& runningDistances, Node* thisLevelNode, Node* loggedNode, int level) -> bool {
		assert((int)thisLevelNode->level == level);
		assert(keepTrackOfNodes.count(loggedNode));

		size_t thisLevelNode_index = nodesBeforeLevels.at(level-start_x) + node2Int(thisLevelNode, level);
		size_t loggedNode_index = keepTrackNode_2_int.at(loggedNode);
		size_t combinedIndex = (thisLevelNode_index * keepTrackNode_2_int.size()) + loggedNode_index;

		assert(combinedIndex >= 0);
		assert(combinedIndex < runningDistances.size());

		std::vector<bool>* assignedNess = getAssignedVector(&runningDistances);
		assert(assignedNess->size() == runningDistances.size());

		return (assignedNess->at(combinedIndex));
	};

	auto setRunningValue = [&](std::vector<backtrackBookkeeper_UGA>& runningDistances, Node* thisLevelNode, Node* loggedNode, int level, backtrackBookkeeper_UGA newValue) -> void {
		assert((int)thisLevelNode->level == level);
		assert(keepTrackOfNodes.count(loggedNode));

		size_t thisLevelNode_index = nodesBeforeLevels.at(level-start_x) + node2Int(thisLevelNode, level);
		size_t loggedNode_index = keepTrackNode_2_int.at(loggedNode);
		size_t combinedIndex = (thisLevelNode_index * keepTrackNode_2_int.size()) + loggedNode_index;

		assert(combinedIndex >= 0);
		assert(combinedIndex < runningDistances.size());

		std::vector<bool>* assignedNess = getAssignedVector(&runningDistances);
		assert(assignedNess->size() == runningDistances.size());

		assignedNess->at(combinedIndex) = true;
		runningDistances.at(combinedIndex) = newValue;
	};

	int lastPrintedDistance;
	for(int lI = start_x; lI <= stop_x; lI++)
	{

		if(printMe)
		{
			if((lI == start_x) || ((lI - lastPrintedDistance) >= 50000))
			{
				std::cout << "\t" << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", level " << lI << "\n" << std::flush;
				lastPrintedDistance = lI;
			}
		}
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{

			Node* n = *nIt;
			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == start_x)
			{
				// jump over node if start_z specified
				bool takeNode = true;
				if(start_z != -1)
				{
					if(z != start_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{

					backtrackBookkeeper_UGA D_notStartInAffineGap_endInAnything;
					D_notStartInAffineGap_endInAnything.S = 0;

					backtrackBookkeeper_UGA D_notStartInAffineGap_endInAffineGap;
					D_notStartInAffineGap_endInAffineGap.S = minusInfinity;

					backtrackBookkeeper_UGA D_startInAffineGap_endInAnything;
					D_startInAffineGap_endInAnything.S = 0;

					backtrackBookkeeper_UGA D_startInAffineGap_endInAffineGap;
					D_startInAffineGap_endInAffineGap.S = (lI == 0) ? minusInfinity : 0;

					setRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, n, n, lI, D_notStartInAffineGap_endInAnything);
					setRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, n, n, lI, D_notStartInAffineGap_endInAffineGap);
					setRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, n, n, lI, D_startInAffineGap_endInAnything);
					setRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, n, n, lI, D_startInAffineGap_endInAffineGap);
				}
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;
					for(std::set<Node*>::iterator nodeIt = keepTrackOfNodes.begin(); nodeIt != keepTrackOfNodes.end(); nodeIt++)
					{
						Node* interestingNode = *nodeIt;

						if(! haveRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, nodeFrom, interestingNode, lI - 1))
						{
							continue;
						}

						const backtrackBookkeeper_UGA* previousDistance_notStartInAffineGap_endInAnything = getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, nodeFrom, interestingNode, lI - 1);
						const backtrackBookkeeper_UGA* previousDistance_notStartInAffineGap_endInAffineGap = getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, nodeFrom, interestingNode, lI - 1);
						const backtrackBookkeeper_UGA* previousDistance_startInAffineGap_endInAnything = getRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, nodeFrom, interestingNode, lI - 1);
						const backtrackBookkeeper_UGA* previousDistance_startInAffineGap_endInAffineGap = getRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, nodeFrom, interestingNode, lI - 1);

						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								backtrackBookkeeper_UGA newDistance_notStartInAffineGap_affine_start; // = previousDistance_notStartInAffineGap_endInAnything;
								backtrackBookkeeper_UGA newDistance_notStartInAffineGap_affine_extend; // = previousDistance_notStartInAffineGap_endInAffineGap;

								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start.S =           minusInfinity;

									newDistance_notStartInAffineGap_affine_extend.S =          previousDistance_notStartInAffineGap_endInAffineGap->S;
									newDistance_notStartInAffineGap_affine_extend.backtrack = (backtrackBookkeeper_UGA*)previousDistance_notStartInAffineGap_endInAffineGap;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_notStartInAffineGap_affine_start.S =           previousDistance_notStartInAffineGap_endInAnything->S + (S_openGap + S_extendGap);
									newDistance_notStartInAffineGap_affine_start.backtrack =   (backtrackBookkeeper_UGA*)previousDistance_notStartInAffineGap_endInAnything;
									newDistance_notStartInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);

									newDistance_notStartInAffineGap_affine_extend.S =          previousDistance_notStartInAffineGap_endInAffineGap->S + S_extendGap;
									newDistance_notStartInAffineGap_affine_extend.backtrack =  (backtrackBookkeeper_UGA*)previousDistance_notStartInAffineGap_endInAffineGap;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}

								backtrackBookkeeper_UGA newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start.S > newDistance_notStartInAffineGap_affine_extend.S) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if (
									(! haveRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI)) ||
									(getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI)->S < newDistance_notStartInAffineGap_affine_maximum.S)
								)
								{
									setRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI, newDistance_notStartInAffineGap_affine_maximum);
								}
							}

							// for the ones that did start with an affine gap
							{
								backtrackBookkeeper_UGA newDistance_startInAffineGap_affine_start; // = previousDistance_startInAffineGap_endInAnything;
								backtrackBookkeeper_UGA newDistance_startInAffineGap_affine_extend; // = previousDistance_startInAffineGap_endInAffineGap;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start.S = minusInfinity;

									newDistance_startInAffineGap_affine_extend.S = previousDistance_startInAffineGap_endInAffineGap->S;
									newDistance_startInAffineGap_affine_extend.backtrack = (backtrackBookkeeper_UGA*)previousDistance_startInAffineGap_endInAffineGap;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_startInAffineGap_affine_start.S = previousDistance_startInAffineGap_endInAnything->S + (S_openGap + S_extendGap);
									newDistance_startInAffineGap_affine_start.backtrack = (backtrackBookkeeper_UGA*)previousDistance_startInAffineGap_endInAnything;
									newDistance_startInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);

									newDistance_startInAffineGap_affine_extend.S = previousDistance_startInAffineGap_endInAffineGap->S + S_extendGap;
									newDistance_startInAffineGap_affine_extend.backtrack = (backtrackBookkeeper_UGA*)previousDistance_startInAffineGap_endInAffineGap;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}

								backtrackBookkeeper_UGA newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start.S > newDistance_startInAffineGap_affine_extend.S) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if (
									(! haveRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI)) ||
									(getRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI)->S < newDistance_startInAffineGap_affine_maximum.S)
								)
								{
									setRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI, newDistance_startInAffineGap_affine_maximum);
								}
							}
						}

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								const backtrackBookkeeper_UGA* alternativeScore_endInAffine = getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI);


								backtrackBookkeeper_UGA followGraphGap; // = previousDistance_notStartInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap.S = previousDistance_notStartInAffineGap_endInAnything;
									followGraphGap.S = previousDistance_notStartInAffineGap_endInAnything->S;
									followGraphGap.backtrack = (backtrackBookkeeper_UGA*)previousDistance_notStartInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper_UGA endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine->S) ? followGraphGap : *alternativeScore_endInAffine;


								if (
									(! haveRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, n, interestingNode, lI)) ||
									(getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, n, interestingNode, lI)->S < endInAnything_maximum.S)
								)
								{
									setRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, n, interestingNode, lI, endInAnything_maximum);
								}

							}

							// for the ones that did start with an affine gap
							{
								const backtrackBookkeeper_UGA* alternativeScore_endInAffine = getRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, n, interestingNode, lI);


								backtrackBookkeeper_UGA followGraphGap;// = previousDistance_startInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap = previousDistance_startInAffineGap_endInAnything;
									followGraphGap.S = previousDistance_startInAffineGap_endInAnything->S;
									followGraphGap.backtrack = (backtrackBookkeeper_UGA*)previousDistance_startInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper_UGA endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine->S) ? followGraphGap : *alternativeScore_endInAffine;


								if (
									(! haveRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, n, interestingNode, lI)) ||
									(getRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, n, interestingNode, lI)->S < endInAnything_maximum.S)
								)
								{
									setRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, n, interestingNode, lI, endInAnything_maximum);
								}
							}
						}
					}
				}
			}

			if(lI == stop_x)
			{
				bool takeNode = true;
				if(stop_z != -1)
				{
					if(z != stop_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{

					auto backtrackThroughGraph = [&](backtrackBookkeeper_UGA start) -> backtrackBookkeeper {

						std::vector<Edge*> usedEdges;
						std::string graphSequence;
						std::vector<int> graphSequence_levels;

						int coveredLevels = stop_x - start_x;
						assert(coveredLevels >= 0);
						usedEdges.reserve(coveredLevels + 1);
						graphSequence.reserve(coveredLevels + 1);
						graphSequence_levels.reserve(coveredLevels + 1);

						backtrackBookkeeper_UGA* current = &start;
						assert(current != 0);

						while(current != 0)
						{
							assert((current->takenNodesAndEdges == 0) || (current->takenNodesAndEdges == 1));
							if(current->takenNodesAndEdges == 1)
							{
//								if(!(current->usedEdges.size() == 1))
//								{
//									std::cerr << "start_x: " << start_x << "\n";
//									std::cerr << "stop_x: " << stop_x << "\n";
//									std::cerr << "current: " << current << "\n";
//									std::cerr << "current->S: " << current->S << "\n";
//
//									std::cerr << "current->usedEdges.size(): " << current->usedEdges.size() << "\n" << std::flush;
//								}
//								assert(current->usedEdges.size() == 1);
//								assert(current->graphSequence.size() == 1);
//								assert(current->graphSequence_levels.size() == 1);

								usedEdges.push_back(current->usedEdges);
								graphSequence.push_back(current->graphSequence);
								graphSequence_levels.push_back(current->graphSequence_levels);
							}
							current = current->backtrack;
						}

						std::reverse(usedEdges.begin(), usedEdges.end());
						std::reverse(graphSequence.begin(), graphSequence.end());
						std::reverse(graphSequence_levels.begin(), graphSequence_levels.end());

						assert(usedEdges.size() == graphSequence.size());
						assert(graphSequence.size() == graphSequence_levels.size());

						int firstEdge_level = usedEdges.front()->From->level;
						int lastEdge_level = usedEdges.back()->To->level;

						assert(firstEdge_level == start_x);
						assert(lastEdge_level == stop_x);

						backtrackBookkeeper forReturn;
						forReturn.S = start.S;
						forReturn.usedEdges = usedEdges;
						forReturn.graphSequence = graphSequence;
						forReturn.graphSequence_levels = graphSequence_levels;

						return forReturn;
					};

					for(std::set<Node*>::iterator startNodeIt = keepTrackOfNodes.begin(); startNodeIt != keepTrackOfNodes.end(); startNodeIt++)
					{
						Node* startNode = *startNodeIt;

						if(! haveRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, n, startNode, lI))
						{
							continue;
						}

						Node* stopNode = n;

						int z_startNode = nodesPerLevel_ordered_rev.at(start_x).at(startNode);
						int z_stopNode = nodesPerLevel_ordered_rev.at(stop_x).at(stopNode);

						{
							// start affinely

							graphPointDistance_withBacktrack distanceSpecifier;

							backtrackBookkeeper_UGA firstBt_endAffinely = *( getRunningValue(v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels, stopNode, startNode, lI) );
							backtrackBookkeeper_UGA firstBt_endInAnything = *( getRunningValue(v_runningNodeDistances_startInAffineGap_endInAnything_allLevels, stopNode, startNode, lI) );

							distanceSpecifier.endAffinely = backtrackThroughGraph(firstBt_endAffinely);
							distanceSpecifier.endInAnything = backtrackThroughGraph(firstBt_endInAnything);
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startAffinely.count(z_startNode) == 0) || (distances_startAffinely.at(z_startNode).count(z_stopNode) == 0));
							distances_startAffinely[z_startNode][z_stopNode] = distanceSpecifier;
							if(printMe)
							{
								std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", distances_startAffinely.size() " << distances_startAffinely.size() << " / " << distances_startAffinely.at(z_startNode).size() << "\n" << std::flush;
							}
						}

						{
							// start arbitrarily

							graphPointDistance_withBacktrack distanceSpecifier;

							backtrackBookkeeper_UGA firstBt_endAffinely = *( getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels, stopNode, startNode, lI) );
							backtrackBookkeeper_UGA firstBt_endInAnything = *( getRunningValue(v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels, stopNode, startNode, lI) );

							distanceSpecifier.endAffinely = backtrackThroughGraph(firstBt_endAffinely);
							distanceSpecifier.endInAnything = backtrackThroughGraph(firstBt_endInAnything);
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startNormally.count(z_startNode) == 0) || (distances_startNormally.at(z_startNode).count(z_stopNode) == 0));
							distances_startNormally[z_startNode][z_stopNode] = distanceSpecifier;
						}
					}
				}
			}
		}
	}

	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine_MTM(..): Now delete!\n" << std::flush;
	}

		v_runningNodeDistances_notStartInAffineGap_endInAnything_allLevels.clear();
		v_runningNodeDistances_notStartInAffineGap_endInAffineGap_allLevels.clear();
		v_runningNodeDistances_startInAffineGap_endInAnything_allLevels.clear();
		v_runningNodeDistances_startInAffineGap_endInAffineGap_allLevels.clear();

	if(printMe)
	{
		std::cout << Utilities::timestamp() << "Thread " << omp_get_thread_num() << ", GraphAlignerUnique::findShortGappedGraphConnection_affine_MTM(..): Done\n" << std::flush;
	}

}



void GraphAlignerUnique::kMerEdgeChains2vNW(VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, int>& currentChains_start, std::map<kMerEdgeChain*, NWPath*>& chains2Paths)
{
	std::set<kMerEdgeChain*> currentChains;

	for(unsigned int positionI = 0; positionI < sequencePositions_covered.size(); positionI++)
	{
		kMerEdgeChain* chain = sequencePositions_covered.at(positionI);
		if(chain != 0)
		{
			currentChains.insert(chain);
			if(currentChains_start.count(chain) == 0)
			{
				currentChains_start[chain] = positionI;
			}
		}
	}

	for(std::set<kMerEdgeChain*>::iterator chainIt = currentChains.begin(); chainIt != currentChains.end(); chainIt++)
	{
		kMerEdgeChain* chain = *chainIt;
		NWPath* chainPath = new NWPath();
		int initialY = currentChains_start.at(chain);
		int traversedEdges_nonGap = 0;
		for(unsigned int eI = 0; eI < chain->traversedEdges.size(); eI++)
		{
			Edge* e = chain->traversedEdges.at(eI);
			std::string e_emission = g->CODE.deCode(e->locus_id, e->emission);

			Node* from = e->From;
			Node* to = e->To;

			int y_from = initialY + traversedEdges_nonGap;
			int y_to = y_from + ((e_emission != "_") ? 1 : 0);

			int x_from = from->level;
			int x_to = to->level;
			assert(x_from < x_to);

			int z_from = nodesPerLevel_ordered_rev.at(x_from).at(from);
			int z_to = nodesPerLevel_ordered_rev.at(x_to).at(to);

			bool firstEdge = (eI == 0);
			bool lastEdge = (eI == (chain->traversedEdges.size() - 1));
			int edgeCode = 0;
			if(firstEdge)
			{
				edgeCode = -1;
			}
			if(lastEdge)
			{
				edgeCode = 1;
			}
			assert(!(firstEdge && lastEdge));

			chainPath->createAndAddEdge(x_from, y_from, z_from, x_to, y_to, z_to, e, edgeCode, edgeCode);

			if(e_emission != "_")
			{
				traversedEdges_nonGap++;
			}
		}

		chains2Paths[chain] = chainPath;

		vNW.addPath(chainPath);
	}
}

void GraphAlignerUnique::printSequenceChainCoverageStats(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered)
{
	// this code has problems with first and last positions, might not catch these gaps!
	// print statistics on chain coverage
	long long positions_coveredByChain = 0;
	std::vector<int> gapLengths_sequence;
	std::vector<int> gapLengths_graph;
	long long lastGapStart_sequence = -1;
	kMerEdgeChain* lastChain = 0;
	for(size_t seqI = 0; seqI < sequence.length(); seqI++)
	{
		if(sequencePositions_covered.at(seqI))
		{
			positions_coveredByChain++;
			if(lastGapStart_sequence != -1)
			{
				long long gapLength_sequence = seqI - lastGapStart_sequence;
				kMerEdgeChain* newChain = sequencePositions_covered.at(seqI);

				long long newChain_firstLevel;
				long long lastChain_lastLevel = 0;

				Edge* newChain_firstEdge = newChain->traversedEdges.front();
				newChain_firstLevel = newChain_firstEdge->From->level;

				if(lastChain != 0)
				{
					Edge* lastChain_lastEdge = lastChain->traversedEdges.back();
					lastChain_lastLevel = lastChain_lastEdge->To->level;
				}

				long long gapLength_graph = newChain_firstLevel - lastChain_lastLevel;
				assert(gapLength_sequence > 0);
				assert(gapLength_graph >= 0);

				gapLengths_sequence.push_back(gapLength_sequence);
				gapLengths_graph.push_back(gapLength_graph);
			}

			lastGapStart_sequence = -1;
			lastChain = sequencePositions_covered.at(seqI);
		}
		else
		{
			if(lastGapStart_sequence == -1)
			{
				lastGapStart_sequence = seqI;
			}
		}
	}

	double gapNumber = gapLengths_sequence.size();
	double gapSequenceLengths_sum = 0;
	double gapGraphLengths_sum = 0;
	double gapComplexities_sum = 0;
	double maxGapComplexity = -1;
	for(size_t gapI = 0; gapI < gapLengths_sequence.size(); gapI++)
	{
		double gapComplexity = gapLengths_sequence.at(gapI) * gapLengths_graph.at(gapI);
		gapComplexities_sum += gapComplexity;
		gapSequenceLengths_sum += gapLengths_sequence.at(gapI);
		gapGraphLengths_sum += gapLengths_graph.at(gapI);
		if(gapComplexity > maxGapComplexity)
		{
			maxGapComplexity = gapComplexity;
		}
	}
	
	if(verbose)
	{
		std::cout << "\t" << "Chains covering " << positions_coveredByChain << " / " << sequence.length() << " positions." << "\n" << std::flush;
		std::cout << "\t" << "Number of gaps: " << gapNumber << "\n" << std::flush;
		std::cout << "\t" << "Average complexity: " << ((gapNumber != 0) ? (gapComplexities_sum/gapNumber) : 0) << "\n" << std::flush;
		std::cout << "\t" << "Maximum complexity: " << ((gapNumber != 0) ? maxGapComplexity : 0) << "\n" << std::flush;
		std::cout << "\t" << "Average sequence length: " << ((gapNumber != 0) ? (gapSequenceLengths_sum/gapNumber) : 0) << "\n" << std::flush;
		std::cout << "\t" << "Average graph length: " << ((gapNumber != 0) ? (gapGraphLengths_sum/gapNumber) : 0) << "\n" << std::flush;
	}
}

void GraphAlignerUnique::fixUniqueChains(std::string& sequence, bool thisIterationRandomization, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered, std::set<kMerEdgeChain*>& selectedChains, 	std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, int>& uniquelyTrimmedChains_doubleUniquekMers, bool rescueNonUnique)
{
	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::fixUniqueChains(..): Fix chains, starting with " << uniquelyTrimmedChains_ordered.size() << " chains. Randomization = " << thisIterationRandomization << ".\n" << std::flush;

	std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>> remainingChains = uniquelyTrimmedChains_ordered;
	sequencePositions_covered.resize(sequence.length(), 0);

	while((remainingChains.size() > 0) && ((uniquelyTrimmedChains_doubleUniquekMers.at(*remainingChains.begin()) >= minimumChainUniqueness) || rescueNonUnique))
	{
		if(thisIterationRandomization)
		{
			std::vector<kMerEdgeChain*> candidateChains;
			std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>::iterator candidateIt = remainingChains.begin();
			assert(candidateIt != remainingChains.end());
			for(int chainI = 1; chainI <= randomizationParameter; chainI++)
			{
				kMerEdgeChain* candidateChain  = *candidateIt;
				if((uniquelyTrimmedChains_doubleUniquekMers.at(*remainingChains.begin()) < minimumChainUniqueness) && (! rescueNonUnique))
				{
					break;
				}
				candidateChains.push_back(candidateChain);
				candidateIt++;
				if(candidateIt == remainingChains.end())
					break;
			}
			assert(candidateChains.size() > 0);

			int selectedCandidate = Utilities::randomNumber_nonCritical(candidateChains.size() - 1, &(rng_seeds.at(omp_get_thread_num())));
			kMerEdgeChain* selectedChain = candidateChains.at(selectedCandidate);

//			std::cerr << "\t" << "select " << selectedCandidate << " of " << candidateChains.size() << ": " << selectedChain << "\n" << std::flush;

			remainingChains.erase(selectedChain);
			selectedChains.insert(selectedChain);

			// std::cerr << "selected chain in sequence from " << selectedChain->sequence_begin << " to " << selectedChain->sequence_end << ", inclusive.\n" << std::flush;

			for(int seqI = selectedChain->sequence_begin; seqI <= selectedChain->sequence_end; seqI++)
			{
				assert(sequencePositions_covered.at(seqI) == 0);
				sequencePositions_covered.at(seqI) = selectedChain;
			}

			cleanChainsAccordingToSelectedChain(selectedChain, remainingChains);
		}
		else
		{

			kMerEdgeChain* selectedChain = *(remainingChains.begin());
			remainingChains.erase(selectedChain);
			selectedChains.insert(selectedChain);

			if(verbose)
			{
				
				// std::cerr << "\t" << "select x of " << remainingChains.size() << ": " << selectedChain << "\n" << std::flush;
				std::cout << "fixUniqueChains(..): selected chain in sequence from " << selectedChain->sequence_begin << " to " << selectedChain->sequence_end << ", inclusive.\n" << std::flush;
			}
			
			for(int seqI = selectedChain->sequence_begin; seqI <= selectedChain->sequence_end; seqI++)
			{
				assert(sequencePositions_covered.at(seqI) == 0);
				sequencePositions_covered.at(seqI) = selectedChain;
			}

			cleanChainsAccordingToSelectedChain(selectedChain, remainingChains);

		}
	}

	selectedChainsConsistencyCheck(selectedChains);

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::fixUniqueChains(..): Done." << "\n" << std::flush;
}


void GraphAlignerUnique::cleanChainsAccordingToSelectedChain(kMerEdgeChain* selectedChain, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& availableChains)
{
	std::set<kMerEdgeChain*> chainsToDelete;

	Edge* selectedChainFirstEdge = selectedChain->traversedEdges.front();
	Edge* selectedChainLastEdge = selectedChain->traversedEdges.back();
	int selectedChainFirstGraphLevel = selectedChainFirstEdge->From->level;
	int selectedChainLastGraphLevel = selectedChainLastEdge->To->level;

	// Part 1: search for easy-to-spot incompatibilities

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::cleanChainsAccordingToSelectedChain(..): Scan for coordinate mismatches.\n" << std::flush;

	for(std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>::iterator chainIt = availableChains.begin(); chainIt != availableChains.end(); chainIt++)
	{
		kMerEdgeChain* potentialChain = *chainIt;
		Edge* potentialChainFirstEdge = potentialChain->traversedEdges.front();
		Edge* potentialChainLastEdge = potentialChain->traversedEdges.back();
		int potentialChainLastGraphLevel = potentialChainLastEdge->To->level;
		int potentialChainFirstGraphLevel = potentialChainFirstEdge->From->level;

		bool deleteChain = false;

		// we need to delete all chains which
		// -- cover this chain
		// -- have overlaps with this chain, in either graph or sequence
		// -- are either before or after this chain, but inconsistently so between graph and sequence

		// cover
		if((potentialChain->sequence_begin <= selectedChain->sequence_begin) && (potentialChain->sequence_end >= selectedChain->sequence_end))
		{
			deleteChain = true;
		}
		if((potentialChainFirstGraphLevel <= selectedChainFirstGraphLevel) && (potentialChainLastGraphLevel >= selectedChainLastGraphLevel))
		{
			deleteChain = true;
		}

		// overlap
		if((potentialChain->sequence_begin >= selectedChain->sequence_begin) && (potentialChain->sequence_begin <= selectedChain->sequence_end))
		{
			deleteChain = true;
		}
		if((potentialChain->sequence_end >= selectedChain->sequence_begin) && (potentialChain->sequence_end <= selectedChain->sequence_end))
		{
			deleteChain = true;
		}
		if((potentialChainFirstGraphLevel > selectedChainFirstGraphLevel) && (potentialChainFirstGraphLevel < selectedChainLastGraphLevel))
		{
			deleteChain = true;
		}
		if((potentialChainLastGraphLevel > selectedChainFirstGraphLevel) && (potentialChainLastGraphLevel < selectedChainLastGraphLevel))
		{
			deleteChain = true;
		}

		// now there is neither overlap nor cover
		if(! deleteChain)
		{
			assert((potentialChain->sequence_end < selectedChain->sequence_begin) || (potentialChain->sequence_begin > selectedChain->sequence_end));
			assert((potentialChainLastGraphLevel <= selectedChainFirstGraphLevel) || (potentialChainFirstGraphLevel >= selectedChainLastGraphLevel));
			if(potentialChain->sequence_end < selectedChain->sequence_begin)
			{
				if(!(potentialChainLastGraphLevel <= selectedChainFirstGraphLevel))
				{
					deleteChain = true;
				}
			}
			else if(potentialChain->sequence_begin > selectedChain->sequence_end)
			{
				if(!(potentialChainFirstGraphLevel >= selectedChainLastGraphLevel))
				{
					deleteChain = true;
				}
			}
			else
			{
				assert( 1 == 0);
			}
		}

		if(deleteChain)
		{
			chainsToDelete.insert(potentialChain);
		}
	}

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::cleanChainsAccordingToSelectedChain(..): Done. Now delete " << chainsToDelete.size() << " overlapping / covering / incompatible chains.\n" << std::flush;

	for(std::set<kMerEdgeChain*>::iterator deleteChainIt = chainsToDelete.begin(); deleteChainIt != chainsToDelete.end(); deleteChainIt++)
	{
		kMerEdgeChain* chainToDelete = *deleteChainIt;
		availableChains.erase(chainToDelete);
	}

	chainsToDelete.clear();
	// Part 2: search for graph-induced incompatibilities
	// all remaining chains are either to the left or to the right of the selected chain

	// we expand from the right side until we have reached a graph region with 1 level - from there onwards we can reach anything.
	// level for level, we keep track of where we can go - we check all other edges, and if one edge is the start or stop
	// edge of a chain, we kick out the chain.

	std::map<Edge*, std::set<kMerEdgeChain*>> remainingChains_starts;
	std::map<Edge*, std::set<kMerEdgeChain*>> remainingChains_stops;
	for(std::set<kMerEdgeChain*>::iterator remainingChainIt = availableChains.begin(); remainingChainIt != availableChains.end(); remainingChainIt++)
	{
		kMerEdgeChain* chain = *remainingChainIt;
		Edge* firstEdge = chain->traversedEdges.front();
		Edge* lastEdge = chain->traversedEdges.back();
		remainingChains_starts[firstEdge].insert(chain);
		remainingChains_stops[lastEdge].insert(chain);
	}

	// scan to the right
	std::map<Node*, bool> runningReachableNodes;
	bool _rightScan_haveInitialNode = false;
	for(int lI = selectedChainLastGraphLevel; lI < (int)g->NodesPerLevel.size(); lI++)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, bool> runningReachableNodes_thisLevel;
		std::set<Edge*> edgesToThisLevel;

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(lI == selectedChainLastGraphLevel)
			{
				if(n == selectedChainLastEdge->To)
				{
					runningReachableNodes_thisLevel[n] = true;
					_rightScan_haveInitialNode = true;
				}
				else
				{
					runningReachableNodes_thisLevel[n] = false;
				}
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);
				runningReachableNodes_thisLevel[n] = false;

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					Node* nodeFrom = e->From;

					assert(runningReachableNodes.count(nodeFrom) > 0);
					runningReachableNodes_thisLevel.at(n) = (runningReachableNodes_thisLevel.at(n) || runningReachableNodes.at(nodeFrom));

					edgesToThisLevel.insert(e);
				}
			}
		}

		for(std::set<Edge*>::iterator edgeIt = edgesToThisLevel.begin(); edgeIt != edgesToThisLevel.end(); edgeIt++)
		{
			Edge* e = *edgeIt;
			Node* fromNode = e->From;
			if(remainingChains_starts.count(e))
			{
				if(runningReachableNodes.at(fromNode))
				{
					// OK, reachable!
				}
				else
				{
					// NOT reachable!
					std::set<kMerEdgeChain*> startingChains = remainingChains_starts.at(e);
					for(std::set<kMerEdgeChain*>::iterator startChainIt = startingChains.begin(); startChainIt != startingChains.end(); startChainIt++)
					{
						kMerEdgeChain* chain = *startChainIt;
						chainsToDelete.insert(chain);
					}
				}
			}
			if(remainingChains_stops.count(e))
			{
				if(runningReachableNodes.at(fromNode))
				{
					// OK, reachable!
				}
				else
				{
					// NOT reachable!
					std::set<kMerEdgeChain*> stoppingChains = remainingChains_stops.at(e);
					for(std::set<kMerEdgeChain*>::iterator stopChainIt = stoppingChains.begin(); stopChainIt != stoppingChains.end(); stopChainIt++)
					{
						kMerEdgeChain* chain = *stopChainIt;
						assert(chainsToDelete.count(chain));
					}
				}
			}
		}

		runningReachableNodes = runningReachableNodes_thisLevel;

		if(g->NodesPerLevel.at(lI).size() == 1)
		{
			break;
		}
	}
	assert(_rightScan_haveInitialNode);

	// scan to the left
	runningReachableNodes.clear();
	bool _leftScan_haveInitialNode = false;
	for(int lI = selectedChainFirstGraphLevel; lI >= 0; lI--)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, bool> runningReachableNodes_thisLevel;
		std::set<Edge*> edgesToThisLevel;

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(lI == selectedChainFirstGraphLevel)
			{
				if(n == selectedChainFirstEdge->From)
				{
					runningReachableNodes_thisLevel[n] = true;
					_leftScan_haveInitialNode = true;
				}
				else
				{
					runningReachableNodes_thisLevel[n] = false;
				}
			}
			else
			{
				std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
				assert(outgoingEdges.size() > 0);
				runningReachableNodes_thisLevel[n] = false;

				for(std::set<Edge*>::iterator eIt = outgoingEdges.begin(); eIt != outgoingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					Node* nodeTo = e->To;

					assert(runningReachableNodes.count(nodeTo) > 0);
					runningReachableNodes_thisLevel.at(n) = (runningReachableNodes_thisLevel.at(n) || runningReachableNodes.at(nodeTo));

					edgesToThisLevel.insert(e);
				}
			}
		}

		for(std::set<Edge*>::iterator edgeIt = edgesToThisLevel.begin(); edgeIt != edgesToThisLevel.end(); edgeIt++)
		{
			Edge* e = *edgeIt;
			Node* toNode = e->To;

			if(remainingChains_starts.count(e))
			{
				if(runningReachableNodes.at(toNode))
				{
					// OK, reachable!
				}
				else
				{
					// NOT reachable!
					std::set<kMerEdgeChain*> startingChains = remainingChains_starts.at(e);
					for(std::set<kMerEdgeChain*>::iterator startChainIt = startingChains.begin(); startChainIt != startingChains.end(); startChainIt++)
					{
						kMerEdgeChain* chain = *startChainIt;
						assert(chainsToDelete.count(chain));
					}
				}
			}
			if(remainingChains_stops.count(e))
			{
				if(runningReachableNodes.at(toNode))
				{
					// OK, reachable!
				}
				else
				{
					// NOT reachable!
					std::set<kMerEdgeChain*> stoppingChains = remainingChains_stops.at(e);
					for(std::set<kMerEdgeChain*>::iterator stopChainIt = stoppingChains.begin(); stopChainIt != stoppingChains.end(); stopChainIt++)
					{
						kMerEdgeChain* chain = *stopChainIt;
						chainsToDelete.insert(chain);
					}
				}
			}

//			if(remainingChains_starts.count(e) || remainingChains_stops.count(e))
//			{
//				if(runningReachableNodes.at(toNode))
//				{
//					// OK, reachable!
//				}
//				else
//				{
//					// NOT reachable!
//					if(!(remainingChains_stops.count(e) || chainsToDelete.count(remainingChains_starts.at(e))))
//					{
//						std::cerr << "! (remainingChains_stops.count(e) || chainsToDelete.count(remainingChains_starts.at(e)))" << "\n";
//						std::cerr << "Level: " << lI << "\n";
//						std::cerr << "selectedChainFirstGraphLevel: " << selectedChainFirstGraphLevel << "\n";
//						std::cerr << "Edge: " << e << "\n";
//						kMerEdgeChain* chain = remainingChains_starts.at(e);
//						std::cerr << "Chain: " << chain << "\n";
//						std::cerr << "chain->traversedEdges.at(0)->From->level: " << chain->traversedEdges.at(0)->From->level << "\n";
//						std::cerr << "chain->traversedEdges.back()->To->level: " << chain->traversedEdges.back()->To->level << "\n" << std::flush;
//
//					}
//					assert(remainingChains_stops.count(e) || chainsToDelete.count(remainingChains_starts.at(e)));
//					kMerEdgeChain* chainToDelete = remainingChains_stops.count(e) ? remainingChains_stops.at(e) : remainingChains_starts.at(e);
//
//					chainsToDelete.insert(chainToDelete);
//					std::cout << "\t\tDELETE\n" << std::flush;
//				}
//			}
		}

		runningReachableNodes = runningReachableNodes_thisLevel;

		if(g->NodesPerLevel.at(lI).size() == 1)
		{
			break;
		}
	}
	assert(_leftScan_haveInitialNode);

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::cleanChainsAccordingToSelectedChain(..): Done. Now delete " << chainsToDelete.size() << " because of graph structure.\n" << std::flush;

	for(std::set<kMerEdgeChain*>::iterator deleteChainIt = chainsToDelete.begin(); deleteChainIt != chainsToDelete.end(); deleteChainIt++)
	{
		kMerEdgeChain* chainToDelete = *deleteChainIt;
		availableChains.erase(chainToDelete);
	}

	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::cleanChainsAccordingToSelectedChain(..): Exit function.\n" << std::flush;

}

void GraphAlignerUnique::selectedChainsConsistencyCheck(std::set<kMerEdgeChain*> selectedChains)
{
	if(verbose) std::cout << Utilities::timestamp() << "GraphAlignerUnique::selectedChainsConsistencyCheck: Enter function on " << selectedChains.size() << " chains.\n" << std::flush;

	std::vector<kMerEdgeChain*> selectedChains_in_order(selectedChains.begin(), selectedChains.end());
	std::sort(selectedChains_in_order.begin(), selectedChains_in_order.end(), [&](kMerEdgeChain* one, kMerEdgeChain* two){
		return (one->sequence_begin < two->sequence_begin);
	});

	if(selectedChains_in_order.size() > 1)
	{
		assert(selectedChains_in_order.at(0)->sequence_begin <= selectedChains_in_order.at(1)->sequence_begin);
	}

	std::set<Node*> interestingNodes;
	for(unsigned int chainI = 0; chainI < selectedChains_in_order.size(); chainI++)
	{
		kMerEdgeChain* chain = selectedChains_in_order.at(chainI);
		assert(chain->sequence_begin <= chain->sequence_end);
		if(chainI < (selectedChains_in_order.size() - 1))
		{
			kMerEdgeChain* nextChain = selectedChains_in_order.at(chainI+1);
			assert(chain->sequence_end < nextChain->sequence_begin);
		}

		Node* chain_from = chain->traversedEdges.front()->From;
		Node* chain_to = chain->traversedEdges.back()->To;

		interestingNodes.insert(chain_from);
		interestingNodes.insert(chain_to);
	}

	int chainToCheckIndex = 0;

	std::map<Node*, std::set<Node*> > runningReachableNodes;
	std::set<Node*> nodes_l0 = g->NodesPerLevel.at(0);
	for(unsigned int lI = 0; lI < g->NodesPerLevel.size(); lI++)
	{
		if((lI % 1000) == 0)
		{
			// std::cout << "\r" << "Level " << lI << " / " << g->NodesPerLevel.size() << std::flush;
		}

		if(chainToCheckIndex >= (int)selectedChains_in_order.size())
		{
			break;
		}

		Node* currentChainEntryNode = selectedChains_in_order.at(chainToCheckIndex)->traversedEdges.at(0)->From;
		Node* currentChainExitNode = selectedChains_in_order.at(chainToCheckIndex)->traversedEdges.back()->To;
		assert(currentChainEntryNode != currentChainExitNode);

		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);
		std::map<Node*, std::set<Node*> > runningReachableNodes_thisLevel;

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(lI == 0)
			{
				runningReachableNodes_thisLevel[n] = std::set<Node*>();
				runningReachableNodes_thisLevel[n].insert(0);
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);
				runningReachableNodes_thisLevel[n] = std::set<Node*>();

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					Node* nodeFrom = e->From;
					assert(runningReachableNodes.count(nodeFrom) > 0);
					runningReachableNodes_thisLevel.at(n).insert(runningReachableNodes.at(nodeFrom).begin(), runningReachableNodes.at(nodeFrom).end());
				}
			}

			if(interestingNodes.count(n))
			{
				runningReachableNodes_thisLevel.at(n).insert(n);
			}

			if(n == currentChainEntryNode)
			{
				if(chainToCheckIndex == 0)
				{
					assert(runningReachableNodes_thisLevel.at(n).count(0));
				}
				else
				{
					Node* previousChainExitNode = selectedChains_in_order.at(chainToCheckIndex-1)->traversedEdges.back()->To;
					assert(runningReachableNodes_thisLevel.at(n).count(previousChainExitNode));
				}
			}
			if(n == currentChainExitNode)
			{
				assert(runningReachableNodes_thisLevel.at(n).count(currentChainEntryNode));
				chainToCheckIndex++;
			}
		}

		runningReachableNodes = runningReachableNodes_thisLevel;

		if(lI == currentChainEntryNode->level)
		{
			for(std::map<Node*, std::set<Node*> >::iterator targetNodeIt = runningReachableNodes.begin(); targetNodeIt != runningReachableNodes.end(); targetNodeIt++)
			{
				std::set<Node*> nodesToRemove;
				for(std::set<Node*>::iterator fromNodeIt = targetNodeIt->second.begin(); fromNodeIt != targetNodeIt->second.end(); fromNodeIt++)
				{
					Node* fromNode = *fromNodeIt;
					if(fromNode != currentChainEntryNode)
					{
						nodesToRemove.insert(fromNode);
					}
				}

				for(std::set<Node*>::iterator deleteNodeIt = nodesToRemove.begin(); deleteNodeIt != nodesToRemove.end(); deleteNodeIt++)
				{
					targetNodeIt->second.erase(*deleteNodeIt);
				}
			}
		}
	}
	assert(chainToCheckIndex == (int)selectedChains_in_order.size());
	if(verbose) std::cout << "\n";
}


std::vector<std::pair<int, Edge*> > GraphAlignerUnique::_graph_get_previous_z_values_and_edges(int x, int z)
{
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x > 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> incomingEdges = thisZ->Incoming_Edges;
	set<Node*> nodesPreviousLevel;
	assert(incomingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* fromNode = e->From;
		int z_for_fromNode = nodesPerLevel_ordered_rev.at(x-1).at(fromNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_fromNode, e));
	}

	return forReturn;
}

std::vector<std::pair<int, Edge*> > GraphAlignerUnique::_graph_get_next_z_values_and_edges(int x, int z)
{
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x >= 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> outgoingEdges = thisZ->Outgoing_Edges;
	set<Node*> nodesNextLevel;
	assert(outgoingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = outgoingEdges.begin(); eIt != outgoingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* toNode = e->To;
		int z_for_toNode = nodesPerLevel_ordered_rev.at(x+1).at(toNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_toNode, e));
	}

	return forReturn;
}

std::vector<int> GraphAlignerUnique::_graph_get_previous_z_values(int x, int z)
{
	std::vector<std::pair<int, Edge*> > previousZs = _graph_get_previous_z_values_and_edges(x, z);

	std::set<int> forReturn;
	for(unsigned int i = 0; i < previousZs.size(); i++)
	{
		int previousZ = previousZs.at(i).first;
		forReturn.insert(previousZ);
	}

	return std::vector<int>(forReturn.begin(), forReturn.end());
}





double GraphAlignerUnique::score(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence)
{
	assert(reconstructed_graph.length() == reconstructed_sequence.length());
	assert(reconstructed_graph_levels.size() == reconstructed_graph.length());
	assert(S_graphGap == 0);

	double score = 0;
	bool inAffineGapInGraph = false;
	bool inAffineGapInSequence = false;
	bool verbose = false;

	std::vector<double> runningScores;

	std::vector<bool> endsStatus_sequence;
	endsStatus_sequence.resize(reconstructed_graph.length(), false);

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::string G = reconstructed_graph.substr(i, 1);

		if((S == "_") && (endsStatus_sequence.at(i) == true))
		{
			score += 0;
			inAffineGapInSequence = false;
			inAffineGapInGraph = false;
		}
		else
		{
			if(S != "_")
			{
				inAffineGapInSequence = false;
			}
			if(G != "_")
			{
				inAffineGapInGraph = false;
			}

			if(G == "_")
			{
				if(reconstructed_graph_levels.at(i) != -1)
				{
					inAffineGapInGraph = false;
					// inAffineGapInSequence = false;

					if(S == "_")
					{
						score += S_graphGap;
					}
					else
					{
						score += S_mismatch;
					}
				}
				else
				{
					assert(S != "_");
					if(inAffineGapInGraph == true)
					{
						score += S_extendGap;
					}
					else
					{
						inAffineGapInGraph = true;
						score += (S_openGap + S_extendGap);
					}
				}
			}
			else
			{
				if(S == "_")
				{
					if(inAffineGapInSequence == true)
					{
						score += S_extendGap;
					}
					else
					{
						inAffineGapInSequence = true;
						score += (S_openGap + S_extendGap);
					}
				}
				else
				{
					if(S == G)
					{
						score += S_match;
					}
					else
					{
						score += S_mismatch;
					}
				}
			}
		}
		runningScores.push_back(score);
	}


	if(verbose)
	{
		for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
		{
			int l = reconstructed_graph_levels.at(i);
			if(l != -1)
				l++;
			std::cout << Utilities::ItoStr(l) << "\t";
		}
		std::cout << "\n";

		for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
		{
			std::string G = reconstructed_graph.substr(i, 1);
			std::cout << G << "\t";
		}
		std::cout << "\n";

		for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
		{
			std::string S = reconstructed_sequence.substr(i, 1);
			std::cout << S << "\t";
		}
		std::cout << "\n";

		for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
		{
			double S = runningScores.at(i);
			std::cout << S << "\t";
		}
		std::cout << "\n" << std::flush;
	}


	return score;
}


std::vector<localExtension_pathDescription> GraphAlignerUnique::fullNeedleman_diagonal_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive, bool returnGlobalScore)
{
	if(directionPositive)
	{
		if(!((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph)))
		{
			std::cerr << "maxLevel_graph: " << maxLevel_graph << "\n";
			std::cerr << "startLevel_graph: " << startLevel_graph << "\n" << std::flush;
		}
		if(!((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence)))
		{
			std::cerr << "maxPosition_sequence: " << maxPosition_sequence << "\n";
			std::cerr << "start_sequence: " << start_sequence << "\n" << std::flush;
		}
		assert((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence));
	}
	else
	{
		assert((maxLevel_graph == -1) || (maxLevel_graph < startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence < start_sequence));
	}

	double minusInfinity = -1 * numeric_limits<double>::max();

	if(returnGlobalScore)
	{
		diagonal_stop_threshold = minusInfinity;
		assert(maxLevel_graph == -1);
	}

	class mScore {
	public:
		double D;
		double GraphGap;
		double SequenceGap;
	};

	class mScore_backtrace {
	public:
		backtraceStep_affine D;
		backtraceStep_affine GraphGap;
		backtraceStep_affine SequenceGap;
	};

	class mScore_alternatives {
	public:
		std::vector<double> D;
		std::vector<double> GraphGap;
		std::vector<double> SequenceGap;
	};

	class mScore_backtrace_alternatives {
	public:
		std::vector<backtraceStep_affine> D;
		std::vector<backtraceStep_affine> GraphGap;
		std::vector<backtraceStep_affine> SequenceGap;
	};

	std::map<int, std::map<int, std::map<int, mScore>> > scores;
	std::map<int, std::map<int, std::map<int, mScore_backtrace>> > scores_backtrace;

	std::vector<std::vector<int> > coordinates_for_backtracking;
	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	if(verbose)
	{
		std::cout << "fullNeedleman_affine_diagonal_extension(..) called, direction " << directionPositive << ".\n" << std::flush;
	}

	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	int diagonals = sequenceLength + levels - 1;
	int max_levelI = levels - 1;
	int max_seqI = sequenceLength;
	int min_levelI = 0;
	int min_seqI = 0;

	if(maxLevel_graph != -1)
	{
		if(directionPositive)
		{
			max_levelI = maxLevel_graph;
		}
		else
		{
			min_levelI = maxLevel_graph;
		}
	}
	if(maxPosition_sequence != -1)
	{
		if(directionPositive)
		{
			max_seqI = maxPosition_sequence;
		}
		else
		{
			min_seqI = maxPosition_sequence;
		}
	}

	assert(startLevel_graph >= min_levelI);
	assert(start_sequence >= min_seqI);
	assert(startLevel_graph <= max_levelI);
	assert(start_sequence <= max_seqI);


	assert(min_levelI >= 0);
	assert(max_levelI <= (levels - 1));
	assert(max_levelI > min_levelI);

	assert(min_seqI >= 0);
	assert(max_seqI <= (int)sequenceLength);
	assert(max_seqI > min_seqI);

	double currentMaximum = 0;
	std::vector<std::vector<int> > currentMaxima_coordinates;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(startLevel_graph).size();
	assert((startZ_graph >= 0) && (startZ_graph < (int)statesPerLevel0));


	// parameters
	// threshold_for_filtering: will remove all cells from NW table in a given diagonal which have value > 15 difference from maximum
	int threshold_for_filtering = 15;
	int maximum_steps_nonIncrease = 40;


	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		if((int)stateI == startZ_graph)
		{
			scores[startLevel_graph][start_sequence][stateI].D = 0;
			scores_backtrace[startLevel_graph][start_sequence][stateI] = mScore_backtrace();
		}
		else
		{
			scores[startLevel_graph][start_sequence][stateI].D = minusInfinity;
		}

		scores[startLevel_graph][start_sequence][stateI].GraphGap = minusInfinity;
		scores[startLevel_graph][start_sequence][stateI].SequenceGap = minusInfinity;

		if((int)stateI == startZ_graph)
		{
			std::vector<int> existingCoordinates;
			existingCoordinates.push_back(startLevel_graph);
			existingCoordinates.push_back(start_sequence);
			existingCoordinates.push_back(stateI);
			m1_diagonal.push_back(existingCoordinates);
			currentMaxima_coordinates.push_back(existingCoordinates);
		}
	}

	std::map<NWPath*, std::pair<double, std::vector<int> > > hit_NW_paths;

//	std::cerr << "fullNeedleman_diagonal_extension:\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tdirectionPositive: " << directionPositive << "\n" << std::flush;

	int lastMaximumIncrease_at_diagonalI = 0;


	for(int diagonalI = 1; diagonalI <= diagonals; diagonalI++)
	{

		if(verbose)
		{
			std::cout << "\t diagonalI " << diagonalI << "/" << diagonals << ".\n" << std::flush;
		}


//		int scores_size = 0;
//		for(std::map<int, std::map<int, std::map<int, mScore>> >::iterator it1 = scores.begin(); it1 != scores.end(); it1++)
//		{
//			for(std::map<int, std::map<int, mScore>>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
//			{
//				for(std::map<int, mScore>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
//				{
//					scores_size++;
//				}
//			}
//		}
		//std::cerr << "\t\tdiagonalI = " << diagonalI << " => scores_size: " << scores_size << "\n" << std::flush;

		if((diagonalI - lastMaximumIncrease_at_diagonalI) > maximum_steps_nonIncrease)
		{
			break;
		}

		std::map<int, std::map<int, std::map<int, mScore_alternatives > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, mScore_backtrace_alternatives > > > thisDiagonal_backtrace;

		if(verbose)
			std::cout << "\t\tfrom m-2 diagonal" << "\n" << std::flush;

		// extend from m-2 diagonal
		for(int m2I = 0; m2I < (int)m2_diagonal.size(); m2I++)
		{
			std::vector<int>& previous_coordinates = m2_diagonal.at(m2I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			int next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);

			if((next_levelI > max_levelI) || (next_seqI > max_seqI))
				continue;

			if((next_levelI < min_levelI) || (next_seqI < min_seqI))
				continue;

			std::string sequenceEmission = (directionPositive ? sequence.substr(previous_seqI, 1) : sequence.substr(previous_seqI-1, 1));

			std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
			assert(nextZs.size() > 0);

			for(unsigned int zI = 0; zI < nextZs.size(); zI++)
			{
				std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
				int next_stateI = thisZjump.first;

				std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
				assert(edgeEmission.size() == 1);

				double score_MatchMismatch = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep_affine backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;
				backtrack_MatchMismatch.sourceMatrix = 0;

				thisDiagonal[next_levelI][next_seqI][next_stateI].D.push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].D.push_back(backtrack_MatchMismatch);
			}
		}

		if(verbose)
			std::cout << "\t\tfrom m-1 diagonal" << "\n" << std::flush;

		// extend from m-1 diagonal
		for(int m1I = 0; m1I < (int)m1_diagonal.size(); m1I++)
		{
			std::vector<int>& previous_coordinates = m1_diagonal.at(m1I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			// gap in graph
			int gapInGraph_next_levelI = previous_coordinates.at(0);
			int gapInGraph_next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);
			if	(
					(directionPositive && (gapInGraph_next_levelI <= max_levelI) && (gapInGraph_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInGraph_next_levelI >= min_levelI) && (gapInGraph_next_seqI >= min_seqI))
				)
			{
				double score_gapInGraph_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_open = minusInfinity;
				}
				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep_affine backtrack_gapInGraph_open;
				backtrack_gapInGraph_open.x = previous_levelI;
				backtrack_gapInGraph_open.y = previous_seqI;
				backtrack_gapInGraph_open.z = previous_stateI;
				backtrack_gapInGraph_open.usedEdge = 0;
				backtrack_gapInGraph_open.sourceMatrix = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_open);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_open);

				double score_gapInGraph_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).GraphGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_extend = minusInfinity;
				}

				backtraceStep_affine backtrack_gapInGraph_extend;
				backtrack_gapInGraph_extend.x = previous_levelI;
				backtrack_gapInGraph_extend.y = previous_seqI;
				backtrack_gapInGraph_extend.z = previous_stateI;
				backtrack_gapInGraph_extend.usedEdge = 0;
				backtrack_gapInGraph_extend.sourceMatrix = 1;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_extend);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_extend);
			}

			// gap in sequence
			int gapInSequence_next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int gapInSequence_next_seqI = previous_coordinates.at(1);

			if	(
					(directionPositive && (gapInSequence_next_levelI <= max_levelI) && (gapInSequence_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInSequence_next_levelI >= min_levelI) && (gapInSequence_next_seqI >= min_seqI))
				)
			{
				std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
				assert(nextZs.size() > 0);

				for(unsigned int zI = 0; zI < nextZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
					int next_stateI = thisZjump.first;

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					// open sequence gap

					double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
					if(edgeEmission == "_")
					{
						score_gapInSequence_open = minusInfinity;
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_open = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_open;
					backtrack_gapInSequence_open.x = previous_levelI;
					backtrack_gapInSequence_open.y = previous_seqI;
					backtrack_gapInSequence_open.z = previous_stateI;
					backtrack_gapInSequence_open.usedEdge = thisZjump.second;
					backtrack_gapInSequence_open.sourceMatrix = 0;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);

					// extend sequence gap

					double score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap + S_extendGap;
//					if(edgeEmission == "_")
//					{
//						score_gapInSequence_extend = minusInfinity;
//					}
					if(edgeEmission == "_")
					{
						if(scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap == minusInfinity)
						{
							score_gapInSequence_extend = minusInfinity;
						}
						else
						{
							score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap +  S_graphGap;
						}
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_extend = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_extend;
					backtrack_gapInSequence_extend.x = previous_levelI;
					backtrack_gapInSequence_extend.y = previous_seqI;
					backtrack_gapInSequence_extend.z = previous_stateI;
					backtrack_gapInSequence_extend.usedEdge = thisZjump.second;
					backtrack_gapInSequence_extend.sourceMatrix = 2;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_extend);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_extend);

					// non-affine sequence gap
					if(edgeEmission == "_")
					{
						double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_graphGap;

						backtraceStep_affine backtrack_gapInSequence_nonAffine;
						backtrack_gapInSequence_nonAffine.x = previous_levelI;
						backtrack_gapInSequence_nonAffine.y = previous_seqI;
						backtrack_gapInSequence_nonAffine.z = previous_stateI;
						backtrack_gapInSequence_nonAffine.usedEdge = thisZjump.second;
						backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(score_gapInSequence_nonAffine);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(backtrack_gapInSequence_nonAffine);
					}
				}
			}
		}

		if(verbose)
			std::cout << "\t\tmaximal" << "\n" << std::flush;

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, mScore_alternatives > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, mScore_alternatives > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, mScore_alternatives >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;

					// call maxima for GraphGap
					std::vector<double>& scores_GraphGap = diagIt2->second.at(stateI).GraphGap;
					std::vector<backtraceStep_affine>& scores_GraphGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).GraphGap;
					double selectedScore_GraphGap;
					backtraceStep_affine selectedStep_GraphGap;

					if(scores_GraphGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_GraphGap = Utilities::findVectorMaxP_nonCritical(scores_GraphGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_GraphGap = maxScore_GraphGap.first;
						selectedStep_GraphGap = scores_GraphGap_bt.at(maxScore_GraphGap.second);
					}
					else
					{
						selectedScore_GraphGap = minusInfinity;
					}

					// call maxima for SequenceGap
					std::vector<double>& scores_SequenceGap = diagIt2->second.at(stateI).SequenceGap;
					std::vector<backtraceStep_affine>& scores_SequenceGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).SequenceGap;

					double selectedScore_SequenceGap;
					backtraceStep_affine selectedStep_SequenceGap;
					if(scores_SequenceGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_SequenceGap = Utilities::findVectorMaxP_nonCritical(scores_SequenceGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_SequenceGap = maxScore_SequenceGap.first;
						selectedStep_SequenceGap = scores_SequenceGap_bt.at(maxScore_SequenceGap.second);
					}
					else
					{
						selectedScore_SequenceGap = minusInfinity;
					}

					// two additional steps for D, jumping from the two gap matrices

					double score_intoD_fromGraphGap = selectedScore_GraphGap;
					backtraceStep_affine step_intoD_fromGraphGap;
					step_intoD_fromGraphGap.x = levelI;
					step_intoD_fromGraphGap.y = seqI;
					step_intoD_fromGraphGap.z = stateI;
					step_intoD_fromGraphGap.sourceMatrix = 1;
					step_intoD_fromGraphGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromGraphGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromGraphGap);

					double score_intoD_fromSequenceGap = selectedScore_SequenceGap;
					backtraceStep_affine step_intoD_fromSequenceGap;
					step_intoD_fromSequenceGap.x = levelI;
					step_intoD_fromSequenceGap.y = seqI;
					step_intoD_fromSequenceGap.z = stateI;
					step_intoD_fromSequenceGap.sourceMatrix = 2;
					step_intoD_fromSequenceGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromSequenceGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromSequenceGap);

					// find final maximum for D
					std::vector<double>& scores_D = diagIt2->second.at(stateI).D;
					std::vector<backtraceStep_affine>& scores_D_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D;
					assert(scores_D.size() == scores_D_bt.size());

					std::pair<double, unsigned int> maxScore_D = Utilities::findVectorMaxP_nonCritical(scores_D, &(rng_seeds.at(omp_get_thread_num())));

					bool blockOutCell = false;
					if(directionPositive)
					{
						std::set<NWEdge*> emanatingNWEdges = blockedPathsTable->getEdgesEmanatingFrom(levelI, seqI, stateI);
						if((emanatingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = emanatingNWEdges.begin(); edgeIt != emanatingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}
					else
					{
						std::set<NWEdge*> incomingNWEdges = blockedPathsTable->getEdgesGoingInto(levelI, seqI, stateI);
						if((incomingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = incomingNWEdges.begin(); edgeIt != incomingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}

					if(blockOutCell)
					{
						double scoreThatWillBeDeleted = maxScore_D.first;
						std::vector<int> thisCoordinates;
						thisCoordinates.push_back(levelI);
						thisCoordinates.push_back(seqI);
						thisCoordinates.push_back(stateI);
						if(scoreThatWillBeDeleted == currentMaximum)
						{
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}
						else if(scoreThatWillBeDeleted > currentMaximum)
						{
							currentMaximum = scoreThatWillBeDeleted;
							currentMaxima_coordinates.clear();
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}

						scores[levelI][seqI][stateI].D = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

						assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

						scores[levelI][seqI][stateI].GraphGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

						scores[levelI][seqI][stateI].SequenceGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;
					}
					else
					{
						if(maxScore_D.first >= diagonal_stop_threshold)
						{
							scores[levelI][seqI][stateI].D = maxScore_D.first;
							scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

							assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

							scores[levelI][seqI][stateI].GraphGap = selectedScore_GraphGap;
							scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

							scores[levelI][seqI][stateI].SequenceGap = selectedScore_SequenceGap;
							scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;

							std::vector<int> thisCoordinates;
							thisCoordinates.push_back(levelI);
							thisCoordinates.push_back(seqI);
							thisCoordinates.push_back(stateI);
							m_thisDiagonal.push_back(thisCoordinates);

							backtraceStep_affine oneRealStepBackwards = scores_backtrace[levelI][seqI][stateI].D;
							while((oneRealStepBackwards.x == levelI) && (oneRealStepBackwards.y == seqI))
							{
								assert(oneRealStepBackwards.sourceMatrix != 0);
								if(oneRealStepBackwards.sourceMatrix == 1)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
								}
								else if(oneRealStepBackwards.sourceMatrix == 2)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
								}
								else
								{
									assert(1 == 0);
								}
							}
							int previousScore;
							if(oneRealStepBackwards.sourceMatrix == 0)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].D;
							}
							else if(oneRealStepBackwards.sourceMatrix == 1)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
							}
							else if(oneRealStepBackwards.sourceMatrix == 2)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
							}
							else
							{
								assert( 1 == 0 );
							}
							int scoreDifference = maxScore_D.first - previousScore;

							if(maxScore_D.first == currentMaximum)
							{
								if(scoreDifference != 0)
								{
									currentMaxima_coordinates.push_back(thisCoordinates);
									lastMaximumIncrease_at_diagonalI = diagonalI;
								}
							}
							else if(maxScore_D.first > currentMaximum)
							{
								currentMaximum = maxScore_D.first;
								currentMaxima_coordinates.clear();
								currentMaxima_coordinates.push_back(thisCoordinates);
								lastMaximumIncrease_at_diagonalI = diagonalI;
							}
						}
						else
						{
							// assert(1 == 0);
						}
					}
				}
			}
		}

		if(verbose)
			std::cout << "\t\tfiltering" << "\n" << std::flush;

		std::vector<std::vector<int> > m_thisDiagonal_filtered;
		if(m_thisDiagonal.size() > 0)
		{
			double max;
			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				if((i == 0) || (max < S))
				{
					max = S;
				}
			}

			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				assert(S <= max);
				if((max - S) <= threshold_for_filtering)
				{
					m_thisDiagonal_filtered.push_back(coordinates);
				}
			}

			m_thisDiagonal = m_thisDiagonal_filtered;
		}

		m2_diagonal = m1_diagonal;
		m1_diagonal = m_thisDiagonal;
	}

	std::vector<localExtension_pathDescription> forReturn;
	auto backtraceFrom = [&](int start_x, int start_y, int start_z, double StartScore) {
		int backtrace_x = start_x;
		int backtrace_y = start_y;
		int backtrace_z = start_z;
		int backtrace_matrix = 0;

		std::string reconstructed_graph;
		std::string reconstructed_sequence;
		std::vector<int> reconstructed_graph_levels;
		std::vector<std::vector<int>> edge_coordinates;
		std::vector<Edge*> utilizedEdges;

		std::vector<int> startCoordinates;
		startCoordinates.push_back(start_x);
		startCoordinates.push_back(start_y);
		startCoordinates.push_back(start_z);
		edge_coordinates.push_back(startCoordinates);

		if(verbose)
		{
			std::cout << "\tbacktraceFrom() called.\n";

		}
		while((backtrace_x != startLevel_graph) || (backtrace_y != start_sequence))
		{
			if(verbose)
			{
				std::cout << "\t\tbacktrace_x: " << backtrace_x << "\n";
				std::cout << "\t\tbacktrace_y: " << backtrace_y << "\n";
				std::cout << "\t\tbacktrace_z: " << backtrace_z << "\n";
				std::cout << "\t\tselected matrix: " << backtrace_matrix << "\n" << std::flush;
			}

			if(directionPositive)
			{
				assert((backtrace_x >= startLevel_graph) && (backtrace_x <= (levels - 1)));
				assert((backtrace_y >= start_sequence) && (backtrace_y <= (int)sequenceLength));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}
			else
			{
				assert((backtrace_x <= startLevel_graph) && (backtrace_x >= 0));
				assert((backtrace_y <= start_sequence) && (backtrace_y >= 0));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}

			double Score;
			backtraceStep_affine step;
			if(backtrace_matrix == 0)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
			}
			else if(backtrace_matrix == 1)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
			}
			else if(backtrace_matrix == 2)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
			}

			if((backtrace_x == start_x) && (backtrace_y == start_y) && (backtrace_z == start_z))
			{
				Score = StartScore;
			}

			if(verbose)
			{
				std::cout << "\t\tScore: " << Score << "\n\n" << std::flush;
			}

			std::string edgeEmission;
			if(step.usedEdge != 0)
			{
				edgeEmission = g->CODE.deCode(step.usedEdge->locus_id, step.usedEdge->emission);
				assert(edgeEmission.size() == 1);
			}
			std::string sequenceEmission;
			if((backtrace_y >= 1) && (directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y - 1, 1);
			}
			if((backtrace_y < max_seqI) && (! directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y, 1);
			}

			int next_x = step.x;
			int next_y = step.y;
			int next_z = step.z;
			int next_matrix = step.sourceMatrix;

			bool dontAddCoordinates = false; // if we jump from one matrix to the other without changing coordinates, we don't store the coordinates...
			if(directionPositive)
			{
				if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);

				}
				else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}
			else
			{
				if((next_x == (backtrace_x + 1)) && (next_y == (backtrace_y + 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y + 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);
				}
				else if((next_x == (backtrace_x + 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}

			if(verbose)
			{
				std::cout << "\t\t\t" << "Emission: " << edgeEmission << "\n";
				std::cout << "\t\t\t" << "next_x: " << next_x << "\n";
				std::cout << "\t\t\t" << "next_y: " << next_y << "\n";
				std::cout << "\t\t\t" << "next_z: " << next_z << "\n";
				std::cout << "\t\t\t" << "next_matrix: " << next_matrix << "\n\n" << std::flush;
			}

			backtrace_x = next_x;
			backtrace_y = next_y;
			backtrace_z = next_z;
			backtrace_matrix = next_matrix;

			if( ! dontAddCoordinates)
			{
				std::vector<int> nextCoordinates;
				nextCoordinates.push_back(next_x);
				nextCoordinates.push_back(next_y);
				nextCoordinates.push_back(next_z);
				edge_coordinates.push_back(nextCoordinates);
			}
		}

		if(directionPositive)
		{
			std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
			std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
			std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());
			std::reverse(utilizedEdges.begin(), utilizedEdges.end());
			std::reverse(edge_coordinates.begin(), edge_coordinates.end());
		}

		assert(reconstructed_graph.size() == reconstructed_sequence.size());
		assert(reconstructed_graph_levels.size() == reconstructed_graph.size());

		localExtension_pathDescription pathReturn;
		pathReturn.Score = StartScore;
		pathReturn.usedEdges = utilizedEdges;
		pathReturn.coordinates = edge_coordinates;
		pathReturn.alignedSequence = reconstructed_sequence;
		pathReturn.alignedGraph = reconstructed_graph;
		pathReturn.alignedGraph_levels = reconstructed_graph_levels;

		forReturn.push_back(pathReturn);
	};

	if(verbose)
	{
		for(unsigned int x = 231412; x <= 231418; x++)
		{
			for(unsigned int y = 0; y <= 5; y++)
			{
				std::cerr << "x = " << x << ", y = " << y << ": ";
				if((scores.count(x) == 0) || (scores.at(x).count(y) == 0) || (scores.at(x).at(y).size() == 0))
				{
					std::cerr << "Undefined.";
				}
				else
				{
					for(unsigned int z = 0; z < scores.at(x).at(y).size(); z++)
					{
						std::cerr << "\n\t z = " << z << ": " << scores.at(x).at(y).at(z).D;
					}
				}
				std::cerr << "\n\n" << std::flush;
			}
		}
	}

	if(! returnGlobalScore)
	{
		if(currentMaximum > 0)
		{
			if(verbose)
			{
				std::cout << "\tMaximum " << currentMaximum << ", achieved at " << currentMaxima_coordinates.size() << " positions!\n" << std::flush;
			}

			for(int maximumI = 0; maximumI < (int)currentMaxima_coordinates.size(); maximumI++)
			{
				std::vector<int> coordinates = currentMaxima_coordinates.at(maximumI);
				if(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D != minusInfinity)
				{
					if(verbose)
						std::cout << " - Start maximum backtrace from " <<coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

					backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D);
				}
			}
		}

		if(verbose)
			std::cout << "Have hit " << hit_NW_paths.size() << " NW paths, backtrace independent of achieved value!\n" << std::flush;

		for(std::map<NWPath*, std::pair<double, std::vector<int> > >::iterator hitPathsIt = hit_NW_paths.begin(); hitPathsIt != hit_NW_paths.end(); hitPathsIt++)
		{
			// NWPath* path = hitPathsIt->first;
			double score = hitPathsIt->second.first;
			std::vector<int> coordinates = hitPathsIt->second.second;
			assert(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D == minusInfinity);

			if(verbose)
				std::cout << " - Start NW path backtrace from " << coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

			backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), score);
		}
	}
	else
	{

		if(directionPositive)
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(levels-1).at(sequenceLength).begin(); altIt != scores.at(levels-1).at(sequenceLength).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(levels - 1, sequenceLength, maxFinalState.second, maxFinalState.first);
		}
		else
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(0).at(0).begin(); altIt != scores.at(0).at(0).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;

				if(verbose)
					std::cout << "Final state z-value " << altIt->first << ", D value " << scores.at(0).at(0).at(altIt->first).D << ", GraphGap: " <<  scores.at(0).at(0).at(altIt->first).GraphGap << ", SequenceGap: " <<  scores.at(0).at(0).at(altIt->first).SequenceGap << "\n" << std::flush;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(0, 0, maxFinalState.second, maxFinalState.first);
		}
	}
	return forReturn;
}


int GraphAlignerUnique::countMatchesInSequence(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence)
{
	assert(reconstructed_graph.length() == reconstructed_sequence.length());
	assert(reconstructed_graph_levels.size() == reconstructed_graph.length());

	int matches = 0;

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::string G = reconstructed_graph.substr(i, 1);

		if(S != "_")
		{
			if(S == G)
			{
				matches++;
			}
		}
	}

	return matches;
}




};

