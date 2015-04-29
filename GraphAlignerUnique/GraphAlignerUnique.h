/*
 * GraphAlignerUnique.h
 *
 *  Created on: 30.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHALIGNERUNIQUE_H_
#define GRAPHALIGNERUNIQUE_H_

#include "GraphAndEdgeIndex.h"
#include "../GraphAligner/GraphAligner.h"
#include "VirtualNWUnique.h"
#include "../NextGen/readSimulator.h"
#include "coveredIntervals.h"


#include <map>
#include <vector>
#include <string>
#include <functional>
#include <utility>

namespace GraphAlignerUnique {


class backtrackBookkeeper_UGA {
public:
	double S;
	Edge* usedEdges;
	char graphSequence;
	int graphSequence_levels;
	backtrackBookkeeper_UGA* backtrack;

	int takenNodesAndEdges;

	void takeNodeAndEdge(Node* nodeFrom, Edge* usedEdge)
	{
		assert(takenNodesAndEdges == 0);
		std::string edgeEmission = usedEdge->From->g->CODE.deCode(usedEdge->locus_id, usedEdge->emission);
		assert(edgeEmission.length() == 1);
		usedEdges = usedEdge;
		graphSequence = edgeEmission.at(0);
		graphSequence_levels = nodeFrom->level;
		assert(usedEdge->From == nodeFrom);
		takenNodesAndEdges++;
	}

	backtrackBookkeeper_UGA()
	{
		takenNodesAndEdges = 0;

		backtrack = 0;
		usedEdges = 0;
		graphSequence = 0;
		graphSequence_levels = 0;
		S = 0;
	}
};


class VirtualNWTable_Unique;

class GraphAlignerUnique {

	Graph* g;
	int kMerSize;
	GraphAndEdgeIndex gI;
	double S_match;
	double S_mismatch;
	double S_gap;
	double S_openGap;
	double S_extendGap;
	double S_graphGap;
	int randomizationParameter;
	int iterationsMainRandomizationLoop;
	int minimumChainUniqueness;
	bool verbose;

	int threads;

	coveredIntervals myGraph_coveredIntervals;

	std::vector<unsigned int> rng_seeds;

	std::vector<std::vector<Node*> > nodesPerLevel_ordered;
	std::vector<std::map<Node*, unsigned int> > nodesPerLevel_ordered_rev;

	void seedAndExtend_init_occurrence_strand_etc(std::string& sequence_nonReverse, bool& useReverse, std::string& sequence, std::vector<std::string>& kMers_sequence, std::map<std::string, int>& kMer_sequence_occurrences);
	bool iskMerDoubleUnique(std::string& kMer, std::map<std::string, int>& occurrencesInSequence);
	std::vector<kMerEdgeChain*> trimChainsForUniqueness(std::vector<kMerEdgeChain*>& inputChains, std::string& sequence, std::map<std::string, int>& occurrencesInSequence);
	void analyzeChainUniqueness(std::string& sequence, std::vector<kMerEdgeChain*>& uniquelyTrimmedChains, std::map<std::string, int>& kMer_sequence_occurrences, std::map<kMerEdgeChain*, int>& uniquelyTrimmedChains_doubleUniquekMers, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered);
	void cleanChainsAccordingToSelectedChain(kMerEdgeChain* selectedChain, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& availableChains);
	void selectedChainsConsistencyCheck(std::set<kMerEdgeChain*> selectedChains);
	void fixUniqueChains(std::string& sequence, bool thisIterationRandomization, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered, std::set<kMerEdgeChain*>& selectedChains, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, int>& uniquelyTrimmedChains_doubleUniquekMers, bool rescueNonUnique);
	void trimChain(kMerEdgeChain* inputChain, std::string& sequence, int removeLeft, int removeRight);

	void fixNonUniqueChains(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered, bool thisIterationRandomization, std::vector<kMerEdgeChain*>& allChains, std::vector<kMerEdgeChain*>& newChains);

	std::vector<std::pair<int, Edge*> > _graph_get_previous_z_values_and_edges(int x, int z);
	std::vector<std::pair<int, Edge*> > _graph_get_next_z_values_and_edges(int x, int z);
	std::vector<int> _graph_get_previous_z_values(int x, int z);

	std::vector<std::pair<int, int> > findGaps_chainCoverage(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered_input);
	void closeOneGap_withNonUniqueChains(std::string& sequence, std::pair<int, int> gapCoordinates, std::vector<kMerEdgeChain*>& sequencePositions_covered, bool thisIterationRandomization, std::vector<kMerEdgeChain*>& chainsToConsider, std::vector<kMerEdgeChain*>& newTrimmedChains, int& assignedChains);

	void cleanChainsAccordingToBoundaries(Node* leftBoundaryNode, Node* rightBoundaryNode, std::vector<kMerEdgeChain*>& chains);
	void analyzeLocalChainUniqueness(std::string& sequence, std::pair<int, int>& spannedRegionCoordinates, std::pair<int, int>& spannedRegionGraphBoundaries, std::vector<kMerEdgeChain*>& availableChains, 	std::map<kMerEdgeChain*, int>& chains_localUniqueness, std::set<kMerEdgeChain*, std::function<bool(kMerEdgeChain*,kMerEdgeChain*)>>& uniquelyTrimmedChains_ordered_output);
	bool iskMerLocallyDoubleUnique(std::string& kMer, std::pair<int, int>& spannedRegionCoordinates, std::pair<int, int>& spannedRegionGraphBoundaries, std::map<std::string, int>& occurrencesInSequence);

	void buildkMerOccurrenceMap(std::string& S, std::map<std::string, int>& occurrencesInSequence);

	void printSequenceChainCoverageStats(std::string& sequence, std::vector<kMerEdgeChain*>& sequencePositions_covered);

	void kMerEdgeChains2vNW(VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, int>& currentChains_start, std::map<kMerEdgeChain*, NWPath*>& chains2Paths);
	void vNW_completeRemainingGaps_and_score(std::string& sequence, VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, NWPath*>& chains2Paths, std::map<kMerEdgeChain*, int>& currentChains_start, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, double& finalScore, int& finalScore_z, NWEdge*& finalScore_backtrack);
	void vNW_completeRemainingGaps_and_score_local(std::string& sequence, VirtualNWTable_Unique& vNW, std::vector<kMerEdgeChain*>& sequencePositions_covered, std::map<kMerEdgeChain*, NWPath*>& chains2Paths, std::map<kMerEdgeChain*, int>& currentChains_start, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, double& finalScore, int& finalScore_z, NWEdge*& finalScore_backtrack, bool greedyLocalExtension = false);

	void findShortGappedGraphConnection_affine(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally);
	void findShortGappedGraphConnection_affine_MTM(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally);

	void seedAndExtend_backtrack(VirtualNWTable_Unique& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap);
	void seedAndExtend_backtrack_local(VirtualNWTable_Unique& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap);
	void affinelyCalculateGraphDistancesForVirtualNW(int fromLevel, int toLevel, std::vector<NWEdge*>& entryEdges_vector, std::vector<NWEdge*>& exitEdges_vector, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal);
	void affinelyCalculateGraphDistancesForVirtualNW_local(int fromLevel, int toLevel, std::vector<NWEdge*>& entryEdges_vector, std::vector<NWEdge*>& exitEdges_vector, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal);

	std::vector<localExtension_pathDescription> fullNeedleman_diagonal_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive,  bool returnGlobalScore);
	std::vector<localExtension_pathDescription> fullNeedleman_diagonal_completeGreedy_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, bool directionPositive);
	
	bool alignmentContainedWithinAreaCoveredByMyGraph(std::string regionID, int start, int stop);

public:
	GraphAlignerUnique(Graph* graph, int k);
	seedAndExtend_return seedAndExtend(std::string sequence);
	seedAndExtend_return_local seedAndExtend_local(std::string sequence, std::vector<seedAndExtend_return_local>& allBacktraces);
	seedAndExtend_return_local seedAndExtend_short(std::string sequence, std::vector<seedAndExtend_return_local>& allBacktraces);

	std::pair<seedAndExtend_return_local, seedAndExtend_return_local> seedAndExtend_local_paired_or_short(oneReadPair readPair, bool usePairing, bool use_short, double insertSize_mean, double insertSize_sd, bool estimateInsertSize, std::map<int, double>& insertSize_posterior_ret);
	std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > seedAndExtend_short_allAlignments(oneReadPair readPair, double insertSize_mean, double insertSize_sd);

	std::vector<seedAndExtend_return_local> seedAndExtend_longlocal_allAlignments(oneRead R);

	double scoreOneAlignment(oneRead& underlyingRead, seedAndExtend_return_local& alignment, int& totalMismatches);


	GraphAndEdgeIndex& getGI()
	{
		return gI;
	}

	double get_S_match()
	{
		return S_match;
	}
	double get_S_mismatch()
	{
		return S_mismatch;
	}
	double get_S_gapOpen()
	{
		return S_openGap;
	}
	double get_S_gapExtend()
	{
		return S_extendGap;
	}

	double score(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);
	int countMatchesInSequence(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);

	int getIterationsMainRandomizationLoop()
	{
		return iterationsMainRandomizationLoop;
	}

	void setIterationsMainRandomizationLoop(int i)
	{
		iterationsMainRandomizationLoop = i;
	}

	void setThreads(int t)
	{
		threads = t;
	}

	int getThreads()
	{
		return threads;
	}

	void printkMerProfile(std::string file);

	friend class VirtualNWTable_Unique;
	friend class NWEdge;
	friend class NWPath;

};

};

#endif /* GRAPHALIGNERUNIQUE_H_ */
