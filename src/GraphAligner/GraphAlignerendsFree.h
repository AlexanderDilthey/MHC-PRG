/*
 * GraphAlignerendsFree.h
 *
 *  Created on: 27.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHALIGNERENDSFREE_H_
#define GRAPHALIGNERENDSFREE_H_

#include "GraphAligner.h"

#include <vector>
#include <utility>
#include <map>
#include <list>
#include <algorithm>

class GraphAligner_endsFree: public GraphAligner {
public:
	GraphAligner_endsFree(Graph* graph, int k);

	void fullNeedleman_endsFree(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned);
	double score_endsFree(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);
	std::vector<localExtension_pathDescription> fullNeedleman_endsFree_diagonal_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int diagonal_stop_threshold, VirtualNWTable* blockedPathsTable, bool directionPositive, bool returnGlobalScore = false);
	seedAndExtend_return seedAndExtend4(std::string sequence, bool noExtension = false);
	void seedAndExtend4_backtrack(VirtualNWTable& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels);
	void findShortGappedGraphConnection_endsFree(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally);
};

extern int threads;

#endif /* GRAPHALIGNERENDSFREE_H_ */
