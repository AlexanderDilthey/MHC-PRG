/*
 * GraphAlignernonAffine.h
 *
 *  Created on: 27.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHALIGNERNONAFFINE_H_
#define GRAPHALIGNERNONAFFINE_H_

#include "GraphAligner.h"

#include <vector>
#include <utility>
#include <map>
#include <list>
#include <algorithm>

class GraphAligner_nonAffine: public GraphAligner {
public:
	GraphAligner_nonAffine(Graph* graph, int k);

	int score_fullNeedleman(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);

	void fullNeedleman(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned);
	std::map<NWEdge*, std::map<NWEdge*, int> > naivelyCalculateGraphDistancesForVirtualNW_nonAffine(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, int> >& lastPositionDistances_perZ);
	void fullNeedleman_diagonal(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned);

	void seedAndExtend(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned, bool testing = false);

	void findShortGappedGraphConnection_nonAffine(int start_x, int start_z, int stop_x, int stop_z, std::vector<Edge*>& usedEdges, std::string& graphSequence, std::vector<int>& graphSequence_levels);

	std::map<NWEdge*, std::map<NWEdge*, int> > fastCalculateGraphDistancesForVirtualNW_nonAffine(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, int> >& lastPositionDistances_perZ);

};

extern int threads;

#endif /* GRAPHALIGNERNONAFFINE_H_ */
