/*
 * GraphAlignerAffine.h
 *
 *  Created on: 27.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHALIGNERAFFINE_H_
#define GRAPHALIGNERAFFINE_H_

#include "GraphAligner.h"

#include <vector>
#include <utility>
#include <map>
#include <list>
#include <algorithm>

class GraphAligner_affine: public GraphAligner {
public:
	GraphAligner_affine(Graph* graph, int k);

	double score_fullNeedleman_affine(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);

	void fullNeedleman_affine(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned);
	void fullNeedleman_affine_diagonal(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned);
	std::vector<localExtension_pathDescription> fullNeedleman_affine_diagonal_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int diagonal_stop_threshold, VirtualNWTable* blockedPathsTable, bool directionPositive, bool returnGlobalScore = false);

	seedAndExtend_return seedAndExtend2(std::string sequence, bool noExtension = false, int distanceCalculationMode = 1);
	seedAndExtend_return seedAndExtend3(std::string sequence, bool noExtension = false);

	void seedAndExtend2_backtrack(VirtualNWTable& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels, 	std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal);
	void seedAndExtend3_backtrack(VirtualNWTable& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels);

	void findShortGappedGraphConnection_affine(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally);

	void affinelyCalculateGraphDistancesForVirtualNW(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal);
	void affinelyCalculateGraphDistancesForVirtualNW_faster(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startInAffineGap,  std::vector<std::map<NWEdge*, graphPointDistance> >& lastPositionDistances_perZ_startNormal, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startAffineGap, std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >& NWedges_graphDist_startNormal);

};

extern int threads;

#endif /* GRAPHALIGNERAFFINE_H_ */
