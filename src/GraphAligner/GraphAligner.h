/*
 * GraphAligner.h
 *
 *  Created on: 18.06.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHALIGNER_H_
#define GRAPHALIGNER_H_

#include "GraphAndIndex.h"
#include "VirtualNW.h"

#include <vector>
#include <utility>
#include <map>
#include <list>
#include <algorithm>
#include <iostream>

class localExtension_pathDescription {
public:
	std::vector<std::vector<int> > coordinates;
	std::vector<Edge*> usedEdges;
	double Score;
	std::string alignedSequence;
	std::string alignedGraph;
	std::vector<int> alignedGraph_levels;

	void _printExtension()
	{
		std::cout << "localExtension_pathDescription object:\n" << std::flush;
		for(unsigned int cI = 0; cI < coordinates.size(); cI++)
		{
			for(unsigned int cII = 0; cII < coordinates.at(cI).size(); cII++)
			{
				std::cout << coordinates.at(cI).at(cII) << " " << std::flush;
			}
			std::cout << "\n" << std::flush;
			if(cI < usedEdges.size())
			{
				std::cout << "\t " << usedEdges.at(cI) << "\t" << std::flush;
				if(usedEdges.at(cI) != 0)
				{
					Edge* e = usedEdges.at(cI);
					Graph* g = usedEdges.at(cI)->From->g;
					assert(g != 0);
					std:cout << g->CODE.deCode(e->locus_id, e->emission) << " " << std::flush;
				}
				std::cout << "\n" << std::flush;
			}
		}
	}
};

class seedAndExtend_return {
public:
	double Score;
	std::string graph_aligned;
	std::vector<int> graph_aligned_levels;
	std::string sequence_aligned;
};

class seedAndExtend_return_local {
public:
	double Score;
	std::string graph_aligned;
	std::vector<int> graph_aligned_levels;
	std::vector<double> certainty_sequence2Graph;
	std::string sequence_aligned;
	size_t kMers_total;
	size_t kMers_unique_total;
	size_t kMers_unique_utilized;
};


class graphPointDistance {
public:
	double Score_endAffinely;
	double Score_endInAnything;
	bool start_in_affine_sequenceGap;
};

class backtrackBookkeeper {
public:
	double S;
	std::vector<Edge*> usedEdges;
	std::string graphSequence;
	std::vector<int> graphSequence_levels;
	backtrackBookkeeper* backtrack;

	void takeNodeAndEdge(Node* nodeFrom, Edge* usedEdge)
	{
		std::string edgeEmission = usedEdge->From->g->CODE.deCode(usedEdge->locus_id, usedEdge->emission);
		assert(edgeEmission.length() == 1);
		usedEdges.push_back(usedEdge);
		graphSequence.append(edgeEmission);
		graphSequence_levels.push_back(nodeFrom->level);
		assert(usedEdge->From == nodeFrom);
	}

	backtrackBookkeeper()
	{
		backtrack = 0;
	}
};



class graphPointDistance_withBacktrack {
public:
	backtrackBookkeeper endAffinely;
	backtrackBookkeeper endInAnything;
	bool start_in_affine_sequenceGap;

	graphPointDistance getGraphPointDistance()
	{
		graphPointDistance forReturn;
		forReturn.Score_endAffinely = endAffinely.S;
		forReturn.Score_endInAnything = endInAnything.S;
		forReturn.start_in_affine_sequenceGap = start_in_affine_sequenceGap;
		return forReturn;
	}
};

class sAE3_nodePositionStorage {
protected:
	bool compare_byY (NWEdge* first, NWEdge* second)
	{
	  int first_y = (first != 0) ? first->to_y : 0;
	  int second_y = (second != 0) ? second->to_y : 0;
	  if(first_y == second_y)
	  {
		  return (first < second);
	  }
	  else
	  {
		  return (first_y < second_y);
	  }
	}

	void bubbleSort_scoreList()
	{
		std::list<NWEdge*>::iterator begin = edges_sortedByScore.begin();
		std::list<NWEdge*>::iterator end = edges_sortedByScore.end();

		if (begin == end)
			return;

		--end;

		bool swapped;
		do
		{
			swapped = false;
			std::list<NWEdge*>::iterator current = begin;
			while (current != end)
			{
				// bubble up
				std::list<NWEdge*>::iterator next = current;
				++next;

				if (edgeScores.at(*current) > edgeScores.at(*next))
				{
					std::iter_swap(current, next);
					swapped = true;
				}
				++current;
			}

			// last element is already sorted now
			--end;
		} while (swapped && begin != end);
	}

	bool paranoid;

public:
	sAE3_nodePositionStorage() : paranoid(true) {}

	std::map<NWEdge*, double> edgeScores;
	std::list<NWEdge*> edges_sortedByScore;
	std::set<NWEdge*> edges_sortedByY;

	void checkConsistency()
	{
		for(std::map<NWEdge*, double>::iterator eScoreIt = edgeScores.begin(); eScoreIt != edgeScores.end(); eScoreIt++)
		{
			NWEdge* e = eScoreIt->first;
			assert(edges_sortedByY.count(e) > 0);
			assert(std::find(edges_sortedByScore.begin(), edges_sortedByScore.end(), e) != edges_sortedByScore.end());
		}

		assert(edgeScores.size() == edges_sortedByScore.size());
		assert(edgeScores.size() == edges_sortedByY.size());

		if(paranoid)
		{
			std::list<NWEdge*> list_2 = edges_sortedByScore;
			list_2.sort([&](NWEdge* first, NWEdge* second){
				return (edgeScores.at(first) < edgeScores.at(second));
			});
			bubbleSort_scoreList();
			assert(list_2.size() == edges_sortedByScore.size());

			std::list<NWEdge*>::iterator iter_edges_sortedByScore = edges_sortedByScore.begin();
			std::list<NWEdge*>::iterator iter_list2 = list_2.begin();
			while(iter_edges_sortedByScore != edges_sortedByScore.end())
			{
				assert(iter_list2 != edges_sortedByScore.end());

				NWEdge* e1 = *iter_edges_sortedByScore;
				NWEdge* e2 = *iter_list2;
				assert(edgeScores.at(e1) == edgeScores.at(e2));

				double score_e1 = edgeScores.at(e1);
				iter_edges_sortedByScore++;
				iter_list2++;

				if(iter_edges_sortedByScore != edges_sortedByScore.end())
				{
					double score_next = edgeScores.at(*iter_edges_sortedByScore);
					assert(score_e1 <= score_next);
				}
			}
		}
	}

	void addNewEdge(NWEdge* e, double initialScore)
	{
		edgeScores[e] = initialScore;

		std::list<NWEdge*>::iterator iter_edges_sortedByScore = edges_sortedByScore.begin();
		while((iter_edges_sortedByScore != edges_sortedByScore.end()) && (edgeScores.at(*iter_edges_sortedByScore) < initialScore))
		{
			iter_edges_sortedByScore++;
		}
		edges_sortedByScore.insert(iter_edges_sortedByScore, e);
		edges_sortedByY.insert(e);
	}

	void updateEdgeValue(NWEdge* e, double newScore)
	{
		edgeScores.at(e) = newScore;
	}

	double getEdgeValue(NWEdge* e)
	{
		return edgeScores.at(e);
	}

	void initFromExistingSources(std::vector<sAE3_nodePositionStorage*> positionStorages)
	{
		for(unsigned int posI = 0; posI < positionStorages.size(); posI++)
		{
			std::list<NWEdge*> edges_to_insert;
			std::map<NWEdge*, double>& posI_edgeScores = positionStorages.at(posI)->edgeScores;
			for(std::map<NWEdge*, double>::iterator eScoreIt = posI_edgeScores.begin(); eScoreIt != posI_edgeScores.end(); eScoreIt++)
			{
				NWEdge* e = eScoreIt->first;
				if(edgeScores.count(e) == 0)
				{
					edgeScores[e] = eScoreIt->second;
					edges_sortedByY.insert(e);
					edges_to_insert.push_back(e);
				}
			}

			edges_sortedByScore.merge(edges_to_insert, [&](NWEdge* first, NWEdge* second){
				return (edgeScores.at(first) < edgeScores.at(second));
			});
		}
	}

	std::vector<NWEdge*> getEdges()
	{
		std::vector<NWEdge*> forReturn;
		forReturn.reserve(edgeScores.size());
		for(std::map<NWEdge*, double>::iterator eIt = edgeScores.begin(); eIt != edgeScores.end(); eIt++)
		{
			forReturn.push_back(eIt->first);
		}
		return forReturn;
	}
};


class backtraceStep {
public:
	int x;
	int y;
	int z;
	Edge* usedEdge;

	backtraceStep() : x(-1), y(-1), z(-1), usedEdge(0)
	{

	}
};

class backtraceStep_affine {
public:
	int x;
	int y;
	int z;
	int sourceMatrix;
	Edge* usedEdge;
	bool _border_lastStep_affine;

	backtraceStep_affine() : x(-1), y(-1), z(-1), sourceMatrix(-1), usedEdge(0), _border_lastStep_affine(false)
	{

	}
};


class GraphAligner {
protected:

	Graph* g;
	int kMerSize;
	GraphAndIndex gI;

	std::vector<std::vector<Node*> > nodesPerLevel_ordered;
	std::vector<std::map<Node*, unsigned int> > nodesPerLevel_ordered_rev;

	double S_match;
	double S_mismatch;
	double S_gap;
	double S_graphGap;

	double S_openGap;
	double S_extendGap;

	double S_longGap_open;
	double S_longGap_perGap;

	std::vector<std::pair<int, Edge*> > _graph_get_previous_z_values_and_edges(int x, int z);
	std::vector<std::pair<int, Edge*> > _graph_get_next_z_values_and_edges(int x, int z);

	std::vector<int> _graph_get_previous_z_values(int x, int z);


public:
	GraphAligner(Graph* graph, int k);

	std::vector<NWPath*> findVirtualNWExactSeeds(std::string& sequence, int sequenceStart, int sequenceStop, int graph_firstLevel, int graph_lastLevel, int minimumInexactMatchLength = -1);
	std::vector<std::pair<std::vector<Edge*>, std::string>> findEdgePathForSequence(std::string sequence, int graph_firstLevel, int graph_lastLevel);

	int countMatchesInSequence(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence);


	friend class VirtualNWTable;
	friend class NWEdge;
	friend class NWPath;

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

	GraphAndIndex* getGI()
	{
		return &gI;
	}
};

extern int threads;

#endif /* GRAPHALIGNER_H_ */
