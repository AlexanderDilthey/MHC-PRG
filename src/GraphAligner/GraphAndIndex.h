/*
 * GraphAndIndex.h
 *
 *  Created on: 18.06.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHANDINDEX_H_
#define GRAPHANDINDEX_H_

#include <map>
#include <vector>

#include "../Graph/Graph.h"

class kMerChain {
public:
	int sequence_begin;
	int sequence_end;
	bool sequence_reverse;

	int graph_firstLevel;
	int graph_lastLevel;

	int lastAttachedkMer_level;
	int lastAttachedkMer_graph_level;

	int matchedkMers;
	bool extended;

	kMerChain() : sequence_begin(-1), sequence_end(-1), sequence_reverse(false), graph_firstLevel(-1), graph_lastLevel(-1), lastAttachedkMer_level(-1), lastAttachedkMer_graph_level(-1), matchedkMers(-1), extended(false)
	{

	}
};

class GraphAndIndex {
	int kMerSize;

	std::map<std::string, std::vector< std::pair<int, int> > > kMerPositions;
	Graph* g;

public:
	GraphAndIndex(Graph* graph, int k);
	void Index();
	void printIndex();
	std::vector< std::pair<int, int> > queryIndex(std::string kMer);

	std::vector<kMerChain> findChainsPlusStrand(std::string sequence);
	std::vector<kMerChain> findChains(std::string sequence);

};

void _printChains(std::vector<kMerChain>& chains);

#endif /* GRAPHANDINDEX_H_ */
