/*
 * GraphAndEdgeIndex.h
 *
 *  Created on: 29.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHANDEDGEINDEX_H_
#define GRAPHANDEDGEINDEX_H_

#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "../Graph/Graph.h"

namespace GraphAlignerUnique {

class kMerInGraphSpec {
public:
	std::vector<Edge*> traversedEdges;
	std::string getID()
	{
		std::string forReturn;
		std::ostringstream idString;
		idString << setw(15) << (void const *) traversedEdges.front();
		idString << " -- ";
		idString << setw(15) << (void const *) traversedEdges.back();
		return idString.str();
	}
};

class kMerEdgeChain {
public:
	int sequence_begin;
	int sequence_end;

	std::vector<Edge*> traversedEdges;
	kMerEdgeChain() : sequence_begin(-1), sequence_end(-1)
	{

	}
};

class GraphAndEdgeIndex {
	int kMerSize;
	Graph* g;

	std::map<std::string, std::vector<kMerInGraphSpec> > kMers;
	std::map<Node*, std::vector<Edge*>> nodes_jumpOverGaps;

	// largely useless!
	LargeGraph* kMerGraph;

	std::set<Node*> generated_nodes;
	std::set<Edge*> generated_edges;


public:
	GraphAndEdgeIndex(Graph* graph, int k);
	~GraphAndEdgeIndex();
	void Index();
	void fillEdgeJumper();
	void printIndex();
	std::vector<kMerEdgeChain*> findChains(std::string sequence);
	std::vector<kMerInGraphSpec> queryIndex(std::string kMer);
	std::vector<std::string> getIndexedkMers();
};
}
#endif /* GRAPHANDEDGEINDEX_H_ */
