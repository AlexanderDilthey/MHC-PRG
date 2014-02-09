/*
 * Graph.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef GRAPH_H_
#define GRAPH_H_

class Graph;
class Node;
class Edge;

#include <vector>
#include <string>
#include <set>
#include "Node.h"
#include "Edge.h"
#include "../MHC-PRG.h"
#include "../LocusCodeAllocation.h"
#include "LargeGraph.h"

using namespace std;

struct levelInfo
{
	int nodes;
	int edges;
	int symbols;
	int symbols_CODE;
};


class Graph {
public:
	Graph();

	set<Node*> Nodes;
	set<Edge*> Edges;
	vector< set<Node*> > NodesPerLevel;

	LocusCodeAllocation CODE;
	void registerNode(Node* n, unsigned int level);
	void registerEdge(Edge* e);
	void unRegisterNode(Node* n);
	void unRegisterEdge(Edge* e);

	void checkConsistency(bool terminalCheck);
	void checkLocusOrderConsistency(vector<string> loci);
	vector<string> getAssignedLoci();
	
	void freeMemory();

	void writeToFile(string filename);
	void readFromFile(string filename);
	void printComplexity (string filename);


	int trimGraph(bool remove2DHLA = false);

	void simulateHaplotypes(int number);
	void makeEdgesGaps(double proportion);

	void removeStarPaths();


	void graphViz(int level_start, int level_stop, std::string output_filename);

	vector<levelInfo> getLevelInfo();
	
	std::string getOneLocusIDforLevel(unsigned int level);

	diploidEdgePointerPath simulateRandomDiploidPath();


};

#endif /* GRAPH_H_ */
