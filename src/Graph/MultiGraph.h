/*
 * Graph.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef MULTIGRAPH_H_
#define  MULTIGRAPH_H_

class MultiGraph;
class Node;
class Edge;

#include <vector>
#include <string>
#include <set>
#include "Node.h"
#include "Edge.h"
#include "LargeGraph.h"
#include "../MHC-PRG.h"
#include "../LargeLocusCodeAllocation.h"


using namespace std;

struct multiHaploLabelPair
{
	vector<string> h1;
	vector<string> h2;
};

struct avg_struct
{
	int count;
	int sum;
};

struct kMerPositionInfo
{
	set<string> kMers;
	map<string, set<int> > kMers2Level;
	map<int, set<string> > Level2Kmers;
	map<string, int > kMerMultiplicity;
	map<string, int> kMerMaxPerLevel;
	map<string, int> kMerEdgeNum;

};


struct kMerUniquenessInfo
{
	map<string, set<int> > kMers2Level;
	set<string> levelUniqueKMers;
};

struct kMerNonUniquenessInfo
{
	map<string, int > kMerMultiplicity;
	map<string, set<int> > kMers2Level;
	set<string> levelNonUniqueKMers;
	map<int, set<string> > level2kMers;
};

class MultiGraph {
public:
	MultiGraph();
	~MultiGraph();

	set<Node*> Nodes;
	set<Edge*> Edges;
	vector< set<Node*> > NodesPerLevel;

	LargeLocusCodeAllocation CODE;

	void registerNode(Node* n, unsigned int level);
	void registerEdge(Edge* e);
	void unRegisterNode(Node* n);
	void unRegisterEdge(Edge* e);

	void checkConsistency(bool terminalCheck);
	vector<string> getAssignedLoci();
	
	void freeMemory();

	void stats();

	int kMerSize;
	LargeGraph* underlyingLargeGraph;

	diploidEdgePointerPath LargeEdgePointerPathToMultiEdgePointerPath(diploidEdgePointerPath largeEdgePointerPath, bool sameLength);
	diploidEdgePointerPath MultiGraphEdgesToLargeGraphEdges(diploidEdgePointerPath& multiEdgesPath);

	multiHaploLabelPair edgePointerPathToLabels(diploidEdgePointerPath& multiEdgesPath);

	void kMerDiagnostics();
	kMerUniquenessInfo kMerUniqueness();
	kMerPositionInfo getkMerPositions();

	kMerNonUniquenessInfo kMerNonUniqueness();

	int maxEdgeNumber();
	int minEdgeNumber();


};

#endif /*  MULTIGRAPH_H_ */
