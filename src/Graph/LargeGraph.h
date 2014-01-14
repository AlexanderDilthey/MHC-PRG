/*
 * Graph.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef LARGEGRAPH_H_
#define LARGEGRAPH_H_

class LargeGraph;
class Node;
class Edge;

struct diploidEdgePath;
struct diploidEdgePointerPath;
struct diploidNucleotidePath;
struct levelInfo;

#include "../NextGen/Validation.h"


#include <vector>
#include <string>
#include <set>
#include "Node.h"
#include "Edge.h"
#include "Graph.h"
#include "../MHC-PRG.h"
#include "../LargeLocusCodeAllocation.h"

using namespace std;

typedef std::vector< std::vector<std::string> > diploidGenomeString;


struct diploidEdgePath
{
	vector<int> h1;
	vector<int> h2;
};

struct diploidEdgePointerPath
{
	vector<Edge*> h1;
	vector<Edge*> h2;
};

struct diploidNucleotidePath
{
	vector<string> h1;
	vector<string> h2;
};

struct readSimulationResults
{
	vector<string> haplotypes;
	vector<string> haplotypes_untrimmed;
	vector< vector<int> > haplotypeStartPositions;
	map<string, long long> simulatedReads;
	map<string, long long> kMerOccurence;
	map<string, set<int> > kMerOccurenceLocalization;
};

class LargeGraph {
public:
	LargeGraph();
	~LargeGraph();
	LargeGraph(LargeGraph& otherGraph);

	// if add new data members, remember copy constructor!

	set<Node*> Nodes;
	set<Edge*> Edges;
	vector< set<Node*> > NodesPerLevel;
	map<Edge*, int> EdgeCounter;
	map<Node*, int> NodeCounter;
	map<int, Node*> idx2Node;
	map<int, Edge*> idx2Edge;
	LargeLocusCodeAllocation CODE;
	int kMerSize;

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

	int trimGraph();

	void simulateHaplotypes(int number);
	void stats();


	vector<string> requiredKMers();

	string newick();

	void assignEdgeNodeIndices();
	diploidEdgePath simulateRandomPath();
	readSimulationResults simulateReadsForPath(diploidEdgePath& path);

	vector<string> haploidPathToNucleotides(vector<int> path);
	diploidNucleotidePath diploidPathToNucleotides(diploidEdgePath& path);

	diploidNucleotidePath diploidPathToAlignedNucleotides(diploidEdgePath& path);
	diploidNucleotidePath diploidPathToAlignedNucleotides(diploidEdgePointerPath& path);

	diploidGenomeString diploidPathToGenomeString(diploidEdgePath& path);
	diploidGenomeString diploidPathToGenomeString(diploidEdgePointerPath& path);

	vector< vector<string> > categorizeEdgeLevels();

	vector<int> edgePointerPathToEdgeIndexPath(vector<Edge*>& edgePointerPath);
	diploidEdgePath edgePointerPathToEdgeIndexPath(diploidEdgePointerPath& edgePointerPath);

	vector<Edge*> edgeIndexPathToEdgePointerPath(vector<int>& edgeIndexPath);
	diploidEdgePointerPath edgeIndexPathToEdgePointerPath(diploidEdgePath& edgeIndexPath);



	vector<levelInfo> getLevelInfo();

	void removeStarPaths();

	void reverseGraph();

};

#endif /* LARGEGRAPH_H_ */
