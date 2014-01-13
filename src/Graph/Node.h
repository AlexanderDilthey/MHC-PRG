/*
 * Node.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef NODE_H_
#define NODE_H_

struct nodeHoldsHaplo {
	int haploID;
	double P;
	int indexInHaplotypeVector;
};

#include <vector>
#include <string>
#include <set>

#include "Graph.h"
#include "LargeGraph.h"
#include "MultiGraph.h"

#include "Edge.h"

struct HaploAttachedToNode {
	Node* n;
	int indexInNodeHaploInfo;
	double P;
};

using namespace std;



class Node {
public:
	Node();

	Graph* g;
	LargeGraph* lG;
	MultiGraph* mG;

	set<Edge*> Incoming_Edges;
	set<Edge*> Outgoing_Edges;

	vector<nodeHoldsHaplo> haplotypes;

	double Sum_Incoming();
	double Sum_Outgoing();

	unsigned int level;
	bool terminal;
};

#endif /* NODE_H_ */

