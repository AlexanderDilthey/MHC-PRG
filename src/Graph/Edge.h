/*
 * Edge.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <vector>
#include <string>
#include <map>

#include "Graph.h"
#include "Node.h"

using namespace std;

class Edge;

class Edge {
public:
	Edge();

	Node* From;
	Node* To;

	double count;
	unsigned char emission;
	int largeEmission;
	int largeEmission2;
	map<int, int> multiEmission;
	vector<int> multiEmission_full;

	Edge* _largeGraphEdge;
	vector<int> levelsNucleotideGraph;

	string locus_id;
	string label;

	bool pgf_protect;
};

#endif /* EDGE_H_ */
