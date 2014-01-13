/*
 * MHC-PRG.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef MHCPRG_H_
#define MHCPRG_H_

#include <string>
#include <vector>
#include <map>
#include "Data/GenotypePanel.h"
#include "Data/HaplotypePanel.h"

#include "Graph/Graph.h"
#include "Graph/Node.h"

struct Config {
	int threads;
	bool quiet;
};

void errEx(std::string message);

extern Config CONFIG;
extern double epsilon;

class positionsSorter {
public:
	map<string, int>* p;
	bool operator() (string i, string j) {
		if(p->count(i) == 0)
			errEx("No position information for locus: "+i);
		if(p->count(j) == 0)
			errEx("No position information for locus: "+j);
		return ((*p)[i] < (*p)[j]);
	}
};

#endif /* MHCPRG_H_ */
