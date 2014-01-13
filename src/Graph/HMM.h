#ifndef HMM_H_
#define HMM_H_

#include "Graph.h"
#include "Edge.h"
#include <vector>
#include <map>
#include "../LocusCodeAllocation.h"
#include "../MHC-PRG.h"

struct stateAlternative
{
	int id;
	double p;
};

struct haplotypePair
{
	basic_string<unsigned char> h1;
	basic_string<unsigned char> h2;
};

struct haploidState
{
	int level;
	int id;
	Edge* e;
	vector<stateAlternative> next;
	vector<stateAlternative> previous;
	vector<MutateAwayTo> mutate;
};

#endif /*HMM_H_*/
