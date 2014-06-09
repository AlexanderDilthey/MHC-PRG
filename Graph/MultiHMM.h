#ifndef MultiHMM_H_
#define MultiHMM_H_

#include "LargeGraph.h"
#include "MultiGraph.h"
#include "Edge.h"
#include "HMM.h"
#include <vector>
#include <map>
#include "../LocusCodeAllocation.h"
#include "../MHC-PRG.h"

struct largeHaplotypePair
{
	basic_string<int> h1;
	basic_string<int> h2;
};


class MultiHMM
{
protected:
	map<Edge*, int> Edge2Int;
	map<Edge*, int> Edge2Level;
	
	vector< map<int, Edge*> > Int2Edge;
	vector< vector<haploidState> > haploidStatesByLevel;
	
	vector<double> haploid_initialProbabilities;
	vector<double> diploid_initialProbabilities;

	vector<int> diploidStatesByLevel;
	vector< vector< vector<stateAlternative> > > diploidStateTransitions;
	vector< vector< vector<stateAlternative> > > diploidStateTransitions_Reverse;

	vector< vector<double> > fw;
	vector< vector<double> > bw;
	
	double fw_underflow_factor;
	int bw_underflow_count;
	
	MultiGraph* myGraph;

	// create haploid states for all edges
	map<int, map< int, set<string> > > SubLevelToUniqueKMers;

	
public:
	MultiHMM(MultiGraph* g, bool verbose=true);
	map<string, long long> estimateEmissions(string kMerCountsSamplePath, string kMerCountsGenomePath);
	double fillForwardBackwardTable(map<string, long long> globalEmission);
	multiHaploLabelPair retrieveProbabilisticSample(bool labelOnly);
};

#endif /*MultiHMM_H_*/
