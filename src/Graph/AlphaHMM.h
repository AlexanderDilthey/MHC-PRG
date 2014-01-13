#ifndef AlphaHMM_H_
#define AlphaHMM_H_

#include "LargeGraph.h"
#include "Edge.h"
#include "HMM.h"
#include "MultiHMM.h"
#include <vector>
#include <map>
#include "../LocusCodeAllocation.h"
#include "../MHC-PRG.h"



struct multiHaplotypePair_and_P
{
	multiHaploLabelPair haploPair;
	double P;
	diploidEdgePointerPath diploidEdgePointers;
};

class AlphaHMM
{
protected:

	vector< vector<double> > populatekMerMatrix(double coverage, map<string, long long>& globalEmission);
	double estimateCoverage(map<string, long long>& globalEmission);

	int genotypingMode;
	
	vector< map<int, Edge*> > Int2Edge;
	
	vector<double> haploid_initialProbabilities;
	vector<double> diploid_initialProbabilities;

	vector<int> diploidStatesByLevel;

	int maxStateKMerCount;
	set<string> interestingKMerStrings;
	vector< vector< map<int, int> > > diploidStateEmissions;
	vector< vector< vector<stateAlternative> > > diploidStateTransitions;

	vector< vector<double> > fw;
	vector< vector<double> > bw;
	
	vector< int > haploid_edge_lengths;
	vector< vector<int> > haploid_edge_lengths_valid;

	double fw_underflow_factor;
	int bw_underflow_count;
	
	MultiGraph* myGraph;

	void fillStateTransitions(int level);

	void fillStateTransitions(int startLevelFill, int stopLevelFill, int previousStepFillStop, int previousStepFillStart);


	// create haploid states for all edges
	map<int, map< int, set<string> > > SubLevelToUniqueKMers;

	int haploidEdgeLength(int level);
	int haploidEdgeValidLength(int level, int diploidIndex, int e12);

	double fillForwardBackwardTable_0(map<string, long long> globalEmission);
	double fillForwardBackwardTable_1(map<string, long long> globalEmission);
	double fillForwardBackwardTable_2(map<string, long long> globalEmission);
	double fillForwardBackwardTable_3(map<string, long long> globalEmission);
	double fillForwardBackwardTable_4(map<string, long long> globalEmission);
	double fillForwardBackwardTable_5(map<string, long long> globalEmission);
	double fillForwardBackwardTable_6(map<string, long long> globalEmission);
	double fillForwardBackwardTable_7(map<string, long long> globalEmission);
	double fillForwardBackwardTable_8(map<string, long long> globalEmission);

	set<string> utilizekMers;
	set<string> ignorekMers;

	kMerPositionInfo kMersInGraphInfo;
	set<string> kMersInGraph;
	kMerUniquenessInfo uniqueMers;
	double error_rate;
	double coverage;

public:

	map<Edge*, int> Edge2Int;
	map<Edge*, int> Edge2Level;
	vector< vector<haploidState> > haploidStatesByLevel;
	vector< vector<double> > fw_Viterbi;
	vector< vector<double> > fw_Viterbi_e1;
	vector< vector<double> > fw_Viterbi_e2;

	vector< vector<double> > fw_Emission;
	vector< vector<double> > fw_Emission_e1;
	vector< vector<double> > fw_Emission_e2;
	vector< vector<double> > fw_adaptedEmission;
	vector< vector<double> > fw_adaptedEmission_e1;
	vector< vector<double> > fw_adaptedEmission_e2;
	vector< vector<double> > fw_expectedMissing_e1;
	vector< vector<double> > fw_expectedMissing_e2;
	vector< vector<double> > fw_Coverage_e1;
	vector< vector<double> > fw_Coverage_e2;

	vector< vector<double> > fw_e1_present;
	vector< vector<int> > fw_e1_length;
	vector< vector<int> > fw_e1_length_valid;

	vector< vector<double> > fw_e2_present;
	vector< vector<int> > fw_e2_length;
	vector< vector<int> > fw_e2_length_valid;


	vector< vector< vector<int> > > fw_Viterbi_backtrack;
	vector< vector< int > > fw_Viterbi_backtrack_int;

	vector< vector< vector<stateAlternative> > > diploidStateTransitions_Reverse;

	void init(MultiGraph* g, int pGenotypingMode, bool verbose);

	AlphaHMM(MultiGraph* g, bool verbose=true);
	AlphaHMM(MultiGraph* g, int pGenotypingMode, bool verbose=true);

	static map<string, long long> estimateEmissions(string kMerCountsSamplePath, string kMerCountsGenomePath, int kMerSize);
	double fillForwardBackwardTable(map<string, long long> globalEmission);

	multiHaplotypePair_and_P retrieveViterbiSample(bool labelOnly);
	diploidEdgePointerPath sampleFromPosterior(map<string, long long>& observedEmissions);

	vector<diploidNucleotidePath> retrieveAllBestNucleotidePaths(map<string, long long>& observedEmissions);
	vector<diploidEdgePointerPath> retrieveAllBestPaths(map<string, long long>& observedEmissions);
	vector<diploidEdgePointerPath> retrieveAllBestPaths(map<string, long long>& observedEmissions, bool getAllPaths);

	vector<diploidEdgePath> retrieveAllBestPathsLargeGraph(map<string, long long>& observedEmissions);
	//vector< vector< string> > retrieveNucleotidePaths(map<string, long long>& observedEmissions);

	vector< vector<int> > retrieveAllBestPaths_rec(int level, int runningState, map<string, long long>& observedEmissions, bool getAllPaths);

	int getGenotypingMode();

	void debug(map<string, long long> globalEmission);

};

#endif /*AlphaHMM_H_*/
