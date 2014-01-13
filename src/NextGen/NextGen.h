
#ifndef NEXTGEN_H_
#define NEXTGEN_H_

#include <string>
#include <map>
#include <vector>
#include <set>
#include "../Graph/Graph.h"
#include "../Graph/Node.h"
#include "../Graph/Edge.h"

using namespace std;

Graph* variationGraph(string input_panel, string positions_file, bool wantPGFprotection = true);

Graph* augmentGraph(Graph* g, string input_panel, string positions_file, string locus, string HLAAlleleAlignmentDir);
void loadAndAugmentGraph(string graph_file, string input_panel, string positions_file, string locus, string HLAAlleleAlignmentDir);

LargeGraph* kMerify(Graph* g, bool quiet = false, int kMerSize = 5, bool wantPGFprotection = true);
MultiGraph* multiBeautify(LargeGraph* g, string kMerCountsGenomePath, bool quiet = false);
MultiGraph* multiBeautifyForAlpha(LargeGraph* g, string kMerCountsGenomePath, bool quiet = false);
MultiGraph* multiBeautifyForAlpha2(LargeGraph* g, string kMerCountsGenomePath, bool quiet = false);

MultiGraph* simplifyAccordingToCoverage(MultiGraph* mG, map<string, long long> estimatedEmissions, set<int> protectLevels);


Graph* HLAalignmentGraph(string hla_allele_dir, string locus);


struct kMerInfo {
	vector<unsigned char> kMer_coded;
	vector<string> kMer_deCoded;
	double p;
	vector<Edge*> traverseEdges;
	string traverseEdges_string;
	bool gapEdge;
	bool allPGF;
};

struct kMerAtNode {
	vector<unsigned char> km1Mer_coded;
	vector<string> km1Mer_deCoded;
	vector<Edge*> traverseEdges;
	string traverseEdges_string;
	Node* lastNewNode;
	bool allPGF;
};

vector<kMerInfo> forwardScan(Node* start, int limit, int firstEdgeGap);


#endif
