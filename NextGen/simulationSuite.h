/*
 * simulationSuite.h
 *
 *  Created on: 17 Apr 2012
 *      Author: dilthey
 */

#ifndef SIMULATIONSUITE_H_
#define SIMULATIONSUITE_H_

#include <string>
#include <map>
#include <vector>
#include <set>
#include "../Graph/Graph.h"
#include "../Graph/Node.h"
#include "../Graph/Edge.h"

using namespace std;

void simulationSuite(string graph_file, string temp_dir, string temp_label, int genotypingMode);
int agreementStats(string estimated_1, string estimated_2, string real_1, string real_2);
int agreementStats(int estimated_1, int estimated_2, int real_1, int real_2);
vector<int> agreementStarStats(string estimated_1, string estimated_2, string real_1, string real_2);
void describeGraph(string graph_file, string temp_dir, string temp_label, bool output_kMer_levels, std::string referenceGenomeCortexGraph);
void describeGraph_25(string graph_file, string temp_dir, string temp_label, bool output_kMer_levels, std::string referenceGenomeCortexGraph);

void describeNucleotideGraph(string graph_file, string temp_dir, string temp_label);

#endif /* SIMULATIONSUITE_H_ */
