/*
 * UniqueAlignerTests.h
 *
 *  Created on: 30.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef UNIQUEALIGNERTESTS_H_
#define UNIQUEALIGNERTESTS_H_

#include <string>
#include "../Graph/Graph.h"
#include <vector>
#include <map>

namespace GraphAlignerUnique {
namespace tests {
	void testChains();

	void sampleExactStringFromGraph(Graph* g, int minLength_string, int maxLength_string, std::string& string_ret, std::vector<Edge*>& traversedEdges_ret);
	void testSeedAndExtend();
	void testSeedAndExtend_local();
	void testSeedAndExtend_local_realGraph(std::string graph_filename, int read_length, double insertSize_mean, double insertSize_sd, std::string qualityMatrixFile);

}
};

#endif /* UNIQUEALIGNERTESTS_H_ */
