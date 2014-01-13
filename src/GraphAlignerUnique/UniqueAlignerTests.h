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
}
};

#endif /* UNIQUEALIGNERTESTS_H_ */
