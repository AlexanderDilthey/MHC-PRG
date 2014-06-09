/*
 * AlignerTests.h
 *
 *  Created on: 27.06.2013
 *      Author: AlexanderDilthey
 */

#ifndef ALIGNERTESTS_H_
#define ALIGNERTESTS_H_

#include "GraphAligner.h"
#include "GraphAndIndex.h"

void testGraphAligner();

diploidGenomeString generateRandomGenome(int minimumLength, bool beginAndEndOneLevel = false);
void _printDiploidGenomeString(diploidGenomeString& gS);
void sampleStringFromGraph(Graph* g, std::string& edgeLabels_ret, int& string_begin_ret, int& string_end_ret, int min_string_length = 10, int max_string_length = 20);
// void sampleStringFromGraph_forSeedAndExtend(Graph* g, std::string& edgeLabels_ret, std::string& underlyingEdges_ret, std::vector<int>& underlyingEdges_levels_ret, int& string_begin_ret, int& string_end_ret, int minTrunkLength, int maxTrunkLength, int minTrunkExtensionsEachSide, int maxSideTrunkExtensions, int minIntermediatePenalty, int maxIntermediatePenalty, bool guaranteeSmoothExtensibility, bool returnWholeGenomeString);
void sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(Graph* g, int kMerSize, int S_match, int S_mismatch, int S_gapOpen, int S_gapExtend, std::string& sequenceLabels_ret, std::string& underlyingEdgeLabels_ret, std::vector<int>& underlyingEdges_levels_ret, int& string_begin_ret, int& string_end_ret, int minTrunkLength, int maxTrunkLength, int minTrunks, int maxTrunks, int minSideTrunkExtensions, int maxSideTrunkExtensions, int minIntermediatePenalty, int maxIntermediatePenalty, int betweenTrunkGapMin, int betweenTrunkGapMax, bool returnWholeGenomeString);

void sampleStringFromGraph_forSeedAndExtend(Graph* g, int kMerSize, int S_match, int S_mismatch, int S_gapOpen, int S_gapExtend, std::string& edgeLabels_ret, std::string& underlyingEdges_ret, std::vector<int>& underlyingEdges_levels_ret, int& string_begin_ret, int& string_end_ret, int minTrunkLength = 10, int maxTrunkLength = 15, int minSideTrunkExtensions = 1, int maxSideTrunkExtensions = 3, int minIntermediatePenalty = 3, int maxIntermediatePenalty = 7, int minExtensionStaticRegionLength = 4, int maxExtensionStaticRegionLength = 4, bool guaranteeSmoothExtensibility = true, bool returnWholeGenomeString = false);

void testSeedAndExtend();
void test_Simple_longRange_SeedAndExtend();
void test_Simple_longRange_SeedAndExtend_2();
void test_Simple_longRange_SeedAndExtend_3();
void test_Simple_longRange_SeedAndExtend_4();


void testSeedAndExtend_Algorithm();

void testDiagonalInverseExtendAlignment();
void testChainFindingAndExtension();
void testExtensionStop();

extern int lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret;
extern int lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret;
extern int lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret;


extern int lastCall_sampleString_forSeedAndExtend_trunk_end_sequence_ret;


#endif /* ALIGNERTESTS_H_ */
