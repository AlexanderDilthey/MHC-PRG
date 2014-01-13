/*
 * test.h
 *
 *  Created on: 05.12.2012
 *      Author: AlexanderDilthey
 */

#ifndef TEST_H_
#define TEST_H_

#include <string>
#include <vector>

class readPair
{
public:
	std::string r1;
	std::string r2;
	int distance;
};

class debugGenomeInfo {
public:
	std::string genomeString;
	std::vector<int> startingReadsPerPosition;
	int maxNoStartingReads;
	int readLength;
};

class test {
public:
	test();
	static void testComplements();
	static void testdeBruijnGraph();
	static void testEdgeFunctions();
	static void testdeBruijnWalker();
	static void testSupernodeFunctions();
	static void testkMerDistance();

	static void testPairedBasicFunctions();
	static void testEntwirrung();
	static void testArcValidation();
	static void testArcCoverage();
	static void testArcPathSearch();

	static void testRefGenomePresence();


	static std::string generateRandomSequence(int length);
	static char randomNucleotide();

	static void pgfSimulation();
	static void testPairedMSSA476();
	static std::vector<readPair> simulateReadsFromString(std::string baseString, int readLength, int coverage,  debugGenomeInfo* dbgInfo  = 0, bool perfectCoverage = false, bool displaySummary = false);

	static void evaluateContigsAgainstBase(std::string baseString, std::vector<std::string> contigs);

	static void testContigEvaluationOnPGF();
	static void testContigEvaluation(std::string baseString);
	static std::vector<std::string>  generateRandomContigs(std::string baseString, double poissonStartPerPosition, double lengthWeibullAlpha, double lengthWeibullBeta, double errorRate);

	static void testRandomWalk();
	static std::string generateRepetitiveSequence(unsigned int numberOfSegments, unsigned int repetitiveElementLength, unsigned int interspersedSequenceLength = 800, unsigned int mediumRepetitiveNumber = 3);

	static void testStraightenAndCompress();

	void testAll();

	static int selectRandomNumber (int min, int max);

	static std::vector<std::string> findCaps(std::string existingSequence, int number, int k);

};

#endif /* TEST_H_ */
