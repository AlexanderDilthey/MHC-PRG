/*
 * filterLongOverlappingReads.h
 *
 *  Created on: 03.03.2015
 *      Author: AlexanderDilthey
 */

#ifndef FILTERLONGOVERLAPPINGREADS_H_
#define FILTERLONGOVERLAPPINGREADS_H_

#include "readFilter.h"
#include "../GraphAlignerUnique/coveredIntervals.h"

class filterLongOverlappingReads {
public:
	std::string input_BAM;
	std::string input_FASTQ;

	std::string output_FASTQ;
	std::string graphDir;

	std::string referenceGenomeFile;
	filterLongOverlappingReads();

	void doFilter();
	void doFilter2();

	int threads;

	int minimumOverlap;
protected:

	static bool alignmentOverlapsRegionCoveredByGraph(GraphAlignerUnique::coveredIntervals& myGraph_coveredInterval, std::string regionID, int start, int stop);

};

#endif /* FILTERLONGOVERLAPPINGREADS_H_ */
