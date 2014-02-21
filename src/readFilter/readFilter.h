/*
 * readFilter.h
 *
 *  Created on: 14.01.2014
 *      Author: AlexanderDilthey
 */

#ifndef READFILTER_H_
#define READFILTER_H_

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <assert.h>
#include <functional>


class readFilter {
public:
	std::string positiveFilter;
	std::string negativeFilter;

	double positiveThreshold;
	double negativeThreshold;

	bool positiveUnique;
	bool negativePreserveUnique;

	int positiveUnique_threshold;
	int negativePreserveUnique_threshold;

	std::string uniqueness_base;
	std::string uniqueness_subtract;

	std::string input_BAM;
	std::string input_FASTQ;

	std::string output_FASTQ;

	int k;

	readFilter();

	void doFilter();
};

class BAMRegionSpecifier {
public:
	std::string ID;
	size_t firstPos;
	size_t lastPos;
};

class BAMalignment {
public:
	std::string readID;
	std::string sequence;
	std::string qualities;
};

class fastq_readPair {
public:
	bool have1;
	bool have2;

	BAMalignment a1;
	BAMalignment a2;

	fastq_readPair() : have1(false), have2(false)
	{

	}

	bool takeAlignment(BAMalignment a, int which)
	{
		bool success = false;
		assert((which == 1) || (which == 2));
		if(which == 1)
		{
			if(have1 == false)
			{
				success = true;
				a1 = a;
				have1 = true;
			}
		}
		else if(which == 2)
		{
			if(have2 == false)
			{
				success = true;
				a2 = a;
				have2 = true;
			}
		}
		return success;
	}

	bool isComplete()
	{
		return (have1 && have2);
	}
};


std::vector<BAMRegionSpecifier> getBAMregions(std::string BAMfile);

void filterBAM(std::string BAMfile, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);
void filterFastQPairs(std::string fastq_basePath, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);
void filterFastQPairs(std::string fastq_1_path, std::string fastq_2_path, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);

#endif /* READFILTER_H_ */
