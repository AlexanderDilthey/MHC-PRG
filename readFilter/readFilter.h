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
#include "../Utilities.h"

#include "api/BamAlignment.h"
#include "../NextGen/readSimulator.h"
#include "../GraphAligner/GraphAligner.h"


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

	std::string referenceGenomeFile;

	int k;

	int threads;

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
	BAMalignment()
	{
		likelihood_from_normalAlignment = -1;
		fromReverse = false;
		fromPosition = -1;
		toPosition = -1;
		likelihood_from_normalAlignment_to = -1;
		likelihood_from_normalAlignment_from = -1;
	}
	std::string readID;
	std::string sequence;
	std::string qualities;
	std::string fromString;

	double likelihood_from_normalAlignment;
	int likelihood_from_normalAlignment_from;
	int likelihood_from_normalAlignment_to;
	
	bool fromReverse;

	std::string fromID;
	int fromPosition;
	int toPosition;
	
	std::string getLikelihoodAndPositionString() const
	{
		assert(likelihood_from_normalAlignment != -1);
		return Utilities::DtoStr(likelihood_from_normalAlignment) + "[" + fromID + ":" + Utilities::ItoStr(likelihood_from_normalAlignment_from) + "-" + Utilities::ItoStr(likelihood_from_normalAlignment_to) + "]";
	}
};

class fastq_readPair;
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

	bool take_another_readPair(fastq_readPair* otherPair)
	{
		assert(! otherPair->isComplete());
		if(isComplete())
		{
			return false;
		}
		else
		{
			if(otherPair->have1)
			{
				return takeAlignment(otherPair->a1, 1);
			}
			if(otherPair->have2)
			{
				return takeAlignment(otherPair->a2, 2);
			}
			assert( 1 == 0 );
			return false;
		}
	}

	std::string getNormalAlignmentString() const
	{
		std::string forReturn = "normalAlignment=";
		std::vector<std::string> fields;
		fields.push_back( (a1.likelihood_from_normalAlignment != -1) ? a1.getLikelihoodAndPositionString() : "NA");
		fields.push_back( (a2.likelihood_from_normalAlignment != -1) ? a2.getLikelihoodAndPositionString() : "NA");

		std::string strandsOK = "NA";
		std::string IS = "NA";
		if(a1.fromID.length() && a2.fromID.length())
		{
			strandsOK = "0";
			if(a1.fromID == a2.fromID)
			{
				if(a1.fromReverse != a2.fromReverse)
				{
					if(a1.fromReverse == false)
					{
						if(a1.fromPosition < a2.fromPosition)
						{
							strandsOK = "1";
							int ISint = a2.fromPosition - a1.toPosition;
							IS = Utilities::ItoStr(ISint);
						}
					}
					else
					{
						if(a1.toPosition > a2.toPosition)
						{
							strandsOK = "1";
							int ISint = a1.fromPosition - a2.toPosition;
							IS = Utilities::ItoStr(ISint);
						}
					}
				}
			}
		}

		fields.push_back(strandsOK);
		fields.push_back(IS);

		forReturn += Utilities::join(fields, ";");

		return forReturn;
	}
};


std::vector<BAMRegionSpecifier> getBAMregions(std::string BAMfile);

void filterBAM(int threads, std::string BAMfile, std::string referenceGenomeFile, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);
void filterFastQPairs(int threads, std::string fastq_basePath, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);
void filterFastQPairs(int threads, std::string fastq_1_path, std::string fastq_2_path, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print);

bool transformBAMreadToInternalAlignment(const std::map<std::string, std::string>& referenceGenome, const std::string& regionID, const BamTools::BamAlignment& al, oneRead& read_forLL, seedAndExtend_return_local& alignment_forLL);

#endif /* READFILTER_H_ */
