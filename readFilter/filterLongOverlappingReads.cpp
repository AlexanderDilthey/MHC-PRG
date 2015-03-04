/*
 * filterLongOverlappingReads.cpp
 *
 *  Created on: 03.03.2015
 *      Author: AlexanderDilthey
 */

#include "filterLongOverlappingReads.h"

filterLongOverlappingReads::filterLongOverlappingReads() {
	// TODO Auto-generated constructor stub
	threads = 2;
	minimumOverlap = 0;
}

void filterLongOverlappingReads::doFilter()
{
	assert(output_FASTQ.length());
	assert(input_BAM.length());
	assert(referenceGenomeFile.length());
	assert(graphDir.length());

	std::string genomicCoverage_file = graphDir + "/genomicMapping.txt";
	if(! Utilities::fileExists(genomicCoverage_file))
	{
		std::cerr << "Error: file " << genomicCoverage_file << " not there!\n" << std::flush;
	}
	assert(Utilities::fileExists(genomicCoverage_file));

	GraphAlignerUnique::coveredIntervals myGraph_coveredIntervals;

	std::ifstream genomicCoverageStream;
	genomicCoverageStream.open(genomicCoverage_file.c_str());
	assert(genomicCoverageStream.is_open());

	std::string line;
	while(genomicCoverageStream.good())
	{
		std::getline(genomicCoverageStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> fields = Utilities::split(line, ",");
			for(unsigned int fI = 0; fI < fields.size(); fI++)
			{
				std::string fP = fields.at(fI);
				std::vector<std::string> fP_parts = Utilities::split(fP, ":");
				assert(fP_parts.size() == 2);

				std::string regionID = fP_parts.at(0);
				if(regionID == "chr6"); // if not, we have to scan the complete BAM and not just chromosome 6 - arguments to filterBAM!
				{
					int part_pos = Utilities::StrtoI(fP_parts.at(1));
					myGraph_coveredIntervals.addPoint(regionID, part_pos);
				}
			}
		}
	}
	genomicCoverageStream.close();

	myGraph_coveredIntervals.validate();

	std::cout << "myGraph_coveredIntervals: Have " << myGraph_coveredIntervals.getNumIntervals() << " intervals.\n" << std::flush;

	myGraph_coveredIntervals.reduceIntervalsBy(minimumOverlap);

	std::string fn = output_FASTQ;

	std::ofstream fastq_output;
	fastq_output.open(output_FASTQ.c_str());
	if(! fastq_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+output_FASTQ);
	}

	std::function<void(const fastq_readPair&)> printFunction = [&](const fastq_readPair& read) -> void {
			std::string normalAlignmentInfoString = ":"+read.getNormalAlignmentString();

			std::string read1_readID = read.a1.readID + normalAlignmentInfoString;		

			assert(read.a2.sequence.length() == 0);
			assert(read.a2.qualities.length() == 0);

			fastq_output << "@" << read1_readID << ":FROM:" << read.a1.fromString << "\n"
					  << read.a1.sequence    << "\n"
					  << "+"         << "\n"
					  << read.a1.qualities   << "\n";
	};

	std::function<bool(const fastq_readPair&, bool)> decisionFunction = [&](const fastq_readPair& read, bool verboseDecisionFunction) -> bool {
		bool forReturn = false;

		assert(read.a2.sequence.length() == 0);
		assert(read.a2.qualities.length() == 0);
		
		if(read.a1.fromID.length())
			forReturn = forReturn || alignmentOverlapsRegionCoveredByGraph(myGraph_coveredIntervals, read.a1.fromID, read.a1.fromPosition, read.a1.toPosition);

		
		//if(read.a2.fromID.length())
		//	forReturn = forReturn || alignmentOverlapsRegionCoveredByGraph(myGraph_coveredIntervals, read.a2.fromID, read.a2.fromPosition, read.a2.toPosition);

		return forReturn;
	};

	filterBAM(threads, input_BAM, referenceGenomeFile, output_FASTQ, &decisionFunction, &printFunction, true);
	
	fastq_output.close();
}



bool filterLongOverlappingReads::alignmentOverlapsRegionCoveredByGraph(GraphAlignerUnique::coveredIntervals& myGraph_coveredIntervals, std::string regionID, int start, int stop)
{
	std::string modifiedRegionID;
	if(regionID.substr(0, 3) == "chr")
	{
		modifiedRegionID = regionID.substr(3);
	}
	else
	{
		modifiedRegionID = "chr"+regionID;
	}

	assert(! ( myGraph_coveredIntervals.knowRegionID(regionID) && myGraph_coveredIntervals.knowRegionID(modifiedRegionID) ) );

	std::string useRegionID = (myGraph_coveredIntervals.knowRegionID(regionID)) ? regionID : modifiedRegionID;

	bool forReturn = myGraph_coveredIntervals.isThere_overlap_with_givenInterval(useRegionID, start, stop);

	return forReturn;
}

