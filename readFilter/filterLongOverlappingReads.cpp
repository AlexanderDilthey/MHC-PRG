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


void filterLongOverlappingReads::doFilter2()
{
	assert(output_FASTQ.length());
	assert(input_FASTQ.length());
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

	std::function<void(const fastq_readPair&)> printFunction = [&](const fastq_readPair& read) -> void {
	};

	std::set<std::string> want_readIDs;
	std::function<bool(const fastq_readPair&, bool)> decisionFunction = [&](const fastq_readPair& read, bool verboseDecisionFunction) -> bool {
		assert(read.a2.sequence.length() == 0);
		assert(read.a2.qualities.length() == 0);
		
		assert(read.a1.readID.length());

		if(alignmentOverlapsRegionCoveredByGraph(myGraph_coveredIntervals, read.a1.fromID, read.a1.fromPosition, read.a1.toPosition))
		{ 
			want_readIDs.insert("@"+read.a1.readID);
		}

		
		//if(read.a2.fromID.length())
		//	forReturn = forReturn || alignmentOverlapsRegionCoveredByGraph(myGraph_coveredIntervals, read.a2.fromID, read.a2.fromPosition, read.a2.toPosition);

		return false;
	};

	filterBAM(threads, input_BAM, referenceGenomeFile, output_FASTQ, &decisionFunction, &printFunction, true);
	
	std::cout << Utilities::timestamp() << "filterLongOverlappingReads::doFilter2(): Now searching " << want_readIDs.size() << " read IDs.\n" << std::flush;
	int have_readIDs = 0;
	for(std::set<std::string>::iterator readIDIt = want_readIDs.begin(); readIDIt != want_readIDs.end(); readIDIt++)
	{
		have_readIDs++;
		if(have_readIDs > 5)
		{
			break;
		}
		std::cout << "\t - " << *readIDIt << "\n" << std::flush;
	}
	
	std::ofstream fastq_output;
	fastq_output.open(output_FASTQ.c_str());
	if(! fastq_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+output_FASTQ);
	}
	
	std::ifstream fastq_input;
	fastq_input.open(input_FASTQ.c_str());
	assert(fastq_input.is_open());
	
	
	auto getLinesFromFastQ = [](std::ifstream& inputStream, unsigned int lines) -> std::vector<std::string> {
		std::vector<std::string> forReturn;
		assert(lines > 0);
		for(unsigned int lI = 0; lI < lines; lI++)
		{
			if(!inputStream.good())
			{
				return forReturn;
			}
			assert(inputStream.good());
			std::string thisLine;
			std::getline(inputStream, thisLine);
			Utilities::eraseNL(thisLine);
			forReturn.push_back(thisLine);
		}
		assert(forReturn.size() == lines);
		return forReturn;
	};

	auto getReadFromFastQ = [&](std::ifstream& inputStream, std::string& ret_readID, std::string& ret_sequence, std::string& ret_qualities) -> void {
		assert(inputStream.good());
		std::vector<std::string> lines = getLinesFromFastQ(inputStream, 4);
		if(lines.size() == 4)
		{
			assert(lines.at(2) == "+");
			ret_readID = lines.at(0);
			ret_sequence = lines.at(1);
			ret_qualities = lines.at(3);
			if(!(ret_sequence.length() == ret_qualities.length()))
			{
				std::cerr << "Problem with input file stream\n";
				std::cerr << "\t" << "ret_sequence.length()" << ": " << ret_sequence.length() << "\n";
				std::cerr << "\t" << "ret_qualities.length()" << ": " << ret_qualities.length() << "\n";
				std::cerr << "\t" << "ret_sequence" << ": " << ret_sequence << "\n";
				std::cerr << "\t" << "ret_qualities" << ": " << ret_qualities << "\n";
				std::cerr << std::flush;
				
			}
			assert(ret_sequence.length() == ret_qualities.length());
		}
		else
		{
			ret_readID.clear();
			ret_sequence.clear();
			ret_qualities.clear();
		}
	};

	std::set<std::string> got_readIDs;
	std::set<std::string> missing_readIDs;
	
	while(fastq_input.good())
	{
		std::string read_ID; std::string read_sequence; std::string read_qualities;
		getReadFromFastQ(fastq_input, read_ID, read_sequence, read_qualities);
		if(want_readIDs.count(read_ID))
		{
			fastq_output << read_ID << "\n";
			fastq_output << read_sequence << "\n";
			fastq_output << "+" << "\n";
			fastq_output << read_qualities << "\n";
			
			got_readIDs.insert(read_ID);
		}
	}
		
	fastq_input.close();
	
	fastq_output.close();
	
	for(std::set<std::string>::iterator readIDIt = want_readIDs.begin(); readIDIt != want_readIDs.end(); readIDIt++)
	{
		if(got_readIDs.count(*readIDIt) == 0)
		{
			missing_readIDs.insert(*readIDIt);
		}
	}
	
	std::cout << Utilities::timestamp() << "filterLongOverlappingReads::doFilter2(): Found " << got_readIDs.size() << " / " << want_readIDs.size() << " read IDs, " << missing_readIDs.size() << " missing.\n" << std::flush;
	
	have_readIDs = 0;
	for(std::set<std::string>::iterator readIDIt = missing_readIDs.begin(); readIDIt != missing_readIDs.end(); readIDIt++)
	{
		have_readIDs++;
		if(have_readIDs > 5)
		{
			break;
		}
		std::cout << "\t - " << *readIDIt << "\n" << std::flush;
	}
	
	
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

