/*
 * Validation.cpp
 *
 *  Created on: 13.06.2013
 *      Author: AlexanderDilthey
 */

#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <assert.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <omp.h>
#include <functional>
#include <cctype>
#include <set>

#include "Validation.h"
#include "../Utilities.h"
#include <boost/math/distributions/poisson.hpp>

#include "../hash/deBruijn/DeBruijnGraph.h"
#include "../GraphAligner/GraphAlignerAffine.h"
#include "../GraphAlignerUnique/GraphAlignerUnique.h"

#include <stdio.h>
#include "dirent.h"

using namespace boost::math::policies;
using namespace boost::math;

typedef boost::math::poisson_distribution< double, policy < discrete_quantile < integer_round_inwards> > > poisson_up;

// predeclarations

template<int m, int k, int colours>
void evaluate_dGS(diploidGenomeString& gS, diploidGenomeString& gS_unresolved, const std::set<std::string>& kMers_reference, DeBruijnGraph<m, k, colours>* graph, std::set<std::string>* kMers_in_dGS, std::set<std::string>* kMers_in_dGS_in_sample, std::map<std::string, double>* kMers_in_dGS_optimality, std::string nameForSummary, ofstream& summaryFileStream, std::string pathForSpatialSummary, std::vector<std::vector<int> > chromotypes_referencePositions);

template<int m, int k, int colours>
std::pair<diploidGenomeString, diploidGenomeString> greedilyResolveDiploidKMerString(diploidGenomeString& original_gS, DeBruijnGraph<m, k, colours>* graph);

void alignedShortReads2SAM(std::ofstream& SAMoutputStream, std::vector<int>& uncompressed_graph_referencePositions, std::string& referenceSequence, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<oneReadPair>& originalReads);
// void read_shortReadAlignments_fromFile (std::string file, std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>>& ret_alignments, std::vector<oneReadPair>& ret_alignments_originalReads, double& ret_IS_mean, double& ret_IS_sd);

// functions

std::vector<std::string> filesInDirectory(std::string path)
{
	std::vector<std::string> forReturn;

    struct dirent *pDirent;
    DIR *pDir;

    pDir = opendir (path.c_str());
    if (pDir == NULL) {
        throw std::runtime_error("Cannot open directory " + path);
    }

    while ((pDirent = readdir(pDir)) != NULL) {
    	std::string fileName(pDirent->d_name);

    	if((fileName == ".") || (fileName == ".."))
    	{
    		continue;
    	}

    	fileName = path + "/" + fileName;

    	// std::cout << "\t" << fileName << "\n";

    	forReturn.push_back(fileName);
    }
    closedir (pDir);

    return forReturn;
}




std::vector<oneRead> getUnpairedReadsFromFastQ(std::string fastq_path)
{
	std::vector<oneRead> forReturn;

	std::ifstream fastQ_stream;
	fastQ_stream.open(fastq_path.c_str());
	assert(fastQ_stream.is_open());


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
			assert(ret_sequence.length() == ret_qualities.length());
		}
		else
		{
			ret_readID.clear();
			ret_sequence.clear();
			ret_qualities.clear();
		}
	};

	while(fastQ_stream.good())
	{
		std::string read1_ID; std::string read1_sequence; std::string read1_qualities;
		getReadFromFastQ(fastQ_stream, read1_ID, read1_sequence, read1_qualities);

		std::string read1_ID_noFrom = Utilities::removeFROM(read1_ID);

		oneRead r1(read1_ID, read1_sequence, read1_qualities);
		forReturn.push_back(r1);
	}

	return forReturn;
}

std::vector<oneReadPair> getReadsFromFastQ(std::string fastq_base_path)
{
	std::string f1 = fastq_base_path + "_1";
	std::string f2 = fastq_base_path + "_2";

	if(! Utilities::fileReadable(f1))
	{
		throw std::runtime_error("Expected file "+f1+" can't be opened.");
	}

	if(! Utilities::fileReadable(f2))
	{
		throw std::runtime_error("Expected file "+f2+" can't be opened.");
	}

	return getReadsFromFastQ(f1, f2);
}


std::vector<oneReadPair> getReadsFromFastQ(std::string fastq_1_path, std::string fastq_2_path)
{
	std::vector<oneReadPair> forReturn;

	std::ifstream fastQ_1_stream;
	fastQ_1_stream.open(fastq_1_path.c_str());
	assert(fastQ_1_stream.is_open());

	std::ifstream fastQ_2_stream;
	fastQ_2_stream.open(fastq_2_path.c_str());
	assert(fastQ_2_stream.is_open());

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
			assert(ret_sequence.length() == ret_qualities.length());
		}
		else
		{
			ret_readID.clear();
			ret_sequence.clear();
			ret_qualities.clear();
		}
	};

	while(fastQ_1_stream.good())
	{
		assert(fastQ_2_stream.good());

		std::string read1_ID; std::string read1_sequence; std::string read1_qualities;
		getReadFromFastQ(fastQ_1_stream, read1_ID, read1_sequence, read1_qualities);

		std::string read2_ID; std::string read2_sequence; std::string read2_qualities;
		getReadFromFastQ(fastQ_2_stream, read2_ID, read2_sequence, read2_qualities);

		assert((read1_ID.length() && read2_ID.length()) || ((!read1_ID.length()) && (!read2_ID.length())));
		if((!read1_ID.length()) && (!read2_ID.length()))
		{
			break;
		}

		std::string read1_ID_noFrom = Utilities::removeFROM(read1_ID);
		std::string read2_ID_noFrom = Utilities::removeFROM(read2_ID);
		assert((read1_ID_noFrom.substr(read1_ID_noFrom.length() - 2, 2) == "/1") || (read1_ID_noFrom.substr(read1_ID_noFrom.length() - 2, 2) == "/2"));
		assert((read2_ID_noFrom.substr(read2_ID_noFrom.length() - 2, 2) == "/1") || (read2_ID_noFrom.substr(read2_ID_noFrom.length() - 2, 2) == "/2"));
		if(!(read1_ID_noFrom.substr(0, read1_ID_noFrom.length() - 2) == read2_ID_noFrom.substr(0, read2_ID_noFrom.length() - 2)))
		{
			std::cerr << "Warning: read IDs don't match! " << read1_ID_noFrom << " vs " << read2_ID_noFrom << "\n";
		}
		assert(read1_ID_noFrom.substr(0, read1_ID_noFrom.length() - 2) == read2_ID_noFrom.substr(0, read2_ID_noFrom.length() - 2));

		oneRead r1(read1_ID, read1_sequence, read1_qualities);
		oneRead r2(read2_ID, read2_sequence, read2_qualities);
		oneReadPair thisPair(r1, r2, 0);
		forReturn.push_back(thisPair);
	}

	return forReturn;
}

void alignedShortReads2SAM(std::ofstream& SAMoutputStream, std::vector<int>& uncompressed_graph_referencePositions, std::string& referenceSequence, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<oneReadPair>& originalReads)
{
	assert(SAMoutputStream.is_open());

	for(unsigned int lI = 0; lI < uncompressed_graph_referencePositions.size(); lI++)
	{
		int referencePosition = uncompressed_graph_referencePositions.at(lI);
		if(!((referencePosition == -1) || ((referencePosition > 0) && (referencePosition <= referenceSequence.length()))))
		{
			std::cerr << "EARLY CHECK!" << "\n";
			std::cerr << "referencePosition" << ": " << referencePosition << "\n";
			std::cerr << "lI" << ": " << lI << "\n";
			std::cerr << "referenceSequence.length()" << ": " << referenceSequence.length() << "\n" << std::flush;
		}
		assert((referencePosition == -1) || ((referencePosition > 0) && (referencePosition <= referenceSequence.length())));
	}

	auto singleAlignment2SAM = [&](seedAndExtend_return_local& alignment, oneRead& originalRead, bool printToStream, int otherReadReferencePosition, bool isFirstRead) -> int {
		std::string alignment_graph = alignment.graph_aligned;
		std::string alignment_sequence = alignment.sequence_aligned;
		std::vector<int> levels_separated = alignment.graph_aligned_levels;

		assert((otherReadReferencePosition == -1) || (otherReadReferencePosition > 0));
		assert(levels_separated.size() == alignment_sequence.length());

		std::vector<int> levels_separated_2_reference;
		levels_separated_2_reference.resize(levels_separated.size(), -1);

		int lastUsedLevel = -1;
		for(unsigned int levelI = 0; levelI < levels_separated.size(); levelI++)
		{
			int level = levels_separated.at(levelI);
			std::string graphCharacter = alignment_graph.substr(levelI, 1);
			std::string sequenceCharacter = alignment_sequence.substr(levelI, 1);

			if(level != -1)
			{
				assert(level < uncompressed_graph_referencePositions.size());
				int referencePosition = uncompressed_graph_referencePositions.at(level);

				if(!((referencePosition == -1) || ((referencePosition > 0) && (referencePosition <= referenceSequence.length()))))
				{
					std::cerr << "referencePosition" << ": " << referencePosition << "\n";
					std::cerr << "levelI" << ": " << levelI << "\n";
					std::cerr << "level" << ": " << level << "\n";
					std::cerr << "referenceSequence.length()" << ": " << referenceSequence.length() << "\n" << std::flush;
				}

				assert((referencePosition == -1) || ((referencePosition > 0) && (referencePosition <= referenceSequence.length())));

				levels_separated_2_reference.at(levelI) = referencePosition;
			}
		}

		std::vector<std::string> levels_separated_2_reference_String;
		std::vector<std::string> levels_referenceCharacters;
		std::vector<std::string> levels_graphCharacters;
		std::vector<std::string> levels_readCharacters;

		int firstReferencePosition = -1;
		int lastReferencePosition = -1;
		for(unsigned int levelI = 0; levelI < levels_separated_2_reference.size(); levelI++)
		{
			levels_separated_2_reference_String.push_back(Utilities::ItoStr(levels_separated_2_reference.at(levelI)));
			int correspondingReferencePosition = levels_separated_2_reference.at(levelI);
			if(correspondingReferencePosition != -1)
			{
				if(! (((correspondingReferencePosition - 1) >= 0) && ((correspondingReferencePosition - 1) < referenceSequence.length())))
				{
					std::cerr << "(correspondingReferencePosition - 1)" << ":" << (correspondingReferencePosition - 1) << "\n";
					std::cerr << "referenceSequence.length()" << ":" << referenceSequence.length() << "\n" << std::flush;
				}

				assert((correspondingReferencePosition - 1) >= 0);
				assert((correspondingReferencePosition - 1) < referenceSequence.length());

				levels_referenceCharacters.push_back(referenceSequence.substr(correspondingReferencePosition - 1, 1));
			}
			else
			{
				levels_referenceCharacters.push_back("_");
			}
			levels_graphCharacters.push_back(alignment_graph.substr(levelI, 1));
			levels_readCharacters.push_back(alignment_sequence.substr(levelI, 1));

			if(firstReferencePosition == -1)
			{
				if((levels_separated_2_reference.at(levelI) != -1) && (alignment_sequence.substr(levelI, 1) != "_"))
				{
					firstReferencePosition = levels_separated_2_reference.at(levelI);
				}
			}

			if((levels_separated_2_reference.at(levelI) != -1) && (alignment_sequence.substr(levelI, 1) != "_"))
			{
				lastReferencePosition = levels_separated_2_reference.at(levelI);
			}

		}

		// assert((firstReferencePosition != -1) && (lastReferencePosition != -1));

		// std::cout << "Contig: " << contigFullID << "\n";
		// std::cout << Utilities::join(levels_separated_2_reference_String, " ") << "\n";
		// std::cout << Utilities::join(levels_referenceCharacters, " ") << "\n";
		// std::cout << Utilities::join(levels_graphCharacters, " ") << "\n";
		// std::cout << Utilities::join(levels_contigCharacters, " ") << "\n\n" << std::flush;

		std::string readIDForSAM = originalRead.name;
		readIDForSAM.erase(std::remove_if(readIDForSAM.begin(), readIDForSAM.end(), [&](char c){return (((int)isalnum(c) == 0) || (c == ' '));}), readIDForSAM.end());

		std::string readCharacters_noGaps = alignment_sequence;
		readCharacters_noGaps.erase(std::remove_if(readCharacters_noGaps.begin(),readCharacters_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), readCharacters_noGaps.end());

		size_t FLAGS_sizeT = 0;
		bool bothReadsProperlyAligned = false;
		//FLAGS_sizeT = (FLAGS_sizeT || 0x2);
//		if((firstReferencePosition != -1) || (lastReferencePosition != -1))
//		{
//			if(isFirstRead)
//			{
//				FLAGS_sizeT = (FLAGS_sizeT || 0x40);
//			}
//			else
//			{
//				FLAGS_sizeT = (FLAGS_sizeT || 0x80);
//			}
//
//			if(otherReadReferencePosition)
//			{
//				FLAGS_sizeT = (FLAGS_sizeT || 0x2);
//				bothReadsProperlyAligned = true;
//			}
//		}
//		else
//		{
//			FLAGS_sizeT = (FLAGS_sizeT || 0x4);
//			bothReadsProperlyAligned = false;
//		}

		if((firstReferencePosition != -1) || (lastReferencePosition != -1))
		{
			FLAGS_sizeT = (FLAGS_sizeT || 0x4);
		}

		int lastPrintedRealReferencePosition = -1;
		std::string CIGAR_uncompressed;
		for(unsigned int levelI = 0; levelI < levels_separated_2_reference.size(); levelI++)
		{
			int reference_position = levels_separated_2_reference.at(levelI);
			std::string alignmentCharacter = alignment_sequence.substr(levelI, 1);

			if(reference_position != -1)
			{
				if(lastPrintedRealReferencePosition != -1)
				{
					if(reference_position != (lastPrintedRealReferencePosition + 1))
					{
						if(!(reference_position > lastPrintedRealReferencePosition))
						{
							std::cerr << "reference_position" << ": " << reference_position << "\n";
							std::cerr << "lastPrintedRealReferencePosition" << ": " << lastPrintedRealReferencePosition << "\n" << std::flush;
						}
						assert(reference_position > lastPrintedRealReferencePosition);
						int missingReferencePositions = reference_position - lastPrintedRealReferencePosition - 1;
						assert(missingReferencePositions > 0);
						for(unsigned int j = 0; (int)j < missingReferencePositions; j++)
						{
							CIGAR_uncompressed.push_back('D');
						}
					}
				}
			}

			if(reference_position == -1)
			{
				if(alignmentCharacter == "_")
				{
					// nothing
				}
				else
				{
					CIGAR_uncompressed.push_back('I');
				}
			}
			else
			{
				if(alignmentCharacter == "_")
				{
					CIGAR_uncompressed.push_back('D');
				}
				else
				{
					CIGAR_uncompressed.push_back('M');
				}
			}

			if(reference_position != -1)
			{
				lastPrintedRealReferencePosition = reference_position;
			}
		}

		std::string CIGAR_compressed;
		size_t CIGAR_sum_M = 0;
		size_t CIGAR_sum_I = 0;
		std::string currentOperation;
		int operationCount = 0;
		for(unsigned int cigarI = 0; cigarI < CIGAR_uncompressed.size(); cigarI++)
		{
			std::string thisOperation = CIGAR_uncompressed.substr(cigarI, 1);
			if(currentOperation.size() == 0)
			{
				currentOperation = thisOperation;
			}
			operationCount++;
			if((cigarI == (CIGAR_uncompressed.size() - 1)) || (CIGAR_uncompressed.substr(cigarI+1, 1) != currentOperation))
			{
				CIGAR_compressed += (Utilities::ItoStr(operationCount) + currentOperation);

				if(currentOperation == "M")
				{
					CIGAR_sum_M += operationCount;
				}
				else if(currentOperation == "I")
				{
					CIGAR_sum_I += operationCount;
				}
				operationCount = 0;
				currentOperation = "";
			}
		}

		assert(readCharacters_noGaps.length() == (CIGAR_sum_M + CIGAR_sum_I));

		if(printToStream)
		{
			double P_mapping_wrong = 0.001;

			std::string QUALforSAM = originalRead.quality;
			if(alignment.reverse)
			{
				std::reverse(QUALforSAM.begin(), QUALforSAM.end());

			}
			std::string QNAME = readIDForSAM;
			std::string FLAG = Utilities::ItoStr(FLAGS_sizeT);
			std::string RNAME = "ref";
			std::string POS = Utilities::ItoStr(firstReferencePosition);
			std::string MAPQ = Utilities::ItoStr(-10.0*log10(P_mapping_wrong)+0.5);
			std::string CIGAR = CIGAR_compressed;
			std::string RNEXT = "*";
			std::string PNEXT = "0";
			std::string TLEN = Utilities::ItoStr(lastReferencePosition - firstReferencePosition + 1);
			std::string SEQ = readCharacters_noGaps;
			std::string QUAL = QUALforSAM;

			SAMoutputStream <<
				QNAME << "\t" <<
				FLAG << "\t" <<
				RNAME << "\t" <<
				POS << "\t" <<
				MAPQ << "\t" <<
				CIGAR << "\t" <<
				RNEXT << "\t" <<
				PNEXT << "\t" <<
				TLEN << "\t" <<
				SEQ << "\t" <<
				QUAL << "\n" << std::flush;

//			double P_mapping_wrong = 1 - alignment.mapQ;
//
//			std::string QUALforSAM = originalRead.quality;
//			if(alignment.reverse)
//			{
//				std::reverse(QUALforSAM.begin(), QUALforSAM.end());
//			}
//			std::string QNAME = readIDForSAM;
//			std::string FLAG = Utilities::ItoStr(FLAGS_sizeT);
//			std::string RNAME = "ref";
//			std::string POS = Utilities::ItoStr(firstReferencePosition);
//			std::string MAPQ = Utilities::ItoStr(-10.0*log10(P_mapping_wrong)+0.5);
//			std::string CIGAR = CIGAR_compressed;
//			std::string RNEXT = (bothReadsProperlyAligned) ? "*" : "=";
//			std::string PNEXT = Utilities::ItoStr(otherReadReferencePosition);
//			std::string TLEN = Utilities::ItoStr(lastReferencePosition - firstReferencePosition + 1);
//			std::string SEQ = readCharacters_noGaps;
//			std::string QUAL = QUALforSAM;
//
//			SAMoutputStream <<
//				QNAME << "\t" <<
//				FLAG << "\t" <<
//				RNAME << "\t" <<
//				POS << "\t" <<
//				MAPQ << "\t" <<
//				CIGAR << "\t" <<
//				RNEXT << "\t" <<
//				PNEXT << "\t" <<
//				TLEN << "\t" <<
//				SEQ << "\t" <<
//				QUAL << "\n" << std::flush;
		}

		return firstReferencePosition;
	};

	for(unsigned int alignmentI = 0; alignmentI < alignments.size(); alignmentI++)
	{

		std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& thisPair_alignment = alignments.at(alignmentI);
		// std::pair<std::string, std::string>& thisPair_readIDs = alignments_readIDs.at(alignmentI);
		
		int firstRead_firstReferencePosition = singleAlignment2SAM(thisPair_alignment.first, originalReads.at(alignmentI).reads.first, false, -1, 1);
		int secondRead_firstReferencePosition = singleAlignment2SAM(thisPair_alignment.second, originalReads.at(alignmentI).reads.second, false, -1, 0);

		singleAlignment2SAM(thisPair_alignment.first, originalReads.at(alignmentI).reads.first, true, secondRead_firstReferencePosition, 1);
		singleAlignment2SAM(thisPair_alignment.second, originalReads.at(alignmentI).reads.second, true, firstRead_firstReferencePosition, 0);
	}

/*
1 QNAME String [!-?A-~]f1,255g Query template NAME
2 FLAG Int [0,216-1] bitwise FLAG
3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
4 POS Int [0,231-1] 1-based leftmost mapping POSition
5 MAPQ Int [0,28-1] MAPping Quality
6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
8 PNEXT Int [0,231-1] Position of the mate/next read
9 TLEN Int [-231+1,231-1] observed Template LENgth
10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
*/

}

void read_shortReadAlignments_fromFile (std::string file, std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>>& ret_alignments, std::vector<oneReadPair>& ret_alignments_originalReads, double& ret_IS_mean, double& ret_IS_sd)
{
	ret_alignments.clear();
	ret_alignments_originalReads.clear();
	
	std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> forReturn;

	auto getLines = [](std::ifstream& inputStream, unsigned int lines) -> std::vector<std::string> {
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


	auto getAlignment = [&](std::ifstream& inputStream, bool& fail, seedAndExtend_return_local& ret_alignment, oneRead& ret_originalRead, std::string insertFirstLine) -> void {
		fail = false;
		int readLines = 9;
		if(insertFirstLine.length())
		{
			readLines = 8;
		}
		std::vector<std::string> lines = getLines(inputStream, readLines);
		
		if(insertFirstLine.length())
		{
			lines.insert(lines.begin(), insertFirstLine);
		}
		
		if(lines.size() != 9)
		{
			fail = true;
		}
		else
		{

			if(!(lines.at(0).substr(0, 5) == "\tRead"))
			{
				std::cerr << "Line 0 should be TABRead, but is not!\n" << lines.at(0) << "\n" << std::flush;
			}
			assert(lines.at(0).substr(0, 6) == "\tRead ");

			std::string str_readID = lines.at(0).substr(6);
			std::string str_score = lines.at(1).substr(2);
			std::string str_reverse = lines.at(2).substr(2);
			std::string str_mapQ = lines.at(3).substr(2);
			std::string str_graph_aligned = lines.at(4).substr(2);
			std::string str_sequence_aligned = lines.at(5).substr(2);
			std::string str_levels = lines.at(6).substr(2);
			std::string str_originalSequence = lines.at(7).substr(2);
			std::string str_qualities = lines.at(8).substr(2);

			std::vector<std::string> mapQs = Utilities::split(str_mapQ, " ");

			seedAndExtend_return_local a;
			a.Score = Utilities::StrtoD(str_score);
			a.reverse = Utilities::StrtoB(str_reverse);
			if(mapQs.size() == 2)
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = Utilities::StrtoD(mapQs.at(1));
			}
			else if(mapQs.size() == 3)
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = Utilities::StrtoD(mapQs.at(1));
				a.mapQ_genomic_perPosition = mapQs.at(2);
				assert(a.mapQ_genomic_perPosition.length() == str_graph_aligned.length());
			}
			else
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = 2;
			}
			a.graph_aligned = str_graph_aligned;
			a.sequence_aligned = str_sequence_aligned;
			a.graph_aligned_levels = Utilities::StrtoI(Utilities::split(str_levels, " "));			
			ret_alignment = a;
			
			oneRead r(str_readID, str_originalSequence, str_qualities);
			ret_originalRead = r;			
		}
	};
	std::ifstream inputStream;
	inputStream.open(file.c_str());
	assert(inputStream.is_open());
	std::string line;
	assert(inputStream.good());
	std::getline(inputStream, line);
	
	std::vector<std::string> firstLine_fields = Utilities::split(line, " ");
	assert(firstLine_fields.size() == 3);
	assert(firstLine_fields.at(0) == "IS");

	ret_IS_mean = Utilities::StrtoD(firstLine_fields.at(1));
	ret_IS_sd = Utilities::StrtoD(firstLine_fields.at(2));
	
	while(inputStream.good())
	{
		std::getline(inputStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{

			std::string lineForInsertion;
			
			assert( (line.substr(0, std::string("Aligned pair").length()) == "Aligned pair") ||
					(line.substr(0, std::string("\tRead").length()) == "\tRead")
			);
			
			if(line.substr(0, std::string("\tRead").length()) == "\tRead")
			{
				lineForInsertion = line;
			}
			
			bool fail_1; bool fail_2;

			seedAndExtend_return_local a1;
			seedAndExtend_return_local a2;
			oneRead r1("", "", "");
			oneRead r2("", "", "");
			
			getAlignment(inputStream, fail_1, a1, r1, lineForInsertion);
			getAlignment(inputStream, fail_2, a2, r2, "");

			oneReadPair rP(r1, r2, 0);			 
			if(fail_1 || fail_2)
			{
				assert(! inputStream.good());
			}

			if((! fail_1) && (! fail_2))
			{
				ret_alignments.push_back(make_pair(a1, a2));
				ret_alignments_originalReads.push_back(rP);
			}
		}
	}
}

void read_longReadAlignments_fromFile (std::string file, std::vector<seedAndExtend_return_local>& ret_alignments, std::vector<oneRead>& ret_alignments_originalReads)
{
	ret_alignments.clear();
	ret_alignments_originalReads.clear();

	std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> forReturn;

	size_t got_lines = 0;
	
	auto getLines = [](std::ifstream& inputStream, unsigned int lines) -> std::vector<std::string> {
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


	auto getAlignment = [&](std::ifstream& inputStream, bool& fail, seedAndExtend_return_local& ret_alignment, oneRead& ret_originalRead) -> void {
		fail = false;
		int readLines = 9;
		std::vector<std::string> lines = getLines(inputStream, readLines);

		got_lines += lines.size();
		
		if(lines.size() != 9)
		{
			fail = true;
		}
		else
		{
			if(!(lines.at(0).substr(0, 5) == "\tRead"))
			{
				std::cerr << "Line 0 should be TABRead, but is not!\n" << lines.at(0) << "\n" << std::flush;
			}
			assert(lines.at(0).substr(0, 6) == "\tRead ");

			std::string str_readID = lines.at(0).substr(6);
			std::string str_score = lines.at(1).substr(2);
			std::string str_reverse = lines.at(2).substr(2);
			std::string str_mapQ = lines.at(3).substr(2);
			std::string str_graph_aligned = lines.at(4).substr(2);
			std::string str_sequence_aligned = lines.at(5).substr(2);
			std::string str_levels = lines.at(6).substr(2);
			std::string str_originalSequence = lines.at(7).substr(2);
			std::string str_qualities = lines.at(8).substr(2);

			std::vector<std::string> mapQs = Utilities::split(str_mapQ, " ");

			seedAndExtend_return_local a;
			a.Score = Utilities::StrtoD(str_score);
			a.reverse = Utilities::StrtoB(str_reverse);
			if(mapQs.size() == 2)
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = Utilities::StrtoD(mapQs.at(1));
			}
			else if(mapQs.size() == 3)
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = Utilities::StrtoD(mapQs.at(1));
				a.mapQ_genomic_perPosition = mapQs.at(2);
				assert(a.mapQ_genomic_perPosition.length() == str_graph_aligned.length());
			}
			else
			{
				a.mapQ = Utilities::StrtoD(mapQs.at(0));
				a.mapQ_genomic = 2;
			}
			a.graph_aligned = str_graph_aligned;
			a.sequence_aligned = str_sequence_aligned;
			a.graph_aligned_levels = Utilities::StrtoI(Utilities::split(str_levels, " "));
			ret_alignment = a;

			oneRead r(str_readID, str_originalSequence, str_qualities);
			ret_originalRead = r;
		}
	};
	std::ifstream inputStream;
	inputStream.open(file.c_str());
	assert(inputStream.is_open());
	std::string line;
	assert(inputStream.good());

	while(inputStream.good())
	{
		std::getline(inputStream, line);
		
		got_lines++;
		
		Utilities::eraseNL(line);
		if(line.length())
		{

			std::string lineForInsertion;

			if(!(line.substr(0, std::string("Aligned unpaired read").length()) == "Aligned unpaired read"))
			{
				std::cerr << "! line.substr(0, std::string('Aligned unpaired read').length()) == 'Aligned unpaired read')" << "\n" << std::flush;
				std::cerr << "'" << line.substr(0, std::string("Aligned unpaired read").length()) << "'" << "\n" << std::flush;
				std::cerr << "got_lines: " <<  got_lines << "\n" << std::flush;
			}
			assert(line.substr(0, std::string("Aligned unpaired read").length()) == "Aligned unpaired read");

			bool fail;

			seedAndExtend_return_local a;

			oneRead r("", "", "");
			getAlignment(inputStream, fail, a, r);

			if(fail)
			{
				assert(! inputStream.good());
			}

			if(! fail)
			{
				ret_alignments.push_back(a);
				ret_alignments_originalReads.push_back(r);
			}
		}
	}
}


void estimateInsertSizeFromGraph(std::string FASTQs, std::string graphDir, std::vector<std::pair<double, double>>& inserSize_mean_sd_perFile, bool MiSeq250bp)
{
	inserSize_mean_sd_perFile.clear();
	
	int aligner_kMerSize = 25;
	int outerThreads = 8;
	int skipPairs_MOD = 1;
	bool useShort = true;

	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	std::cout << Utilities::timestamp() << "estimateInsertSizeFromGraph(..): Loading graph.\n" << std::flush;
	Graph* g = new Graph();
	g->readFromFile(graph);

	std::cout << Utilities::timestamp() << "estimateInsertSizeFromGraph(..): Create GraphAlignerUnique(s).\n" << std::flush;
	GraphAlignerUnique::GraphAlignerUnique* gA = new GraphAlignerUnique::GraphAlignerUnique(g, aligner_kMerSize);
	gA->setIterationsMainRandomizationLoop(4);
	gA->setThreads(1);

	omp_set_num_threads(1);

	double targetReadPairs = 2000;
	
	std::vector<std::string> FASTQ_files = Utilities::split(FASTQs, ",");
	for(unsigned int fI = 0; fI < FASTQ_files.size(); fI++)
	{

		std::string FASTQ = FASTQ_files.at(fI);

		std::cout << "File: " << FASTQ << "\n";
				
		std::cout << Utilities::timestamp() << "estimateInsertSizeFromGraph(..): Loading reads from " << FASTQ << ".\n" << std::flush;
		std::vector<oneReadPair> combinedPairs_for_alignment = getReadsFromFastQ(FASTQ);

		int MOD_for_filtering = 0;

		std::cout << "Read pairs before filtering: " << combinedPairs_for_alignment.size() << "\n";
		
		if(combinedPairs_for_alignment.size() > targetReadPairs)
		{
			MOD_for_filtering = combinedPairs_for_alignment.size() / targetReadPairs;
		}

		if(MOD_for_filtering != 0)
		{
			std::cout << "\tMOD: " << MOD_for_filtering << "\n";
			
			std::vector<oneReadPair> combinedPairs_for_alignment_F;
			for(unsigned int i = 0; i < combinedPairs_for_alignment.size(); i++)
			{
				if((i % MOD_for_filtering) == 0)
				{
					combinedPairs_for_alignment_F.push_back(combinedPairs_for_alignment.at(i));
				}
			}

			assert(combinedPairs_for_alignment_F.size() >= targetReadPairs);

			combinedPairs_for_alignment = combinedPairs_for_alignment_F;
		}
		
		std::cout << "Read pairs after filtering: " << combinedPairs_for_alignment.size() << "\n";

		double IS_total_size = 0;
		std::map<int, double> IS_combined_counts;
		std::set<int> IS_keys;

		for(unsigned int pI = 0; pI < combinedPairs_for_alignment.size(); pI++)
		{
			std::map<int, double> IS_p;
			std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_pair = gA->seedAndExtend_local_paired_or_short(combinedPairs_for_alignment.at(pI), true, useShort, 1, 1, true, IS_p, MiSeq250bp);

			for(std::map<int, double>::iterator ISit = IS_p.begin(); ISit != IS_p.end(); ISit++)
			{
				int d = ISit->first;
				double p = ISit->second;

				if(IS_combined_counts.count(d) == 0)
				{
					IS_combined_counts[d] = 0;
					IS_keys.insert(d);
				}

				IS_combined_counts.at(d) += p;
				IS_total_size += p;
			}

		}

		double cumulative_sum = 0;
		double weighted_median = 0;
		double weighted_20 = 0;
		double weighted_80 = 0;
		bool set_median = false;
		bool set_weighted_20 = false;
		bool set_weighted_80 = false;

		std::cout << "\n\nIS histogram over " << IS_total_size << " read pairs:\n";
		for(std::set<int>::iterator ISit = IS_keys.begin(); ISit != IS_keys.end(); ISit++)
		{
			int d = *ISit;
			double c = IS_combined_counts.at(d);
			std::cout << "\t" << d << ": " << c << "\n";
			cumulative_sum += IS_combined_counts.at(d);
			if((set_median == false) && (cumulative_sum >= (IS_total_size * 0.5)))
			{
				weighted_median = d;
				set_median = true;
			}
			if((set_weighted_20 == false) && (cumulative_sum >= (IS_total_size * 0.2)))
			{
				weighted_20 = d;
				set_weighted_20 = true;
			}
			if((set_weighted_80 == false) && (cumulative_sum >= (IS_total_size * 0.8)))
			{
				weighted_80 = d;
				set_weighted_80 = true;
			}
		}
		std::cout << "\n" << std::flush;

		std::cout << "Summary statistics:\n";
		std::cout << "\t" << "Median: " << weighted_median << "\n";
		std::cout << "\t" << "20%: " << weighted_20 << "\n";
		std::cout << "\t" << "80%: " << weighted_80 << "\n";

		std::cout << "\n" << std::flush;

		double f_mean_ret = weighted_median;
		weighted_20 = abs(weighted_median - weighted_20);
		weighted_80 = abs(weighted_median - weighted_80);
		double f_sd_ret = (weighted_20 > weighted_80) ? weighted_20 : weighted_80;
		
		std::cout << "\n\nAssumed values for " << FASTQ << ":\n";
		std::cout << "\t" << "Mean: " << f_mean_ret << "\n";
		std::cout << "\t" << "SD: " << f_sd_ret << "\n";
		std::cout << "\n" << std::flush;		
		
		inserSize_mean_sd_perFile.push_back(make_pair(f_mean_ret, f_sd_ret));
	}
	
	delete(gA);
	delete(g);
}

void alignLongUnpairedReadsToHLAGraph(std::string FASTQs, std::string graphDir, std::string referenceGenomeFile)
{
	int aligner_kMerSize = 15;
	int outerThreads = 8;
	int skipPairs_MOD = 1;
	bool useShort = true;

	unsigned int print_max_alignments = 1;

	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	std::map<std::string, std::string> referenceChromosomes = Utilities::readFASTA(referenceGenomeFile);
	assert(referenceChromosomes.count("ref"));

	std::cout << Utilities::timestamp() << "alignLongUnpairedReadsToHLAGraph(..): Loading graph.\n" << std::flush;

	std::cout << Utilities::timestamp() << "alignLongUnpairedReadsToHLAGraph(..): Create GraphAlignerUnique(s).\n" << std::flush;

	omp_set_num_threads(outerThreads);

	std::vector<std::map<int, int>> printedAlignments_perRead_perThread;
	std::vector<std::map<int, int>> printedAlignments_perRead_combinedLL_perThread;
	std::vector<std::map<int, int>> printedAlignments_perRead_combinedGenomicLL_perThread;

	std::vector<Graph*> graphs;
	std::vector<GraphAlignerUnique::GraphAlignerUnique*> graphAligners;
	graphs.resize(outerThreads);
	graphAligners.resize(outerThreads);

	#pragma omp parallel for
	for(int tI = 0; tI < outerThreads; tI++)
	{
		// std::cout << "Thread " << tI << "\n" << std::flush;

		Graph* g = new Graph();
		g->readFromFile(graph);

		graphAligners.at(tI) = new GraphAlignerUnique::GraphAlignerUnique(g, aligner_kMerSize);
		graphAligners.at(tI)->setIterationsMainRandomizationLoop(4);
		graphAligners.at(tI)->setThreads(1);

		graphs.at(tI) = g;
	}

	printedAlignments_perRead_perThread.resize(outerThreads);
	printedAlignments_perRead_combinedLL_perThread.resize(outerThreads);
	printedAlignments_perRead_combinedGenomicLL_perThread.resize(outerThreads);

	auto alignReads = [&](std::vector<oneRead>& reads, std::vector< std::vector< std::vector<seedAndExtend_return_local> > >& alignments_perThread, std::vector< std::vector<int> >& alignments_readPairI_perThread) -> void
	{
		unsigned int readI = 0;
		unsigned int readIMax = reads.size();
		#pragma omp parallel for schedule(dynamic)
		for(readI = 0; readI < readIMax; readI++)
		{
			int tI = omp_get_thread_num();
			assert(omp_get_num_threads() == outerThreads);
			assert((tI >= 0) && (tI < outerThreads));
				
			assert((readI >= 0) && (readI < reads.size()));
			oneRead R= reads.at(readI);

			assert((tI >= 0) && (tI < graphAligners.size()));

			std::map<int, double> _IS_ignore;
						
			std::vector<seedAndExtend_return_local> alignments = graphAligners.at(tI)->seedAndExtend_longlocal_allAlignments(R);


			if(alignments.size() > 1)
			{
				assert(alignments.at(0).mapQ_genomic >= alignments.at(1).mapQ_genomic);
				assert(alignments.at(0).mapQ >= alignments.at(1).mapQ);
			}

			std::vector<seedAndExtend_return_local> alignments_forPrint;
			unsigned int max_print_index = (alignments.size() > print_max_alignments) ? print_max_alignments : alignments.size();
			double accumulated_LL = 0;
			double accumulated_genomic_LL = 0;   
			assert(alignments.size() > 0);
			assert(max_print_index <= alignments.size());
			for(unsigned int i = 0; i < max_print_index; i++)
			{
				alignments_forPrint.push_back(alignments.at(i));
				accumulated_LL += alignments.at(i).mapQ;
				accumulated_genomic_LL += alignments.at(i).mapQ_genomic;
			}

			if(printedAlignments_perRead_perThread.at(tI).count(max_print_index) == 0)
			{
				printedAlignments_perRead_perThread.at(tI)[max_print_index] = 0;
			}
			printedAlignments_perRead_perThread.at(tI).at(max_print_index)++;

			int accumulated_LL_asInt = int(accumulated_LL * 10.0);
			if(printedAlignments_perRead_combinedLL_perThread.at(tI).count(accumulated_LL_asInt) == 0)
			{
				printedAlignments_perRead_combinedLL_perThread.at(tI)[accumulated_LL_asInt] = 0;
			}
			printedAlignments_perRead_combinedLL_perThread.at(tI).at(accumulated_LL_asInt)++;

			int accumulated_genomic_LL_asInt = int(accumulated_genomic_LL * 10.0);
			if(printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).count(accumulated_genomic_LL_asInt) == 0)
			{
				printedAlignments_perRead_combinedGenomicLL_perThread.at(tI)[accumulated_genomic_LL_asInt] = 0;
			}
			printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).at(accumulated_genomic_LL_asInt)++;

			alignments_perThread.at(tI).push_back(alignments_forPrint);
			alignments_readPairI_perThread.at(tI).push_back(readI);

			if(tI == 0)
			{
				std::cout  << Utilities::timestamp() << "\t\t" << "Thread " << tI << ": align read " << readI << "\n" << std::flush;
			}			
		}
	};


	int totalPrintedBases = 0;
	int totalPrintedBases_recovered = 0;

	// auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<std::pair<std::string, std::string> > alignments_readIDs) -> void {
	auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::vector<seedAndExtend_return_local>>& alignments, std::vector<oneRead>& originalReads) -> void {

		double printingThreshold = 0.01;

		std::ofstream outputStream;
		outputStream.open(outputFilename.c_str());
		assert(outputStream.is_open());

		auto printOneRead = [&](seedAndExtend_return_local alignment, oneRead originalRead, int aI) -> void {
			outputStream << "\t" << "Read " << originalRead.name << "A" << aI << "\n";
			outputStream << "\t\t" << alignment.Score << "\n";
			outputStream << "\t\t" << alignment.reverse << "\n";
			outputStream << "\t\t" << alignment.mapQ << " " << alignment.mapQ_genomic << " " << alignment.mapQ_genomic_perPosition << "\n";
			outputStream << "\t\t" << alignment.graph_aligned << "\n";
			outputStream << "\t\t" << alignment.sequence_aligned << "\n";
			outputStream << "\t\t" << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), " ") << "\n";
			outputStream << "\t\t" << originalRead.sequence << "\n";
			outputStream << "\t\t" << originalRead.quality << "\n";

			totalPrintedBases += alignment.graph_aligned.length();
			if((alignment.mapQ_genomic < 0.9) && (alignment.mapQ_genomic_perPosition.length()))
			{
				assert(alignment.mapQ_genomic_perPosition.length() == alignment.graph_aligned.length());
				for(unsigned int i = 0; i < alignment.mapQ_genomic_perPosition.length(); i++)
				{
					char c = alignment.mapQ_genomic_perPosition.at(i);
					if(Utilities::PhredToPCorrect(c) > 0.9)
					{
						totalPrintedBases_recovered++;
					}
				}
			}

		};

		for(unsigned int readI = 0; readI < alignments.size(); readI++)
		{
			bool printAtLeastOneAlignment = false;
			for(unsigned int aI = 0; aI < alignments.at(readI).size(); aI++)
			{
				if(alignments.at(readI).at(aI).mapQ_genomic > printingThreshold)
				{
					printAtLeastOneAlignment = true;
					break;
				}
			}

			if(printAtLeastOneAlignment)
			{
				outputStream << "Aligned unpaired read " << readI << "\n";
				for(unsigned int aI = 0; aI < alignments.at(readI).size(); aI++)
				{
					if(alignments.at(readI).at(aI).mapQ_genomic > printingThreshold)
					{
						printOneRead(alignments.at(readI).at(aI), originalReads.at(readI), aI);
					}
				}
			}
		}
		outputStream.close();
	};


	std::vector<std::string> FASTQ_files = Utilities::split(FASTQs, ",");
	for(unsigned int fI = 0; fI < FASTQ_files.size(); fI++)
	{
		for(int tI = 0; tI < outerThreads; tI++)
		{
			printedAlignments_perRead_perThread.at(tI).clear();
			printedAlignments_perRead_combinedLL_perThread.at(tI).clear();
			printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).clear();
		}

		std::string FASTQ = FASTQ_files.at(fI);

		std::cout << Utilities::timestamp() << "alignLongUnpairedReadsToHLAGraph(..): Loading reads from " << FASTQ << ".\n" << std::flush;
		std::vector<oneRead> combinedReads_for_alignment = getUnpairedReadsFromFastQ(FASTQ);

		if(skipPairs_MOD != 1)
		{
			std::vector<oneRead> reads_after_filtering;
			for(unsigned int readI = 0; readI < combinedReads_for_alignment.size(); readI++)
			{
				if((readI % skipPairs_MOD) == 0)
				{
					reads_after_filtering.push_back(combinedReads_for_alignment.at(readI));
				}
			}
			combinedReads_for_alignment = reads_after_filtering;
		}

		std::cout << "\t" << "Now align " << combinedReads_for_alignment.size() << " unpaired reads." << "\n" << std::flush;

		std::vector< std::vector< std::vector<seedAndExtend_return_local> > > alignments_perThread;
		std::vector< std::vector<int> > alignments_readI_perThread;

		alignments_perThread.resize(outerThreads);
		alignments_readI_perThread.resize(outerThreads);

		alignReads(combinedReads_for_alignment, alignments_perThread, alignments_readI_perThread);

		// merge

		std::cout  << Utilities::timestamp() << "\t\t" << "All reads aligned - merge.\n" << std::flush;

		std::vector< std::vector< seedAndExtend_return_local > > alignments;
		std::vector< int > alignments_readI;
		for(unsigned int tI = 0; (int)tI < outerThreads; tI++)
		{
			alignments.insert(alignments.end(), alignments_perThread.at(tI).begin(), alignments_perThread.at(tI).end());
			alignments_readI.insert(alignments_readI.end(), alignments_readI_perThread.at(tI).begin(), alignments_readI_perThread.at(tI).end());
		}

		// Produce normal output file
		std::string alignments_output_file = FASTQ + ".aligned";
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Produce alignments file " << alignments_output_file << ".\n" << std::flush;

		// std::vector<std::pair<std::string, std::string> > alignments_readIDs;
		std::vector<oneRead> combinedReads_for_alignment_inAlignmentOrder;
		for(unsigned int readII = 0; readII < alignments_readI.size(); readII++)
		{
			int readI = alignments_readI.at(readII);
			// std::string ID1 = combinedPairs_for_alignment.at(readPairI).reads.first.name;
			// std::string ID2 = combinedPairs_for_alignment.at(readPairI).reads.second.name;
			// alignments_readIDs.push_back(make_pair(ID1, ID2));
			combinedReads_for_alignment_inAlignmentOrder.push_back(combinedReads_for_alignment.at(readI));
		}

		assert(combinedReads_for_alignment_inAlignmentOrder.size() == alignments.size());

		printAlignmentsToFile(alignments_output_file, alignments, combinedReads_for_alignment_inAlignmentOrder);

		std::map<int, int> printedAlignments_perRead;
		std::map<int, int> printedAlignments_perRead_combinedLL;
		std::map<int, int> printedAlignments_perRead_combinedGenomicLL;

		for(int tI = 0; tI < outerThreads; tI++)
		{
			for(std::map<int, int>::iterator printedAlignmentsIt = printedAlignments_perRead_perThread.at(tI).begin(); printedAlignmentsIt != printedAlignments_perRead_perThread.at(tI).end(); printedAlignmentsIt++)
			{
				if(printedAlignments_perRead.count(printedAlignmentsIt->first) == 0)
				{
					printedAlignments_perRead[printedAlignmentsIt->first] = 0;
				}
				printedAlignments_perRead.at(printedAlignmentsIt->first) += printedAlignmentsIt->second;
			}

			for(std::map<int, int>::iterator alignmentPrintedProbabilityIt = printedAlignments_perRead_combinedLL_perThread.at(tI).begin(); alignmentPrintedProbabilityIt != printedAlignments_perRead_combinedLL_perThread.at(tI).end(); alignmentPrintedProbabilityIt++)
			{
				if(printedAlignments_perRead_combinedLL.count(alignmentPrintedProbabilityIt->first) == 0)
				{
					printedAlignments_perRead_combinedLL[alignmentPrintedProbabilityIt->first] = 0;
				}
				printedAlignments_perRead_combinedLL.at(alignmentPrintedProbabilityIt->first) += alignmentPrintedProbabilityIt->second;
			}

			for(std::map<int, int>::iterator alignmentPrintedGenomicProbabilityIt = printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).begin(); alignmentPrintedGenomicProbabilityIt != printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).end(); alignmentPrintedGenomicProbabilityIt++)
			{
				if(printedAlignments_perRead_combinedGenomicLL.count(alignmentPrintedGenomicProbabilityIt->first) == 0)
				{
					printedAlignments_perRead_combinedGenomicLL[alignmentPrintedGenomicProbabilityIt->first] = 0;
				}
				printedAlignments_perRead_combinedGenomicLL.at(alignmentPrintedGenomicProbabilityIt->first) += alignmentPrintedGenomicProbabilityIt->second;
			}
		}

		std::cout  << Utilities::timestamp() << "\t\t\t" << "Print multi-alignment statistics\n" << std::flush;
		std::cout << "\t\t\t\t" << "Printed alingments per read:\n";
		for(std::map<int, int>::iterator printedAlignmentsIt = printedAlignments_perRead.begin(); printedAlignmentsIt != printedAlignments_perRead.end(); printedAlignmentsIt++)
		{
			std::cout << "\t\t\t\t\t" << printedAlignmentsIt->first << " => " << (int)printedAlignmentsIt->second << "\n";
		}

		std::cout << "\t\t\t\t" << "Accumulated mapping quality per read (x 10):\n";
		for(std::map<int, int>::iterator alignmentPrintedProbabilityIt = printedAlignments_perRead_combinedLL.begin(); alignmentPrintedProbabilityIt != printedAlignments_perRead_combinedLL.end(); alignmentPrintedProbabilityIt++)
		{
			std::cout << "\t\t\t\t\t" << alignmentPrintedProbabilityIt->first << " => " << (int)alignmentPrintedProbabilityIt->second << "\n";
		}

		std::cout << "\t\t\t\t" << "Accumulated genomic mapping quality per read (x 10):\n";
		for(std::map<int, int>::iterator alignmentPrintedGenomicProbabilityIt = printedAlignments_perRead_combinedGenomicLL.begin(); alignmentPrintedGenomicProbabilityIt != printedAlignments_perRead_combinedGenomicLL.end(); alignmentPrintedGenomicProbabilityIt++)
		{
			std::cout << "\t\t\t\t\t" << alignmentPrintedGenomicProbabilityIt->first << " => " << (int)alignmentPrintedGenomicProbabilityIt->second << "\n";
		}
	}

	std::cout << "\ntotalPrintedBases: " << totalPrintedBases << "\n";
	std::cout << "totalPrintedBases_recovered: " << totalPrintedBases_recovered << "\n\n" << std::flush;

	std::cout  << Utilities::timestamp() << "\t\t\t" << "All alignments done, free memory.\n" << std::flush;
	for(int tI = 0; tI < outerThreads; tI++)
	{
		delete(graphAligners.at(tI));
		delete(graphs.at(tI));
	}
}

void alignShortReadsToHLAGraph_multipleAlignments(std::string FASTQs, std::string graphDir, std::string referenceGenomeFile, std::vector<std::pair<double, double>> inserSize_mean_sd_perFile, bool debug, bool MiSeq250bp)
{
	int aligner_kMerSize = 25;
	int outerThreads = (debug ? 1: 40);
	//int outerThreads = 1;
	int skipPairs_MOD = 1;
	bool useShort = true;

	unsigned int print_max_alignments = 1;

	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	std::map<std::string, std::string> referenceChromosomes = Utilities::readFASTA(referenceGenomeFile);
	assert(referenceChromosomes.count("ref"));

	std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph_multipleAlignments(..): Loading graph.\n" << std::flush;

	std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph_multipleAlignments(..): Create GraphAlignerUnique(s).\n" << std::flush;

	omp_set_num_threads(outerThreads);

	std::vector<std::map<int, int>> printedAlignments_perRead_perThread;
	std::vector<std::map<int, int>> printedAlignments_perRead_combinedLL_perThread;
	std::vector<std::map<int, int>> printedAlignments_perRead_combinedGenomicLL_perThread;

	std::vector<Graph*> graphs;
	std::vector<GraphAlignerUnique::GraphAlignerUnique*> graphAligners;
	graphs.resize(outerThreads);
	graphAligners.resize(outerThreads);

	#pragma omp parallel for
	for(int tI = 0; tI < outerThreads; tI++)
	{
		// std::cout << "Thread " << tI << "\n" << std::flush;

		Graph* g = new Graph();
		g->readFromFile(graph);

		graphAligners.at(tI) = new GraphAlignerUnique::GraphAlignerUnique(g, aligner_kMerSize);
		graphAligners.at(tI)->setIterationsMainRandomizationLoop(4);
		graphAligners.at(tI)->setThreads(1);

		graphs.at(tI) = g;
	}

	printedAlignments_perRead_perThread.resize(outerThreads);
	printedAlignments_perRead_combinedLL_perThread.resize(outerThreads);
	printedAlignments_perRead_combinedGenomicLL_perThread.resize(outerThreads);

	auto alignReadPairs = [&](std::vector<oneReadPair>& readPairs, std::vector< std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> > >& alignments_perThread, std::vector< std::vector<int> >& alignments_readPairI_perThread, bool usePairing, double insertSize_mean, double insertSize_sd) -> void
	{
		unsigned int pairI = 0;
		unsigned int pairMax = readPairs.size();
		#pragma omp parallel for schedule(dynamic)
		for(pairI = 0; pairI < pairMax; pairI++)
		{
			int tI = omp_get_thread_num();
			assert(omp_get_num_threads() == outerThreads);
			assert((tI >= 0) && (tI < outerThreads));

			assert((pairI >= 0) && (pairI < readPairs.size()));
			oneReadPair rP = readPairs.at(pairI);

			assert((tI >= 0) && (tI < graphAligners.size()));

			std::map<int, double> _IS_ignore;
			std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> alignment_pairs = graphAligners.at(tI)->seedAndExtend_short_allAlignments(rP, insertSize_mean, insertSize_sd, true, debug, MiSeq250bp);

			if(alignment_pairs.size() > 1)
			{
				assert(alignment_pairs.at(0).first.mapQ_genomic >= alignment_pairs.at(1).first.mapQ_genomic);
				assert(alignment_pairs.at(0).first.mapQ >= alignment_pairs.at(1).first.mapQ);
			}
			
			std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> alignment_pairs_forPrint;
			unsigned int max_print_index = (alignment_pairs.size() > print_max_alignments) ? print_max_alignments : alignment_pairs.size();
			double accumulated_LL = 0;
			double accumulated_genomic_LL = 0;
			for(unsigned int i = 0; i < max_print_index; i++)
			{
				alignment_pairs_forPrint.push_back(alignment_pairs.at(i));
				accumulated_LL += alignment_pairs.at(i).first.mapQ;
				accumulated_genomic_LL += alignment_pairs.at(i).first.mapQ_genomic;
			}

			if(printedAlignments_perRead_perThread.at(tI).count(max_print_index) == 0)
			{
				printedAlignments_perRead_perThread.at(tI)[max_print_index] = 0;
			}
			printedAlignments_perRead_perThread.at(tI).at(max_print_index)++;

			int accumulated_LL_asInt = int(accumulated_LL * 10.0);
			if(printedAlignments_perRead_combinedLL_perThread.at(tI).count(accumulated_LL_asInt) == 0)
			{
				printedAlignments_perRead_combinedLL_perThread.at(tI)[accumulated_LL_asInt] = 0;
			}
			printedAlignments_perRead_combinedLL_perThread.at(tI).at(accumulated_LL_asInt)++;
			
			int accumulated_genomic_LL_asInt = int(accumulated_genomic_LL * 10.0);
			if(printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).count(accumulated_genomic_LL_asInt) == 0)
			{
				printedAlignments_perRead_combinedGenomicLL_perThread.at(tI)[accumulated_genomic_LL_asInt] = 0;
			}
			printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).at(accumulated_genomic_LL_asInt)++;


			alignments_perThread.at(tI).push_back(alignment_pairs_forPrint);
			alignments_readPairI_perThread.at(tI).push_back(pairI);

			if(tI == 0)
			{
				std::cout  << Utilities::timestamp() << "\t\t" << "Thread " << tI << ": align pair " << pairI << "\n" << std::flush;
			}
		}
	};


	int totalPrintedBases = 0;
	int totalPrintedBases_recovered = 0;

	// auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<std::pair<std::string, std::string> > alignments_readIDs) -> void {
	auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > >& alignments, std::vector<oneReadPair>& originalReads, double insertSize_mean, double insertSize_sd) -> void {
	
		double printingThreshold = 0.01;
		
		std::ofstream outputStream;
		outputStream.open(outputFilename.c_str());
		assert(outputStream.is_open());

		outputStream << "IS " << insertSize_mean << " " << insertSize_sd << "\n";
		
		auto printOneRead = [&](seedAndExtend_return_local alignment, oneRead originalRead, int aI) -> void {
			outputStream << "\t" << "Read " << originalRead.name << "A" << aI << "\n";
			outputStream << "\t\t" << alignment.Score << "\n";
			outputStream << "\t\t" << alignment.reverse << "\n";
			outputStream << "\t\t" << alignment.mapQ << " " << alignment.mapQ_genomic << " " << alignment.mapQ_genomic_perPosition << "\n";
			outputStream << "\t\t" << alignment.graph_aligned << "\n";
			outputStream << "\t\t" << alignment.sequence_aligned << "\n";
			outputStream << "\t\t" << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), " ") << "\n";
			outputStream << "\t\t" << originalRead.sequence << "\n";
			outputStream << "\t\t" << originalRead.quality << "\n";

			totalPrintedBases += alignment.graph_aligned.length();
			if((alignment.mapQ_genomic < 0.9) && (alignment.mapQ_genomic_perPosition.length()))
			{
				assert(alignment.mapQ_genomic_perPosition.length() == alignment.graph_aligned.length());
				for(unsigned int i = 0; i < alignment.mapQ_genomic_perPosition.length(); i++)
				{
					char c = alignment.mapQ_genomic_perPosition.at(i);
					if(Utilities::PhredToPCorrect(c) > 0.9)
					{
						totalPrintedBases_recovered++;
					}
				}
			}

		};

		for(unsigned int pairI = 0; pairI < alignments.size(); pairI++)
		{
			bool printAtLeastOneAlignment = false;
			for(unsigned int aI = 0; aI < alignments.at(pairI).size(); aI++)
			{
				if(alignments.at(pairI).at(aI).first.mapQ_genomic > printingThreshold)
				{
					printAtLeastOneAlignment = true;
					break;
				}					
			}			
			
			if(printAtLeastOneAlignment)
			{
				outputStream << "Aligned pair " << pairI << "\n";   
				for(unsigned int aI = 0; aI < alignments.at(pairI).size(); aI++)
				{
					if(alignments.at(pairI).at(aI).first.mapQ_genomic > printingThreshold)
					{
						printOneRead(alignments.at(pairI).at(aI).first, originalReads.at(pairI).reads.first, aI);
						printOneRead(alignments.at(pairI).at(aI).second, originalReads.at(pairI).reads.second, aI);
					}
				}
			}
		}
		outputStream.close();
	};       


	std::vector<std::string> FASTQ_files = Utilities::split(FASTQs, ",");
	assert(inserSize_mean_sd_perFile.size() == FASTQ_files.size());
	
	for(unsigned int fI = 0; fI < FASTQ_files.size(); fI++)
	{
		for(int tI = 0; tI < outerThreads; tI++)
		{
			printedAlignments_perRead_perThread.at(tI).clear();
			printedAlignments_perRead_combinedLL_perThread.at(tI).clear();
			printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).clear();
		}

		std::string FASTQ = FASTQ_files.at(fI);

		std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph_multipleAlignments(..): Loading reads from " << FASTQ << ".\n" << std::flush;
		std::vector<oneReadPair> combinedPairs_for_alignment = getReadsFromFastQ(FASTQ);

		if(skipPairs_MOD != 1)
		{
			std::vector<oneReadPair> pairs_after_filtering;
			for(unsigned int pairI = 0; pairI < combinedPairs_for_alignment.size(); pairI++)
			{
				if((pairI % skipPairs_MOD) == 0)
				{
					pairs_after_filtering.push_back(combinedPairs_for_alignment.at(pairI));
				}
			}
			combinedPairs_for_alignment = pairs_after_filtering;
		}

		std::cout << "\t" << "Now align " << combinedPairs_for_alignment.size() << " read pairs." << "\n" << std::flush;

		std::vector< std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > > > withPairing_alignments_perThread;
		std::vector< std::vector<int> > withPairing_alignments_readPairI_perThread;

		withPairing_alignments_perThread.resize(outerThreads);
		withPairing_alignments_readPairI_perThread.resize(outerThreads);

		alignReadPairs(combinedPairs_for_alignment, withPairing_alignments_perThread, withPairing_alignments_readPairI_perThread, true, inserSize_mean_sd_perFile.at(fI).first, inserSize_mean_sd_perFile.at(fI).second);

		// merge

		std::cout  << Utilities::timestamp() << "\t\t" << "All pairs aligned - merge.\n" << std::flush;

		std::vector< std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > > withPairing_alignments;
		std::vector< int > withPairing_alignments_readPairI;
		for(unsigned int tI = 0; (int)tI < outerThreads; tI++)
		{
			withPairing_alignments.insert(withPairing_alignments.end(), withPairing_alignments_perThread.at(tI).begin(), withPairing_alignments_perThread.at(tI).end());
			withPairing_alignments_readPairI.insert(withPairing_alignments_readPairI.end(), withPairing_alignments_readPairI_perThread.at(tI).begin(), withPairing_alignments_readPairI_perThread.at(tI).end());
		}

		// Produce normal output file
		std::string alignments_output_file = FASTQ + ".aligned";
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Produce alignments file " << alignments_output_file << ".\n" << std::flush;

		// std::vector<std::pair<std::string, std::string> > alignments_readIDs;
		std::vector<oneReadPair> combinedPairs_for_alignment_inAlignmentOrder;
		for(unsigned int pairII = 0; pairII < withPairing_alignments_readPairI.size(); pairII++)
		{
			int readPairI = withPairing_alignments_readPairI.at(pairII);
			// std::string ID1 = combinedPairs_for_alignment.at(readPairI).reads.first.name;
			// std::string ID2 = combinedPairs_for_alignment.at(readPairI).reads.second.name;
			// alignments_readIDs.push_back(make_pair(ID1, ID2));
			combinedPairs_for_alignment_inAlignmentOrder.push_back(combinedPairs_for_alignment.at(readPairI));
		}

		assert(combinedPairs_for_alignment_inAlignmentOrder.size() == withPairing_alignments.size());


		printAlignmentsToFile(alignments_output_file, withPairing_alignments, combinedPairs_for_alignment_inAlignmentOrder, inserSize_mean_sd_perFile.at(fI).first, inserSize_mean_sd_perFile.at(fI).second);

		

		std::map<int, int> printedAlignments_perRead;
		std::map<int, int> printedAlignments_perRead_combinedLL;
		std::map<int, int> printedAlignments_perRead_combinedGenomicLL;

	
		for(int tI = 0; tI < outerThreads; tI++)
		{
			for(std::map<int, int>::iterator printedAlignmentsIt = printedAlignments_perRead_perThread.at(tI).begin(); printedAlignmentsIt != printedAlignments_perRead_perThread.at(tI).end(); printedAlignmentsIt++)
			{
				if(printedAlignments_perRead.count(printedAlignmentsIt->first) == 0)
				{
					printedAlignments_perRead[printedAlignmentsIt->first] = 0;
				}
				printedAlignments_perRead.at(printedAlignmentsIt->first) += printedAlignmentsIt->second;
			}
			
			for(std::map<int, int>::iterator alignmentPrintedProbabilityIt = printedAlignments_perRead_combinedLL_perThread.at(tI).begin(); alignmentPrintedProbabilityIt != printedAlignments_perRead_combinedLL_perThread.at(tI).end(); alignmentPrintedProbabilityIt++)
			{
				if(printedAlignments_perRead_combinedLL.count(alignmentPrintedProbabilityIt->first) == 0)
				{
					printedAlignments_perRead_combinedLL[alignmentPrintedProbabilityIt->first] = 0;
				}
				printedAlignments_perRead_combinedLL.at(alignmentPrintedProbabilityIt->first) += alignmentPrintedProbabilityIt->second;
			}			

			for(std::map<int, int>::iterator alignmentPrintedGenomicProbabilityIt = printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).begin(); alignmentPrintedGenomicProbabilityIt != printedAlignments_perRead_combinedGenomicLL_perThread.at(tI).end(); alignmentPrintedGenomicProbabilityIt++)
			{
				if(printedAlignments_perRead_combinedGenomicLL.count(alignmentPrintedGenomicProbabilityIt->first) == 0)
				{
					printedAlignments_perRead_combinedGenomicLL[alignmentPrintedGenomicProbabilityIt->first] = 0;
				}
				printedAlignments_perRead_combinedGenomicLL.at(alignmentPrintedGenomicProbabilityIt->first) += alignmentPrintedGenomicProbabilityIt->second;
			}	
		}

		std::cout  << Utilities::timestamp() << "\t\t\t" << "Print multi-alignment statistics\n" << std::flush;
		std::cout << "\t\t\t\t" << "Printed alingments per read:\n";
		for(std::map<int, int>::iterator printedAlignmentsIt = printedAlignments_perRead.begin(); printedAlignmentsIt != printedAlignments_perRead.end(); printedAlignmentsIt++)
		{
			std::cout << "\t\t\t\t\t" << printedAlignmentsIt->first << " => " << (int)printedAlignmentsIt->second << "\n";
		}		

		std::cout << "\t\t\t\t" << "Accumulated mapping quality per read (x 10):\n";
		for(std::map<int, int>::iterator alignmentPrintedProbabilityIt = printedAlignments_perRead_combinedLL.begin(); alignmentPrintedProbabilityIt != printedAlignments_perRead_combinedLL.end(); alignmentPrintedProbabilityIt++)
		{
			std::cout << "\t\t\t\t\t" << alignmentPrintedProbabilityIt->first << " => " << (int)alignmentPrintedProbabilityIt->second << "\n";
		}

		std::cout << "\t\t\t\t" << "Accumulated genomic mapping quality per read (x 10):\n";
		for(std::map<int, int>::iterator alignmentPrintedGenomicProbabilityIt = printedAlignments_perRead_combinedGenomicLL.begin(); alignmentPrintedGenomicProbabilityIt != printedAlignments_perRead_combinedGenomicLL.end(); alignmentPrintedGenomicProbabilityIt++)
		{
			std::cout << "\t\t\t\t\t" << alignmentPrintedGenomicProbabilityIt->first << " => " << (int)alignmentPrintedGenomicProbabilityIt->second << "\n";
		}		

		// Produce SAM -- deactivated
		continue;

		std::string SAM_output_file = FASTQ + ".sam";
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Produce SAM " << SAM_output_file << ".\n" << std::flush;


		// Get graph loci and reference positions
		std::ofstream SAM_output_stream;
		SAM_output_stream.open(SAM_output_file.c_str());
		assert(SAM_output_stream.is_open());

		std::vector<int> uncompressed_graph_referencePositions;
		std::vector<std::string> graphLoci = readGraphLoci(graphDir);

		int lastReferencePosition = -1;
		for(unsigned int i = 0; i < graphLoci.size(); i++)
		{
			std::string locusID = graphLoci.at(i);
			std::vector<std::string> locusParts = Utilities::split(locusID, "_");
			if(locusParts.size() != 3)
			{
				throw std::runtime_error("alignShortReadsToHLAGraph(..): Cannot decompose locus ID " +locusID);
			}
			int thisLocus_refPos = Utilities::StrtoI(locusParts.at(2));
			if(thisLocus_refPos != -1)
			{
				thisLocus_refPos = thisLocus_refPos + 1;
			}
			if((i == 0) || (lastReferencePosition != thisLocus_refPos))
			{
				uncompressed_graph_referencePositions.push_back(thisLocus_refPos);
				lastReferencePosition = thisLocus_refPos;
			}
			else
			{
				uncompressed_graph_referencePositions.push_back(-1);
			}
		}

		// alignedShortReads2SAM(SAM_output_stream, uncompressed_graph_referencePositions, referenceChromosomes.at("ref"), withPairing_alignments, combinedPairs_for_alignment_inAlignmentOrder);

		std::cout  << Utilities::timestamp() << "\t\t\t" << "Done. Output in " << SAM_output_file << ".\n" << std::flush;
	}

	std::cout << "\ntotalPrintedBases: " << totalPrintedBases << "\n";
	std::cout << "totalPrintedBases_recovered: " << totalPrintedBases_recovered << "\n\n" << std::flush;

	std::cout  << Utilities::timestamp() << "\t\t\t" << "All alignments done, free memory.\n" << std::flush;
	for(int tI = 0; tI < outerThreads; tI++)
	{
		delete(graphAligners.at(tI));
		delete(graphs.at(tI));
	}
}

void alignShortReadsToHLAGraph(std::string FASTQs, std::string graphDir, std::string referenceGenomeFile, std::vector<std::pair<double, double>> inserSize_mean_sd_perFile)
{
	int aligner_kMerSize = 25;
	int outerThreads = 8;
	int skipPairs_MOD = 1;
	bool useShort = true;

	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	std::map<std::string, std::string> referenceChromosomes = Utilities::readFASTA(referenceGenomeFile);
	assert(referenceChromosomes.count("ref"));

	std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph(..): Loading graph.\n" << std::flush;

	std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph(..): Create GraphAlignerUnique(s).\n" << std::flush;

	omp_set_num_threads(outerThreads);

	std::vector<Graph*> graphs;
	std::vector<GraphAlignerUnique::GraphAlignerUnique*> graphAligners;
	graphs.resize(outerThreads);
	graphAligners.resize(outerThreads);

	#pragma omp parallel for
	for(int tI = 0; tI < outerThreads; tI++)
	{
		// std::cout << "Thread " << tI << "\n" << std::flush;

		Graph* g = new Graph();
		g->readFromFile(graph);

		graphAligners.at(tI) = new GraphAlignerUnique::GraphAlignerUnique(g, aligner_kMerSize);
		graphAligners.at(tI)->setIterationsMainRandomizationLoop(4);
		graphAligners.at(tI)->setThreads(1);

		graphs.at(tI) = g;
	}


	auto alignReadPairs = [&](std::vector<oneReadPair>& readPairs, std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> >& alignments_perThread, std::vector< std::vector<int> >& alignments_readPairI_perThread, bool usePairing, double insertSize_mean, double insertSize_sd) -> void
	{
		unsigned int pairI = 0;
		unsigned int pairMax = readPairs.size();
		#pragma omp parallel for schedule(dynamic)
		for(pairI = 0; pairI < pairMax; pairI++)
		{
			int tI = omp_get_thread_num();
			assert(omp_get_num_threads() == outerThreads);
			assert((tI >= 0) && (tI < outerThreads));

			assert((pairI >= 0) && (pairI < readPairs.size()));
			oneReadPair rP = readPairs.at(pairI);

			assert((tI >= 0) && (tI < graphAligners.size()));

			std::map<int, double> _IS_ignore;
			std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignment_pair = graphAligners.at(tI)->seedAndExtend_local_paired_or_short(rP, usePairing, useShort, insertSize_mean, insertSize_sd, false, _IS_ignore);

			alignments_perThread.at(tI).push_back(alignment_pair);
			alignments_readPairI_perThread.at(tI).push_back(pairI);

			if(tI == 0)
			{
				std::cout  << Utilities::timestamp() << "\t\t" << "Thread " << tI << ": align pair " << pairI << "\n" << std::flush;
			}
		}
	};



	// auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<std::pair<std::string, std::string> > alignments_readIDs) -> void {
	auto printAlignmentsToFile = [&](std::string outputFilename, std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> >& alignments, std::vector<oneReadPair>& originalReads, double insertSize_mean, double insertSize_sd) -> void {
		std::ofstream outputStream;
		outputStream.open(outputFilename.c_str());
		assert(outputStream.is_open());

		outputStream << "IS " << insertSize_mean << " " << insertSize_sd << "\n";

		auto printOneRead = [&](seedAndExtend_return_local alignment, oneRead originalRead) -> void {
			outputStream << "\t" << "Read " << originalRead.name << "\n";
			outputStream << "\t\t" << alignment.Score << "\n";
			outputStream << "\t\t" << alignment.reverse << "\n";
			outputStream << "\t\t" << alignment.mapQ << " " << alignment.mapQ_genomic << "\n";
			outputStream << "\t\t" << alignment.graph_aligned << "\n";
			outputStream << "\t\t" << alignment.sequence_aligned << "\n";
			outputStream << "\t\t" << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), " ") << "\n";
			outputStream << "\t\t" << originalRead.sequence << "\n";
			outputStream << "\t\t" << originalRead.quality << "\n";
		};

		for(unsigned int pairI = 0; pairI < alignments.size(); pairI++)
		{
			outputStream << "Aligned pair " << pairI << "\n";
			printOneRead(alignments.at(pairI).first, originalReads.at(pairI).reads.first);
			printOneRead(alignments.at(pairI).second, originalReads.at(pairI).reads.second);
		}
		outputStream.close();
	};


	std::vector<std::string> FASTQ_files = Utilities::split(FASTQs, ",");
	assert(inserSize_mean_sd_perFile.size() == FASTQ_files.size());

	for(unsigned int fI = 0; fI < FASTQ_files.size(); fI++)
	{
		std::string FASTQ = FASTQ_files.at(fI);

		std::cout << Utilities::timestamp() << "alignShortReadsToHLAGraph(..): Loading reads from " << FASTQ << ".\n" << std::flush;
		std::vector<oneReadPair> combinedPairs_for_alignment = getReadsFromFastQ(FASTQ);

		if(skipPairs_MOD != 1)
		{
			std::vector<oneReadPair> pairs_after_filtering;
			for(unsigned int pairI = 0; pairI < combinedPairs_for_alignment.size(); pairI++)
			{
				if((pairI % skipPairs_MOD) == 0)
				{
					pairs_after_filtering.push_back(combinedPairs_for_alignment.at(pairI));
				}
			}
			combinedPairs_for_alignment = pairs_after_filtering;
		}

		std::cout << "\t" << "Now align " << combinedPairs_for_alignment.size() << " read pairs." << "\n" << std::flush;

		std::vector< std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> > withPairing_alignments_perThread;
		std::vector< std::vector<int> > withPairing_alignments_readPairI_perThread;

		withPairing_alignments_perThread.resize(outerThreads);
		withPairing_alignments_readPairI_perThread.resize(outerThreads);

		alignReadPairs(combinedPairs_for_alignment, withPairing_alignments_perThread, withPairing_alignments_readPairI_perThread, true, inserSize_mean_sd_perFile.at(fI).first, inserSize_mean_sd_perFile.at(fI).second);

		// merge

		std::cout  << Utilities::timestamp() << "\t\t" << "All pairs aligned - merge.\n" << std::flush;

		std::vector< std::pair<seedAndExtend_return_local, seedAndExtend_return_local> > withPairing_alignments;
		std::vector< int > withPairing_alignments_readPairI;
		for(unsigned int tI = 0; (int)tI < outerThreads; tI++)
		{
			withPairing_alignments.insert(withPairing_alignments.end(), withPairing_alignments_perThread.at(tI).begin(), withPairing_alignments_perThread.at(tI).end());
			withPairing_alignments_readPairI.insert(withPairing_alignments_readPairI.end(), withPairing_alignments_readPairI_perThread.at(tI).begin(), withPairing_alignments_readPairI_perThread.at(tI).end());
		}

		// Produce normal output file
		std::string alignments_output_file = FASTQ + ".aligned";
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Produce alignments file " << alignments_output_file << ".\n" << std::flush;

		// std::vector<std::pair<std::string, std::string> > alignments_readIDs;
		std::vector<oneReadPair> combinedPairs_for_alignment_inAlignmentOrder;
		for(unsigned int pairII = 0; pairII < withPairing_alignments_readPairI.size(); pairII++)
		{
			int readPairI = withPairing_alignments_readPairI.at(pairII);
			// std::string ID1 = combinedPairs_for_alignment.at(readPairI).reads.first.name;
			// std::string ID2 = combinedPairs_for_alignment.at(readPairI).reads.second.name;
			// alignments_readIDs.push_back(make_pair(ID1, ID2));
			combinedPairs_for_alignment_inAlignmentOrder.push_back(combinedPairs_for_alignment.at(readPairI));
		}

		assert(combinedPairs_for_alignment_inAlignmentOrder.size() == withPairing_alignments.size());


		printAlignmentsToFile(alignments_output_file, withPairing_alignments, combinedPairs_for_alignment_inAlignmentOrder, inserSize_mean_sd_perFile.at(fI).first, inserSize_mean_sd_perFile.at(fI).second);


		// Produce SAM
		std::string SAM_output_file = FASTQ + ".sam";
		std::cout  << Utilities::timestamp() << "\t\t\t" << "Produce SAM " << SAM_output_file << ".\n" << std::flush;


		// Get graph loci and reference positions
		std::ofstream SAM_output_stream;
		SAM_output_stream.open(SAM_output_file.c_str());
		assert(SAM_output_stream.is_open());

		std::vector<int> uncompressed_graph_referencePositions;
		std::vector<std::string> graphLoci = readGraphLoci(graphDir);

		int lastReferencePosition = -1;
		for(unsigned int i = 0; i < graphLoci.size(); i++)
		{
			std::string locusID = graphLoci.at(i);
			std::vector<std::string> locusParts = Utilities::split(locusID, "_");
			if(locusParts.size() != 3)
			{
				throw std::runtime_error("alignShortReadsToHLAGraph(..): Cannot decompose locus ID " +locusID);
			}
			int thisLocus_refPos = Utilities::StrtoI(locusParts.at(2));
			if(thisLocus_refPos != -1)
			{
				thisLocus_refPos = thisLocus_refPos + 1;
			}
			if((i == 0) || (lastReferencePosition != thisLocus_refPos))
			{
				uncompressed_graph_referencePositions.push_back(thisLocus_refPos);
				lastReferencePosition = thisLocus_refPos;
			}
			else
			{
				uncompressed_graph_referencePositions.push_back(-1);
			}
		}

		alignedShortReads2SAM(SAM_output_stream, uncompressed_graph_referencePositions, referenceChromosomes.at("ref"), withPairing_alignments, combinedPairs_for_alignment_inAlignmentOrder);

		std::cout  << Utilities::timestamp() << "\t\t\t" << "Done. Output in " << SAM_output_file << ".\n" << std::flush;
	}

	std::cout  << Utilities::timestamp() << "\t\t\t" << "All alignments done, free memory.\n" << std::flush;
	for(int tI = 0; tI < outerThreads; tI++)
	{
		delete(graphAligners.at(tI));
		delete(graphs.at(tI));
	}
}




void validateChromotypesVsVCF(std::string chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width)
{
	assert(1 == 0);
	assert(kMer_size == 31);

	/*
	 *
	 * Will not compile, missing arguments for dGS_evaluate...
	 *
	std::cout << Utilities::timestamp() << "Load PRG-Viterbi chromotypes...\n" << std::flush;
	diploidGenomeString chromotypes = readGenomeStringFromFile(chromotypes_file);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << Utilities::timestamp() << "Compress PRG-Viterbi chromotypes...\n" << std::flush;
	diploidGenomeString chromotypes_compressed = compressGenomeString(chromotypes);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << "Convert VCF to diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress VCF diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
	std::cout << "\tdone\n" << std::flush;
		
	std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
	diploidGenomeString gS_ref = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, true);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress reference diploid GS..\n" << std::flush;
	diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;	
	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;
	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";
		
	std::cout << "Resolve chromotypes from VCF..\n" << std::flush;
	diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Resolve PRG-Viterbi chromotypes..\n" << std::flush;
	diploidGenomeString chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << "\nEvaluate reference-only chromotypes\n" << std::flush;
	evaluate_dGS(gS_ref_compressed, &myGraph);
	std::cout << "\n\n" << std::flush;
	
	std::cout << "\nEvaluate VCF chromotypes\n" << std::flush;
	evaluate_dGS(VCF_chromotypes_resolved, &myGraph);
	std::cout << "\n\n" << std::flush;

	std::cout << "\nEvaluate PRG-Viterbi chromotypes\n" << std::flush;
	evaluate_dGS(chromotypes_resolved, &myGraph);
	std::cout << "\n\n" << std::flush;
	*/
}


void alignContigsToAllChromotypes(std::string chromotypes_file, std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width, std::string outputDir_contigs, std::string contigsFile_Fasta, std::string graphDir)
{
	std::string chromosomeID = "6";

//
//	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;
//	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
//	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;
//	myGraph.loadMultiColourBinary(deBruijnGraph);
//	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
//	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";
//
//	std::cout << "Resolve chromotypes from VCF..\n" << std::flush;
//	diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_chromotypes_compressed, &myGraph).second;
//	std::cout << "\tdone\n" << std::flush;
//
//	std::cout << "Resolve PRG-Viterbi chromotypes..\n" << std::flush;
//	diploidGenomeString chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(chromotypes_compressed, &myGraph).second;
//	std::cout << "\tdone\n" << std::flush;
//
//	std::cout << "Resolve amended PRG-Viterbi chromotypes..\n" << std::flush;
//	diploidGenomeString amended_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(amended_chromotypes_compressed, &myGraph).second;
//	std::cout << "\tdone\n" << std::flush;

//	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;
//	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
//	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;
//	myGraph.loadMultiColourBinary(deBruijnGraph);
//	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
//	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";
//
//	std::cout << "Resolve PRG-Viterbi chromotypes..\n" << std::flush;
//	diploidGenomeString chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(chromotypes_compressed, &myGraph).second;
//	std::cout << "\tdone\n" << std::flush;

	std::map<std::string, std::string> fasta_contigs = Utilities::readFASTA(contigsFile_Fasta);

	auto alignAllContigsToGraph = [&](Graph* graph, std::string specificOutputDir) -> void {

		std::string outputDir = outputDir_contigs + "/" + specificOutputDir;
		std::string statusFile = outputDir + "/all.status";

		std::string outputFile_indexedkMerUniqueness = outputDir_contigs + "/indexedkMers_" + specificOutputDir + ".txt";

		if(Utilities::readStatus(statusFile) == 1)
		{
			std::cout << "All contigs in " << outputDir << " aligned already - continue!\n";
			return;
		}
		std::cout << Utilities::timestamp() << "Create GraphAligner...\n" << std::flush;
		int aligner_kMerSize = 31;
		GraphAlignerUnique::GraphAlignerUnique gA(graph, aligner_kMerSize);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		gA.printkMerProfile(outputFile_indexedkMerUniqueness);



		for(std::map<std::string, std::string>::iterator contigIt = fasta_contigs.begin(); contigIt != fasta_contigs.end(); contigIt++)
		{
			std::string normalizedContigName = contigIt->first;
			
			std::function<bool(char)> nAlNum = [&](char c) {
				return (! std::isalnum(c));
			};	
			normalizedContigName.erase(std::remove_if(normalizedContigName.begin(), normalizedContigName.end(), nAlNum), normalizedContigName.end());

			std::string output_file = outputDir + "/" + normalizedContigName + ".alignment";
			std::string output_file_status = output_file + ".status";

			if(Utilities::readStatus(output_file_status) == 1)
			{
				continue;
			}

			Utilities::writeStatus(output_file_status, 0);

			ofstream outputContigsStream;
			outputContigsStream.open (output_file.c_str(), ios::out);
			if(! outputContigsStream.is_open())
			{
				throw std::runtime_error("Cannot open file "+output_file+" for writing.\n");
			}
			assert(outputContigsStream.is_open());

			std::string contigID = contigIt->first;
			std::string sequence = contigIt->second;

			std::cout << Utilities::timestamp() << "Now aligning: " << contigID << " ...\n" << std::flush;

			gA.setIterationsMainRandomizationLoop(11);
			gA.setThreads(6);
			std::vector<seedAndExtend_return_local> allBacktraces;
			seedAndExtend_return_local wholeString_alignment = gA.seedAndExtend_local(sequence, allBacktraces);

			assert(wholeString_alignment.graph_aligned_levels.size() == wholeString_alignment.graph_aligned.length());
			std::vector<std::string> graph_aligned_levels_Strings;
			for(unsigned int i = 0; i < wholeString_alignment.graph_aligned_levels.size(); i++)
			{
				graph_aligned_levels_Strings.push_back(Utilities::ItoStr(wholeString_alignment.graph_aligned_levels.at(i)));
			}

			std::vector<std::string> sequence2graph_certainty_Strings;
			for(unsigned int i = 0; i < wholeString_alignment.certainty_sequence2Graph.size(); i++)
			{
				sequence2graph_certainty_Strings.push_back(Utilities::DtoStr(wholeString_alignment.certainty_sequence2Graph.at(i)));
			}

			// std::cout << "\n - done - score: " << wholeString_alignment.Score << "\n";
			// std::cout << wholeString_alignment.graph_aligned << "\n";
			// std::cout << wholeString_alignment.sequence_aligned << "\n" << std::flush;

			outputContigsStream << ">" << contigID << " - graph - score " << wholeString_alignment.Score << "\n";
			outputContigsStream << wholeString_alignment.graph_aligned << "\n";
			outputContigsStream << ">" << contigID << " - graph - levels" << "\n";
			outputContigsStream << Utilities::join(graph_aligned_levels_Strings, " ") << "\n";
			outputContigsStream << ">" << contigID << " - sequence" << "\n";
			outputContigsStream << wholeString_alignment.sequence_aligned << "\n";
			outputContigsStream << ">" << contigID << " - sequence2graph_certainty" << "\n";
			outputContigsStream << Utilities::join(sequence2graph_certainty_Strings, " ") << "\n";
			outputContigsStream << ">" << contigID << " - uniqueness" << "\n";
			outputContigsStream << Utilities::ItoStr(wholeString_alignment.kMers_total) << " " << Utilities::ItoStr(wholeString_alignment.kMers_unique_total) << " " << Utilities::ItoStr(wholeString_alignment.kMers_unique_utilized) << "\n";

			outputContigsStream.close();

			Utilities::writeStatus(output_file_status, 1);
		}
		
		Utilities::writeStatus(statusFile, 1);
	};
	

	{
		std::vector<std::vector<int> > VCF_chromotypes_referencePositions;
		std::cout << "[VCF] Convert VCF to diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes = VCF2GenomeString(chromosomeID, VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, VCF_chromotypes_referencePositions);
		std::cout << "\tdone\n" << std::flush;

		std::cout << "[VCF] Compress VCF diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
		std::cout << "\tdone\n" << std::flush;

		std::cout << "[VCF] Generate graph from VCF...\n" << std::flush;
		Graph* gS_graph_VCF = genomeString2Graph(VCF_chromotypes_compressed, true);

		std::cout << "[VCF] Align contigs to VCF chromotypes...\n" << std::flush;
		alignAllContigsToGraph(gS_graph_VCF, "toVCF");

		gS_graph_VCF->freeMemory();
		delete(gS_graph_VCF);
	}

	{
		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Load PRG-Viterbi chromotypes...\n" << std::flush;
		diploidGenomeString Viterbi_chromotypes = readGenomeStringFromFile(chromotypes_file, true);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Compress PRG-Viterbi chromotypes...\n" << std::flush;
		diploidGenomeString Viterbi_chromotypes_compressed = compressGenomeString(Viterbi_chromotypes);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		_addPufferChromotypes(Viterbi_chromotypes_compressed);

		std::cout << "[PRG-Viterbi] Generate graph from compressed chromotypes...\n" << std::flush;
		Graph* gS_graph_Viterbi = genomeString2Graph(Viterbi_chromotypes_compressed, true);

		std::cout << "[PRG-Viterbi] Align contigs to Viterbi chromotypes...\n" << std::flush;
		alignAllContigsToGraph(gS_graph_Viterbi, "toViterbiChromotypes");

		gS_graph_Viterbi->freeMemory();
		delete(gS_graph_Viterbi);
	}

	{
		std::vector<std::vector<int> > reference_chromotype_referencePositions;
		std::cout << "[Reference] Get non-modified diploidGS for reference genome...\n" << std::flush;
		diploidGenomeString gS_ref = VCF2GenomeString(chromosomeID, VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, reference_chromotype_referencePositions, true);
		std::cout << "\tdone\n" << std::flush;

		std::cout << "[Reference] Compress reference diploid GS..\n" << std::flush;
		diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
		std::cout << "\tdone\n" << std::flush;

		std::cout << "[Reference] Generate graph from reference chromotypes...\n" << std::flush;
		Graph* gS_graph_referenceCompressed = genomeString2Graph(gS_ref_compressed, true);

		std::cout << "[Reference] Align contigs to reference chromotypes...\n" << std::flush;
		alignAllContigsToGraph(gS_graph_referenceCompressed, "toReference");

		gS_graph_referenceCompressed->freeMemory();
		delete(gS_graph_referenceCompressed);
	}

	{
		std::vector<int> amendedChromotypes_genomicGraphLoci;
		std::cout << Utilities::timestamp() << "Load PRG-amended chromotypes from " << amended_chromotypes_file << "\n" << std::flush;
		diploidGenomeString amended_chromotypes = readGenomeStringFromChromotypesFile(amended_chromotypes_file, kMer_size, graphDir, VCF_minRange, amendedChromotypes_genomicGraphLoci);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "Compress PRG-amended chromotypes...\n" << std::flush;
		diploidGenomeString amended_chromotypes_compressed = compressGenomeString(amended_chromotypes);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
		
		std::cout << "[PRG-Amended] Generate graph from amended chromotypes...\n" << std::flush;
		Graph* gS_graph_amended_chromotypes_compressedd = genomeString2Graph(amended_chromotypes_compressed, true);

		std::cout << "[PRG-Amended] Align contigs to amended chromotypes...\n" << std::flush;
		alignAllContigsToGraph(gS_graph_amended_chromotypes_compressedd, "toAmendedChromotypes");

		gS_graph_amended_chromotypes_compressedd->freeMemory();
		delete(gS_graph_amended_chromotypes_compressedd);
	}

	auto evaluateSingleAlignment = [](std::string alignment_graph, std::string alignment_sequence, std::vector<int> graph_levels, int& ret_graphIntrinsicGap_sequenceGap, int& ret_graphIntrinsicGap_sequenceCharacter, int& ret_matches, int& ret_mismatches, int& ret_graphNovelGap, int& ret_graphNonGap_sequenceGap) -> void {

		ret_graphIntrinsicGap_sequenceGap = 0;
		ret_graphIntrinsicGap_sequenceCharacter = 0;
		ret_matches = 0;
		ret_mismatches = 0;
		ret_graphNovelGap = 0;
		ret_graphNonGap_sequenceGap = 0;
		
		assert(alignment_graph.length() == alignment_sequence.length());
		assert(graph_levels.size() == alignment_graph.length());
		
		for(unsigned int posI = 0; posI < alignment_graph.length(); posI++)
		{
			std::string C_graph = alignment_graph.substr(posI, 1);
			std::string C_sequence = alignment_sequence.substr(posI, 1);
			int L = graph_levels.at(posI);
			
			if(L == -1)
			{
				assert(C_graph == "_");
				assert(C_sequence != "_");
				
				ret_graphNovelGap++;				
			}
			else
			{
				if(C_graph == "_")
				{
					assert(L != -1);
					if(C_sequence == "_")
					{
						ret_graphIntrinsicGap_sequenceGap++;
					}
					else
					{	
						ret_graphIntrinsicGap_sequenceCharacter++;
					}
				}
				else if(C_sequence == "_")
				{
					assert(C_graph != "_");
					ret_graphNonGap_sequenceGap++;					
				}
				else
				{
					assert(C_graph != "_");
					assert(C_sequence != "_");
					if(C_graph == C_sequence)
					{
						ret_matches++;
					}
					else
					{
						ret_mismatches++;
					}
				}
			}
		}
	};

	// unfiltered

	std::string summaryOutputFile = outputDir_contigs + "/contigEvaluation.txt";
	ofstream summaryOutputStream;
	summaryOutputStream.open(summaryOutputFile.c_str());

	if(! summaryOutputStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ summaryOutputFile);
	}

	std::string detailsOutputFile = outputDir_contigs + "/contigEvaluation_details.txt";
	ofstream detailsOutputStream;
	detailsOutputStream.open(detailsOutputFile.c_str());

	if(! detailsOutputStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ detailsOutputFile);
	}

	
	std::string chromotypeSupport_summary_file = outputDir_contigs + "/chromotypeSupport_summary.txt";
	ofstream chromotypeSupportSummaryStream;
	chromotypeSupportSummaryStream.open(chromotypeSupport_summary_file.c_str());
	assert(chromotypeSupportSummaryStream.is_open());

	// filtered

	std::string summaryOutputFile_filtered = outputDir_contigs + "/contigEvaluation_FILTERED.txt";
	ofstream summaryOutputStream_filtered;
	summaryOutputStream_filtered.open(summaryOutputFile_filtered.c_str());

	if(! summaryOutputStream_filtered.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ summaryOutputFile_filtered);
	}

	std::string detailsOutputFile_filtered = outputDir_contigs + "/contigEvaluation_details_FILTERED.txt";
	ofstream detailsOutputStream_filtered;
	detailsOutputStream_filtered.open(detailsOutputFile_filtered.c_str());

	if(! detailsOutputStream_filtered.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ detailsOutputFile_filtered);
	}

	std::string dismissedContigsFile_filtered = outputDir_contigs + "/contigEvaluation_dismissed_FILTERED.txt";
	ofstream dismissedContigsStream_filtered;
	dismissedContigsStream_filtered.open(dismissedContigsFile_filtered.c_str());
	if(! dismissedContigsStream_filtered.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ dismissedContigsFile_filtered);
	}


	std::string chromotypeSupport_summary_file_filtered = outputDir_contigs + "/chromotypeSupport_summary_FILTERED.txt";
	ofstream chromotypeSupportSummaryStream_filtered;
	chromotypeSupportSummaryStream_filtered.open(chromotypeSupport_summary_file_filtered.c_str());
	assert(chromotypeSupportSummaryStream_filtered.is_open());

	// non-filtered first lines

	summaryOutputStream <<
			"Method" << "\t" <<		
			"TotalAlignmentLength" << "\t" <<
			"AverageAlignmentSequenceCertainty" << "\t" <<
			"AveragekMerUniqueness" << "\t" <<
			"AverageUtilizedkMerUniqueness" << "\t" <<
			"graphIntrinsicGap_sequenceGap" << "\t" <<
			"graphIntrinsicGap_sequenceCharacter" << "\t" <<
			"Matches" << "\t" <<
			"Mismatches" << "\t" <<
			"graphNovelGap" << "\t" <<
			"graphNonGap_sequenceGap" << "\n";
			
	detailsOutputStream <<
			"Method" << "\t" <<
			"ID" << "\t" <<
			"AlignmentLength" << "\t" <<
			"AverageAlignmentSequenceCertainty" << "\t" <<
			"kMer_sequence" << "\t" <<
			"kMer_sequence_unique" << "\t" <<
			"kMer_sequence_unique_utilized" << "\t" <<
			"graphIntrinsicGap_sequenceGap" << "\t" <<
			"graphIntrinsicGap_sequenceCharacter" << "\t" <<
			"Matches" << "\t" <<
			"Mismatches" << "\t" <<
			"graphNovelGap" << "\t" <<
			"graphNonGap_sequenceGap" << "\n";

	chromotypeSupportSummaryStream <<
			"Method" << "\t" <<
			"AverageCoverage" << "\t" <<
			"AverageSupport" << "\t" <<
			"Length" << "\t" <<
			"Threshold_Support_perAllele" << "\t" <<
			"Positions_Support_0" << "\t" <<
			"Positions_Support_1" << "\t" <<
			"Positions_Support_2" << "\t" <<
			"Positions_Support_Fraction" << "\n";

	// filtered first lines

	summaryOutputStream_filtered <<
			"Method" << "\t" <<
			"TotalAlignmentLength" << "\t" <<
			"AverageAlignmentSequenceCertainty" << "\t" <<
			"AveragekMerUniqueness" << "\t" <<
			"AverageUtilizedkMerUniqueness" << "\t" <<
			"graphIntrinsicGap_sequenceGap" << "\t" <<
			"graphIntrinsicGap_sequenceCharacter" << "\t" <<
			"Matches" << "\t" <<
			"Mismatches" << "\t" <<
			"graphNovelGap" << "\t" <<
			"graphNonGap_sequenceGap" << "\n";

	detailsOutputStream_filtered <<
			"Method" << "\t" <<
			"ID" << "\t" <<
			"AlignmentLength" << "\t" <<
			"AverageAlignmentSequenceCertainty" << "\t" <<
			"kMer_sequence" << "\t" <<
			"kMer_sequence_unique" << "\t" <<
			"kMer_sequence_unique_utilized" << "\t" <<
			"graphIntrinsicGap_sequenceGap" << "\t" <<
			"graphIntrinsicGap_sequenceCharacter" << "\t" <<
			"Matches" << "\t" <<
			"Mismatches" << "\t" <<
			"graphNovelGap" << "\t" <<
			"graphNonGap_sequenceGap" << "\n";


	dismissedContigsStream_filtered <<
			"Method" << "\t" <<
			"ID" << "\t" <<
			"Dismissed" << "\n";

	chromotypeSupportSummaryStream_filtered <<
			"Method" << "\t" <<
			"AverageCoverage" << "\t" <<
			"AverageSupport" << "\t" <<
			"Length" << "\t" <<
			"Threshold_Support_perAllele" << "\t" <<
			"Positions_Support_0" << "\t" <<
			"Positions_Support_1" << "\t" <<
			"Positions_Support_2" << "\t" <<
			"Positions_Support_Fraction" << "\n";
			

	std::map<std::string, std::string> genomeReference = Utilities::readFASTA(referenceGenome);

	auto contigAlignment2SAM = [&](std::string specificOutputDir, ofstream& SAMoutputStream, diploidGenomeString& chromotypes, std::vector<std::vector<int> > chromotypes_referencePositions) -> void {
		assert(SAMoutputStream.is_open());
		assert(genomeReference.count(chromosomeID));

		std::vector<std::vector<std::string> > uncompressed_chromotypes;
		std::vector<int> uncompressed_chromotypes_referencePositions;
		for(unsigned int outerI = 0; outerI < chromotypes.size(); outerI++)
		{
			std::vector<std::string>& segment = chromotypes.at(outerI);
			assert((segment.size() == 1) || (segment.size() == 2));
			if(segment.size() == 2)
			{
				assert(segment.at(0).size() == segment.at(1).size());
			}

			unsigned int requiredSpace = segment.at(0).size();

			for(unsigned int posInSegment = 0; posInSegment < requiredSpace; posInSegment++)
			{
				std::vector<std::string> posVector;
				std::vector<double> posCoverageVector;
				std::vector<double> posSupportVector;
				std::vector<double> posInsertionsVector;

				if(segment.size() == 2)
				{
					posVector.push_back(segment.at(0).substr(posInSegment, 1));
					posVector.push_back(segment.at(1).substr(posInSegment, 1));
					posCoverageVector.push_back(0);
					posCoverageVector.push_back(0);
					posSupportVector.push_back(0);
					posSupportVector.push_back(0);
					posInsertionsVector.push_back(0);
					posInsertionsVector.push_back(0);
				}
				else
				{
					posVector.push_back(segment.at(0).substr(posInSegment, 1));
					posCoverageVector.push_back(0);
					posSupportVector.push_back(0);
					posInsertionsVector.push_back(0);
				}

				uncompressed_chromotypes.push_back(posVector);
				
				if(chromotypes_referencePositions.size() > 1)
				{
					uncompressed_chromotypes_referencePositions.push_back(chromotypes_referencePositions.at(outerI).at(posInSegment));
				}
			}
		}
		
		if(chromotypes_referencePositions.size() == 1)
		{
			uncompressed_chromotypes_referencePositions = chromotypes_referencePositions.at(0);
		}
		for(unsigned int lI = 0; lI < uncompressed_chromotypes_referencePositions.size(); lI++)
		{
			int referencePosition = uncompressed_chromotypes_referencePositions.at(lI);
			if(!((referencePosition == -1) || ((referencePosition >= 0) && (referencePosition < (int)genomeReference.at(chromosomeID).length()))))
			{
				std::cerr << "EARLY CHECK!" << "\n";
				std::cerr << "referencePosition" << ": " << referencePosition << "\n";
				std::cerr << "lI" << ": " << lI << "\n";
				std::cerr << "genomeReference.at(chromosomeID).length()" << ": " << genomeReference.at(chromosomeID).length() << "\n" << std::flush;
			}			
			assert((referencePosition == -1) || ((referencePosition >= 0) && (referencePosition < (int)genomeReference.at(chromosomeID).length())));
		}

		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size());

		std::cerr << "\n\nUncompressed chromotypes: " << uncompressed_chromotypes.size() << " levels!\n\n" << std::flush;

		std::string outputDir = outputDir_contigs + "/" + specificOutputDir;
		std::string statusFile = outputDir + "/all.status";
		if(Utilities::readStatus(statusFile) != 1)
		{
			throw std::runtime_error("Not all contigs in " + outputDir + " aligned -- abort!\n");
		}

		std::vector<std::string> files = filesInDirectory(outputDir);

		size_t contigNumber = 0;
		
		for(int fileI = 0; fileI < files.size(); fileI++)
		{
			std::string file = files.at(fileI);

			std::string statusSuffix = ".status";

			if(file.length() > statusSuffix.length())
			{
				std::string file_potentialSuffix = file.substr(file.length() - statusSuffix.length(), statusSuffix.length());
				assert(file_potentialSuffix.length() == statusSuffix.length());
				if(file_potentialSuffix == statusSuffix)
				{
					continue;
				}
			}

			std::map<std::string, std::string> aligned_contigs = Utilities::readFASTA(file, true);


			for(std::map<std::string, std::string>::iterator contigIt = aligned_contigs.begin(); contigIt != aligned_contigs.end(); contigIt++)
			{
				contigNumber++;
								
				std::string contigFullID = contigIt->first;
				std::string alignment_graph = contigIt->second;

				if(contigFullID.find("- graph - score ") != std::string::npos)
				{
				
					// todo remove later   
					if(contigNumber > 10)
					{
						//  return;
					}
				
					size_t position_graphString = contigFullID.find("- graph -");
					std::string string_without_graphString(contigFullID.begin(), contigFullID.begin() + position_graphString);
					std::string string_with_sequenceString = string_without_graphString + "- sequence";
					if(aligned_contigs.count(string_with_sequenceString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding sequence alignment //" << string_with_sequenceString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_sequenceString));

					std::string string_with_levelsString = string_without_graphString + "- graph - levels";
					if(aligned_contigs.count(string_with_levelsString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding levels string //" << string_with_levelsString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_levelsString));

					std::string string_with_certaintyString = string_without_graphString + "- sequence2graph_certainty";
					if(aligned_contigs.count(string_with_certaintyString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding certainty string //" << string_with_certaintyString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_certaintyString));

					std::string string_with_uniquenessString = string_without_graphString + "- uniqueness";
					if(aligned_contigs.count(string_with_uniquenessString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding certainty string //" << string_with_uniquenessString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_uniquenessString));



					std::string alignment_sequence = aligned_contigs.at(string_with_sequenceString);
					std::string alignment_levels = aligned_contigs.at(string_with_levelsString);
					std::string sequence2Graph_certainties = aligned_contigs.at(string_with_certaintyString);
					std::string sequence2Graph_uniquenesses = aligned_contigs.at(string_with_uniquenessString);

					// get levels for alignment
					std::vector<std::string> levels_separated_strings = Utilities::split(alignment_levels, " ");
					std::vector<int> levels_separated;
					for(unsigned int i = 0; i < levels_separated_strings.size(); i++)
					{
						std::string levelString = levels_separated_strings.at(i);
						int L = Utilities::StrtoI(levelString);
						levels_separated.push_back(L);
					}

					assert(levels_separated.size() == alignment_sequence.length());

					// get character certainties
					double certainties_sum = 0;
					std::vector<std::string> certainties_separated_strings = Utilities::split(sequence2Graph_certainties, " ");
					std::vector<int> certainties_separated;
					for(unsigned int i = 0; i < certainties_separated_strings.size(); i++)
					{
						std::string certaintyString = certainties_separated_strings.at(i);
						double C = Utilities::StrtoD(certaintyString);
						certainties_separated.push_back(C);
						certainties_sum += C;
					}
					assert(certainties_separated_strings.size() > 0);	
					double average_certainty = 1;
					if(certainties_separated_strings.size() > 0)
					{
						average_certainty = certainties_sum / (double)certainties_separated_strings.size();
					}
					
					// now go to levels and see how they support our chromotypes!

					std::vector<int> levels_separated_2_reference;
					levels_separated_2_reference.resize(levels_separated.size(), -1);

					int lastUsedLevel = -1;
					for(unsigned int levelI = 0; levelI < levels_separated.size(); levelI++)
					{
						int level = levels_separated.at(levelI);
						std::string graphCharacter = alignment_graph.substr(levelI, 1);
						std::string sequenceCharacter = alignment_sequence.substr(levelI, 1);

						if(level != -1)
						{
							assert(level < uncompressed_chromotypes.size());
							int referencePosition = uncompressed_chromotypes_referencePositions.at(level);
							
							if(!((referencePosition == -1) || ((referencePosition >= 0) && (referencePosition < genomeReference.at(chromosomeID).length()))))
							{
								std::cerr << "referencePosition" << ": " << referencePosition << "\n";
								std::cerr << "levelI" << ": " << levelI << "\n";
								std::cerr << "level" << ": " << level << "\n";
								std::cerr << "genomeReference.at(chromosomeID).length()" << ": " << genomeReference.at(chromosomeID).length() << "\n" << std::flush;
							}
							
							assert((referencePosition == -1) || ((referencePosition >= 0) && (referencePosition < genomeReference.at(chromosomeID).length())));
							
							levels_separated_2_reference.at(levelI) = referencePosition;
						}
					}

					std::vector<std::string> levels_separated_2_reference_String;
					std::vector<std::string> levels_referenceCharacters;
					std::vector<std::string> levels_graphCharacters;
					std::vector<std::string> levels_contigCharacters;

					int firstReferencePosition = -1;
					int lastReferencePosition = -1;
					for(unsigned int levelI = 0; levelI < levels_separated_2_reference.size(); levelI++)
					{
						levels_separated_2_reference_String.push_back(Utilities::ItoStr(levels_separated_2_reference.at(levelI)));
						int correspondingReferencePosition = levels_separated_2_reference.at(levelI);
						if(correspondingReferencePosition != -1)
						{
							if(! (((correspondingReferencePosition - 1) >= 0) && ((correspondingReferencePosition - 1) < (int)genomeReference.at(chromosomeID).length())))
							{
								std::cerr << "chromosomeID" << ":" << chromosomeID << "\n";
								std::cerr << "(correspondingReferencePosition - 1)" << ":" << (correspondingReferencePosition - 1) << "\n";
								std::cerr << "genomeReference.at(chromosomeID).length()" << ":" << genomeReference.at(chromosomeID).length() << "\n" << std::flush;
							}
							
							assert((correspondingReferencePosition - 1) >= 0);
							assert((correspondingReferencePosition - 1) < (int)genomeReference.at(chromosomeID).length());
							
							levels_referenceCharacters.push_back(genomeReference.at(chromosomeID).substr(correspondingReferencePosition - 1, 1));
						}
						else
						{
							levels_referenceCharacters.push_back("_");
						}
						levels_graphCharacters.push_back(alignment_graph.substr(levelI, 1));
						levels_contigCharacters.push_back(alignment_sequence.substr(levelI, 1));
						
						if(firstReferencePosition == -1)
						{
							if((levels_separated_2_reference.at(levelI) != -1) && (alignment_sequence.substr(levelI, 1) != "_"))
							{
								firstReferencePosition = levels_separated_2_reference.at(levelI);
							}
						}
						
						if((levels_separated_2_reference.at(levelI) != -1) && (alignment_sequence.substr(levelI, 1) != "_"))
						{
							lastReferencePosition = levels_separated_2_reference.at(levelI);
						}						
						
					}
					
					// assert((firstReferencePosition != -1) && (lastReferencePosition != -1));

					// std::cout << "Contig: " << contigFullID << "\n";
					// std::cout << Utilities::join(levels_separated_2_reference_String, " ") << "\n";
					// std::cout << Utilities::join(levels_referenceCharacters, " ") << "\n";
					// std::cout << Utilities::join(levels_graphCharacters, " ") << "\n";
					// std::cout << Utilities::join(levels_contigCharacters, " ") << "\n\n" << std::flush;
					
					std::string contigIDforSAM = string_without_graphString;
					contigIDforSAM.erase(std::remove_if(contigIDforSAM.begin(), contigIDforSAM.end(), [&](char c){return (((int)isalnum(c) == 0) || (c == ' '));}), contigIDforSAM.end());					
					
					std::string contigCharacters_noGaps = alignment_sequence;
					contigCharacters_noGaps.erase(std::remove_if(contigCharacters_noGaps.begin(),contigCharacters_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), contigCharacters_noGaps.end());
					
					size_t FLAGS_sizeT = 0;
					//FLAGS_sizeT = (FLAGS_sizeT || 0x2);
					if((firstReferencePosition != -1) || (lastReferencePosition != -1))
					{	
						FLAGS_sizeT = (FLAGS_sizeT || 0x4);
					}					
					
					int lastPrintedRealReferencePosition = -1;
					std::string CIGAR_uncompressed;
					for(unsigned int levelI = 0; levelI < levels_separated_2_reference.size(); levelI++)
					{
						int reference_position = levels_separated_2_reference.at(levelI);
						std::string alignmentCharacter = alignment_sequence.substr(levelI, 1);
						
						if(reference_position != -1)
						{
							if(lastPrintedRealReferencePosition != -1)
							{
								if(reference_position != (lastPrintedRealReferencePosition + 1))
								{
									if(!(reference_position > lastPrintedRealReferencePosition))
									{
										std::cerr << "reference_position" << ": " << reference_position << "\n";
										std::cerr << "lastPrintedRealReferencePosition" << ": " << lastPrintedRealReferencePosition << "\n" << std::flush;
									}
									assert(reference_position > lastPrintedRealReferencePosition);
									int missingReferencePositions = reference_position - lastPrintedRealReferencePosition - 1;
									assert(missingReferencePositions > 0);
									for(unsigned int j = 0; (int)j < missingReferencePositions; j++)
									{
										CIGAR_uncompressed.push_back('D');
									}
								}
							}	
						}
						
						if(reference_position == -1)
						{
							if(alignmentCharacter == "_")
							{
								// nothing
							}
							else
							{
								CIGAR_uncompressed.push_back('I');								
							}
						}
						else
						{
							if(alignmentCharacter == "_")
							{
								CIGAR_uncompressed.push_back('D');
							}
							else
							{
								CIGAR_uncompressed.push_back('M');							
							}
						}
						
						if(reference_position != -1)
						{
							lastPrintedRealReferencePosition = reference_position;
						}
					}

					std::string CIGAR_compressed;
					size_t CIGAR_sum_M = 0;
					size_t CIGAR_sum_I = 0;
					std::string currentOperation;
					int operationCount = 0;
					for(unsigned int cigarI = 0; cigarI < CIGAR_uncompressed.size(); cigarI++)
					{
						std::string thisOperation = CIGAR_uncompressed.substr(cigarI, 1);
						if(currentOperation.size() == 0)
						{
							currentOperation = thisOperation;
						}
						operationCount++;
						if((cigarI == (CIGAR_uncompressed.size() - 1)) || (CIGAR_uncompressed.substr(cigarI+1, 1) != currentOperation))
						{
							CIGAR_compressed += (Utilities::ItoStr(operationCount) + currentOperation);
							
							if(currentOperation == "M")
							{
								CIGAR_sum_M += operationCount;
							}
							else if(currentOperation == "I")
							{
								CIGAR_sum_I += operationCount;
							}
							operationCount = 0;
							currentOperation = "";
						}
					}
					
					assert(contigCharacters_noGaps.length() == (CIGAR_sum_M + CIGAR_sum_I));
					double P_mapping_wrong = 1 - average_certainty;
					if(P_mapping_wrong == 0)
					{
						P_mapping_wrong = 0.0001;
					}
					
					std::string QNAME = contigIDforSAM;
					std::string FLAG = Utilities::ItoStr(FLAGS_sizeT);
					std::string RNAME = chromosomeID;
					std::string POS = Utilities::ItoStr(firstReferencePosition);
					std::string MAPQ = Utilities::ItoStr(-10.0*log10(P_mapping_wrong)+0.5);
					std::string CIGAR = CIGAR_compressed;
					std::string RNEXT = "*";
					std::string PNEXT = "0";
					std::string TLEN = Utilities::ItoStr(lastReferencePosition - firstReferencePosition + 1);
					std::string SEQ = contigCharacters_noGaps;
					std::string QUAL = "*";
					
		
					SAMoutputStream <<
						QNAME << "\t" << 
						FLAG << "\t" << 
						RNAME << "\t" << 
						POS << "\t" << 
						MAPQ << "\t" << 
						CIGAR << "\t" << 
						RNEXT << "\t" << 
						PNEXT << "\t" << 
						TLEN << "\t" << 
						SEQ << "\t" << 
						QUAL << "\n" << std::flush;
					
/*			
1 QNAME String [!-?A-~]f1,255g Query template NAME
2 FLAG Int [0,216-1] bitwise FLAG
3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
4 POS Int [0,231-1] 1-based leftmost mapping POSition
5 MAPQ Int [0,28-1] MAPping Quality
6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
8 PNEXT Int [0,231-1] Position of the mate/next read
9 TLEN Int [-231+1,231-1] observed Template LENgth
10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
*/					
				}
			}
		}
	};

	auto evaluateContigAlignment = [&](std::string specificOutputDir, diploidGenomeString& chromotypes, std::vector<std::vector<int> > chromotypes_referencePositions, bool filter) -> void {

		std::string chromotypeSupport_output_file = outputDir_contigs + "/chromotypeSupport_by_"+specificOutputDir+".txt";
		if(filter)
		{
			chromotypeSupport_output_file = outputDir_contigs + "/chromotypeSupport_by_"+specificOutputDir+"_FILTERED.txt";
		}

		ofstream chromotypeSupportOutputStream;
		chromotypeSupportOutputStream.open(chromotypeSupport_output_file.c_str());

		std::vector<std::vector<std::string> > uncompressed_chromotypes;
		std::vector<std::vector<double> > uncompressed_chromotypes_coverage;		
		std::vector<std::vector<double> > uncompressed_chromotypes_support;
		std::vector<std::vector<double> > uncompressed_chromotypes_beforeInsertions;
		std::vector<std::vector<double> > uncompressed_chromotypes_beforeInsertions_II;
		std::vector<std::vector<double> > uncompressed_chromotypes_deletions;

		std::vector<std::vector<double> > uncompressed_chromotypes_mismatches;
		std::vector<std::vector<double> > uncompressed_chromotypes_matches;
		std::vector<std::vector<double> > uncompressed_chromotypes_graphIntrinsicGap_sequenceGap;
		std::vector<std::vector<double> > uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter;

		std::vector<int> uncompressed_chromotypes_referencePositions;		
		
		std::vector<int> uncompressed_chromotypes_nGaps;
		std::vector<int> uncompressed_chromotypes_isGap_1;
		std::vector<int> uncompressed_chromotypes_isGap_2;


		std::cerr << "chromotypes.size(): " << chromotypes.size() << "\n";
		std::cerr << "chromotypes_referencePositions.size(): " << chromotypes_referencePositions.size() << "\n" << std::flush;
		
		for(unsigned int outerI = 0; outerI < chromotypes.size(); outerI++)
		{
			std::vector<std::string>& segment = chromotypes.at(outerI);
			assert((segment.size() == 1) || (segment.size() == 2));
			if(segment.size() == 2)
			{
				assert(segment.at(0).size() == segment.at(1).size());
			}

			unsigned int requiredSpace = segment.at(0).size();

			for(unsigned int posInSegment = 0; posInSegment < requiredSpace; posInSegment++)
			{
				std::vector<std::string> posVector;
				std::vector<double> posCoverageVector;
				std::vector<double> posSupportVector;
				std::vector<double> posInsertionsVector;

				int isGap_1 = 0;
				int isGap_2 = 0;
				if(segment.size() == 2)
				{
					posVector.push_back(segment.at(0).substr(posInSegment, 1));
					posVector.push_back(segment.at(1).substr(posInSegment, 1));
					posCoverageVector.push_back(0);
					posCoverageVector.push_back(0);
					posSupportVector.push_back(0);
					posSupportVector.push_back(0);
					posInsertionsVector.push_back(0);
					posInsertionsVector.push_back(0);

					if(segment.at(0).substr(posInSegment, 1) == "_")
					{
						isGap_1 = 1;
					}

					if(segment.at(1).substr(posInSegment, 1) == "_")
					{
						isGap_2 = 1;
					}

				}
				else
				{
					posVector.push_back(segment.at(0).substr(posInSegment, 1));
					posCoverageVector.push_back(0);					
					posSupportVector.push_back(0);
					posInsertionsVector.push_back(0);

					if(segment.at(0).substr(posInSegment, 1) == "_")
					{
						isGap_1 = 1;
						isGap_2 = 1;
					}
				}

				int nGaps = isGap_1 + isGap_2;

				uncompressed_chromotypes.push_back(posVector);
				uncompressed_chromotypes_coverage.push_back(posCoverageVector);				
				uncompressed_chromotypes_support.push_back(posSupportVector);
				uncompressed_chromotypes_beforeInsertions.push_back(posInsertionsVector);
				uncompressed_chromotypes_beforeInsertions_II.push_back(posInsertionsVector);
				uncompressed_chromotypes_deletions.push_back(posInsertionsVector);

				uncompressed_chromotypes_mismatches.push_back(posInsertionsVector);
				uncompressed_chromotypes_matches.push_back(posInsertionsVector);
				uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.push_back(posInsertionsVector);
				uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.push_back(posInsertionsVector);

				
				uncompressed_chromotypes_nGaps.push_back(nGaps);
				uncompressed_chromotypes_isGap_1.push_back(isGap_1);
				uncompressed_chromotypes_isGap_2.push_back(isGap_2);

				if(chromotypes_referencePositions.size() > 1)
				{
					uncompressed_chromotypes_referencePositions.push_back(chromotypes_referencePositions.at(outerI).at(posInSegment));
				}
			}
		}  
		
		if(chromotypes_referencePositions.size() == 1)
		{
			uncompressed_chromotypes_referencePositions = chromotypes_referencePositions.at(0);
		}
		
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_nGaps.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_isGap_1.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_isGap_2.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_mismatches.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_matches.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.size());
		assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.size());

		std::vector<int> uncompressed_chromotypes_diploid;
		std::vector<int> uncompressed_chromotypes_losePhasing;
		for(unsigned int outerI = 0; outerI < chromotypes.size(); outerI++)
		{
			std::vector<std::string>& segment = chromotypes.at(outerI);
			assert((segment.size() == 1) || (segment.size() == 2));
			if(segment.size() == 2)
			{
				assert(segment.at(0).size() == segment.at(1).size());
			}

			for(unsigned int posInString = 0; posInString < segment.at(0).length(); posInString++)
			{
				uncompressed_chromotypes_diploid.push_back(segment.size() - 1);
				if((posInString == 0) && (segment.size() == 2))
				{
					uncompressed_chromotypes_losePhasing.push_back(1);
				}
				else
				{
					uncompressed_chromotypes_losePhasing.push_back(0);
				}
			}
		}

		assert(uncompressed_chromotypes_diploid.size() == uncompressed_chromotypes.size());
		assert(uncompressed_chromotypes_losePhasing.size() == uncompressed_chromotypes.size());

		std::vector<std::vector<std::map<std::string, double> >> uncompressed_chromotypes_contigInducedGenotypes;
		uncompressed_chromotypes_contigInducedGenotypes.resize(uncompressed_chromotypes.size());
		std::vector<std::vector<std::map<std::string, double> >> uncompressed_chromotypes_contigInducedInsertions;
		uncompressed_chromotypes_contigInducedInsertions.resize(uncompressed_chromotypes.size());

		for(unsigned int levelI = 0; levelI < uncompressed_chromotypes.size(); levelI++)
		{
			if(uncompressed_chromotypes_diploid.at(levelI) == 1)
			{
				uncompressed_chromotypes_contigInducedGenotypes.at(levelI).resize(2);
				uncompressed_chromotypes_contigInducedInsertions.at(levelI).resize(2);
			}
			else
			{
				uncompressed_chromotypes_contigInducedGenotypes.at(levelI).resize(1);
				uncompressed_chromotypes_contigInducedInsertions.at(levelI).resize(1);
			}
		}

		std::cerr << "\n\nUncompressed chromotypes: " << uncompressed_chromotypes.size() << " levels!\n\n" << std::flush;

		int allAlignment_length = 0;
		int allAlignments_graphIntrinsicGap_sequenceGap = 0;
		int allAlignments_graphIntrinsicGap_sequenceCharacter = 0;
		int allAlignments_matches = 0;
		int allAlignments_mismatches = 0;
		int allAlignments_graphNovelGap = 0;
		int allAlignments_graphNonGap_sequenceGap = 0;
		
		double allAlignments_certainties_sum = 0;
		double allAlignments_kMer_uniqueness_sum = 0;
		double allAlignments_kMer_utilized_uniqueness_sum = 0;

		size_t allAlignments_number = 0;

		std::string outputDir = outputDir_contigs + "/" + specificOutputDir;
		std::string statusFile = outputDir + "/all.status";

		if(Utilities::readStatus(statusFile) != 1)
		{
			throw std::runtime_error("Not all contigs in " + outputDir + " aligned -- abort!\n");
		}

		std::vector<std::string> files = filesInDirectory(outputDir);

		for(int fileI = 0; fileI < files.size(); fileI++)
		{
			std::string file = files.at(fileI);

			std::string statusSuffix = ".status";

			if(file.length() > statusSuffix.length())
			{
				std::string file_potentialSuffix = file.substr(file.length() - statusSuffix.length(), statusSuffix.length());
				assert(file_potentialSuffix.length() == statusSuffix.length());
				if(file_potentialSuffix == statusSuffix)
				{
					continue;
				}
			}

			// std::cerr << "\t Reading contig file " << file << "\n" << std::flush;

			std::map<std::string, std::string> aligned_contigs = Utilities::readFASTA(file, true);

			// std::cerr << "Contigs: " << aligned_contigs.size() << "\n" << std::flush;
			
			for(std::map<std::string, std::string>::iterator contigIt = aligned_contigs.begin(); contigIt != aligned_contigs.end(); contigIt++)
			{
				std::string contigFullID = contigIt->first;
				std::string alignment_graph = contigIt->second;

				if(contigFullID.find("- graph - score ") != std::string::npos)
				{
					size_t position_graphString = contigFullID.find("- graph -");
					std::string string_without_graphString(contigFullID.begin(), contigFullID.begin() + position_graphString);
					std::string string_with_sequenceString = string_without_graphString + "- sequence";
					if(aligned_contigs.count(string_with_sequenceString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding sequence alignment //" << string_with_sequenceString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_sequenceString));

					std::string string_with_levelsString = string_without_graphString + "- graph - levels";
					if(aligned_contigs.count(string_with_levelsString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding levels string //" << string_with_levelsString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_levelsString));

					std::string string_with_certaintyString = string_without_graphString + "- sequence2graph_certainty";
					if(aligned_contigs.count(string_with_certaintyString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding certainty string //" << string_with_certaintyString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_certaintyString));

					std::string string_with_uniquenessString = string_without_graphString + "- uniqueness";
					if(aligned_contigs.count(string_with_uniquenessString) == 0)
					{
						std::cerr << "For contig ID //" << contigFullID << "//, cannot find the corresponding certainty string //" << string_with_uniquenessString << "//\n" << std::flush;
					}
					assert(aligned_contigs.count(string_with_uniquenessString));



					std::string alignment_sequence = aligned_contigs.at(string_with_sequenceString);
					std::string alignment_levels = aligned_contigs.at(string_with_levelsString);
					std::string sequence2Graph_certainties = aligned_contigs.at(string_with_certaintyString);
					std::string sequence2Graph_uniquenesses = aligned_contigs.at(string_with_uniquenessString);

					// get levels for alignment
					std::vector<std::string> levels_separated_strings = Utilities::split(alignment_levels, " ");
					std::vector<int> levels_separated;
					for(unsigned int i = 0; i < levels_separated_strings.size(); i++)
					{
						std::string levelString = levels_separated_strings.at(i);
						int L = Utilities::StrtoI(levelString);
						levels_separated.push_back(L);
					}
					
					assert(levels_separated.size() == alignment_sequence.length());

					// get character certainties
					double certainties_sum = 0;
					std::vector<std::string> certainties_separated_strings = Utilities::split(sequence2Graph_certainties, " ");
					std::vector<int> certainties_separated;
					for(unsigned int i = 0; i < certainties_separated_strings.size(); i++)
					{
						std::string certaintyString = certainties_separated_strings.at(i);
						double C = Utilities::StrtoD(certaintyString);
						certainties_separated.push_back(C);
						certainties_sum += C;
					}
					assert(certainties_separated_strings.size() > 0);
					double certainties_average = certainties_sum / certainties_separated_strings.size();

					// get uniquenesses
					std::vector<std::string> uniquenesses_separated_strings = Utilities::split(sequence2Graph_uniquenesses, " ");
					int alignedSequence_kMers_total = Utilities::StrtoI(uniquenesses_separated_strings.at(0));
					int thisAlignment_kMers_unique = Utilities::StrtoI(uniquenesses_separated_strings.at(1));
					int thisAlignment_kMers_unique_utilized = Utilities::StrtoI(uniquenesses_separated_strings.at(2));

					double kMers_unique_average = 0;
					double kMers_unique_utilized_average = 0;

					if(alignedSequence_kMers_total > 0)
					{
						assert(thisAlignment_kMers_unique <= alignedSequence_kMers_total);
						assert(thisAlignment_kMers_unique_utilized <= alignedSequence_kMers_total);
						kMers_unique_average = (double)thisAlignment_kMers_unique / (double)alignedSequence_kMers_total;
						kMers_unique_utilized_average = (double)thisAlignment_kMers_unique_utilized / (double)alignedSequence_kMers_total;
					}

					if(filter)
					{
						double alignedSequenceLength = 0;
						double totalSequenceGaps = 0;
						double totalGraphGaps = 0;

						for(unsigned int levelI = 0; levelI < levels_separated.size(); levelI++)
						{
							int level = levels_separated.at(levelI);
							std::string graphCharacter = alignment_graph.substr(levelI, 1);
							std::string sequenceCharacter = alignment_sequence.substr(levelI, 1);

							if((graphCharacter == "_") && (sequenceCharacter == "_"))
							{
								assert(level != -1);
							}
							else if(graphCharacter == "_")
							{
								assert(sequenceCharacter != "_");
								if(level != -1)
								{
									// graph intrinsic gap, sequence non-gap -> "mismatch"
								}
								else
								{
									assert(level == -1);
									// proper graph gap
									totalGraphGaps++;
								}
							}
							else if(sequenceCharacter == "_")
							{
								assert(graphCharacter != "_");
								assert(level != -1);

								// proper sequence gap
								totalSequenceGaps++;
							}

							if(sequenceCharacter != "_")
							{
								alignedSequenceLength++;
							}
						}

						bool dismiss = false;
						if(totalGraphGaps >= (0.5* alignedSequenceLength))
						{
							dismiss = true;
						}
						if(totalSequenceGaps >= 150000)
						{
							dismiss = true;
						}

						if(dismiss)
						{
							dismissedContigsStream_filtered <<
									specificOutputDir << "\t" <<
									string_without_graphString << "\t" <<
									1 << "\n";

							continue;
						}
					}

					allAlignments_kMer_uniqueness_sum += kMers_unique_average;
					allAlignments_kMer_utilized_uniqueness_sum += kMers_unique_utilized_average;

					int thisAlignment_graphIntrinsicGap_sequenceGap = 0;
					int thisAlignment_graphIntrinsicGap_sequenceCharacter = 0;
					int thisAlignment_matches = 0;
					int thisAlignment_mismatches = 0;
					int thisAlignment_graphNovelGap = 0;
					int thisAlignment_graphNonGap_sequenceGap = 0;
		
					evaluateSingleAlignment(
							alignment_graph,
							alignment_sequence,
							levels_separated,
							thisAlignment_graphIntrinsicGap_sequenceGap,
							thisAlignment_graphIntrinsicGap_sequenceCharacter,
							thisAlignment_matches,
							thisAlignment_mismatches,
							thisAlignment_graphNovelGap,
							thisAlignment_graphNonGap_sequenceGap
					);

					// std::cerr << "\tContig alignment " << contigFullID << "\n";
					// std::cerr << "\t\t" << "thisAlignment_matches" << ": " << thisAlignment_matches << "\n";
					// std::cerr << "\t\t\t" << "thisAlignment_matches_doubleGap" << ": " << thisAlignment_matches_doubleGap << "\n";
					// std::cerr << "\t\t" << "thisAlignment_mismatches" << ": " << thisAlignment_mismatches << "\n";
					// std::cerr << "\t\t" << "thisAlignment_gaps_total" << ": " << thisAlignment_gaps_total << "\n";
					// std::cerr << "\t\t\t" << "thisAlignment_gaps_graph" << ": " << thisAlignment_gaps_graph << "\n";
					// std::cerr << "\t\t\t" << "thisAlignment_gaps_sequence" << ": " << thisAlignment_gaps_sequence << "\n" << std::flush;

					allAlignments_graphIntrinsicGap_sequenceGap += thisAlignment_graphIntrinsicGap_sequenceGap;
					allAlignments_graphIntrinsicGap_sequenceCharacter += thisAlignment_graphIntrinsicGap_sequenceCharacter;
					allAlignments_matches += thisAlignment_matches;
					allAlignments_mismatches += thisAlignment_mismatches;
					allAlignments_graphNovelGap += thisAlignment_graphNovelGap;
					allAlignments_graphNonGap_sequenceGap += thisAlignment_graphNonGap_sequenceGap;
					
					allAlignment_length += alignment_sequence.length();
					allAlignments_certainties_sum += certainties_average;
					allAlignments_number++;

					if(filter)
					{
						detailsOutputStream_filtered <<
								specificOutputDir << "\t" <<
								string_without_graphString << "\t" <<
								alignment_sequence.length() << "\t" <<
								certainties_average << "\t" <<
								alignedSequence_kMers_total << "\t" <<
								thisAlignment_kMers_unique << "\t" <<
								thisAlignment_kMers_unique_utilized << "\t" <<
								thisAlignment_graphIntrinsicGap_sequenceGap << "\t" <<
								thisAlignment_graphIntrinsicGap_sequenceCharacter << "\t" <<
								thisAlignment_matches << "\t" <<
								thisAlignment_mismatches << "\t" <<
								thisAlignment_graphNovelGap << "\t" <<
								thisAlignment_graphNonGap_sequenceGap << "\n";
					}
					else
					{
						detailsOutputStream <<
								specificOutputDir << "\t" <<
								string_without_graphString << "\t" <<
								alignment_sequence.length() << "\t" <<
								certainties_average << "\t" <<
								alignedSequence_kMers_total << "\t" <<
								thisAlignment_kMers_unique << "\t" <<
								thisAlignment_kMers_unique_utilized << "\t" <<
								thisAlignment_graphIntrinsicGap_sequenceGap << "\t" <<
								thisAlignment_graphIntrinsicGap_sequenceCharacter << "\t" <<
								thisAlignment_matches << "\t" <<
								thisAlignment_mismatches << "\t" <<
								thisAlignment_graphNovelGap << "\t" <<
								thisAlignment_graphNonGap_sequenceGap << "\n";
					}
					// now go to levels and see how they support our chromotypes!

					std::string runningInsertion;

					int lastUsedLevel = -1;
					for(unsigned int levelI = 0; levelI < levels_separated.size(); levelI++)
					{
						int level = levels_separated.at(levelI);
						std::string graphCharacter = alignment_graph.substr(levelI, 1);
						std::string sequenceCharacter = alignment_sequence.substr(levelI, 1);

						if(level != -1)
						{
							lastUsedLevel = level;
							
							assert(level < uncompressed_chromotypes.size());
							bool matchOneChromotype = false;
							unsigned int matchingChromotype;
							
							double supportMatchingChromotype = ((sequenceCharacter == graphCharacter) ? 1 : 0);
							
							if(uncompressed_chromotypes.at(level).size() == 1)
							{
								matchOneChromotype = (uncompressed_chromotypes.at(level).at(0) == graphCharacter);
								matchingChromotype = 0;
							}
							else
							{
								matchOneChromotype = ( 	(uncompressed_chromotypes.at(level).at(0) == graphCharacter) ||
														(uncompressed_chromotypes.at(level).at(1) == graphCharacter));

								if((uncompressed_chromotypes.at(level).at(0) == graphCharacter) && (uncompressed_chromotypes.at(level).at(1) == graphCharacter))
								{
									matchingChromotype = 2;
								}
								else if(uncompressed_chromotypes.at(level).at(0) == graphCharacter)
								{
									matchingChromotype = 0;
								}
								else if(uncompressed_chromotypes.at(level).at(1) == graphCharacter)
								{
									matchingChromotype = 1;
								}
							}

							assert(matchOneChromotype);

							if(matchingChromotype == 2)
							{
								uncompressed_chromotypes_coverage.at(level).at(0) += 0.5;
								uncompressed_chromotypes_coverage.at(level).at(1) += 0.5;
								
								uncompressed_chromotypes_support.at(level).at(0) += 0.5*supportMatchingChromotype;
								uncompressed_chromotypes_support.at(level).at(1) += 0.5*supportMatchingChromotype;

								if(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(0).count(sequenceCharacter) == 0)
								{
									uncompressed_chromotypes_contigInducedGenotypes.at(level).at(0)[sequenceCharacter] = 0;
								}
								uncompressed_chromotypes_contigInducedGenotypes.at(level).at(0).at(sequenceCharacter) += 0.5;
								if(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(1).count(sequenceCharacter) == 0)
								{
									uncompressed_chromotypes_contigInducedGenotypes.at(level).at(1)[sequenceCharacter] = 0;
								}
								uncompressed_chromotypes_contigInducedGenotypes.at(level).at(1).at(sequenceCharacter) += 0.5;

								if(runningInsertion.length() > 0)
								{
									if(uncompressed_chromotypes_contigInducedInsertions.at(level).at(0).count(runningInsertion) == 0)
									{
										uncompressed_chromotypes_contigInducedInsertions.at(level).at(0)[runningInsertion] = 0;
									}
									uncompressed_chromotypes_contigInducedInsertions.at(level).at(0).at(runningInsertion) += 0.5;
									uncompressed_chromotypes_beforeInsertions_II.at(level).at(0) += 0.5;

									if(uncompressed_chromotypes_contigInducedInsertions.at(level).at(1).count(runningInsertion) == 0)
									{
										uncompressed_chromotypes_contigInducedInsertions.at(level).at(1)[runningInsertion] = 0;
									}
									uncompressed_chromotypes_contigInducedInsertions.at(level).at(1).at(runningInsertion) += 0.5;
									uncompressed_chromotypes_beforeInsertions_II.at(level).at(1) += 0.5;
								}

								if((graphCharacter != "_") && (sequenceCharacter == "_"))
								{
									uncompressed_chromotypes_deletions.at(level).at(0) += 0.5;
									uncompressed_chromotypes_deletions.at(level).at(1) += 0.5;
								}

								if((graphCharacter != "_") && (sequenceCharacter != "_") && (sequenceCharacter == graphCharacter))
								{
									uncompressed_chromotypes_matches.at(level).at(0) += 0.5;
									uncompressed_chromotypes_matches.at(level).at(1) += 0.5;
									assert(supportMatchingChromotype == 1);
								}
								if((graphCharacter != "_") && (sequenceCharacter != "_") && (sequenceCharacter != graphCharacter))
								{
									uncompressed_chromotypes_mismatches.at(level).at(0) += 0.5;
									uncompressed_chromotypes_mismatches.at(level).at(1) += 0.5;
									assert(supportMatchingChromotype == 0);
								}

								if((graphCharacter == "_") && (sequenceCharacter == "_"))
								{
									uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(0) += 0.5;
									uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(1) += 0.5;
									assert(level != -1);
								}

								if((graphCharacter == "_") && (sequenceCharacter != "_"))
								{
									uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(0) += 0.5;
									uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(1) += 0.5;
									assert(level != -1);
								}

							}
							else
							{
								uncompressed_chromotypes_coverage.at(level).at(matchingChromotype)++;
								uncompressed_chromotypes_support.at(level).at(matchingChromotype) += supportMatchingChromotype;

								if(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(matchingChromotype).count(sequenceCharacter) == 0)
								{
									uncompressed_chromotypes_contigInducedGenotypes.at(level).at(matchingChromotype)[sequenceCharacter] = 0;
								}
								uncompressed_chromotypes_contigInducedGenotypes.at(level).at(matchingChromotype).at(sequenceCharacter)++;


								if(runningInsertion.length() > 0)
								{
									if(uncompressed_chromotypes_contigInducedInsertions.at(level).at(matchingChromotype).count(runningInsertion) == 0)
									{
										uncompressed_chromotypes_contigInducedInsertions.at(level).at(matchingChromotype)[runningInsertion] = 0;
									}
									uncompressed_chromotypes_contigInducedInsertions.at(level).at(matchingChromotype).at(runningInsertion)++;
									uncompressed_chromotypes_beforeInsertions_II.at(level).at(matchingChromotype)++;
								}

								if((graphCharacter != "_") && (sequenceCharacter == "_"))
								{
									uncompressed_chromotypes_deletions.at(level).at(matchingChromotype)++;
								}


								if((graphCharacter != "_") && (sequenceCharacter != "_") && (sequenceCharacter == graphCharacter))
								{
									uncompressed_chromotypes_matches.at(level).at(matchingChromotype)++;
									assert(supportMatchingChromotype == 1);
								}
								if((graphCharacter != "_") && (sequenceCharacter != "_") && (sequenceCharacter != graphCharacter))
								{
									uncompressed_chromotypes_mismatches.at(level).at(matchingChromotype)++;
									assert(supportMatchingChromotype == 0);
								}

								if((graphCharacter == "_") && (sequenceCharacter == "_"))
								{
									uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(matchingChromotype)++;
									assert(level != -1);
								}

								if((graphCharacter == "_") && (sequenceCharacter != "_"))
								{
									uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(matchingChromotype)++;
									assert(level != -1);
								}


							}

							runningInsertion = "";
						}
						else
						{
							runningInsertion.append(sequenceCharacter);

							int nextLevel = -1;
							std::string nextGraphCharacter;
							
							for(unsigned int levelI2 = levelI+1; levelI2 < levels_separated.size(); levelI2++)
							{
								if(levels_separated.at(levelI2) != -1)
								{
									nextLevel = levels_separated.at(levelI2);
									if(lastUsedLevel != -1)
									{
										if(!(nextLevel == (lastUsedLevel + 1)))
										{
											std::cerr << " !(nextLevel == (lastUsedLevel + 1)) \n";
											std::cerr << "nextLevel = " << nextLevel << "\n";
											std::cerr << "(lastUsedLevel + 1) = " << (lastUsedLevel + 1) << "\n";
										}
										assert(nextLevel == (lastUsedLevel + 1));
									}
									nextGraphCharacter = alignment_graph.substr(levelI2, 1);
									break;
								}
							}

							if(nextLevel == -1)
							{
								if((lastUsedLevel != -1) && ((lastUsedLevel + 1) < uncompressed_chromotypes.size()))
								{
									nextLevel = lastUsedLevel + 1;
								}
							}

							// make memory allocation for runningInsertion a bit more efficient...
							int wantCapacity = nextLevel - level + 1;
							if(wantCapacity > (int)runningInsertion.capacity())
							{
								runningInsertion.reserve(wantCapacity);
							}

							if(nextLevel != -1)
							{
								unsigned int matchingChromotype;
								if(nextGraphCharacter.length() == 0)
								{
									if(uncompressed_chromotypes.at(nextLevel).size() == 1)
									{
										matchingChromotype = 0;
									}
									else
									{
										matchingChromotype = 2;
									}
								}
								else
								{
									matchingChromotype = 5; // needs to be != 5 afterwards!

									if(uncompressed_chromotypes.at(nextLevel).size() == 2)
									{
										if((uncompressed_chromotypes.at(nextLevel).at(0) == nextGraphCharacter) && (uncompressed_chromotypes.at(nextLevel).at(1) == nextGraphCharacter))
										{
											matchingChromotype = 2;
										}
										else if(uncompressed_chromotypes.at(nextLevel).at(0) == nextGraphCharacter)
										{
											matchingChromotype = 0;
										}
										else if(uncompressed_chromotypes.at(nextLevel).at(1) == nextGraphCharacter)
										{
											matchingChromotype = 1;
										}
									}
									else
									{
										if(uncompressed_chromotypes.at(nextLevel).at(0) == nextGraphCharacter)
										{
											matchingChromotype = 0;
										}
									}

									if(matchingChromotype == 5)
									{
										std::cerr << "!!! matchingChromotype == 5 \n";
										std::cerr << "nextGraphCharacter = " << nextGraphCharacter << "\n";
										std::cerr << "nextLevel = " << nextLevel << "\n";
																		
										for(unsigned int j = 0; j < uncompressed_chromotypes.at(nextLevel).size(); j++)
										{
											std::cerr << "uncompressed_chromotypes.at(nextLevel).at(" << j << "): " << uncompressed_chromotypes.at(nextLevel).at(j) << "\n";
										}
										std::cerr << "\n" << std::flush;
									}
									assert(matchingChromotype != 5);
								}

								if(matchingChromotype == 2)
								{
									uncompressed_chromotypes_beforeInsertions.at(nextLevel).at(0) += 0.5;
									uncompressed_chromotypes_beforeInsertions.at(nextLevel).at(1) += 0.5;
								}
								else
								{
									uncompressed_chromotypes_beforeInsertions.at(nextLevel).at(matchingChromotype)++;
								}
							}
						}
					}

					if(runningInsertion.length() > 0)
					{
						if(lastUsedLevel != -1)
						{
							int insertionToLevel = lastUsedLevel + 1;
							unsigned int matchingChromotype;

							if(uncompressed_chromotypes.at(insertionToLevel).size() == 1)
							{
								matchingChromotype = 0;
							}
							else
							{
								matchingChromotype = 2;
							}

							if(matchingChromotype == 2)
							{

								if(runningInsertion.length() > 0)
								{
									if(uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(0).count(runningInsertion) == 0)
									{
										uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(0)[runningInsertion] = 0;
									}
									uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(0).at(runningInsertion) += 0.5;
									uncompressed_chromotypes_beforeInsertions_II.at(insertionToLevel).at(0) += 0.5;

									if(uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(1).count(runningInsertion) == 0)
									{
										uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(1)[runningInsertion] = 0;
									}
									uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(1).at(runningInsertion) += 0.5;
									uncompressed_chromotypes_beforeInsertions_II.at(insertionToLevel).at(1) += 0.5;

								}
							}
							else
							{
								if(uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(matchingChromotype).count(runningInsertion) == 0)
								{
									uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(matchingChromotype)[runningInsertion] = 0;
								}
								uncompressed_chromotypes_contigInducedInsertions.at(insertionToLevel).at(matchingChromotype).at(runningInsertion)++;
								uncompressed_chromotypes_beforeInsertions_II.at(insertionToLevel).at(matchingChromotype)++;
							}
						}
					}
				}
			}
		}
		

		std::cerr << "\n\nSUMMARY for " << specificOutputDir << ((filter) ? " -- filtered " : " -- not filtered") << "\n";
		std::cerr << "\t" << "allAlignments_graphIntrinsicGap_sequenceGap" << ": " << allAlignments_graphIntrinsicGap_sequenceGap << "\n";
		std::cerr << "\t\t" << "allAlignments_graphIntrinsicGap_sequenceCharacter" << ": " << allAlignments_graphIntrinsicGap_sequenceCharacter << "\n";
		std::cerr << "\t" << "allAlignments_matches" << ": " << allAlignments_matches << "\n";
		std::cerr << "\t" << "allAlignments_mismatches" << ": " << allAlignments_mismatches << "\n";
		std::cerr << "\t\t" << "allAlignments_graphNovelGap" << ": " << allAlignments_graphNovelGap << "\n";
		std::cerr << "\t\t" << "allAlignments_graphNonGap_sequenceGap" << ": " << allAlignments_graphNonGap_sequenceGap << "\n" << std::flush;

		double allAlignments_certainty_average = -1;
		double allAlignments_kMer_uniqueness_average = -1;
		double allAlignments_kMer_utilized_uniqueness_average = -1;
		if(allAlignments_number > 0)
		{
			allAlignments_certainty_average = allAlignments_certainties_sum / (double)allAlignments_number;
			allAlignments_kMer_uniqueness_average = allAlignments_kMer_uniqueness_sum / (double)allAlignments_number;
			allAlignments_kMer_utilized_uniqueness_average = allAlignments_kMer_utilized_uniqueness_sum / (double)allAlignments_number;
		} 

		if(filter)
		{
			summaryOutputStream_filtered <<
					specificOutputDir << "\t" <<
					allAlignment_length << "\t" <<
					allAlignments_certainty_average << "\t" <<
					allAlignments_kMer_uniqueness_average << "\t" <<
					allAlignments_kMer_utilized_uniqueness_average << "\t" <<
					allAlignments_graphIntrinsicGap_sequenceGap << "\t" <<
					allAlignments_graphIntrinsicGap_sequenceCharacter << "\t" <<
					allAlignments_matches << "\t" <<
					allAlignments_mismatches << "\t" <<
					allAlignments_graphNovelGap << "\t" <<
					allAlignments_graphNonGap_sequenceGap << "\n";
		}
		else
		{
			summaryOutputStream <<
					specificOutputDir << "\t" <<
					allAlignment_length << "\t" <<
					allAlignments_certainty_average << "\t" <<
					allAlignments_kMer_uniqueness_average << "\t" <<
					allAlignments_kMer_utilized_uniqueness_average << "\t" <<
					allAlignments_graphIntrinsicGap_sequenceGap << "\t" <<
					allAlignments_graphIntrinsicGap_sequenceCharacter << "\t" <<
					allAlignments_matches << "\t" <<
					allAlignments_mismatches << "\t" <<
					allAlignments_graphNovelGap << "\t" <<
					allAlignments_graphNonGap_sequenceGap << "\n";
		}

		double totalCoverage = 0;
		double totalSupport = 0;
		for(unsigned int level = 0; level < uncompressed_chromotypes_coverage.size(); level++)
		{
			for(unsigned int haplotypeI = 0; haplotypeI < uncompressed_chromotypes_coverage.at(level).size(); haplotypeI++)
			{
				totalCoverage += uncompressed_chromotypes_coverage.at(level).at(haplotypeI);
				totalSupport += uncompressed_chromotypes_support.at(level).at(haplotypeI);
			}
		}
		double averageCoverage = totalCoverage / (double) uncompressed_chromotypes_coverage.size();
		double averageSupport = totalSupport / (double) uncompressed_chromotypes_coverage.size();

		std::cerr << "[" << specificOutputDir << "] Coverage statistics: " << totalCoverage << " total, length " << uncompressed_chromotypes_support.size() << ", average: " << averageCoverage << "\n" << std::flush;
		std::cerr << "[" << specificOutputDir << "] Support statistics: " << totalSupport << " total, length " << uncompressed_chromotypes_support.size() << ", average: " << averageSupport << "\n" << std::flush;

		double supportThreshold = averageSupport / 8.0;
		std::cerr << "[" << specificOutputDir << "] Set support threshold " << supportThreshold << "\n\n" << std::flush;


		chromotypeSupportOutputStream <<
					"Level" << "\t" <<
					"ReferenceCoordinate" << "\t" <<
					"DiploidChromotype" << "\t" <<
					"ChromotypeLostPhase" << "\t" <<
					"ChromotypeCharacter_1" << "\t" <<
					"ChromotypeCharacter_2" << "\t" <<

					"Coverage_1" << "\t" <<
					"Coverage_2" << "\t" <<
					"Support_1" << "\t" <<
					"Support_2" << "\t" <<
					"graphNovelGap_I_1" << "\t"<<
					"graphNovelGap_I_2" << "\t" <<
					"graphNovelGap_II_1" << "\t"<<
					"graphNovelGap_II_2" << "\t" <<
					"graphNonGap_sequenceGap_1" << "\t"<<
					"graphNonGap_sequenceGap_2" << "\t"<<

					"Mismatches_1" << "\t" <<
					"Mismatches_2" << "\t" <<
					"Matches_1" << "\t" <<
					"Matches_2" << "\t" <<
					"graphIntrinsicGap_sequenceGap_1" << "\t" <<
					"graphIntrinsicGap_sequenceGap_2" << "\t" <<
					"graphIntrinsicGap_sequenceCharacter_1" << "\t" <<
					"graphIntrinsicGap_sequenceCharacter_2" << "\t" <<


					"graphNovelGap_I" << "\t" <<
					"graphNovelGap_II" << "\t" <<
					"graphNonGap_sequenceGap" << "\t" <<
					"Mismatches" << "\t" <<
					"Matches" << "\t" <<
					"graphIntrinsicGap_sequenceGap" << "\t" <<
					"graphIntrinsicGap_sequenceCharacter" << "\t" <<

					"chromotypeIsGap_1" << "\t" <<
					"chromotypeisGap_2" << "\t" <<
					"chromotypeGaps" << "\t" <<
					"chromotypeGaps_Support" << "\t" <<
					"ContigInducedGenotypes" << "\t" <<
					"ContigInducedInsertions" <<
					"\n";

		int totalSupport_0 = 0;
		int totalSupport_1 = 0;
		int totalSupport_2 = 0;

		for(unsigned int level = 0; level < uncompressed_chromotypes_coverage.size(); level++)
		{


			double Coverage_1;
			double Coverage_2;
			double Support_1;
			double Support_2;

			double graphNovelGap_I_1;
			double graphNovelGap_I_2;
			double graphNovelGap_II_1;
			double graphNovelGap_II_2;
			double graphNonGap_sequenceGap_1;
			double graphNonGap_sequenceGap_2;
			double Mismatches_1;
			double Mismatches_2;
			double Matches_1;
			double Matches_2;
			double graphIntrinsicGap_sequenceGap_1;
			double graphIntrinsicGap_sequenceGap_2;
			double graphIntrinsicGap_sequenceCharacter_1;
			double graphIntrinsicGap_sequenceCharacter_2;

			double gaps_Support = 0;

			std::string ChromotypeCharacter_1;
			std::string ChromotypeCharacter_2;
			std::string thisPosition_contigsInducedGenotype;
			std::string thisPosition_contigsInduducedInsertions;

			auto reduceHash = [](std::map<std::string, double> H) -> std::string {
				std::vector<std::string> genotypes;
				for(std::map<std::string, double>::iterator gtIt = H.begin(); gtIt != H.end(); gtIt++)
				{
					std::string genotype = gtIt->first;
					genotypes.push_back(genotype);
				}
				std::sort(genotypes.begin(), genotypes.end(), [&](std::string a, std::string b){
					return(H.at(b) < H.at(a));
				});
				if(genotypes.size() > 1)
				{
					assert(H.at(genotypes.at(0)) >= H.at(genotypes.at(1)));
				}

				assert(genotypes.size() == H.size());

				std::vector<std::string> forReturn_parts;
				for(unsigned int i = 0; i < genotypes.size(); i++)
				{
					std::string p = genotypes.at(i) + ":" + Utilities::DtoStr(H.at(genotypes.at(i)));
					forReturn_parts.push_back(p);
				}
				std::string forReturn = Utilities::join(forReturn_parts, ",");
				return forReturn;
			};

			if(uncompressed_chromotypes_coverage.at(level).size() == 1)
			{
				Coverage_1 = uncompressed_chromotypes_coverage.at(level).at(0) / 2.0;
				Coverage_2 = uncompressed_chromotypes_coverage.at(level).at(0) / 2.0;
				Support_1 = uncompressed_chromotypes_support.at(level).at(0) / 2.0;
				Support_2 = uncompressed_chromotypes_support.at(level).at(0) / 2.0;				
				graphNovelGap_I_1 = uncompressed_chromotypes_beforeInsertions.at(level).at(0) / 2.0;
				graphNovelGap_I_2 = uncompressed_chromotypes_beforeInsertions.at(level).at(0) / 2.0;
				graphNovelGap_II_1 = uncompressed_chromotypes_beforeInsertions_II.at(level).at(0) / 2.0;
				graphNovelGap_II_2 = uncompressed_chromotypes_beforeInsertions_II.at(level).at(0) / 2.0;
				graphNonGap_sequenceGap_1 = uncompressed_chromotypes_deletions.at(level).at(0) / 2.0;
				graphNonGap_sequenceGap_2 = uncompressed_chromotypes_deletions.at(level).at(0) / 2.0;

				Mismatches_1 = uncompressed_chromotypes_mismatches.at(level).at(0) / 2.0;
				Mismatches_2 =  uncompressed_chromotypes_mismatches.at(level).at(0) / 2.0;
				Matches_1 = uncompressed_chromotypes_matches.at(level).at(0) / 2.0;
				Matches_2 = uncompressed_chromotypes_matches.at(level).at(0) / 2.0;
				graphIntrinsicGap_sequenceGap_1 = uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(0) / 2.0;
				graphIntrinsicGap_sequenceGap_2 = uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(0) / 2.0;
				graphIntrinsicGap_sequenceCharacter_1 = uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(0) / 2.0;
				graphIntrinsicGap_sequenceCharacter_2 = uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(0) / 2.0;

				if(uncompressed_chromotypes_nGaps.at(level) > 0)
				{
					assert(uncompressed_chromotypes_nGaps.at(level) == 2);
					gaps_Support = uncompressed_chromotypes_support.at(level).at(0);
				}

				thisPosition_contigsInducedGenotype = reduceHash(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(0));
				thisPosition_contigsInduducedInsertions = reduceHash(uncompressed_chromotypes_contigInducedInsertions.at(level).at(0));

				ChromotypeCharacter_1 = uncompressed_chromotypes.at(level).at(0).substr(0, 1);
				ChromotypeCharacter_2 = uncompressed_chromotypes.at(level).at(0).substr(0, 1);
			}
			else
			{
				Coverage_1 = uncompressed_chromotypes_coverage.at(level).at(0);
				Coverage_2 = uncompressed_chromotypes_coverage.at(level).at(1);
				Support_1 = uncompressed_chromotypes_support.at(level).at(0);
				Support_2 = uncompressed_chromotypes_support.at(level).at(1);								
				graphNovelGap_I_1 = uncompressed_chromotypes_beforeInsertions.at(level).at(0);
				graphNovelGap_I_2 = uncompressed_chromotypes_beforeInsertions.at(level).at(1);
				graphNovelGap_II_1 = uncompressed_chromotypes_beforeInsertions_II.at(level).at(0);
				graphNovelGap_II_2 = uncompressed_chromotypes_beforeInsertions_II.at(level).at(1);
				graphNonGap_sequenceGap_1 = uncompressed_chromotypes_deletions.at(level).at(0);
				graphNonGap_sequenceGap_2 = uncompressed_chromotypes_deletions.at(level).at(1);

				Mismatches_1 = uncompressed_chromotypes_mismatches.at(level).at(0);
				Mismatches_2 =  uncompressed_chromotypes_mismatches.at(level).at(1);
				Matches_1 = uncompressed_chromotypes_matches.at(level).at(0);
				Matches_2 = uncompressed_chromotypes_matches.at(level).at(1);
				graphIntrinsicGap_sequenceGap_1 = uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(0);
				graphIntrinsicGap_sequenceGap_2 = uncompressed_chromotypes_graphIntrinsicGap_sequenceGap.at(level).at(1);
				graphIntrinsicGap_sequenceCharacter_1 = uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(0);
				graphIntrinsicGap_sequenceCharacter_2 = uncompressed_chromotypes_graphIntrinsicGap_sequenceCharacter.at(level).at(1);

				if(uncompressed_chromotypes_isGap_1.at(level) > 0)
				{
					assert(uncompressed_chromotypes_isGap_1.at(level) == 1);
					gaps_Support += uncompressed_chromotypes_support.at(level).at(0);
				}


				if(uncompressed_chromotypes_isGap_2.at(level) > 0)
				{
					assert(uncompressed_chromotypes_isGap_2.at(level) == 1);
					gaps_Support += uncompressed_chromotypes_support.at(level).at(1);
				}

				thisPosition_contigsInducedGenotype = "C1="+reduceHash(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(0))+"; C2="+reduceHash(uncompressed_chromotypes_contigInducedGenotypes.at(level).at(1));
				if((uncompressed_chromotypes_contigInducedInsertions.at(level).at(0).size() > 0) || (uncompressed_chromotypes_contigInducedInsertions.at(level).at(1).size() > 0))
				{
					thisPosition_contigsInduducedInsertions = "C1="+reduceHash(uncompressed_chromotypes_contigInducedInsertions.at(level).at(0))+"; C2="+reduceHash(uncompressed_chromotypes_contigInducedInsertions.at(level).at(1));
				}

				ChromotypeCharacter_1 = uncompressed_chromotypes.at(level).at(0).substr(0, 1);
				ChromotypeCharacter_2 = uncompressed_chromotypes.at(level).at(1).substr(0, 1);
			}


			std::vector<std::string> printFields;

			assert((uncompressed_chromotypes_isGap_1.at(level) + uncompressed_chromotypes_isGap_2.at(level)) == uncompressed_chromotypes_nGaps.at(level));
			printFields.push_back(Utilities::ItoStr(level));
			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_referencePositions.at(level)));
			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_diploid.at(level)));
			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_losePhasing.at(level)));
			printFields.push_back(ChromotypeCharacter_1);
			printFields.push_back(ChromotypeCharacter_2);

			printFields.push_back(Utilities::DtoStr(Coverage_1));
			printFields.push_back(Utilities::DtoStr(Coverage_2));
			printFields.push_back(Utilities::DtoStr(Support_1));
			printFields.push_back(Utilities::DtoStr(Support_2));

			printFields.push_back(Utilities::DtoStr(graphNovelGap_I_1));
			printFields.push_back(Utilities::DtoStr(graphNovelGap_I_2));
			printFields.push_back(Utilities::DtoStr(graphNovelGap_II_1));
			printFields.push_back(Utilities::DtoStr(graphNovelGap_II_2));
			printFields.push_back(Utilities::DtoStr(graphNonGap_sequenceGap_1));
			printFields.push_back(Utilities::DtoStr(graphNonGap_sequenceGap_2));
			printFields.push_back(Utilities::DtoStr(Mismatches_1));
			printFields.push_back(Utilities::DtoStr(Mismatches_2));
			printFields.push_back(Utilities::DtoStr(Matches_1));
			printFields.push_back(Utilities::DtoStr(Matches_2));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceGap_1));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceGap_2));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceCharacter_1));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceCharacter_2));

			printFields.push_back(Utilities::DtoStr(graphNovelGap_I_1+graphNovelGap_I_2));
			printFields.push_back(Utilities::DtoStr(graphNovelGap_II_1+graphNovelGap_II_2));
			printFields.push_back(Utilities::DtoStr(graphNonGap_sequenceGap_1+graphNonGap_sequenceGap_2));
			printFields.push_back(Utilities::DtoStr(Mismatches_1+Mismatches_2));
			printFields.push_back(Utilities::DtoStr(Matches_1+Matches_2));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceGap_1+graphIntrinsicGap_sequenceGap_2));
			printFields.push_back(Utilities::DtoStr(graphIntrinsicGap_sequenceCharacter_1+graphIntrinsicGap_sequenceCharacter_2));


			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_isGap_1.at(level)));
			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_isGap_2.at(level)));
			printFields.push_back(Utilities::ItoStr(uncompressed_chromotypes_nGaps.at(level)));
			printFields.push_back(Utilities::DtoStr(gaps_Support));

			printFields.push_back(thisPosition_contigsInducedGenotype);
			printFields.push_back(thisPosition_contigsInduducedInsertions);

			chromotypeSupportOutputStream << Utilities::join(printFields, "\t") << "\n";
			
			if((Support_1 >= supportThreshold) && (Support_2 >= supportThreshold))
			{
				totalSupport_2++;
			}
			else if((Support_1 >= supportThreshold) || (Support_2 >= supportThreshold))
			{
				totalSupport_1++;			
			}
			else
			{
				totalSupport_0++;
			}
		}
		
		std::cerr << "[" << specificOutputDir << "] Positions with: " << "\n";
		std::cerr << "\t\tSupport for 0 alleles: " << totalSupport_0 << "\n";
		std::cerr << "\t\tSupport for 1 alleles: " << totalSupport_1 << "\n";
		std::cerr << "\t\tSupport for 2 alleles: " << totalSupport_2 << "\n\n";
		
		int totalEvaluatedAlleles = 2 * (totalSupport_0 + totalSupport_1 + totalSupport_2);
		int totalAllelesOK = 2*totalSupport_2 + totalSupport_1;
		double alleleOKRate = (double)totalAllelesOK/(double)totalEvaluatedAlleles;
		
		std::cerr << "\t\tEvaluated alleles: " << totalEvaluatedAlleles << "\n\n";
		std::cerr << "\t\tEvaluated alleles OK: " << totalAllelesOK << "\n\n";
		std::cerr << "\t\tFraction OK: " << alleleOKRate << "\n\n" << std::flush;
		
		chromotypeSupportOutputStream.close();
			
		if(filter)
		{
			chromotypeSupportSummaryStream_filtered <<
					specificOutputDir << "\t" <<
					averageCoverage << "\t" <<
					averageSupport << "\t" <<
					uncompressed_chromotypes_support.size() << "\t" <<
					supportThreshold << "\t" <<
					totalSupport_0 << "\t" <<
					totalSupport_1 << "\t" <<
					totalSupport_2 << "\t" <<
					alleleOKRate << "\n" << std::flush;
		}
		else
		{
			chromotypeSupportSummaryStream <<
					specificOutputDir << "\t" <<
					averageCoverage << "\t" <<
					averageSupport << "\t" <<
					uncompressed_chromotypes_support.size() << "\t" <<
					supportThreshold << "\t" <<
					totalSupport_0 << "\t" <<
					totalSupport_1 << "\t" <<
					totalSupport_2 << "\t" <<
					alleleOKRate << "\n" << std::flush;
		}
	};

	
	auto initSAMOutputStream = [&](std::string specificOutputDir, ofstream& SAMoutputStream) -> void {
		std::string SAM_output_file = outputDir_contigs + "/SAM_"+specificOutputDir+".sam";
		SAMoutputStream.open(SAM_output_file.c_str());
		assert(SAMoutputStream.is_open());
		std::cerr << "Now writing: " << SAM_output_file << "\n";
		SAMoutputStream << "@HD VN:1.5" << "\t" << "SO:unknown" << "\n";
		SAMoutputStream << "@SQ SN:" << chromosomeID << "\t" << "LN:" << genomeReference.at(chromosomeID).length() << "\n";
	};
	
	{
		std::vector<std::vector<int> > VCF_chromotypes_referencePositions;
		std::cout << Utilities::timestamp() << "[VCF] Convert VCF to diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes = VCF2GenomeString(chromosomeID, VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, VCF_chromotypes_referencePositions);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[VCF] Compress VCF diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		// todo later make SAM production first step


		std::cout << Utilities::timestamp() << "[VCF] Evaluate contig alignment, filter out unreasonable alignments ...\n" << std::flush;

		evaluateContigAlignment("toVCF", VCF_chromotypes, VCF_chromotypes_referencePositions, true);

		std::cout << Utilities::timestamp() << "[VCF] Evaluation done.\n\n" << std::flush;


		// todo end

		std::cout << Utilities::timestamp() << "[VCF] Produce SAM file ...\n" << std::flush;

		ofstream SAMoutputStream;
		initSAMOutputStream("toVCF", SAMoutputStream);
		contigAlignment2SAM("toVCF", SAMoutputStream, VCF_chromotypes, VCF_chromotypes_referencePositions);
		SAMoutputStream.close();

		std::cout << std::flush;

		std::cout << Utilities::timestamp() << "[VCF] Evaluate contig alignment, don't filter alignments ...\n" << std::flush;

		evaluateContigAlignment("toVCF", VCF_chromotypes, VCF_chromotypes_referencePositions, false);

		std::cout << Utilities::timestamp() << "[VCF] Evaluation done.\n\n" << std::flush;



		std::cout << std::flush;

	}

	{
		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Load PRG-Viterbi chromotypes...\n" << std::flush;
		diploidGenomeString Viterbi_chromotypes = readGenomeStringFromFile(chromotypes_file, true);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Compress PRG-Viterbi chromotypes...\n" << std::flush;
		diploidGenomeString Viterbi_chromotypes_compressed = compressGenomeString(Viterbi_chromotypes);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Produce SAM file ...\n" << std::flush;
		
		ofstream SAMoutputStream;
		
		// get reference positions for Viterbi chromotypes

		std::vector<int> graphLoci_2_genomicPositions = getGenomicGraphLoci(graphDir, VCF_minRange);
		// for(unsigned int i = 0; i < 100; i++)
		// {
			// std::cerr << i << ": " << graphLoci_2_genomicPositions.at(i) << "\n" << std::flush;
		// }
		// assert(2 == 4);
		
		std::vector<std::vector<int>> graphLoci_2_genomicPositions_envelope;
		graphLoci_2_genomicPositions_envelope.push_back(graphLoci_2_genomicPositions);
		
		initSAMOutputStream("toViterbiChromotypes", SAMoutputStream);
		contigAlignment2SAM("toViterbiChromotypes", SAMoutputStream, Viterbi_chromotypes_compressed, graphLoci_2_genomicPositions_envelope);
		SAMoutputStream.close();
		
		
		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Evaluate contig alignment, don't filter alignments ...\n" << std::flush;
		
		evaluateContigAlignment("toViterbiChromotypes", Viterbi_chromotypes_compressed, graphLoci_2_genomicPositions_envelope, false);

		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Evaluation done.\n\n" << std::flush;


		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Evaluate contig alignment, filter out unreasonable alignments ...\n" << std::flush;

		evaluateContigAlignment("toViterbiChromotypes", Viterbi_chromotypes_compressed, graphLoci_2_genomicPositions_envelope, true);
				
		std::cout << Utilities::timestamp() << "[PRG-Viterbi] Evaluation done.\n\n" << std::flush;

		
	}

	{
		std::vector<std::vector<int> > reference_chromotypes_referencePositions;
		std::cout << Utilities::timestamp()  << "[Reference] Get non-modified diploidGS for reference genome...\n" << std::flush;
		diploidGenomeString gS_ref = VCF2GenomeString(chromosomeID, VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, reference_chromotypes_referencePositions, true);
		std::cout << Utilities::timestamp()  << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp()  << "[Reference] Compress reference diploid GS..\n" << std::flush;
		diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
		std::cout << Utilities::timestamp()  << "\tdone\n" << std::flush;
		
		std::cout << Utilities::timestamp() << "[Reference] Produce SAM file ...\n" << std::flush;
				
		ofstream SAMoutputStream;
		initSAMOutputStream("toReference", SAMoutputStream);
		contigAlignment2SAM("toReference", SAMoutputStream, gS_ref, reference_chromotypes_referencePositions);
		SAMoutputStream.close();
		

		std::cout << Utilities::timestamp() << "[Reference] Evaluate contig alignment, don't filter alignments ...\n" << std::flush;

		evaluateContigAlignment("toReference", gS_ref, reference_chromotypes_referencePositions, false);
				
		std::cout << Utilities::timestamp() << "[Reference] Evaluation done.\n\n" << std::flush;


		std::cout << Utilities::timestamp() << "[Reference] Evaluate contig alignment, filter out unreasonable alignments ...\n" << std::flush;

		evaluateContigAlignment("toReference", gS_ref, reference_chromotypes_referencePositions, true);

		std::cout << Utilities::timestamp() << "[Reference] Evaluation done.\n\n" << std::flush;
	}

	{
		std::vector<int> amendedChromotypes_genomicGraphLoci;
		std::cout << Utilities::timestamp() << "[PRG-Amended] Load PRG-amended chromotypes from " << amended_chromotypes_file << "\n" << std::flush;
		diploidGenomeString amended_chromotypes = readGenomeStringFromChromotypesFile(amended_chromotypes_file, kMer_size, graphDir, VCF_minRange, amendedChromotypes_genomicGraphLoci);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[PRG-Amended] Compress PRG-amended chromotypes...\n" << std::flush;
		diploidGenomeString amended_chromotypes_compressed = compressGenomeString(amended_chromotypes);
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << "[PRG-Amended] Produce SAM file ...\n" << std::flush;
		
		ofstream SAMoutputStream;
		std::vector<std::vector<int>> amendedChromotypes_genomicGraphLoci_envelope;
		amendedChromotypes_genomicGraphLoci_envelope.push_back(amendedChromotypes_genomicGraphLoci);		
		initSAMOutputStream("toAmendedChromotypes", SAMoutputStream);
		contigAlignment2SAM("toAmendedChromotypes", SAMoutputStream, amended_chromotypes_compressed, amendedChromotypes_genomicGraphLoci_envelope);
		SAMoutputStream.close();

				
		std::cout << Utilities::timestamp() << "[PRG-Amended] Evaluate contig alignment, don't filter alignments ...\n" << std::flush;
				
		evaluateContigAlignment("toAmendedChromotypes", amended_chromotypes_compressed, amendedChromotypes_genomicGraphLoci_envelope, false);
				
		std::cout << Utilities::timestamp() << "[PRG-Amended] Evaluation done.\n\n" << std::flush;


		std::cout << Utilities::timestamp() << "[PRG-Amended] Evaluate contig alignment, filter out unreasonable alignments ...\n" << std::flush;

		evaluateContigAlignment("toAmendedChromotypes", amended_chromotypes_compressed, amendedChromotypes_genomicGraphLoci_envelope, true);

		std::cout << Utilities::timestamp() << "[PRG-Amended] Evaluation done.\n\n" << std::flush;
	}	

  
	summaryOutputStream.close();
	chromotypeSupportSummaryStream.close();


	summaryOutputStream_filtered.close();
	detailsOutputStream_filtered.close();
	dismissedContigsStream_filtered.close();
	chromotypeSupportSummaryStream_filtered.close();
	detailsOutputStream.close();


}


void validateCompleteVCF(std::string VCFfile, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width, std::string outputDirectory)
{
	std::string configFileOutputPath= outputDirectory + "/validationDetails.txt";

	ofstream configOutputStream;

	configOutputStream.open(configFileOutputPath.c_str());
	if(! configOutputStream.is_open())
	{
		throw std::runtime_error("validateCompleteVCF(..): Want top open config summary file for writing, but can't! Path:\n"+configFileOutputPath);
	}
	configOutputStream << "k = " << kMer_size << "\n";
	configOutputStream << "deBruijnGraph" << " = " << deBruijnGraph << "\n";
	configOutputStream << "VCF" << " = " << VCFfile << "\n";
	configOutputStream << "referenceGenome" << " = " << referenceGenome << "\n";
	configOutputStream.close();

	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;
	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;

	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";

	std::vector<std::string> chromosomes;
	for(int i = 1; i <= 22; i++)
	{
		chromosomes.push_back(Utilities::ItoStr(i));
	}
	chromosomes.push_back("X");
	chromosomes.push_back("Y");

	
	// todo remove
	// chromosomes.clear();   
	// chromosomes.push_back("22");
	
	// chromosomes.clear();
	// chromosomes.push_back("6");
	
	std::string summaryFilePath = outputDirectory + "/summaryPerChromosome.txt";
	ofstream summaryFileStream;
	summaryFileStream.open(summaryFilePath.c_str());
	assert(summaryFileStream.is_open());

	std::vector<std::string> headerFields;
	headerFields.push_back("Method");
	headerFields.push_back("Total characters");
	headerFields.push_back("Total non-gap characters");
	headerFields.push_back("# kMers");
	headerFields.push_back("# kMers invalid");
	headerFields.push_back("# kMers present");
	headerFields.push_back("Unweighted optimality");
	headerFields.push_back("Coverage-weighted optimality");

	summaryFileStream << Utilities::join(headerFields, "\t") << "\n";

	for(unsigned int chromosomeI = 0; chromosomeI < chromosomes.size(); chromosomeI++)
	{
		std::string chromosome = chromosomes.at(chromosomeI);

		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", convert VCF to diploidGS...\n" << std::flush;
		std::vector<std::vector<int> > VCF_chromotypes_referencePositions;
		diploidGenomeString VCF_chromotypes = VCF2GenomeString(chromosome, -1, -1, VCFfile, referenceGenome, VCF_chromotypes_referencePositions);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", compress VCF diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", resolve chromotypes from VCF..\n" << std::flush;
		diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_chromotypes_compressed, &myGraph).second;
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

	
		
		std::set<std::string> set_kMers_reference;
		std::set<std::string> set_kMers_VCF;
		std::set<std::string> set_kMers_VCF_present;
		std::map<std::string, double> set_kMers_VCF_optimalities;


		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", evaluate VCF chromotypes\n" << std::flush;
		evaluate_dGS(VCF_chromotypes_resolved, VCF_chromotypes_compressed, set_kMers_reference, &myGraph, 0, 0, 0, "VCF"+chromosome, summaryFileStream, outputDirectory, VCF_chromotypes_referencePositions);
   
		std::cout << Utilities::timestamp() << "\n\tdone" <<  "\n\n" << std::flush;
				
		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", convert VCF (reference-only) to diploidGS...\n" << std::flush;
		
			
		std::vector<std::vector<int> > reference_chromotypes_referencePositions;
		std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
		diploidGenomeString gS_ref = VCF2GenomeString(chromosome, -1, -1, VCFfile, referenceGenome, reference_chromotypes_referencePositions, true);
		std::cout << "\tdone\n" << std::flush;		
		
		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", compress VCF (reference-only) diploidGS...\n" << std::flush;
		diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		set_kMers_reference.clear();
		set_kMers_VCF.clear();
		set_kMers_VCF_present.clear();
		set_kMers_VCF_optimalities.clear();
		
		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", evaluate VCF (reference-only) chromotypes\n" << std::flush;
		evaluate_dGS(gS_ref_compressed, gS_ref_compressed, set_kMers_reference, &myGraph, 0, 0, 0, "toReference"+chromosome, summaryFileStream, outputDirectory, reference_chromotypes_referencePositions);
   
		std::cout << Utilities::timestamp() << "\n\tdone" <<  "\n\n" << std::flush;
	}
}




void validateAllChromotypesVsVCF(std::string chromotypes_file, std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width, std::string outputDirectory, std::string graphDir)
{
	std::string configFileOutputPath= outputDirectory + "/validationDetails.txt";
	ofstream configOutputStream;
	configOutputStream.open(configFileOutputPath.c_str());
	if(! configOutputStream.is_open())
	{
		throw std::runtime_error("validateAllChromotypesVsVCF(..): Want top open config summary file for writing, but can't! Path:\n"+configFileOutputPath);
	}
	configOutputStream << "k = " << kMer_size << "\n";
	configOutputStream << "deBruijnGraph" << " = " << deBruijnGraph << "\n";
	configOutputStream << "PRG-Viterbi chromotypes" << " = " << chromotypes_file << "\n";
	configOutputStream << "PRG-amended chromotypes" << " = " << amended_chromotypes_file << "\n";
	configOutputStream << "VCF" << " = " << VCFfile << "\n";
	configOutputStream << "referenceGenome" << " = " << referenceGenome << "\n";
	configOutputStream << "PASS" << "\n";
	
	configOutputStream.close();

	assert(kMer_size == 31);

	std::cout << Utilities::timestamp() << "Load PRG-Viterbi chromotypes...\n" << std::flush;
	diploidGenomeString chromotypes = readGenomeStringFromFile(chromotypes_file, true);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << Utilities::timestamp() << "Compress PRG-Viterbi chromotypes...\n" << std::flush;
	diploidGenomeString chromotypes_compressed = compressGenomeString(chromotypes);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::vector<int> amendedChromotypes_genomicGraphLoci;
	std::cout << Utilities::timestamp() << "Load PRG-amended chromotypes from " << amended_chromotypes_file << "\n" << std::flush;
	diploidGenomeString amended_chromotypes = readGenomeStringFromChromotypesFile(amended_chromotypes_file, kMer_size, graphDir, VCF_minRange, amendedChromotypes_genomicGraphLoci);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << Utilities::timestamp() << "Compress PRG-amended chromotypes...\n" << std::flush;
	diploidGenomeString amended_chromotypes_compressed = compressGenomeString(amended_chromotypes);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	
	std::vector<std::vector<int> > VCF_chromotypes_referencePositions;
	std::cout << "Convert VCF to diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, VCF_chromotypes_referencePositions, false);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress VCF diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
	std::cout << "\tdone\n" << std::flush;
		
	std::vector<std::vector<int> > VCF_filtered_chromotypes_referencePositions;
	std::cout << "Convert VCF (filtered) to diploidGS...\n" << std::flush;
	diploidGenomeString VCF_filtered_chromotypes = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, VCF_filtered_chromotypes_referencePositions, false, true);
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << "Compress VCF (filtered) diploidGS...\n" << std::flush;
	diploidGenomeString VCF_filtered_chromotypes_compressed = compressGenomeString(VCF_filtered_chromotypes);
	std::cout << "\tdone\n" << std::flush;
	
	std::vector<std::vector<int> > reference_chromotypes_referencePositions;
	std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
	diploidGenomeString gS_ref = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, reference_chromotypes_referencePositions, true);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress reference diploid GS..\n" << std::flush;
	diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;	
	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;
	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";
		
	std::cout << "Resolve chromotypes from VCF..\n" << std::flush;
	diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Resolve chromotypes from VCF (filtered)..\n" << std::flush;
	diploidGenomeString VCF_filtered_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_filtered_chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << "Resolve PRG-Viterbi chromotypes..\n" << std::flush;
	diploidGenomeString chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Resolve amended PRG-Viterbi chromotypes..\n" << std::flush;
	diploidGenomeString amended_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(amended_chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;
	
	std::vector<std::string> set_names;
	std::vector<std::set<std::string>* > sets_kMers;
	std::vector<std::set<std::string>* > sets_kMers_present;
	std::vector<std::map<std::string, double>* > sets_kMers_optimalities;

	std::string summaryFilePath = outputDirectory + "/summaryPerMethod.txt";
	ofstream summaryFileStream;
	summaryFileStream.open(summaryFilePath.c_str());
	assert(summaryFileStream.is_open());

	std::vector<std::string> headerFields;
	headerFields.push_back("Method");
	headerFields.push_back("Total characters");
	headerFields.push_back("Total non-gap characters");
	headerFields.push_back("# kMers");
	headerFields.push_back("# kMers invalid");
	headerFields.push_back("# kMers present");
	headerFields.push_back("Unweighted optimality");
	headerFields.push_back("Coverage-weighted optimality");

	summaryFileStream << Utilities::join(headerFields, "\t") << "\n";

	
	
	std::set<std::string> set_kMers_reference;
	std::set<std::string> set_kMers_reference_present;
	std::map<std::string, double> set_kMers_reference_optimalities;
	std::cout << "\nEvaluate reference-only chromotypes\n" << std::flush;
	evaluate_dGS(gS_ref_compressed, gS_ref_compressed, set_kMers_reference, &myGraph, &set_kMers_reference, &set_kMers_reference_present, &set_kMers_reference_optimalities, "toReference", summaryFileStream, outputDirectory, reference_chromotypes_referencePositions);
	std::cout << "\n\n" << std::flush;
	set_names.push_back("Reference");
	sets_kMers.push_back(&set_kMers_reference);
	sets_kMers_present.push_back(&set_kMers_reference_present);
	sets_kMers_optimalities.push_back(&set_kMers_reference_optimalities);


	std::set<std::string> set_kMers_VCF;
	std::set<std::string> set_kMers_VCF_present;
	std::map<std::string, double> set_kMers_VCF_optimalities;
	std::cout << "\nEvaluate VCF chromotypes\n" << std::flush;
	evaluate_dGS(VCF_chromotypes_resolved, VCF_chromotypes_compressed, set_kMers_reference, &myGraph, &set_kMers_VCF, &set_kMers_VCF_present, &set_kMers_VCF_optimalities, "toVCF", summaryFileStream, outputDirectory, VCF_chromotypes_referencePositions);
	std::cout << "\n\n" << std::flush;
	set_names.push_back("VCF");
	sets_kMers.push_back(&set_kMers_VCF);
	sets_kMers_present.push_back(&set_kMers_VCF_present);
	sets_kMers_optimalities.push_back(&set_kMers_VCF_optimalities);


	std::set<std::string> set_kMers_VCF_filtered;
	std::set<std::string> set_kMers_VCF_filtered_present;
	std::map<std::string, double> set_kMers_VCF_filtered_optimalities;
	std::cout << "\nEvaluate VCF (filtered) chromotypes\n" << std::flush;
	evaluate_dGS(VCF_filtered_chromotypes_resolved, VCF_filtered_chromotypes_compressed, set_kMers_reference, &myGraph, &set_kMers_VCF_filtered, &set_kMers_VCF_filtered_present, &set_kMers_VCF_filtered_optimalities, "toFilteredVCF", summaryFileStream, outputDirectory, VCF_filtered_chromotypes_referencePositions);
	std::cout << "\n\n" << std::flush;
	// set_names.push_back("FilteredVCF");
	// sets_kMers.push_back(&set_kMers_VCF_filtered);
	// sets_kMers_present.push_back(&set_kMers_VCF_filtered_present);
	// sets_kMers_optimalities.push_back(&set_kMers_VCF_filtered_optimalities);
	
	std::set<std::string> set_kMers_viterbiPRG;
	std::set<std::string> set_kMers_viterbiPRG_present;
	std::map<std::string, double> set_kMers_viterbiPRG_optimalities;
	std::cout << "\nEvaluate PRG-Viterbi chromotypes\n" << std::flush;
	std::vector<int> Viterbi_graphLoci_2_genomicPositions = getGenomicGraphLoci(graphDir, VCF_minRange);
	std::vector<std::vector<int>> Viterbi_graphLoci_2_genomicPositions_envelope;
	Viterbi_graphLoci_2_genomicPositions_envelope.push_back(Viterbi_graphLoci_2_genomicPositions);	
	evaluate_dGS(chromotypes_resolved, chromotypes_compressed, set_kMers_reference, &myGraph, &set_kMers_viterbiPRG, &set_kMers_viterbiPRG_present, &set_kMers_viterbiPRG_optimalities, "toViterbiChromotypes", summaryFileStream, outputDirectory, Viterbi_graphLoci_2_genomicPositions_envelope);
	std::cout << "\n\n" << std::flush;
	set_names.push_back("PRG-Viterbi");
	sets_kMers.push_back(&set_kMers_viterbiPRG);
	sets_kMers_present.push_back(&set_kMers_viterbiPRG_present);
	sets_kMers_optimalities.push_back(&set_kMers_viterbiPRG_optimalities);


	std::set<std::string> set_kMers_amendedPRG;
	std::set<std::string> set_kMers_amendedPRG_present;
	std::map<std::string, double> set_kMers_amendedPRG_optimalities;
	std::cout << "\nEvaluate amended PRG-Viterbi chromotypes\n" << std::flush;
	std::vector<std::vector<int>> amendedChromotypes_genomicGraphLoci_envelope;
	amendedChromotypes_genomicGraphLoci_envelope.push_back(amendedChromotypes_genomicGraphLoci);		  
	evaluate_dGS(amended_chromotypes_resolved, amended_chromotypes_compressed, set_kMers_reference, &myGraph, &set_kMers_amendedPRG, &set_kMers_amendedPRG_present, &set_kMers_amendedPRG_optimalities, "toAmendedChromotypes", summaryFileStream, outputDirectory, amendedChromotypes_genomicGraphLoci_envelope);
	std::cout << "\n\n" << std::flush;
	set_names.push_back("PRG-Amended");
	sets_kMers.push_back(&set_kMers_amendedPRG);
	sets_kMers_present.push_back(&set_kMers_amendedPRG_present);
	sets_kMers_optimalities.push_back(&set_kMers_amendedPRG_optimalities);

	summaryFileStream.close();

		
	std::string vennDiagramsOutputFile = outputDirectory + "/vennDiagrams.txt";

	vennDiagrams(set_names, sets_kMers, sets_kMers_present, sets_kMers_optimalities, vennDiagramsOutputFile);

}

void validateAmendedChromotypesVsVCF(std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width)
{
	assert(kMer_size == 31);

	assert(1 == 0);
	/*
	 *
	 * won't compile
	std::cout << Utilities::timestamp() << "Load PRG-amended chromotypes from " << amended_chromotypes_file << "\n" << std::flush;
	diploidGenomeString chromotypes = readGenomeStringFromChromotypesFile(amended_chromotypes_file, kMer_size);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << Utilities::timestamp() << "Compress PRG-Viterbi chromotypes...\n" << std::flush;
	diploidGenomeString chromotypes_compressed = compressGenomeString(chromotypes);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

	std::cout << "Convert VCF to diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress VCF diploidGS...\n" << std::flush;
	diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
	diploidGenomeString gS_ref = VCF2GenomeString("6", VCF_minRange, VCF_maxRange, VCFfile, referenceGenome, true);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress reference diploid GS..\n" << std::flush;
	diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
	std::cout << "\tdone\n" << std::flush;

	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;
	DeBruijnGraph<1, 31, 1> myGraph(cortex_height, cortex_width);
	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;
	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";

	std::cout << "Resolve chromotypes from VCF..\n" << std::flush;
	diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(VCF_chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Resolve PRG-Viterbi chromotypes..\n" << std::flush;
	diploidGenomeString chromotypes_resolved = greedilyResolveDiploidKMerString<1, 31, 1>(chromotypes_compressed, &myGraph).second;
	std::cout << "\tdone\n" << std::flush;

	std::vector<std::set<std::string>* > sets_kMers;
	std::vector<std::set<std::string>* > sets_kMers_present;
	std::vector<std::map<std::string, double>* > sets_kMers_optimalities;


	std::cout << "\nEvaluate reference-only chromotypes\n" << std::flush;
	evaluate_dGS(gS_ref_compressed, &myGraph);

	std::cout << "\nEvaluate VCF-induced chromotypes\n" << std::flush;
	evaluate_dGS(VCF_chromotypes_resolved, &myGraph);
	std::cout << "\n\n" << std::flush;

	std::cout << "\nEvaluate amended PRG-Viterbi chromotypes\n" << std::flush;
	evaluate_dGS(chromotypes_resolved, &myGraph);
	*/
}


void testValidation(std::string viterbi_diploid_gS, std::string deBruijnGraph, int k, std::string VCFfile, std::string referenceGenomeFile)
{
	assert(k == 31);
	assert(1 == 0);

	/*
	 *
	 * won't compile
	DeBruijnGraph<1, 31, 1> myGraph(7, 100);
	std::cout << "Graph allocated, loading binary...\n" << std::flush;
	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";

	std::cout << "Convert VCF to diploidGS...\n" << std::flush;
	diploidGenomeString gS = VCF2GenomeString("1", 1, 10000, VCFfile, referenceGenomeFile);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
	diploidGenomeString gS_ref = VCF2GenomeString("1", 1, 10000, VCFfile, referenceGenomeFile, true);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "Compress reference diploid GS..\n" << std::flush;
	diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "(Unnecessarily) resolve reference diploid GS..\n" << std::flush;
	greedilyResolveDiploidKMerString<1, 31, 1>(gS_ref, &myGraph);
	std::cout << "\tdone\n" << std::flush;
	
	std::cout << "Resolve diploid GS...\n" << std::flush;
	std::pair<diploidGenomeString, diploidGenomeString> resolved_gS = greedilyResolveDiploidKMerString<1, 31, 1>(gS, &myGraph);
	std::cout << "\tdone\n" << std::flush;

	std::cout << "\nTest evaluation of reference gS\n" << std::flush;
	evaluate_dGS(gS_ref_compressed, &myGraph);
	std::cout << "\n\n" << std::flush;

//	std::cout << "Test evaluation of original string\n" << std::flush;
//	evaluate_dGS(gS, &myGraph);
//	std::cout << "\n\n" << std::flush;

	std::cout << "Test evaluation of called string\n" << std::flush;
	evaluate_dGS(resolved_gS.second, &myGraph);
	std::cout << "\n\n" << std::flush;
	*/
}

std::vector<double> kMer_PP(size_t observedCoverage, size_t upperBound, double coverage)
{
	std::vector<double> forReturn;
	forReturn.resize(upperBound + 1, 0);

	assert(upperBound < 1000);
	assert(coverage > 0);
	
	double sum_to_upper = 0;
	if(observedCoverage >= (upperBound * coverage))
	{
		forReturn.clear();
		forReturn.resize(upperBound+1, 0);
		forReturn.at(upperBound) = 1;
		sum_to_upper = 1;	
	}
	else
	{
		for(unsigned int underlyingCount = 0; underlyingCount <= (upperBound+1); underlyingCount++)
		{
			double lambda;
			if(underlyingCount == 0)
			{
				lambda = 0.01;
			}
			else
			{
				lambda = underlyingCount * coverage;
			}

			poisson_up poisson(lambda);
			double likelihood = pdf(poisson, observedCoverage);

			if(underlyingCount <= upperBound)
			{
				forReturn.at(underlyingCount) = likelihood;
				sum_to_upper += likelihood;
			}
			else
			{
				if(likelihood > forReturn.at(upperBound))
				{
					forReturn.clear();
					forReturn.resize(upperBound+1, 0);
					forReturn.at(upperBound) = 1;
					sum_to_upper = 1;
				}
			}
		}
	}
	
	if(sum_to_upper == 0)
	{
		std::cerr << "sum_to_upper = 0!!!\n";
		std::cerr << "observedCoverage: " << observedCoverage << "\n";
		std::cerr << "upperBound: " << upperBound << "\n";
		std::cerr << "coverage: " << coverage << "\n";
			
		for(unsigned int underlyingCount = 0; underlyingCount <= (upperBound+1); underlyingCount++)
		{
			double lambda;
			if(underlyingCount == 0)
			{
				lambda = 0.01;
			}
			else
			{
				lambda = underlyingCount * coverage;
			}

			poisson_up poisson(lambda);
			double likelihood = pdf(poisson, observedCoverage);

			std::cerr << "\t underlyingCount = " << underlyingCount << " => likelihood " << likelihood << "\n";
		}
	}
	
	assert(sum_to_upper != 0);
	

	for(unsigned int underlyingCount = 0; underlyingCount <= upperBound; underlyingCount++)
	{
		forReturn.at(underlyingCount) = forReturn.at(underlyingCount) / sum_to_upper;
	}

	double p_sum = 0;
	for(unsigned int underlyingCount = 0; underlyingCount <= upperBound; underlyingCount++)
	{
		double thisElement = forReturn.at(underlyingCount);
		assert(thisElement >= 0);
		assert(thisElement <= 1);
		p_sum += thisElement;
	}
	assert(abs(p_sum - 1) < 1e-5);
	
	return forReturn;
}

template<int m, int k, int colours>
void evaluate_dGS(diploidGenomeString& gS, diploidGenomeString& gS_unresolved, const std::set<std::string>& kMers_reference, DeBruijnGraph<m, k, colours>* graph, std::set<std::string>* kMers_in_dGS, std::set<std::string>* kMers_in_dGS_in_sample, std::map<std::string, double>* kMers_in_dGS_optimality, std::string nameForSummary, ofstream& summaryFileStream, std::string pathForSpatialSummary, std::vector<std::vector<int> > chromotypes_referencePositions)
{
	// We only want to deal with completely resolved genomestrings up to level k
	// That is, we want genomestrings which
	// - have all subsequent homozygous stretches connected
	// - have at least k homozygous characters between any two heterozygous positions
	// - so that if we connect them together, all k-mers are totally determined

	size_t totalLevels = 0;
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "A" << "\n" << std::flush; // todo remove
	
	for(unsigned int l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			if(!((l == 0) || (gS.at(l-1).size() != 1)))
			{
				std::cerr << "!((l == 0) || (gS.at(l-1).size() != 1)), l = " << l << ", gS.at(l-1).size() = " << gS.at(l-1).size() << ", gS.size()= " << gS.size() << "\n" << std::flush;
			}
			if(!((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1)))
			{
				std::cerr << "!((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1)), l = " << l << ", gS.at(l+1).size() = " << gS.at(l+1).size() << ", gS.size()= " << gS.size() << "\n" << std::flush;
			}
			assert((l == 0) || (gS.at(l-1).size() != 1));
			assert((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1));
		}
		

		unsigned int requiredSpace = gS.at(l).at(0).size();

		totalLevels += requiredSpace;
	}
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "B" << "\n" << std::flush; // todo remove
	
	std::vector<std::vector<std::string> > uncompressed_chromotypes;
	std::vector<int> uncompressed_chromotypes_nGaps;
	
	uncompressed_chromotypes.reserve(totalLevels);
	uncompressed_chromotypes_nGaps.reserve(totalLevels);
	
	for(unsigned int outerI = 0; outerI < gS.size(); outerI++)
	{
		std::vector<std::string>& segment = gS.at(outerI);
		assert((segment.size() == 1) || (segment.size() == 2));
		if(segment.size() == 2)
		{
			assert(segment.at(0).size() == segment.at(1).size());
		}

		unsigned int requiredSpace = segment.at(0).size();

		for(unsigned int posInSegment = 0; posInSegment < requiredSpace; posInSegment++)
		{
			std::vector<std::string> posVector;
			std::vector<double> posCoverageVector;
			std::vector<double> posSupportVector;
			std::vector<double> posInsertionsVector;

			if(segment.size() == 2)
			{
				posVector.push_back(segment.at(0).substr(posInSegment, 1));
				posVector.push_back(segment.at(1).substr(posInSegment, 1));
				posCoverageVector.push_back(0);
				posCoverageVector.push_back(0);
				posSupportVector.push_back(0);
				posSupportVector.push_back(0);
				posInsertionsVector.push_back(0);
				posInsertionsVector.push_back(0);
			}
			else
			{
				posVector.push_back(segment.at(0).substr(posInSegment, 1));
				posCoverageVector.push_back(0);
				posSupportVector.push_back(0);
				posInsertionsVector.push_back(0);
			}

			uncompressed_chromotypes.push_back(posVector);

			int numGaps = 0;
			if(segment.size() == 2)
			{
				if(segment.at(0).substr(posInSegment, 1) == "_")
				{
					numGaps++;
				}
				if(segment.at(1).substr(posInSegment, 1) == "_")
				{
					numGaps++;
				}
			}
			else
			{
				if(segment.at(0).substr(posInSegment, 1) == "_")
				{
					numGaps = 2;
				}
			}
			uncompressed_chromotypes_nGaps.push_back(numGaps);
		}
	}
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_nGaps.size());
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "C" << "\n" << std::flush; // todo remove
	
	std::vector<int> uncompressed_chromotypes_referencePositions;
	uncompressed_chromotypes_referencePositions.reserve(totalLevels);
	
	for(unsigned int outerI = 0; outerI < chromotypes_referencePositions.size(); outerI++)
	{
		for(unsigned int posInSegment = 0; posInSegment < chromotypes_referencePositions.at(outerI).size(); posInSegment++)
		{
			uncompressed_chromotypes_referencePositions.push_back(chromotypes_referencePositions.at(outerI).at(posInSegment));
		}
	}
	
	if(!(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size()))
	{
		std::cerr << "! (uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size())" << "\n";
		std::cerr << "uncompressed_chromotypes.size()" << ": " << uncompressed_chromotypes.size() << "\n";
		std::cerr << "uncompressed_chromotypes_referencePositions.size()" << ": " << uncompressed_chromotypes_referencePositions.size() << "\n\n" << std::flush;
		std::cerr << "gS.size()" << ": " << gS.size() << "\n";
		std::cerr << "chromotypes_referencePositions.size()" << ": " << chromotypes_referencePositions.size() << "\n\n";

		std::cerr << std::flush;
	}
	
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size());
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "D" << "\n" << std::flush; // todo remove
	
	std::vector<int> uncompressed_chromotypes_diploid;
	std::vector<int> uncompressed_chromotypes_losePhasing;	
	uncompressed_chromotypes_diploid.reserve(totalLevels);
	uncompressed_chromotypes_losePhasing.reserve(totalLevels);
	
	for(unsigned int outerI = 0; outerI < gS_unresolved.size(); outerI++)
	{
		std::vector<std::string>& segment = gS_unresolved.at(outerI);
		assert((segment.size() == 1) || (segment.size() == 2));
		if(segment.size() == 2)
		{
			assert(segment.at(0).size() == segment.at(1).size());
		}

		for(unsigned int posInString = 0; posInString < segment.at(0).length(); posInString++)
		{
			uncompressed_chromotypes_diploid.push_back(segment.size() - 1);
			if((posInString == 0) && (segment.size() == 2))
			{
				uncompressed_chromotypes_losePhasing.push_back(1);
			}
			else
			{
				uncompressed_chromotypes_losePhasing.push_back(0);
			}
		}
	}



	if(!(uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size()))
	{
		std::cerr << "! uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size()" << "\n";
		std::cerr << "uncompressed_chromotypes.size()" << ": " << uncompressed_chromotypes.size() << "\n";
		std::cerr << "uncompressed_chromotypes_diploid.size()" << ": " << uncompressed_chromotypes_diploid.size() << "\n\n" << std::flush;
		std::cerr << "gS.size()" << ": " << gS.size() << "\n";
		std::cerr << "uncompressed_chromotypes_referencePositions.size()" << ": " << uncompressed_chromotypes_referencePositions.size() << "\n\n";

		std::cerr << std::flush;
	}	

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "E" << "\n" << std::flush; // todo remove
	
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size());
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_losePhasing.size());

	std::string s1; std::string s2;
	s1.reserve(totalLevels);
	s2.reserve(totalLevels);
	for(unsigned int l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			s1.append(gS.at(l).at(0));
			s2.append(gS.at(l).at(0));

			unsigned int nonGap_characters = 0;
			for(unsigned int i = 0; i < gS.at(l).at(0).size(); i++)
			{
				if(gS.at(l).at(0).at(i) != '_')
				{
					nonGap_characters++;
				}
			}

			if((l != 0) && (l != (gS.size() - 1)))
			{
				assert(nonGap_characters >= k);
			}
		}
		else
		{
			assert(gS.at(l).size() == 2);
			assert((l == 0) || (gS.at(l-1).size() == 1));
			assert((l == (gS.size() - 1)) || (gS.at(l+1).size() == 1));

			s1.append(gS.at(l).at(0));
			s2.append(gS.at(l).at(1));
		}
	}
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "F" << "\n" << std::flush; // todo remove
	

	// now count k-Mers
	size_t total_characters = 0;
	size_t total_nonGap_characters = 0;

	std::map<std::string, int> kMer_counts;
	size_t kMersCannotEvaluate = 0;	

	std::vector<int> inducedkMers;
	std::vector<int> inducedkMers_OK;
	std::vector<int> inducedkMers_invalid;
	std::vector<int> inducedkMers_inReference;

	std::string filenameForSpatialSummary = pathForSpatialSummary + "/spatialSummary_" + nameForSummary + ".txt";
	ofstream spatialSummaryStream;
	spatialSummaryStream.open(filenameForSpatialSummary.c_str());
	assert(spatialSummaryStream.is_open());

	spatialSummaryStream <<
			"Level" << "\t" <<
			"ReferenceCoordinate" << "\t" << 
			"kMersInduced" << "\t" <<
			"kMersInvalid" << "\t" <<
			"kMersOK" << "\t" <<
			"InducedkMersInReference" << "\t" <<
			"DiploidChromotype" << "\t" <<
			"ChromotypeLostPhase" << "\t" <<
			"GapsAtLevel" << "\n";

	auto partitionStringIntokMers_overGaps = [](std::string s, int kMerSize) -> std::vector<std::string> {
		std::vector<std::string> forReturn;

		forReturn.reserve((int)s.length() - kMerSize);
		
		for(int i = 0; i <= ((int)s.length() - kMerSize); i++)
		{
			std::string C = s.substr(i, 1);
			if(C == "_")
			{
				forReturn.push_back(std::string(""));
			}
			else
			{
				int foundMoreNonGaps = 0;
				int currentStop = i;
				while((currentStop < ((int)s.length()-1)) && (foundMoreNonGaps < (kMerSize-1)))
				{
					currentStop++;
					std::string justAddedC = s.substr(currentStop,1);
					if(justAddedC != "_")
					{
						foundMoreNonGaps++;
					}
				}

				if(foundMoreNonGaps == (kMerSize-1))
				{
					std::string kMer = C + s.substr(i+1, (currentStop-(i+1)+1));
					kMer.erase(std::remove_if(kMer.begin(),kMer.end(), [&](char c){return ((c == '_') ? true : false);}), kMer.end());
					assert((int)kMer.length() == kMerSize);

					forReturn.push_back(kMer);
				}
				else
				{
					forReturn.push_back(std::string(""));
				}

			}
		}

		assert(forReturn.size() == (s.length() - kMerSize + 1));

		return forReturn;
	};

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "G" << "\n" << std::flush; // todo remove
	
	std::string testString = "AHCCHDHRHJSKSJSKKSKLSSOMCSNCSDKDKDKDKDKCMCMCKWL";
	std::vector<std::string> kMers_1 = partitionStringIntokMers(testString, 5);
	std::vector<std::string> kMers_2 = partitionStringIntokMers_overGaps(testString, 5);
	assert(kMers_1.size() == kMers_2.size());
	for(unsigned int i = 0; i < kMers_1.size(); i++)
	{
		assert(kMers_1.at(i) == kMers_2.at(i));
	}

	auto processString = [&](std::string s) {
		for(size_t i = 0; i < s.length(); i++)
		{
			total_characters++;
			if(s.at(i) != '_')
			{
				total_nonGap_characters++;
			}
		}
		//s.erase(std::remove_if(s.begin(),s.end(), [&](char c){return ((c == '_') ? true : false);}), s.end());

		std::vector<std::string> kMers = partitionStringIntokMers_overGaps(s, k);

		for(size_t i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);
			
			bool containsGap = (kMer.length() == 0);
			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMer.length(); kI++)
			{
				char kMerC = kMer.at(kI);
				assert(kMerC != '_');
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}
		
			if(containsGap)
			{
				// nothing
			}
			else
			{
				inducedkMers.at(i)++;

				if(kMerOK)
				{
					if(kMer_counts.count(kMer) == 0)
						kMer_counts[kMer] = 0;
					kMer_counts[kMer]++;

					if(graph->kMerinGraph(kMer))
					{
						inducedkMers_OK.at(i)++;
					}
				}
				else
				{
					kMersCannotEvaluate++;
					inducedkMers_invalid.at(i)++;
				}

				if(kMers_reference.size() == 0)
				{

				}
				else
				{
					if(kMers_reference.count(kMer))
					{
						inducedkMers_inReference.at(i)++;
					}
				}
			}
		}
	};
	
	assert(s1.length() > k);
	inducedkMers.resize(s1.length()-k+1,0);
	inducedkMers_invalid.resize(s1.length()-k+1,0);
	inducedkMers_OK.resize(s1.length()-k+1,0);
	inducedkMers_inReference.resize(s1.length()-k+1,0);
	assert(s1.length() == s2.length());
	processString(s1);
	processString(s2);

	for(unsigned int level = 0; level < inducedkMers.size(); level++)
	{
		spatialSummaryStream <<
				level << "\t" <<
				uncompressed_chromotypes_referencePositions.at(level) << "\t" <<
				inducedkMers.at(level) << "\t" <<
				inducedkMers_invalid.at(level) << "\t" <<
				inducedkMers_OK.at(level) << "\t" <<
				inducedkMers_inReference.at(level) << "\t" <<
				uncompressed_chromotypes_diploid.at(level) << "\t" <<
				uncompressed_chromotypes_losePhasing.at(level) << "\t" <<
				uncompressed_chromotypes_nGaps.at(level) << "\n";
	}

	spatialSummaryStream.close();



	size_t total_kMers_present = 0;
	size_t maximumkMerCount = 0;
	size_t maximumAllowedkMerCount = 100;
	
	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "H" << "\n" << std::flush; // todo remove
	
	for(std::map<std::string, int>::iterator kMerIt = kMer_counts.begin(); kMerIt != kMer_counts.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		int kMerCount = kMerIt->second;

		if(kMers_in_dGS != 0)
		{
			kMers_in_dGS->insert(kMer);
		}

		if(graph->kMerinGraph(kMer))
		{
			total_kMers_present++;
			if(kMers_in_dGS_in_sample != 0)
			{
				kMers_in_dGS_in_sample->insert(kMer);
			}
		}

		if(kMerCount > maximumkMerCount)
			maximumkMerCount = kMerCount;		
	}

	double haploidGenomeSize = 3.2e9;
	size_t graph_totalCoverage = graph->totalCoverage();
	double averageCoverage = double(graph_totalCoverage)/(2.0 * haploidGenomeSize);
	std::cout << "Estimate coverage posterior probabilties for max count of " << maximumkMerCount << ", assuming WG average coverage of " << averageCoverage << " [total coverage on graph: " <<  graph_totalCoverage << "]\n";

	if(maximumkMerCount > maximumAllowedkMerCount)
	{	
		maximumkMerCount = maximumAllowedkMerCount;
		std::cout << "Restrict maximum kMer count (underlying) to " << maximumAllowedkMerCount << " - all kMers with higher underlying numbers will be ignored!" << "\n" << std::flush;
	}
	
	if(averageCoverage == 0)
	{
		averageCoverage = 1;
	}
	assert(averageCoverage != 0);
	std::map<size_t, std::vector<double> > _pp_cache;
	auto getPPs = [&](size_t observedCoverage) -> std::vector<double> {
		if(_pp_cache.count(observedCoverage) > 0)
		{
			return _pp_cache.at(observedCoverage);
		}
		else
		{
			std::vector<double> forReturn = kMer_PP(observedCoverage, maximumkMerCount, averageCoverage);
			if(observedCoverage < (2 * averageCoverage))
			{
				_pp_cache[observedCoverage] = forReturn;
			}
			return forReturn;
		}
	};

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "I" << "\n" << std::flush; // todo remove
	
	double expected_kMer_presence = 0;
	double expected_kMer_presence_sum = 0;
	for(std::map<std::string, int>::iterator kMerIt = kMer_counts.begin(); kMerIt != kMer_counts.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		int kMerCount = kMerIt->second;

		int kMerCoverage = graph->kMer_getCoverage(kMer);
		std::vector<double> pps;
		if(kMerCount <= maximumAllowedkMerCount)
		{
			pps = getPPs(kMerCoverage);
			
			double expected_missing = 0;
			for(int i = 0; i < kMerCount; i++)
			{
				int e_missing = kMerCount - i;
				double e_p = pps.at(i);
				expected_missing += e_missing * e_p;
			}
			
			assert(expected_missing >= 0);
			assert(expected_missing <= kMerCount);
			
			double this_kMer_presence = (kMerCount - expected_missing);
			expected_kMer_presence += this_kMer_presence;
			expected_kMer_presence_sum += kMerCount;

			if(kMers_in_dGS_optimality != 0)
			{
				(*kMers_in_dGS_optimality)[kMer] = (double)this_kMer_presence / (double)kMerCount;
			}			
		}
	}

	assert(expected_kMer_presence >= 0);
	assert(expected_kMer_presence <= expected_kMer_presence_sum);
	assert(expected_kMer_presence_sum > 0);
	
	double unweighted_optimality = (kMer_counts.size() != 0) ? (double(total_kMers_present)/double(kMer_counts.size())) : 1;
	double weighted_probabilistic_optimality = expected_kMer_presence / expected_kMer_presence_sum;

	std::cout << "evaluate_dGS(..):\n";
	std::cout << "\tLength of gS: " << gS.size() << "\n";	
	std::cout << "\tTotal characters: " << total_characters << ", of which non-gap: " << total_nonGap_characters << "\n";
	std::cout << "\tTotal kMers: " << kMer_counts.size() << "\n";
	std::cout << "\tTotal invalid kMers (driven by non-A/C/G/T): " << kMersCannotEvaluate << "\n";
	
	std::cout << "\tTotal kMers present: " << total_kMers_present << "\n";
	std::cout << "\tUnweighted optimality: " << unweighted_optimality << "\n";
	std::cout << "\tWeighted probabilistic optimality: " << weighted_probabilistic_optimality << "\n";

	std::vector<std::string> fieldsForPrinting;

	fieldsForPrinting.push_back(nameForSummary);
	fieldsForPrinting.push_back(Utilities::ItoStr(total_characters));
	fieldsForPrinting.push_back(Utilities::ItoStr(total_nonGap_characters));
	fieldsForPrinting.push_back(Utilities::ItoStr(kMer_counts.size()));
	fieldsForPrinting.push_back(Utilities::ItoStr(kMersCannotEvaluate));
	fieldsForPrinting.push_back(Utilities::ItoStr(total_kMers_present));
	fieldsForPrinting.push_back(Utilities::DtoStr(unweighted_optimality));
	fieldsForPrinting.push_back(Utilities::DtoStr(weighted_probabilistic_optimality));

	summaryFileStream << Utilities::join(fieldsForPrinting, "\t") << "\n";
}

template<int m, int k, int colours>
std::pair<diploidGenomeString, diploidGenomeString> greedilyResolveDiploidKMerString(diploidGenomeString& original_gS, DeBruijnGraph<m, k, colours>* graph)
{
	int threads = 40;
	omp_set_num_threads(threads);
	unsigned int _compartmentI;
	
	size_t gS_length_before = original_gS.size();
	
	std::vector<std::pair<std::string, std::string> > openPairs;

	std::pair<std::string, std::string> resolved;
	std::pair<std::string, std::string> resolved_noGaps;
	
	diploidGenomeString resolved_gS;

	int addedHomozygousSeq = 0;

	auto evaluateSequence = [&](std::string seq) -> double {
		seq.erase(std::remove_if(seq.begin(),seq.end(), [&](char c){return ((c == '_') ? true : false);}), seq.end());
		std::vector<std::string> kMers = partitionStringIntokMers(seq, k);
		std::map<std::string, int> kMers_table;
		for(unsigned int i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);
			if(kMers_table.count(kMer) == 0)
			{
				kMers_table[kMer] = 0;
			}
			kMers_table[kMer]++;
		}

		double optimality_sum = 0;
		double optimality_score = 0;

		for(std::map<std::string, int>::iterator kMerIt = kMers_table.begin(); kMerIt != kMers_table.end(); kMerIt++)
		{
			std::string kMerSeq = kMerIt->first;
			
			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMerSeq.length(); kI++)
			{
				char kMerC = kMerSeq.at(kI);
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}
			
			if(kMerOK)
			{
				int count = kMerIt->second;
				
				bool present = graph->kMerinGraph(kMerSeq);
				
				if(present)
					optimality_score += count;

				optimality_sum += count;
			}
		}

		if(optimality_sum == 0)
		{
			return 1;
		}
		else
		{
			return optimality_score/optimality_sum;
		}
	};

	auto evaluateSequence2 = [&](std::string seq, int& optimality_sum, int& optimality_score) -> void {
		seq.erase(std::remove_if(seq.begin(),seq.end(), [&](char c){return ((c == '_') ? true : false);}), seq.end());
		std::vector<std::string> kMers = partitionStringIntokMers(seq, k);
		std::map<std::string, int> kMers_table;
		for(unsigned int i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);
			if(kMers_table.count(kMer) == 0)
			{
				kMers_table[kMer] = 0;
			}
			kMers_table[kMer]++;
		}

		optimality_sum = 0;
		optimality_score = 0;

		for(std::map<std::string, int>::iterator kMerIt = kMers_table.begin(); kMerIt != kMers_table.end(); kMerIt++)
		{
			std::string kMerSeq = kMerIt->first;
			
			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMerSeq.length(); kI++)
			{
				char kMerC = kMerSeq.at(kI);
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}
			
			if(kMerOK)
			{
				// int count = kMerIt->second;
				
				bool present = graph->kMerinGraph(kMerSeq);
				
				if(present)
					optimality_score += 1;

				optimality_sum += 1;
			}
		}
	};

	
	auto evaluateAndReduce = [&]() {
		assert(resolved.first.length() == resolved.second.length());

		size_t size_before = openPairs.size();
		if(openPairs.size() > 0)
		{
			std::vector<double> optimality;
			optimality.resize(openPairs.size());

			long long max_i = optimality.size() - 1;
			long long chunk_size = max_i / threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == threads);
				long long thisThread = omp_get_thread_num();
				long long firstPair = thisThread * chunk_size;
				long long lastPair = (thisThread+1) * chunk_size - 1;
				if((thisThread == (threads-1)) && (lastPair < max_i))
				{
					lastPair = max_i;
				}

				for(long long i = firstPair; i <= lastPair; i++)
				{
					int copyCharactersFromPrevious = k;
					if((copyCharactersFromPrevious > resolved_noGaps.first.length()) || (copyCharactersFromPrevious > resolved_noGaps.second.length()))
					{
						copyCharactersFromPrevious = (resolved_noGaps.first.length() < resolved_noGaps.second.length()) ? resolved_noGaps.first.length() : resolved_noGaps.second.length() ;
					}

					std::string previous_first = resolved_noGaps.first.substr(resolved_noGaps.first.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);
					std::string previous_second = resolved_noGaps.second.substr(resolved_noGaps.second.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);

					assert(previous_first.length() == copyCharactersFromPrevious);
					assert(previous_second.length() == copyCharactersFromPrevious);
					assert(previous_first == previous_second);

					std::string seq_evaluate_1 = previous_first + openPairs.at(i).first;
					std::string seq_evaluate_2 = previous_second + openPairs.at(i).second;

					int optimality_sum_seq1 = 0;
					int optimality_score_seq1 = 0;
					int optimality_sum_seq2 = 0;
					int optimality_score_seq2 = 0;					
					evaluateSequence2(seq_evaluate_1, optimality_sum_seq1, optimality_score_seq1);
					evaluateSequence2(seq_evaluate_2, optimality_sum_seq2, optimality_score_seq2);

					double O = 1;
					if((optimality_sum_seq1 + optimality_sum_seq2) > 0)
					{
						O = (double)(optimality_score_seq1 + optimality_score_seq2)/(double)(optimality_sum_seq1 + optimality_sum_seq2);
					}
					
					optimality.at(i) = O;
				}
			}

			double maxOptimality; size_t whereOptimality;
			for(size_t i = 0; i < optimality.size(); i++)
			{
				if((i == 0) || (maxOptimality < optimality.at(i)))
				{
					maxOptimality = optimality.at(i);
					whereOptimality = i;
				}
			}
			
			if(maxOptimality == 0)
			{
				std::vector<std::pair<std::string, std::string> > new_openPairs;			
				new_openPairs.push_back(openPairs.at(0));
				openPairs = new_openPairs;				
			}
			else
			{
				std::vector<std::pair<std::string, std::string> > new_openPairs;
				if(optimality.size() <= 1000)
				{
					for(size_t i = 0; i < optimality.size(); i++)
					{
						if((optimality.at(i)/maxOptimality) > 0.9)
						{
							new_openPairs.push_back(openPairs.at(i));
						}
					}
				}
				else
				{
					std::vector<unsigned int> optimality_sort_indices;
					optimality_sort_indices.resize(optimality.size());
					for(unsigned int i = 0; i < optimality_sort_indices.size(); i++)
					{
						optimality_sort_indices.at(i) = i;
					}
					std::sort(optimality_sort_indices.begin(), optimality_sort_indices.end(), [&](unsigned int a, unsigned int b) -> bool {
							return (optimality.at(b) < optimality.at(a));
						}
					);

					assert(optimality_sort_indices.size() == optimality.size());
					unsigned int maxIndex = 1000;
					for(unsigned int optimalityI = 0; optimalityI <= maxIndex; optimalityI++)
					{
						unsigned int realOptimalityI = optimality_sort_indices.at(optimalityI);
						double realOptimality = optimality.at(realOptimalityI);
						if(optimalityI > 0)
						{
							unsigned int previous_OptimalityI = optimality_sort_indices.at(optimalityI - 1);
							double previous_Optimality = optimality.at(previous_OptimalityI);
							assert(previous_Optimality >= realOptimality);
						}

						new_openPairs.push_back(openPairs.at(realOptimalityI));
					}
				}
				openPairs = new_openPairs;				
			}
			
			size_t size_after = openPairs.size();
			if(((double)size_after/(double)size_before) > 0.5)
			{
				std::cerr << "\nevaluateAndReduce(..): Not very effective. Reduction from " << size_before << " to " << size_after << ", maxOptimality = " << maxOptimality << " -- but does this really matter??\n" << std::flush;
//				std::cerr << "\t" << "resolved.first.length(): " << resolved.first.length() << "\n";
//				std::cerr << "\t" << "openPairs.at(0).first.length(): " << openPairs.at(0).first.length() << " // " << openPairs.at(0).first << "\n" << std::flush;
//				assert(size_after > 100);
//				std::vector<std::pair<std::string, std::string> > new_openPairs;
//				for(size_t i = 0; i < openPairs.size(); i++)
//				{
//					if(Utilities::randomDouble() < 0.1)
//					{
//						new_openPairs.push_back(openPairs.at(i));
//					}
//				}
//				openPairs = new_openPairs;
//				std::cerr << "Random reduction to " << openPairs.size() << " elements. " << "\n" << std::flush;
			}			
		}
	};
	
	auto evaluateAndResolve = [&]() {
		assert(resolved.first.length() == resolved.second.length());

		if(openPairs.size() > 0)
		{
			std::vector<double> optimality;
			optimality.resize(openPairs.size());

			long long max_i = optimality.size() - 1;
			long long chunk_size = max_i / threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == threads);
				long long thisThread = omp_get_thread_num();
				long long firstPair = thisThread * chunk_size;
				long long lastPair = (thisThread+1) * chunk_size - 1;
				if((thisThread == (threads-1)) && (lastPair < max_i))
				{
					lastPair = max_i;
				}

				for(long long i = firstPair; i <= lastPair; i++)
				{
					int copyCharactersFromPrevious = k;
					if((copyCharactersFromPrevious > resolved_noGaps.first.length()) || (copyCharactersFromPrevious > resolved_noGaps.second.length()))
					{
						copyCharactersFromPrevious = (resolved_noGaps.first.length() < resolved_noGaps.second.length()) ? resolved_noGaps.first.length() : resolved_noGaps.second.length() ;
					}

					std::string previous_first = resolved_noGaps.first.substr(resolved_noGaps.first.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);
					std::string previous_second = resolved_noGaps.second.substr(resolved_noGaps.second.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);

					assert(previous_first.length() == copyCharactersFromPrevious);
					assert(previous_second.length() == copyCharactersFromPrevious);
					assert(previous_first == previous_second);

					std::string seq_evaluate_1 = previous_first + openPairs.at(i).first;
					std::string seq_evaluate_2 = previous_second + openPairs.at(i).second;

					double O = evaluateSequence(seq_evaluate_1) + evaluateSequence(seq_evaluate_2);

					optimality.at(i) = O;
				}
			}

			double maxOptimality; size_t whereOptimality;
			for(size_t i = 0; i < optimality.size(); i++)
			{
				if((i == 0) || (maxOptimality < optimality.at(i)))
				{
					maxOptimality = optimality.at(i);
					whereOptimality = i;
				}
			}

			resolved.first.append(openPairs.at(whereOptimality).first);
			resolved.second.append(openPairs.at(whereOptimality).second);

			
			std::string chosen_s1 = openPairs.at(whereOptimality).first;
			std::string chosen_s2 = openPairs.at(whereOptimality).second;
			
			std::string chosen_s1_noGaps = chosen_s1;
			std::string chosen_s2_noGaps = chosen_s2;
			chosen_s1_noGaps.erase(std::remove_if(chosen_s1_noGaps.begin(),chosen_s1_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), chosen_s1_noGaps.end());
			chosen_s2_noGaps.erase(std::remove_if(chosen_s2_noGaps.begin(),chosen_s2_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), chosen_s2_noGaps.end());
			resolved_noGaps.first.append(chosen_s1_noGaps);
			resolved_noGaps.second.append(chosen_s2_noGaps);
			
			
			std::string chosen_hom;
			
			size_t last_character_s1 = chosen_s1.size();			
			size_t last_character_s2 = chosen_s2.size();
			size_t have_hom_s1 = 0;
			size_t have_hom_s2 = 0;
			
			if(addedHomozygousSeq > 0)
			{
				do{
					last_character_s1--;
					if(chosen_s1.at(last_character_s1) != '_')
					{
						have_hom_s1++;
					}
				} while (have_hom_s1 != addedHomozygousSeq);
				
				do{
					last_character_s2--;
					if(chosen_s2.at(last_character_s2) != '_')
					{
						have_hom_s2++;
					}
				} while (have_hom_s2 != addedHomozygousSeq);
				
				std::string suffix_s1 = chosen_s1.substr(last_character_s1);
				std::string suffix_s2 = chosen_s2.substr(last_character_s2);
				
				assert(suffix_s1 == suffix_s2);
				std::string original_chosen_1 = chosen_s1;
				chosen_hom = suffix_s1;
				
				chosen_s1 = chosen_s1.substr(0, last_character_s1);
				chosen_s2 = chosen_s2.substr(0, last_character_s2);
				assert((chosen_s1 + chosen_hom) == original_chosen_1);				
				assert(chosen_hom.length() >= addedHomozygousSeq);
			}
			
			std::vector<std::string> newCompartment;
			newCompartment.push_back(chosen_s1);
			newCompartment.push_back(chosen_s2);
			resolved_gS.push_back(newCompartment);
			
			if(chosen_hom.size() > 0)
			{
				std::vector<std::string> newHomCompartment;
				newHomCompartment.push_back(chosen_hom);
				resolved_gS.push_back(newHomCompartment);		
			}

			openPairs.clear();
		}
		addedHomozygousSeq = 0;
	};

	auto addOpenPair = [&](std::pair<std::string, std::string> p) {
		assert(p.first.length() == p.second.length());

		if(openPairs.size() == 0)
		{
			openPairs.push_back(p);
		}
		else
		{
			std::vector<std::pair<std::string, std::string> > new_openPairs;

			for(unsigned int existingI = 0; existingI < openPairs.size(); existingI++)
			{
				std::pair<std::string, std::string> newpair_1 = openPairs.at(existingI);
				std::pair<std::string, std::string> newpair_2 = openPairs.at(existingI);

				newpair_1.first.append(p.first);
				newpair_1.second.append(p.second);

				newpair_2.first.append(p.second);
				newpair_2.second.append(p.first);

				new_openPairs.push_back(newpair_1);
				new_openPairs.push_back(newpair_2);
			}

			openPairs = new_openPairs;
			
			if(openPairs.size() > 100000)
			{
				evaluateAndReduce();
			}
			
			assert(openPairs.size() < 10000000);
		}

		addedHomozygousSeq = 0;
	};

	auto addOpenSeq = [&](std::string seq) {
		for(unsigned int existingI = 0; existingI < openPairs.size(); existingI++)
		{
			openPairs.at(existingI).first.append(seq);
			openPairs.at(existingI).second.append(seq);
		}

		for(int i = 0; i < seq.length(); i++)
		{
			if(seq.at(i) != '_')
				addedHomozygousSeq++;
		}

		if(addedHomozygousSeq >= k)
		{
			evaluateAndResolve();
		}
	};

	size_t expectedResolvedLength = 0;
	std::cout << "\n";
	for(unsigned int compartmentI = 0; compartmentI < original_gS.size(); compartmentI++)
	{
		_compartmentI = compartmentI;
		
		if((compartmentI % 100) == 0)
		{
			std::cout << "\r" << "compartmentI " << compartmentI << " / " << original_gS.size() << "   openPairs.size(): " << openPairs.size() << "     " << std::flush;
		}

		assert(original_gS.at(compartmentI).size() >= 1);
		assert(original_gS.at(compartmentI).size() <= 2);

		if(original_gS.at(compartmentI).size() == 2)
		{
			assert(original_gS.at(compartmentI).at(0).size() == original_gS.at(compartmentI).at(1).size());
		}

		unsigned int compartmentLength = original_gS.at(compartmentI).at(0).size();
		bool diploid = (original_gS.at(compartmentI).size() == 2);

		if(diploid)
		{
			std::pair<std::string, std::string> oP;
			oP.first = original_gS.at(compartmentI).at(0);
			oP.second = original_gS.at(compartmentI).at(1);

			addOpenPair(oP);
		}
		else
		{
			std::string homozygousSeq = original_gS.at(compartmentI).at(0);
			if(openPairs.size() > 0)
			{
				addOpenSeq(homozygousSeq);
			}
			else
			{
				if(resolved_gS.size() == 0)
				{
					std::vector<std::string> firstCompartment;
					firstCompartment.push_back(homozygousSeq);
					resolved_gS.push_back(firstCompartment);
				}
				else
				{
					if(resolved_gS.at(resolved_gS.size() - 1).size() == 1)
					{
						resolved_gS.at(resolved_gS.size() - 1).at(0).append(homozygousSeq);
					}
					else
					{
						std::vector<std::string> newCompartment;
						newCompartment.push_back(homozygousSeq);
						resolved_gS.push_back(newCompartment);				
					}
				}
		
				resolved.first.append(homozygousSeq);
				resolved.second.append(homozygousSeq);
				
				std::string homozygousSeq_noGaps = homozygousSeq;
				homozygousSeq_noGaps.erase(std::remove_if(homozygousSeq_noGaps.begin(),homozygousSeq_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), homozygousSeq_noGaps.end());
				
				resolved_noGaps.first.append(homozygousSeq_noGaps);
				resolved_noGaps.second.append(homozygousSeq_noGaps);
			}
		}

		expectedResolvedLength += compartmentLength;
		
		// std::cout << "greedilyResolveDiploidKMerString(..): Compartment " << compartmentI << " of " << original_gS.size() << "\n";
		// std::cout << "\topenPairs.size(): " << openPairs.size() << "\n";
		// std::cout << "\taddedHomozygousSeq: " << addedHomozygousSeq << "\n";
		// std::cout << "\texpectedResolvedLength: " << expectedResolvedLength << "\n";
		// std::cout << std::flush;
	}

	std::cout << "\n";

	evaluateAndResolve();
	assert(resolved.first.length() == expectedResolvedLength);
	assert(resolved.second.length() == expectedResolvedLength);

	std::pair<diploidGenomeString, diploidGenomeString> forReturn;
	std::vector<string> forReturn_oneCompartment;
	forReturn_oneCompartment.push_back(resolved.first);
	forReturn_oneCompartment.push_back(resolved.second);
	forReturn.first.push_back(forReturn_oneCompartment);
	forReturn.second = resolved_gS;

	std::cout << "\tgreedilyResolveDiploidKMerString(..): length of gS before: " << gS_length_before << ", and after: " << resolved_gS.size() << "\n" << std::flush;

	
	return forReturn;
}

void validateGenomeString(diploidGenomeString& gS, std::string deBruijnGraph_fileName, int kMer_size)
{
	assert(kMer_size == 31);

	DeBruijnGraph<2, 55, 1> myGraph(25, 100);
	std::cout << "Graph allocated, loading binary...\n" << std::flush;
	myGraph.loadMultiColourBinary("/ddn/projects3/mcvean_res/alex/CortexGraphs/AA02O9Q_A1_55_cleaned_2.ctx");
}

Graph* VCF2Graph(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::string referenceGenomePath)
{
	std::vector<std::vector<int> > VCF_chromotype_referencePositions;
	diploidGenomeString gS = VCF2GenomeString(chromosomeID, positionStart, positionStop, VCFpath, referenceGenomePath, VCF_chromotype_referencePositions);
	return genomeString2Graph(gS);
}

Graph* genomeString2Graph(diploidGenomeString gS, bool verbose)
{
	std::vector< std::vector<std::string> >& chars_2_graph = gS;

	Graph* forReturn = new Graph();

	Node* n0 = new Node();
	n0->level = 0;
	n0->terminal = false;
	forReturn->registerNode(n0, 0);

//	std::cout << "A\n" << std::flush;

	int runningGraphLevel = 0;
	for(unsigned int lI = 0; lI < chars_2_graph.size(); lI++)
	{
		if(verbose)
		{
			if((lI % 1000) == 0)
			{
				std::cerr << "\r" << lI << " / " << chars_2_graph.size();
			}
		}
//		std::cout << "B " << lI << "\n" << std::flush;

		if(!((forReturn->NodesPerLevel.at(runningGraphLevel).size() == 1)))
		{
			std::cout << "!(forReturn->NodesPerLevel.at(lI).size() == 1)\n";
			std::cout << "lI: " << lI << "\n";
			std::cout << "chars_2_graph.size(): " << chars_2_graph.size() << "\n";
			std::cout << "forReturn->NodesPerLevel.at(lI).size(): " << forReturn->NodesPerLevel.at(lI).size() << "\n" << std::flush;
		}
		assert(forReturn->NodesPerLevel.at(runningGraphLevel).size() == 1);
		Node* previousNode = *(forReturn->NodesPerLevel.at(runningGraphLevel).begin());
		if(chars_2_graph.at(lI).size() == 1)
		{
//			std::cout << "C " << lI << "\n" << std::flush;

			Node* local_previousNode = previousNode;
			for(unsigned int cI = 0; cI < chars_2_graph.at(lI).at(0).size(); cI++)
			{
				string c = chars_2_graph.at(lI).at(0).substr(cI, 1);
				Node* n = new Node();
				int level = runningGraphLevel + cI + 1;
				n->level = level;
				n->terminal = false;
				forReturn->registerNode(n, level);

				Edge* newE = new Edge();

				std::string levelID = "L"+Utilities::ItoStr(level);

				newE->count = 1;
				newE->emission = forReturn->CODE.doCode(levelID, c);
				newE->locus_id = levelID;

				forReturn->registerEdge(newE);
				newE->From = local_previousNode;
				local_previousNode->Outgoing_Edges.insert(newE);

				newE->To = n;
				n->Incoming_Edges.insert(newE);

				local_previousNode = n;
			}
		}
		else
		{

//			std::cout << "D " << lI << "\n" << std::flush;

			assert(chars_2_graph.at(lI).size() == 2);
			for(int haploI = 0; haploI < (int)chars_2_graph.at(lI).size(); haploI++)
			{
				Node* local_previousNode = previousNode;
				assert(chars_2_graph.at(lI).at(haploI).size() == chars_2_graph.at(lI).at(0).size());
				for(int cI = 0; cI < (int)chars_2_graph.at(lI).at(haploI).size(); cI++)
				{
					string c = chars_2_graph.at(lI).at(haploI).substr(cI, 1);
					Node* n;

					int level = runningGraphLevel + cI + 1;

					if((cI < (int(chars_2_graph.at(lI).at(haploI).size()) - 1)) || (haploI == 0))
					{
						n = new Node();

						n->level = level;
						n->terminal = false;
						forReturn->registerNode(n, level);
					}
					else
					{
						assert(forReturn->NodesPerLevel.at(level).size() == 1);
						n = *(forReturn->NodesPerLevel.at(level).begin());
					}

					Edge* newE = new Edge();

					std::string levelID = "L"+Utilities::ItoStr(level);

					newE->count = 1;
					newE->emission = forReturn->CODE.doCode(levelID, c);
					newE->locus_id = levelID;

					forReturn->registerEdge(newE);
					newE->From = local_previousNode;
					local_previousNode->Outgoing_Edges.insert(newE);

					newE->To = n;
					n->Incoming_Edges.insert(newE);

					local_previousNode = n;
				}
			}
		}

		runningGraphLevel += chars_2_graph.at(lI).at(0).size();
	}

	assert(forReturn->NodesPerLevel.at(runningGraphLevel).size() == 1);
	(*forReturn->NodesPerLevel.at(runningGraphLevel).begin())->terminal = true;

	forReturn->checkConsistency(true);

	if(verbose)
	{
		std::cout << "\n" << std::flush;
	}

	return forReturn;
}

void vennDiagrams(std::vector<std::string> setNames, std::vector<std::set<std::string>*> kMers, std::vector<std::set<std::string>*> kMers_present, std::vector<std::map<std::string, double>* > kMer_optimalities, std::string outputFile)
{
	ofstream outputStream;
	outputStream.open(outputFile.c_str());
	if(! outputStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ outputFile);
	}

	assert(setNames.size() == kMers.size());
	assert(kMers.size() == kMers_present.size());
	assert(kMers_present.size() == kMer_optimalities.size());

	unsigned int sets = kMers.size();
	assert(sets >= 1);

	unsigned int highest_pseudoNumber = pow(2, sets) - 1;

	std::cout << "vennDiagrams(..): Compute " << highest_pseudoNumber << " constellations.\n" << std::flush;

	std::set<std::string> kMers_combined;
	for(unsigned int i = 0; i < kMers.size(); i++)
	{
		kMers_combined.insert(kMers.at(i)->begin(), kMers.at(i)->end());
	}

	std::vector<std::string> header_fields;
	for(unsigned int sI = 0; sI < setNames.size(); sI++)
	{
		header_fields.push_back(setNames.at(sI));
	}

	header_fields.push_back("Name");
	header_fields.push_back("TOTAL");
	header_fields.push_back("N_THISFIELD");
	header_fields.push_back("Presence");
	header_fields.push_back("Optimality");

	outputStream << Utilities::join(header_fields, "\t") << "\n" << std::flush;
	for(unsigned int constellationI = 1; constellationI <= highest_pseudoNumber; constellationI++)
	{
		std::set<std::string> thisConstellation_kMers;
		std::vector<std::string> thisLine_fields;

		std::string nameForConstellation;
		std::set<std::string> presentSets;
		for(unsigned int sI = 0; sI < setNames.size(); sI++)
		{
			bool want_this_field = Utilities::extractBit(constellationI, sI);
			thisLine_fields.push_back(Utilities::ItoStr(want_this_field));
			if(want_this_field)
			{
				presentSets.insert(setNames.at(sI));
			}
		}
		if(presentSets.size() == 1)
		{
			nameForConstellation = "Only " + *(presentSets.begin());
		}
		else
		{
			std::vector<std::string> _sets(presentSets.begin(), presentSets.end());
			nameForConstellation = Utilities::join(_sets, ", ");
		}

		thisLine_fields.push_back(nameForConstellation);


		for(std::set<std::string>::iterator kMerIt = kMers_combined.begin(); kMerIt != kMers_combined.end(); kMerIt++)
		{
			std::string kMer = *kMerIt;
			bool keepkMer = true;
			for(unsigned int sI = 0; sI < setNames.size(); sI++)
			{
				bool want_this_field = Utilities::extractBit(constellationI, sI);
				if(want_this_field)
				{
					if(kMers.at(sI)->count(kMer) == 0)
					{
						keepkMer = false;
						break;
					}
				}
				else
				{
					if(kMers.at(sI)->count(kMer))
					{
						keepkMer = false;
						break;
					}
				}
			}

			if(keepkMer)
			{
				thisConstellation_kMers.insert(kMer);
			}
		}

		thisLine_fields.push_back(Utilities::ItoStr(kMers_combined.size()));
		thisLine_fields.push_back(Utilities::ItoStr(thisConstellation_kMers.size()));

		double presence_sum = 0;
		double presence_N = 0;
		double optimality_sum = 0;
		double optimality_N = 0;
		for(std::set<std::string>::iterator kMerIt = thisConstellation_kMers.begin(); kMerIt != thisConstellation_kMers.end(); kMerIt++)
		{
			std::string kMer = *kMerIt;
			for(unsigned int sI = 0; sI < setNames.size(); sI++)
			{
				bool want_this_field = Utilities::extractBit(constellationI, sI);
				if(want_this_field)
				{
					if(kMers_present.at(sI)->count(kMer))
					{
						presence_sum++;
					}
					optimality_sum += kMer_optimalities.at(sI)->at(kMer);
					break;
				}
			}
			presence_N++;
			optimality_N++;
		}

		thisLine_fields.push_back(Utilities::DtoStr(presence_sum/presence_N));
		thisLine_fields.push_back(Utilities::DtoStr(optimality_sum/optimality_N));

		outputStream << Utilities::join(thisLine_fields, "\t") << "\n" << std::flush;
	}

	outputStream.close();
}
diploidGenomeString compressGenomeString(diploidGenomeString gS)
{
	size_t gS_length_before = gS.size();
	
	diploidGenomeString forReturn;
	for(size_t l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			if((forReturn.size() == 0) || (forReturn.at(forReturn.size() - 1).size() != 1))
			{
				forReturn.push_back(gS.at(l));
			}
			else
			{
				forReturn.at(forReturn.size() - 1).at(0).append(gS.at(l).at(0));
			}
		}
		else
		{
			forReturn.push_back(gS.at(l));
		}
	}
	
	std::cout << "\tcompressGenomeString(..): length of gS before: " << gS_length_before << ", and after: " << forReturn.size() << "\n" << std::flush;
	
	return forReturn;
}

diploidGenomeString VCF2GenomeString(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::string referenceGenomePath, std::vector<std::vector<int> >& ret_graph_referencePositions, bool ignoreVCF, bool onlyPASS)
{

	map<string, string> referenceGenome = Utilities::readFASTA(referenceGenomePath);
	if(! referenceGenome.count(chromosomeID))
	{
		std::cerr << "Error: reference genome " << referenceGenomePath << " does not seem to contain chromosome " << chromosomeID << "\n";
		std::cerr << "Available chromosomes:\n";
		for(map<string, string>::iterator chrIt = referenceGenome.begin(); chrIt != referenceGenome.end(); chrIt++)
		{
			std::cerr << " - " << chrIt->first << "\n";
		}
		std::cerr << std::flush;
	}
	assert(referenceGenome.count(chromosomeID));

	ret_graph_referencePositions.clear();

	std::string referenceChromosome = referenceGenome.at(chromosomeID);

	ifstream VCFstream;

	VCFstream.open(VCFpath.c_str());

	if(! VCFstream.is_open())
	{
		throw std::runtime_error("Cannot open VCF file: "+ VCFpath);
	}

	bool ignoreBoundaries = ((positionStart == -1) && (positionStop == -1));

	if(!ignoreBoundaries)
		assert(positionStop > positionStart);

	if(ignoreBoundaries)
	{
		positionStart = 1;
		positionStop = referenceGenome.at(chromosomeID).length();
	}

	std::vector< std::vector<std::string> > chars_2_graph;

	std::string line;
	size_t lineCounter = 0;

	std::vector<std::string> header_fields;
	std::string sampleID;
	size_t fieldIndex_sample_genotypes;
	int lastExtractedPosition = 0;

	if(! ignoreVCF)
	{
		while(VCFstream.good())
		{
			std::getline(VCFstream, line);
			lineCounter++;

			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			if((line.length() > 1) && (line.substr(0, 2) == "##"))
				continue;

			if(line.substr(0, 1) == "#")
			{
				header_fields = Utilities::split(line, "\t");
				if(!(header_fields.size() >= 5))
				{
					std::cerr << "Not enough fields in header line of file " << VCFpath << "; expect at least 5, have " << header_fields.size() << "\n" << std::flush;
				}
				assert(header_fields.size() >= 5);

				for(unsigned int fI = 0; fI < header_fields.size(); fI++)
				{
					if(header_fields.at(fI) == "FORMAT")
					{
						assert(fI == (header_fields.size() - 2));
						fieldIndex_sample_genotypes = fI + 1;
						sampleID = header_fields.at(fieldIndex_sample_genotypes);
					}
				}
				assert(header_fields.at(3) == "REF");
				assert(header_fields.at(4) == "ALT");
				assert(header_fields.at(6) == "FILTER");

				continue;
			}

			assert(sampleID != "");

			std::vector<std::string> line_fields = Utilities::split(line, "\t");
			assert(line_fields.size() == header_fields.size());

			std::string chromosome = line_fields.at(0);
			int position = Utilities::StrtoI(line_fields.at(1));
			assert(position > 0);

			if((position % 1000) == 0)
				std::cerr << "\r" << "Chromosome " << chromosomeID << " position " << position << std::flush;
				
			if(chromosome != chromosomeID)
				continue;

			if(position < positionStart)
				continue;

			if(position <= lastExtractedPosition)
				continue;

			if(onlyPASS)
			{
				if(line_fields.at(6) != "PASS")
				{
					continue;
				}
			}
			if(position > positionStop)
			{
				break;
			}

			if(lastExtractedPosition != (position - 1))
			{
				if(lastExtractedPosition == 0)
				{
					for(int pos = positionStart; pos < position; pos++)
					{
						std::vector<std::string> thisPosVec;
						thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
						chars_2_graph.push_back(thisPosVec);

						std::vector<int> thisPosVec_referencePositions;
						thisPosVec_referencePositions.push_back(pos);
						ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
					}
				}
				else
				{
					for(int pos = lastExtractedPosition + 1; pos < position; pos++)
					{
						std::vector<std::string> thisPosVec;
						thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
						chars_2_graph.push_back(thisPosVec);

						std::vector<int> thisPosVec_referencePositions;
						thisPosVec_referencePositions.push_back(pos);
						ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
					}
				}
			}

			std::string reference_allele = line_fields.at(3);

			std::vector<int> thisPosVec_referencePositions;

			for(size_t pos = 0; pos < reference_allele.size(); pos++)
			{
				if(referenceChromosome.at(position + pos - 1) != reference_allele.at(pos))
				{
					std::cerr << "Reference allele error at position " << position << ", file " << VCFpath << "[in-allele position " << pos << "]\n";
					std::cerr << "\tVCF says it is " << reference_allele.at(pos) << ", but reference genome says " << referenceChromosome.at(position + pos - 1) << "\n" << std::flush;
				}
				assert(referenceChromosome.at(position + pos - 1) == reference_allele.at(pos));
			}

			std::vector<std::string> alleles;
			alleles.push_back(reference_allele);

			std::vector<std::string> alternative_alleles = Utilities::split(line_fields.at(4), ",");
			alleles.insert(alleles.end(), alternative_alleles.begin(), alternative_alleles.end());

			std::string formatString = line_fields.at(fieldIndex_sample_genotypes - 1);
			std::vector<std::string> formatString_elements = Utilities::split(formatString, ":");
			assert(formatString_elements.at(0) == "GT");

			std::string dataString = line_fields.at(fieldIndex_sample_genotypes);
			std::vector<std::string> dataString_elements = Utilities::split(dataString, ":");

			std::string genotypesString = dataString_elements.at(0);
			std::vector<std::string> genotypesString_elements = Utilities::split(genotypesString, "/");
			assert(genotypesString_elements.size() == 2);

			std::vector<std::string> called_alleles;
			unsigned int alleles_maxLength = 0;
			for(unsigned int eI = 0; eI < genotypesString_elements.size(); eI++)
			{
				int genotypeIndex = Utilities::StrtoI(genotypesString_elements.at(eI));
				std::string allele = alleles.at(genotypeIndex);
				if((eI > 0) && (allele == called_alleles.at(0)))
				{
					continue;
				}
				called_alleles.push_back(allele);				
				alleles_maxLength = (allele.length() > alleles_maxLength) ? allele.length() : alleles_maxLength;
			}
			assert(alleles_maxLength != 0);

			for(unsigned int aI = 0; aI < called_alleles.size(); aI++)
			{
				int missingGaps = alleles_maxLength - called_alleles.at(aI).length();
				if(missingGaps > 0)
				{
					for(int mI = 0; mI < missingGaps; mI++)
					{
						called_alleles.at(aI).push_back('_');
					}
				}
			}
			
			for(unsigned int aI = 0; aI < called_alleles.size(); aI++)
			{
				assert(called_alleles.at(aI).size() == called_alleles.at(0).size());
			}
			
			for(unsigned int pI = 0; pI < alleles_maxLength; pI++)
			{
				if(pI < reference_allele.size())
				{
					thisPosVec_referencePositions.push_back(position + pI);
				}
				else
				{
					thisPosVec_referencePositions.push_back(-1);
				}
			}
			assert(thisPosVec_referencePositions.size() == alleles_maxLength);

			chars_2_graph.push_back(called_alleles);

			lastExtractedPosition = position + reference_allele.length() - 1;

			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}

	if(chars_2_graph.size() > 0)
	{
		for(int pos = lastExtractedPosition + 1; pos <= positionStop; pos++)
		{
			if(pos >= (int)referenceChromosome.length())
				break;

			std::vector<std::string> thisPosVec;
			thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
			chars_2_graph.push_back(thisPosVec);

			std::vector<int> thisPosVec_referencePositions;
			thisPosVec_referencePositions.push_back(pos);
			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}
	else
	{
		for(int pos = positionStart; pos <= positionStop; pos++)
		{
			if(pos >= (int)referenceChromosome.length())
				break;

			std::vector<std::string> thisPosVec;
			thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
			chars_2_graph.push_back(thisPosVec);

			std::vector<int> thisPosVec_referencePositions;
			thisPosVec_referencePositions.push_back(pos);
			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}

	VCFstream.close();
	
	assert(ret_graph_referencePositions.size() == chars_2_graph.size());
	for(unsigned int segmentI = 0; segmentI < chars_2_graph.size(); segmentI++)
	{
		for(unsigned int haploI = 0; haploI < chars_2_graph.at(segmentI).size(); haploI++)
		{
			assert(chars_2_graph.at(segmentI).at(haploI).length() == ret_graph_referencePositions.at(segmentI).size());
		}
	}
	return chars_2_graph;
}

void _addPufferChromotypes(diploidGenomeString& gS)
{
	std::vector<std::string> puffer;
	puffer.push_back("**");
	gS.push_back(puffer);
}

void storeGenomeStringInFile(diploidGenomeString& gS, std::string filename)
{
	ofstream genomeStringStream;
	genomeStringStream.open(filename.c_str());

	if(! genomeStringStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for writing: "+ filename);
	}

	for(size_t gSl = 0; gSl < gS.size(); gSl++)
	{
		std::vector<std::string> gS_compartment = gS.at(gSl);
		assert((gS_compartment.size() == 1) || (gS_compartment.size() == 2));
		if(gS_compartment.size() == 1)
		{
			gS_compartment.push_back("");
		}

		genomeStringStream << Utilities::join(gS_compartment, ",") << "\n";
	}

	genomeStringStream.close();
}

diploidGenomeString readGenomeStringFromChromotypesFile(std::string filename, int kMer_size, std::string graphDir, int chromotypes_startCoordinate, std::vector<int>& amendedChromotypes_genomicGraphLoci)
{
	std::vector<int> genomicGraphLoci_originalChromotypes = getGenomicGraphLoci(graphDir, chromotypes_startCoordinate);
	amendedChromotypes_genomicGraphLoci.clear();
	
	diploidGenomeString gS;

	ifstream genomeStringStream;
	genomeStringStream.open(filename.c_str());

	if(! genomeStringStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for reading: "+ filename);
	}

	std::string line;
	size_t lineCounter = 0;

	std::vector<std::string> chromotype_1;
	std::vector<std::string> chromotype_2;
	
	while(genomeStringStream.good())
	{
		std::getline(genomeStringStream, line);

		Utilities::eraseNL(line);

		if(line.length() == 0)
			continue;
			
		lineCounter++;
					
		std::vector<std::string> lineFields = Utilities::split(line, " ");

		if(lineCounter == 1)
		{
			chromotype_1 = std::vector<std::string>(lineFields.begin() + 2, lineFields.end());
		}
		else if(lineCounter == 2)
		{
			chromotype_2 = std::vector<std::string>(lineFields.begin() + 2, lineFields.end());
		}
	}

	assert(lineCounter == 2);
	genomeStringStream.close();
	
	assert(chromotype_1.size() == chromotype_2.size());
	assert(chromotype_1.size() == genomicGraphLoci_originalChromotypes.size());
	
	std::vector<int> nextKidentical;
	nextKidentical.resize(chromotype_1.size(), false);
	int kIdentical = 0;
	for(int i = chromotype_1.size() - 1; i >= 0; i--)
	{
		if(chromotype_1.at(i) == chromotype_2.at(i))
		{
			kIdentical++;
		}
		else
		{
			kIdentical = 0;
		}
		nextKidentical.at(i) = kIdentical;
	}
	
	class genoStringBlock {
	public:
		bool diploid;
		int start;
		int stop;
		genoStringBlock()
		{
			start = -1;
			stop = -1;
		}
	};
	
	// we will only use the values for this first block if the first segment is non-diploid
	genoStringBlock currentBlock;
	// currentBlock.start = 0;
	// currentBlock.diploid = true;
	// currentBlock.stop = 0;
	
	std::vector<genoStringBlock> foundBlocks;
	for(int i = 0; i < (int)chromotype_1.size(); i++)
	{
		if(nextKidentical.at(i) >= kMer_size)
		{
			if(currentBlock.start != -1)
			{
				foundBlocks.push_back(currentBlock);
			}
			
			currentBlock = genoStringBlock();
			currentBlock.start = i;
			currentBlock.stop = i + nextKidentical.at(i) - 1;
			currentBlock.diploid = false;
			foundBlocks.push_back(currentBlock);
			
			i = i + nextKidentical.at(i) - 1;
			
			currentBlock = genoStringBlock();
		}
		else
		{
			if(currentBlock.start == -1)
			{
				currentBlock.start = i;
				currentBlock.stop = i;
				currentBlock.diploid = true;
			}
			else
			{
				currentBlock.stop++;
			}
		}
	}
	
	if(currentBlock.start != -1)
	{
		foundBlocks.push_back(currentBlock);
	}	

	int total_block_characters = 0;
	for(unsigned int blockI = 0; blockI < foundBlocks.size(); blockI++)
	{
		genoStringBlock thisBlock = foundBlocks.at(blockI);
		
		if(blockI > 0)
		{
			if(thisBlock.diploid)
			{
				assert(foundBlocks.at(blockI - 1).diploid == false);
			}
			else
			{
				assert(foundBlocks.at(blockI - 1).diploid == true);
			}
		}
		
		total_block_characters += (thisBlock.stop - thisBlock.start + 1);
		
		std::vector<std::string> gS_segment;
		if(thisBlock.diploid)
		{
			std::string s1;
			std::string s2;
			
			for(int cI = thisBlock.start; cI <= thisBlock.stop; cI++)
			{
				assert(cI >= 0);
				if(chromotype_1.at(cI).length() != chromotype_2.at(cI).length())
				{
					if(chromotype_1.at(cI).length() < chromotype_2.at(cI).length())
					{
						int missing_gaps_c1 = chromotype_2.at(cI).length() - chromotype_1.at(cI).length();
						for(int gapI = 1; gapI <= missing_gaps_c1; gapI++)
						{
							chromotype_1.at(cI).push_back('_');
						}
					}
					else
					{
						int missing_gaps_c2 = chromotype_1.at(cI).length() - chromotype_2.at(cI).length();
						assert(missing_gaps_c2 > 0);
						for(int gapI = 1; gapI <= missing_gaps_c2; gapI++)
						{
							chromotype_2.at(cI).push_back('_');
						}				
					}
				}
				s1.append(chromotype_1.at(cI));
				s2.append(chromotype_2.at(cI));
				
				amendedChromotypes_genomicGraphLoci.push_back(genomicGraphLoci_originalChromotypes.at(cI));
				assert(chromotype_1.at(cI).size() == chromotype_2.at(cI).size());
				for(unsigned int j = 1; j < chromotype_1.at(cI).size(); j++)
				{
					amendedChromotypes_genomicGraphLoci.push_back(-1);						
				}
			}
			gS_segment.push_back(s1);
			gS_segment.push_back(s2);
		}
		else
		{
			std::string s;
			for(unsigned int cI = thisBlock.start; cI <= thisBlock.stop; cI++)
			{
				// if(chromotype_1.at(cI).length() != 1)
				// {
					// std::cerr << "cI: " << cI << "\n";
					// std::cerr << "chromotype_1.at(cI).length(): " << chromotype_1.at(cI).length() << "\n";
					// std::cerr << "chromotype_1.at(cI): " << chromotype_1.at(cI) << "\n";
					// std::cerr << "chromotype_2.at(cI): " << chromotype_2.at(cI) << "\n";
				// }
				// assert(chromotype_1.at(cI).length() == 1);			
				assert(chromotype_1.at(cI) == chromotype_2.at(cI));
				s.append(chromotype_1.at(cI));
				
				amendedChromotypes_genomicGraphLoci.push_back(genomicGraphLoci_originalChromotypes.at(cI));				
				for(unsigned int j = 1; j < chromotype_1.at(cI).size(); j++)
				{
					amendedChromotypes_genomicGraphLoci.push_back(-1);						
				}
			}	
			gS_segment.push_back(s);			
		}

		gS.push_back(gS_segment);
	}
	
	assert(total_block_characters == chromotype_1.size());
	
	return gS;
}

diploidGenomeString readGenomeStringFromFile(std::string filename, bool ignorePuffer)
{
	diploidGenomeString gS;

	ifstream genomeStringStream;
	genomeStringStream.open(filename.c_str());

	if(! genomeStringStream.is_open())
	{
		throw std::runtime_error("Cannot open genome string file for reading: "+ filename);
	}

	std::string line;
	size_t lineCounter = 0;

	while(genomeStringStream.good())
	{
		std::getline(genomeStringStream, line);
		lineCounter++;

		Utilities::eraseNL(line);

		if(line.length() == 0)
			continue;

		std::vector<std::string> lineFields = Utilities::split(line, ",");
		assert(lineFields.size() == 2);

		assert(lineFields.at(0).size() > 0);

		std::vector<std::string> gS_compartment;
		gS_compartment.push_back(lineFields.at(0));

		if(lineFields.at(1).size() > 0)
		{
			gS_compartment.push_back(lineFields.at(1));
		}

		gS.push_back(gS_compartment);
	}
	
	if(ignorePuffer)
	{
		size_t lastCompartmentI = gS.size() - 1;
		assert(gS.at(lastCompartmentI).at(0).length() > 2);
		
		gS.at(lastCompartmentI).at(0) = gS.at(lastCompartmentI).at(0).substr(0, gS.at(lastCompartmentI).at(0).length()-2);
		if(gS.at(lastCompartmentI).size() > 1)
		{
			gS.at(lastCompartmentI).at(1) = gS.at(lastCompartmentI).at(1).substr(0, gS.at(lastCompartmentI).at(1).length()-2);	
			assert(gS.at(lastCompartmentI).at(0).length() == gS.at(lastCompartmentI).at(1).length());
		}
	}
	genomeStringStream.close();

	return gS;

}

std::vector<std::string> readGraphLoci(std::string graphDir)
{
	std::vector<std::string> forReturn;

	std::string segmentsFile = graphDir + "/segments.txt";

	ifstream segmentsStream;
	segmentsStream.open(segmentsFile.c_str());
	if(! segmentsStream.is_open())
	{
		throw std::runtime_error("Cannot open segments file: "+segmentsFile);
	}

	std::vector<std::string> individualSegmentFiles;
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			individualSegmentFiles.push_back(graphDir+"/"+line);
		}
	}
	segmentsStream.close();
	
	for(unsigned int segmentI = 0; segmentI < individualSegmentFiles.size(); segmentI++)
	{
		std::string F = individualSegmentFiles.at(segmentI);
		ifstream thisSegmentStream;
		thisSegmentStream.open(F.c_str());
		if(! thisSegmentStream.is_open())
		{
			throw std::runtime_error("Cannot open one segment file: "+F);
		}
		assert(thisSegmentStream.good());
		std::string firstLine;
		std::getline(thisSegmentStream, firstLine);
		Utilities::eraseNL(firstLine);

		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
			
		for(unsigned int fI = 1; fI < firstLine_fields.size(); fI++)
		{
			forReturn.push_back(firstLine_fields.at(fI));
		}

		thisSegmentStream.close();
	}

	return forReturn;
}

std::vector<int> graphLoci_2_PGFpositions(std::vector<std::string> graphLoci)
{
	std::vector<int> PGFpositions;
	int lastPGFposition;
	for(unsigned int i = 0; i < graphLoci.size(); i++)
	{
		std::string locusID = graphLoci.at(i);
		std::vector<std::string> locusParts = Utilities::split(locusID, "_");
		if(locusParts.size() != 3)
		{
			throw std::runtime_error("graphLoci_2_PGFpositions(..): Cannot decompose locus ID " +locusID);
		}
		int thisLocus_PGF = Utilities::StrtoI(locusParts.at(2));
		if((i == 0) || (lastPGFposition != thisLocus_PGF))
		{
			PGFpositions.push_back(thisLocus_PGF);
			lastPGFposition = thisLocus_PGF;
		}
		else
		{
			PGFpositions.push_back(-1);
		}
	}
	assert(PGFpositions.size() == graphLoci.size());
	return PGFpositions;
}

std::vector<int> getGenomicGraphLoci(std::string graphDir, int chromotypes_startCoordinate)
{
	std::vector<int> graphLoci_2_genomicPositions;
	std::vector<std::string> lociInGraph = readGraphLoci(graphDir);
	std::vector<int> lociInGraph_PGFpositions = graphLoci_2_PGFpositions(lociInGraph);
	int firstOffSet = lociInGraph_PGFpositions.at(0);
	for(unsigned int i = 0; i < lociInGraph_PGFpositions.size(); i++)
	{
		if(lociInGraph_PGFpositions.at(i) == -1)
		{
			graphLoci_2_genomicPositions.push_back(-1);
		}
		else
		{
			graphLoci_2_genomicPositions.push_back(lociInGraph_PGFpositions.at(i) - firstOffSet + chromotypes_startCoordinate);
		}
	}
	assert(graphLoci_2_genomicPositions.size() == lociInGraph.size());
	return graphLoci_2_genomicPositions;
}		
