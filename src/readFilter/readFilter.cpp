/*
 * readFilter.cpp
 *
 *  Created on: 14.01.2014
 *      Author: AlexanderDilthey
 */

#include "readFilter.h"

#include <assert.h>
#include <exception>
#include <stdexcept>
#include <functional>
#include <fstream>
#include <ostream>
#include <istream>
#include <set>

#include "../Utilities.h"

#include "../hash/deBruijn/DeBruijnGraph.h"
#include "../hash/sequence/basic.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"
#include "utils/bamtools_utilities.h"

readFilter::readFilter() {
	positiveThreshold = -1;
	negativeThreshold = -1;
	k = -1;
	positiveUnique = false;
	negativePreserveUnique = false;

	positiveUnique_threshold = 1;
	negativePreserveUnique_threshold = 1;
}


void readFilter::doFilter()
{
	if(!(positiveFilter.length() || negativeFilter.length()))
	{
		throw std::runtime_error("Please specify either positive filter or negative filter.");
	}

	if(!(input_BAM.length() || input_FASTQ.length()))
	{
		throw std::runtime_error("Please specify either input BAM or input FASTQ.");
	}

	if(input_BAM.length() && input_FASTQ.length())
	{
		throw std::runtime_error("Both input BAM and input FASTQ set - can't deal with that.");
	}

	assert((positiveThreshold >= 0) && (positiveThreshold <= 1));
	assert((negativeThreshold >= 0) && (negativeThreshold <= 1));

	if(positiveUnique || negativePreserveUnique)
	{
		assert(uniqueness_base.length());
		assert(uniqueness_subtract.length());
	}

	std::string fn_1 = output_FASTQ + "_1";
	std::string fn_2 = output_FASTQ + "_2";

	std::ofstream fastq_1_output;
	fastq_1_output.open(fn_1.c_str());
	if(! fastq_1_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+fn_1);
	}

	std::ofstream fastq_2_output;
	fastq_2_output.open(fn_2.c_str());
	if(! fastq_2_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+fn_2);
	}

	std::set<std::string> positive_kMers;

	bool apply_filter_positive = (positiveFilter.length() > 0);
	bool apply_filter_negative = (negativeFilter.length() > 0);

	auto load_positive_kMers_file = [&](std::string file) -> std::set<std::string> {
		std::set<std::string> forReturn;

		std::ifstream positive_kMers_stream;
		positive_kMers_stream.open(file.c_str());
		if(! positive_kMers_stream.is_open())
		{
			throw std::runtime_error("readFilter::doFilter(): Cannot open kMers file containing the positive-filter kMers: "+file);
		}
		std::string line;
		size_t line_number = 0;
		while(positive_kMers_stream.good())
		{
			line_number++;

			getline (positive_kMers_stream, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			std::string kMer = line;

			if(kMer.length() != k)
			{
				throw std::runtime_error("readFilter::doFilter(): Expect kMers of length " + Utilities::ItoStr(k) + ", but " + positiveFilter + " contains one of length " + Utilities::ItoStr(kMer.length()) + " (line " + Utilities::ItoStr(line_number) + " ).");
			}

			forReturn.insert(kMer);
		}

		positive_kMers_stream.close();

		return forReturn;
	};

	if(apply_filter_positive)
	{
		positive_kMers = load_positive_kMers_file(positiveFilter);
	}


	int cortex_height = 26;
	int cortex_width = 50;

	std::set<std::string> unique_kMers;
	if(positiveUnique || negativePreserveUnique)
	{
		unique_kMers = load_positive_kMers_file(uniqueness_base);
		DeBruijnGraph<1, 25, 1> subtract_kMers_graph(cortex_height, cortex_width);
		for(std::set<std::string>::iterator kMerIt = unique_kMers.begin(); kMerIt != unique_kMers.end(); kMerIt++)
		{
			std::string kMer = *kMerIt;
			if(subtract_kMers_graph.kMerinGraph(kMer))
			{
				unique_kMers.erase(kMer);
			}
		}
	}

	DeBruijnGraph<1, 25, 1>* negative_kMers;
	if(apply_filter_negative)
	{

		std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;

		assert(k == 25);
		negative_kMers = new DeBruijnGraph<1, 25, 1>(cortex_height, cortex_width);

		std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;

		negative_kMers->loadMultiColourBinary(negativeFilter);

		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		std::cout << "\tTotal coverage: " << negative_kMers->totalCoverage() << "\n";
	}

	std::function<bool(const fastq_readPair&)> decisionFunction = [&](const fastq_readPair& read) -> bool {

		std::vector<std::string> kMers_1 = partitionStringIntokMers(read.a1.sequence, k);
		std::vector<std::string> kMers_2 = partitionStringIntokMers(read.a2.sequence, k);

		bool pass_positive = true;
		if(apply_filter_positive)
		{
			// forward check

			double kMers_1_forward_OK = 0;
			double kMers_2_forward_OK = 0;
			double kMers_1_forward_TOTAL = 0;
			double kMers_2_forward_TOTAL = 0;

			int kMers_1_forward_unique = 0;
			int kMers_2_forward_unique = 0;

			kMers_1_forward_TOTAL += kMers_1.size();
			kMers_2_forward_TOTAL += kMers_2.size();

			for(unsigned int kI = 0; kI < kMers_1.size(); kI++)
			{
				std::string kMer = kMers_1.at(kI);
				if(positive_kMers.count(kMer))
				{
					kMers_1_forward_OK++;
				}

				if(unique_kMers.count(kMer))
				{
					kMers_1_forward_unique++;
				}
			}

			for(unsigned int kI = 0; kI < kMers_2.size(); kI++)
			{
				std::string kMer = kMers_2.at(kI);
				if(positive_kMers.count(kMer))
				{
					kMers_2_forward_OK++;
				}

				if(unique_kMers.count(kMer))
				{
					kMers_2_forward_unique++;
				}
			}

			double forward_1_optim = (kMers_1_forward_TOTAL == 0) ? 0 : (kMers_1_forward_OK / kMers_1_forward_TOTAL);
			double forward_2_optim = (kMers_2_forward_TOTAL == 0) ? 0 : (kMers_2_forward_OK / kMers_2_forward_TOTAL);
			double forward_combined_optim = ((kMers_1_forward_TOTAL + kMers_2_forward_TOTAL) == 0) ? 0 : ((kMers_1_forward_OK + kMers_2_forward_OK) / (kMers_1_forward_TOTAL + kMers_2_forward_TOTAL));

			// reverse check

			double kMers_1_reverse_OK = 0;
			double kMers_2_reverse_OK = 0;
			double kMers_1_reverse_TOTAL = 0;
			double kMers_2_reverse_TOTAL = 0;

			int kMers_1_reverse_unique = 0;
			int kMers_2_reverse_unique = 0;

			kMers_1_reverse_TOTAL += kMers_1.size();
			kMers_2_reverse_TOTAL += kMers_2.size();


			for(unsigned int kI = 0; kI < kMers_1.size(); kI++)
			{
				std::string kMer = seq_reverse_complement(kMers_1.at(kI));
				if(positive_kMers.count(kMer))
				{
					kMers_1_reverse_OK++;
				}
				if(unique_kMers.count(kMer))
				{
					kMers_1_reverse_unique++;
				}
			}

			for(unsigned int kI = 0; kI < kMers_2.size(); kI++)
			{
				std::string kMer = seq_reverse_complement(kMers_2.at(kI));
				if(positive_kMers.count(kMer))
				{
					kMers_2_reverse_OK++;
				}
				if(unique_kMers.count(kMer))
				{
					kMers_2_reverse_unique++;
				}
			}

//			double reverse_1_optim = (kMers_1_reverse_TOTAL == 0) ? 0 : (kMers_1_reverse_OK / kMers_1_reverse_TOTAL);
//			double reverse_2_optim = (kMers_2_reverse_TOTAL == 0) ? 0 : (kMers_2_reverse_OK / kMers_2_reverse_TOTAL);
			double reverse_combined_optim = ((kMers_1_reverse_TOTAL + kMers_2_reverse_TOTAL) == 0) ? 0 : ((kMers_1_reverse_OK + kMers_2_reverse_OK) / (kMers_1_reverse_TOTAL + kMers_2_reverse_TOTAL));

			int forward_combined_unique = kMers_1_forward_unique + kMers_2_forward_unique;
			int reverse_combined_unique = kMers_1_reverse_unique + kMers_2_reverse_unique;

			// std::cout << read.a1.sequence << " " << read.a2.sequence << "\n";
			// std::cout << forward_combined_optim << " " << reverse_combined_optim << "\n\n"; 
			
			// std::cout << read.a1.readID << " // " << read.a2.readID << ": " << forward_combined_optim << " / " << reverse_combined_optim << "\n";

			pass_positive = ((forward_combined_optim >= positiveThreshold) || (reverse_combined_optim >= positiveThreshold));

			if(positiveUnique)
			{
				pass_positive = ( pass_positive || ((forward_combined_unique >= positiveUnique_threshold) || (reverse_combined_unique >= positiveThreshold)) );
			}
		}

		bool pass_negative = false;
		if(pass_positive && apply_filter_negative)
		{
			// forward check

			double kMers_1_notOK = 0;
			double kMers_2_notOK = 0;
			double kMers_1_TOTAL = 0;
			double kMers_2_TOTAL = 0;

			kMers_1_TOTAL += kMers_1.size();
			kMers_2_TOTAL += kMers_2.size();

			for(unsigned int kI = 0; kI < kMers_1.size(); kI++)
			{
				std::string kMer = kMers_1.at(kI);
				if(negative_kMers->kMerinGraph(kMer))
				{
					kMers_1_notOK++;
				}
			}

			for(unsigned int kI = 0; kI < kMers_2.size(); kI++)
			{
				std::string kMer = kMers_2.at(kI);
				if(negative_kMers->kMerinGraph(kMer))
				{
					kMers_2_notOK++;
				}
			}

//			double negativity_1 = (kMers_1_TOTAL == 0) ? 1 : (kMers_1_notOK / kMers_1_TOTAL);
//			double negativity_2 = (kMers_2_TOTAL == 0) ? 1 : (kMers_2_notOK / kMers_2_TOTAL);
			double combined_negativity = ((kMers_1_TOTAL + kMers_2_TOTAL) == 0) ? 1 : ((kMers_1_notOK + kMers_2_notOK) / (kMers_1_TOTAL + kMers_2_TOTAL));

			// std::cout << read.a1.readID << " // " << read.a2.readID << ": " << combined_negativity << "\n";
			
			pass_negative = (combined_negativity <= negativeThreshold);

			return (pass_positive && pass_negative);
		}
		else
		{
			return pass_positive;
		}
	};

	std::function<void(const fastq_readPair&)> printFunction = [&](const fastq_readPair& read) -> void {

		fastq_1_output << "@" << read.a1.readID << "\n"
				  << read.a1.sequence    << "\n"
				  << "+"         << "\n"
				  << read.a1.qualities   << "\n";

		fastq_2_output << "@" << read.a2.readID << "\n"
						  << read.a2.sequence    << "\n"
						  << "+"         << "\n"
						  << read.a2.qualities   << "\n";
	};

	if(input_BAM.length())
	{
		filterBAM(input_BAM, output_FASTQ, &decisionFunction, &printFunction);
	}
	else
	{
		filterFastQPairs(input_FASTQ, output_FASTQ, &decisionFunction, &printFunction);
	}

	fastq_1_output.close();
	fastq_1_output.close();

	if(apply_filter_negative)
	{
		delete(negative_kMers);
	}
}

void filterFastQPairs(std::string fastq_basePath, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print)
{
	std::string file_1 = fastq_basePath + "_1";
	std::string file_2 = fastq_basePath + "_2";

	if(! Utilities::fileReadable(file_1))
	{
		throw std::runtime_error("Expected file "+file_1+" can't be opened.");
	}
	if(! Utilities::fileReadable(file_2))
	{
		throw std::runtime_error("Expected file "+file_2+" can't be opened.");
	}

	filterFastQPairs(file_1, file_2, outputFile, decide, print);
}

void filterFastQPairs(std::string fastq_1_path, std::string fastq_2_path, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print)
{
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
		
		BAMalignment simpleAlignment_1;
		simpleAlignment_1.readID = read1_ID;
		simpleAlignment_1.qualities = read1_qualities;
		simpleAlignment_1.sequence = read1_sequence;

		BAMalignment simpleAlignment_2;
		simpleAlignment_2.readID = read2_ID;
		simpleAlignment_2.qualities = read2_qualities;
		simpleAlignment_2.sequence = read2_sequence;

		fastq_readPair thisPair;
		bool success_1 = thisPair.takeAlignment(simpleAlignment_1, 1);
		assert(success_1);
		bool success_2 = thisPair.takeAlignment(simpleAlignment_2, 2);
		assert(success_2);
		assert(thisPair.isComplete());

		assert((read1_ID.substr(read1_ID.length() - 2, 2) == "/1") || (read1_ID.substr(read1_ID.length() - 2, 2) == "/2"));
		assert((read2_ID.substr(read2_ID.length() - 2, 2) == "/1") || (read2_ID.substr(read2_ID.length() - 2, 2) == "/2"));
		
		if(!(read1_ID.substr(0, read1_ID.length() - 2) == read2_ID.substr(0, read2_ID.length() - 2)))
		{
			std::cerr << "Warning: read IDs don't match! " << read1_ID << " vs " << read2_ID << "\n";
		}
		assert(read1_ID.substr(0, read1_ID.length() - 2) == read2_ID.substr(0, read2_ID.length() - 2));		
		
		if((*decide)(thisPair))
		{
			(*print)(thisPair);
		}
	}
}

void filterBAM(std::string BAMfile, std::string outputFile, std::function<bool(const fastq_readPair&)>* decide, std::function<void(const fastq_readPair&)>* print)
{
	BamTools::BamReader reader;
	reader.Open(BAMfile);

	reader.LocateIndex();
    if ( ! reader.HasIndex() )
    {
		throw std::runtime_error("File "+BAMfile+" does not seem to be indexed - please specify indexed BAM!");
    }

    std::map<std::string, fastq_readPair> reads;

	std::vector<BAMRegionSpecifier> BAM_regions = getBAMregions(BAMfile);
	for(unsigned int rI = 0; rI < BAM_regions.size(); rI++)
	{
		const BAMRegionSpecifier& thisStretch = BAM_regions.at(rI);

		int refIDidx = reader.GetReferenceID(thisStretch.ID);
		assert(refIDidx != -1);

		const BamTools::RefData& stretchSpec_BAMTools = reader.GetReferenceData().at(refIDidx);
		assert(thisStretch.lastPos < stretchSpec_BAMTools.RefLength);

		std::cout << "\t" << Utilities::timestamp() << " read " << thisStretch.ID << " from " << thisStretch.firstPos << " to " << thisStretch.lastPos + 1 << "\n" << std::flush;

		BamTools::BamRegion stretch_region_BAMTools;
		stretch_region_BAMTools.LeftRefID = refIDidx;
		stretch_region_BAMTools.LeftPosition = thisStretch.firstPos;
		stretch_region_BAMTools.RightRefID = refIDidx;;
		stretch_region_BAMTools.RightPosition =  thisStretch.lastPos + 1;

		reader.SetRegion(stretch_region_BAMTools);

		BamTools::BamAlignment al;
		#pragma omp critical
		{
			while(reader.GetNextAlignment(al))
			{
				std::string name = al.Name;
			   std::string nameWithPairID = name;

			   int whichMate = 0;
			   assert(al.IsPaired());
			   if ( al.IsPaired() )
			   {
				   nameWithPairID.append( (al.IsFirstMate() ? "/1" : "/2") );
				   whichMate =  (al.IsFirstMate()) ? 1 : 2;
			   }

				// handle reverse strand alignment - bases & qualities
				std::string qualities = al.Qualities;
				std::string sequence  = al.QueryBases;
				if ( al.IsReverseStrand() ) {
					std::reverse(qualities.begin(), qualities.end());
					Utilities::seq_reverse_complement(sequence);
				}

				BAMalignment simpleAlignment;
				simpleAlignment.readID = nameWithPairID;
				simpleAlignment.qualities = qualities;
				simpleAlignment.sequence = sequence;

				// std::cout << name << " " << sequence << "\n";
				
				// std::cout << name << " " << reads.count(name) << "\n";
				if(reads.count(name) == 0)
				{
					fastq_readPair p;
					bool success = p.takeAlignment(simpleAlignment, whichMate);
					assert(success);
					reads[name] = p;
				}
				else
				{
					fastq_readPair& thisPair = reads.at(name);
					bool success = thisPair.takeAlignment(simpleAlignment, whichMate);
					if(! success)
					{					
						std::cerr << "There is a problem with the read IDs in this BAM.\n";
						std::cerr << "Read ID: " << name << " / " << nameWithPairID << "\n";
						std::cerr << "whichMate: " << whichMate << "\n";
						std::cerr << "thisPair.have1: " << thisPair.have1 << " with ID " << thisPair.a1.readID << "\n";
						std::cerr << "thisPair.have2: " << thisPair.have2 << " with ID " << thisPair.a2.readID << "\n" << std::flush;
					}
					assert(success);
					if(thisPair.isComplete())
					{
						// process
						if((*decide)(thisPair))
						{
							(*print)(thisPair);
						}

						reads.erase(name);
					}
				}

			}
		}
	}

	if(reads.size() > 0)
	{
		std::cerr << "\n\n!!!!!!!!!!!!!!!!!!!!!!!\n\nAfter processing " << BAMfile << ", have " << reads.size() << " dangling reads.\n\n!!!!!!!!!!!!!!!!!!!!!!!\n\n";
	}

}

std::vector<BAMRegionSpecifier> getBAMregions(std::string BAMfile)
{
	std::vector<BAMRegionSpecifier> forReturn;

	BamTools::BamReader metaReader;
	metaReader.Open(BAMfile);

	metaReader.LocateIndex();
    if ( ! metaReader.HasIndex() )
    {
		throw std::runtime_error("File "+BAMfile+" does not seem to be indexed - please specify indexed BAM!");
    }

    BamTools::RefVector availableRegions = metaReader.GetReferenceData();
    for(unsigned int i = 0; i < availableRegions.size(); i++)
    {
    	BamTools::RefData thisRegion = availableRegions.at(i);

		size_t firstPos = 0;
		size_t lastPos =  thisRegion.RefLength - 1;

		BAMRegionSpecifier thisStretch;
		thisStretch.ID = thisRegion.RefName;
		thisStretch.firstPos = firstPos;
		thisStretch.lastPos = lastPos;

		forReturn.push_back(thisStretch);
    }

    return forReturn;
}
