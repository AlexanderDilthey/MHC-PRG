/*
 * readSimulator.cpp
 *
 *  Created on: 21.05.2013
 *      Author: AlexanderDilthey
 */

#include <omp.h>
#include <algorithm>
#include <utility>
#include <fstream>
#include <cmath>

#include "readSimulator.h"
#include "../Utilities.h"


#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>

std::string readName_field_separator = "|||";

double readSimulator::averageErrorRate(std::vector<std::map<char, double> > q_freq, std::vector<std::map<char, double> > error_freq_conditional_q)
{
	unsigned int readLength = q_freq.size();
	assert(error_freq_conditional_q.size() == readLength);

	double forReturn = 0;

	for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
	{
		double thisPos_errorRate = 0;
		assert(q_freq.at(positionInRead).size() == error_freq_conditional_q.at(positionInRead).size());
		Utilities::check_map_is_normalized(q_freq.at(positionInRead));
		for(std::map<char, double>::iterator qIt = q_freq.at(positionInRead).begin(); qIt != q_freq.at(positionInRead).end(); qIt++)
		{
			char quality = qIt->first;
			double p_quality = qIt->second;
			double p_error = 1 - error_freq_conditional_q.at(positionInRead).at(quality);
			assert(p_error >= 0);
			assert(p_error <= 1);
			thisPos_errorRate += p_quality * p_error;
		}

		forReturn += thisPos_errorRate;
	}

	forReturn = forReturn / (double) readLength;

	return forReturn;
}

readSimulator::readSimulator(std::string qualityMatrixFile, unsigned int readLength, bool interpolateLength_, char removeUpperBaseQualityIndices, char additional_2ndRead_removeUpperBaseQualityIndices) {

	// attention - these are chars - make sure that there are no overflows
	assert(removeUpperBaseQualityIndices >= 0);
	assert(removeUpperBaseQualityIndices <= 20);
	assert(additional_2ndRead_removeUpperBaseQualityIndices >= 0);
	assert(additional_2ndRead_removeUpperBaseQualityIndices <= 20);
	
	std::ifstream matrixStream;
	matrixStream.open (qualityMatrixFile.c_str(), std::ios::in);

	interpolateLength = interpolateLength_;


	read_quality_frequencies.resize(readLength);
	read_quality_correctness.resize(readLength);
	read_INDEL_freq.resize(readLength, 0);

	std::map< int, std::vector<std::map<char, double> > > read_quality_frequencies_perLength;
	std::map< int, std::vector<std::map<char, double> > > read_quality_correctness_perLength;
	std::map< int, std::vector<double> > read_INDEL_freq_perLength;
	std::map< int, std::vector<double> > read_nonINDEL_freq_perLength;

	if(matrixStream.is_open())
	{
		std::string line;
		std::string currentIdentifier;
		std::string currentSequence;

		size_t lineCounter = 0;

		while(matrixStream.good())
		{
			std::getline(matrixStream, line);
			lineCounter++;

			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			std::vector<std::string> line_fields = Utilities::split(line, "\t");

			if(lineCounter == 1)
			{
				if(!(line_fields.size() == 6))
				{
					throw std::runtime_error("readSimulator::readSimulator(): file "+qualityMatrixFile+", expect 6 fields, but have "+Utilities::ItoStr(line_fields.size())+"; line "+Utilities::ItoStr(lineCounter)+".\n"+line+"\n");
				}
				assert(line_fields.at(0) == "readLength");
				assert(line_fields.at(1) == "qualityScore");
				assert(line_fields.at(2) == "positionInRead");
				assert(line_fields.at(3) == "N");
				assert(line_fields.at(4) == "ExpectedCorrect");
				assert(line_fields.at(5) == "EmpiricalCorrect");
			}
			else
			{
				int thisLine_readLength = Utilities::StrtoI(line_fields.at(0));

				if(read_quality_frequencies_perLength.count(thisLine_readLength) == 0)
				{
					read_quality_frequencies_perLength[thisLine_readLength].resize(thisLine_readLength);
					read_quality_correctness_perLength[thisLine_readLength].resize(thisLine_readLength);
					read_INDEL_freq_perLength[thisLine_readLength].resize(thisLine_readLength, 0);
					read_nonINDEL_freq_perLength[thisLine_readLength].resize(thisLine_readLength, 0);
				}

				char quality = line_fields.at(1).at(0);
				unsigned int positionInRead = Utilities::StrtoI(line_fields.at(2));
				int datapoints = Utilities::StrtoI(line_fields.at(3));

				assert(positionInRead < readLength);


				if(quality == 0)
				{
					read_INDEL_freq_perLength.at(thisLine_readLength).at(positionInRead) += datapoints;
				}
				else
				{
					read_nonINDEL_freq_perLength.at(thisLine_readLength).at(positionInRead) += datapoints;

					if(read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0)
					{
						read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead)[quality] = 0;
					}
					read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead)[quality] += datapoints;

					if(!(read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0))
					{
						throw std::runtime_error("readSimulator::readSimulator(): file "+qualityMatrixFile+", double-defined quality. Line "+Utilities::ItoStr(lineCounter)+".\n"+line+"\n");
					}
					assert(read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0);
					read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead)[quality] += Utilities::StrtoD(line_fields.at(5));
				}
			}
		}
		matrixStream.close();

		std::vector<double> read_nonINDEL_freq;


		if(read_quality_correctness_perLength.count(readLength))
		{
			read_quality_frequencies = read_quality_frequencies_perLength.at(readLength);
			read_quality_correctness = read_quality_correctness_perLength.at(readLength);
			read_INDEL_freq = read_INDEL_freq_perLength.at(readLength);
			read_nonINDEL_freq = read_nonINDEL_freq_perLength.at(readLength);
		}
		else
		{
			assert(interpolateLength);

			assert(read_quality_frequencies_perLength.size() > 0);
			int closestNeighbour_which = -1;
			int closestNeighbour_distance = -1;

			for(std::map< int, std::vector<std::map<char, double> > >::iterator correctnessIt = read_quality_correctness_perLength.begin(); correctnessIt != read_quality_correctness_perLength.end(); correctnessIt++ )
			{
				int thisDist = abs(correctnessIt->first - (int)readLength);
				if((correctnessIt == read_quality_correctness_perLength.begin()) || (thisDist < closestNeighbour_distance))
				{
					closestNeighbour_which = correctnessIt->first;
					closestNeighbour_distance = thisDist;
				}
			}
			assert(closestNeighbour_distance >= 0);

			read_quality_frequencies.resize(readLength);
			read_quality_correctness.resize(readLength);
			read_INDEL_freq.resize(readLength);
			read_nonINDEL_freq.resize(readLength);

			for(unsigned int posInRead = 0; posInRead < readLength; posInRead++)
			{
				double fraction = (double)posInRead / (double) readLength;
				double idx_target = fraction * (closestNeighbour_which - 1);
				int idx_target_int = round(idx_target);
				assert(idx_target_int >= 0);
				assert(idx_target_int < (int)read_quality_correctness_perLength.at(closestNeighbour_which).size());

				read_quality_frequencies.at(posInRead) = read_quality_frequencies_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_quality_correctness.at(posInRead) = read_quality_correctness_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_INDEL_freq.at(posInRead) = read_INDEL_freq_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_nonINDEL_freq.at(posInRead) = read_nonINDEL_freq_perLength.at(closestNeighbour_which).at(idx_target_int);
			}
		}

		for(unsigned int i = 0; i < readLength; i++)
		{
			assert(read_quality_correctness.at(i).size() > 0);
			assert(read_quality_frequencies.at(i).size() > 0);
			read_quality_frequencies.at(i) = Utilities::normalize_map(read_quality_frequencies.at(i));

			double total_data_position = read_INDEL_freq.at(i) + read_nonINDEL_freq.at(i);
			assert(total_data_position != 0);
			// std::cout << i << "\t" << read_INDEL_freq.at(i) << "\t" << total_data_position << "\n" << std::flush;
			read_INDEL_freq.at(i) = read_INDEL_freq.at(i) / total_data_position;
			if(read_INDEL_freq.at(i) == 0)
			{
				read_INDEL_freq.at(i) = 1e-4;
			}
		}
		
		read_quality_frequencies_2nd = read_quality_frequencies;
		read_quality_correctness_2nd = read_quality_correctness;
		read_INDEL_freq_2nd = read_INDEL_freq;		
		

		if(removeUpperBaseQualityIndices || additional_2ndRead_removeUpperBaseQualityIndices)
		{
			std::cout << "Before adjusting base qualities by " << (int)removeUpperBaseQualityIndices << " and " << (int)additional_2ndRead_removeUpperBaseQualityIndices << ", have average error:\n";
			std::cout << "\tR1: " << averageErrorRate(read_quality_frequencies, read_quality_correctness) << "\n";
			std::cout << "\tR2: " << averageErrorRate(read_quality_frequencies_2nd, read_quality_correctness_2nd) << "\n";
			std::cout << std::flush;
		}

		if(removeUpperBaseQualityIndices)
		{
			for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
			{
				std::map<char, double> new_qualityFrequencies;
				std::map<char, double> new_correctness;
				std::vector<char> existingQualities;
				for(std::map<char, double>::iterator existingFrequencyIt = read_quality_frequencies.at(positionInRead).begin(); existingFrequencyIt != read_quality_frequencies.at(positionInRead).end(); existingFrequencyIt++ )
				{
					existingQualities.push_back(existingFrequencyIt->first);
				}
				std::sort(existingQualities.begin(), existingQualities.end());
				if(existingQualities.size() > 1)
				{
					assert(existingQualities.at(0) <= existingQualities.at(1));
				}
				assert(removeUpperBaseQualityIndices >= 0);
				assert(removeUpperBaseQualityIndices < (int)(existingQualities.size()-1));
				for(int baseQualityI = 0; baseQualityI < (int)(existingQualities.size() - removeUpperBaseQualityIndices); baseQualityI++)
				{
					char qualityCharacter = existingQualities.at(baseQualityI);
					double f = read_quality_frequencies.at(positionInRead).at(qualityCharacter);
					new_qualityFrequencies[qualityCharacter] = f;
					new_correctness[qualityCharacter] =  read_quality_correctness.at(positionInRead).at(qualityCharacter);
				}
				new_qualityFrequencies = Utilities::normalize_map(new_qualityFrequencies);
				read_quality_frequencies.at(positionInRead) = new_qualityFrequencies;
				read_quality_correctness.at(positionInRead) = new_correctness;
			}
		}


		if(removeUpperBaseQualityIndices || additional_2ndRead_removeUpperBaseQualityIndices)
		{
			char shift_read_2 = removeUpperBaseQualityIndices + additional_2ndRead_removeUpperBaseQualityIndices;
			
			for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
			{
				std::map<char, double> new_qualityFrequencies;
				std::map<char, double> new_correctness;

				std::vector<char> existingQualities;
				for(std::map<char, double>::iterator existingFrequencyIt = read_quality_frequencies_2nd.at(positionInRead).begin(); existingFrequencyIt != read_quality_frequencies_2nd.at(positionInRead).end(); existingFrequencyIt++ )
				{
					existingQualities.push_back(existingFrequencyIt->first);
				}
				std::sort(existingQualities.begin(), existingQualities.end());
				if(existingQualities.size() > 1)
				{
					assert(existingQualities.at(0) <= existingQualities.at(1));
				}
				assert(shift_read_2 >= 0);
				assert(shift_read_2 < (int)(existingQualities.size()-1));
				for(int baseQualityI = 0; baseQualityI < (int)(existingQualities.size() - shift_read_2); baseQualityI++)
				{
					char qualityCharacter = existingQualities.at(baseQualityI);
					double f = read_quality_frequencies_2nd.at(positionInRead).at(qualityCharacter);
					new_qualityFrequencies[qualityCharacter] = f;
					new_correctness[qualityCharacter] =  read_quality_correctness_2nd.at(positionInRead).at(qualityCharacter);
				}

				new_qualityFrequencies = Utilities::normalize_map(new_qualityFrequencies);
				read_quality_frequencies_2nd.at(positionInRead) = new_qualityFrequencies;
				read_quality_correctness_2nd.at(positionInRead) = new_correctness;

			}
		}

	}
	else
	{
		throw std::runtime_error("readSimulator::readSimulator(): Cannot open file "+qualityMatrixFile);
	}

	for(unsigned i = 0; i < readLength; i++)
	{
		Utilities::check_map_is_normalized(read_quality_frequencies.at(i));
		Utilities::check_map_is_normalized(read_quality_frequencies_2nd.at(i));
	}
	
	if(1 == 1)
	{
		std::cout << "Summary error probabilities:\n";
		std::cout << "After adjusting base qualities by " << (int)removeUpperBaseQualityIndices << " and " << (int)additional_2ndRead_removeUpperBaseQualityIndices << ", have average error:\n";
		std::cout << "\tR1: " << averageErrorRate(read_quality_frequencies, read_quality_correctness) << "\n";
		std::cout << "\tR2: " << averageErrorRate(read_quality_frequencies_2nd, read_quality_correctness_2nd) << "\n";
		std::cout << std::flush;	
	}

	read_length = readLength;
	threads = 40;
	paranoid = true;
}


std::vector<oneReadPair> readSimulator::simulate_paired_reads_from_string(std::string S, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd, bool perfectly)
{
	std::vector<oneReadPair> forReturn;


	std::string edgePath_string = S;

	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_paired_reads_from_edgePath(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;

	{
		boost::mt19937 rnd_gen;

		auto seed = boost::random::random_device()();
		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);

			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);

				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}


			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}

			assert((! success) || (coordinates_string.size() == read.size()));

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		long long rPI = 0;
		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));

				rPI++;

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);
				sampleRead(i + this->read_length + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					std::reverse(read2_qualities.begin(), read2_qualities.end());

					// std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);
					// std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);

					std::string read1_name = "p" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);
					std::string read2_name = read1_name;
					
					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					std::vector<int> read1_coordinates_edgePath;
					std::vector<int> read2_coordinates_edgePath;
//					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
//					{
//						int c = read1_coordinates_string.at(cI);
//					}
//					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
//					{
//						int c = read2_coordinates_string.at(cI);
//					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;

					std::reverse(read2_coordinates_string.begin(), read2_coordinates_string.end());
					std::reverse(read2_coordinates_edgePath.begin(), read2_coordinates_edgePath.end());

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());
					assert(r2.coordinates_string.size() == r2.sequence.size());

					oneReadPair rP(r1, r2, jumpSize);

					if(Utilities::oneBernoulliTrial(0.5, rnd_gen))
					{
						rP.invert();
					}

					forReturn.push_back(rP);
				}
			}
		}

		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_string(..): Simulated " << forReturn.size() << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	return forReturn;
}


std::vector<oneRead> readSimulator::simulate_unpaired_reads_from_string(std::string S, double expected_haploid_coverage, bool perfectly)
{
	std::vector<oneRead> forReturn;

	std::string edgePath_string = S;

	double poissonStartRate = expected_haploid_coverage / (double)read_length; // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_unpaired_reads_from_string(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;

	{
		boost::mt19937 rnd_gen;

		auto seed = boost::random::random_device()();
		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);

			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);

				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}


			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}

			assert((! success) || (coordinates_string.size() == read.size()));

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		long long rPI = 0;
		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				rPI++;

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);

				if(read1_success)
				{
					thread_read_pairs++;

					// std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);
					// std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);

					std::string read1_name = "p" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);

					oneRead r1(read1_name, read1, read1_qualities);

					std::vector<int> read1_coordinates_edgePath;
//					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
//					{
//						int c = read1_coordinates_string.at(cI);
//					}
//					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
//					{
//						int c = read2_coordinates_string.at(cI);
//					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());

					forReturn.push_back(r1);
				}
			}
		}

		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_string(..): Simulated " << forReturn.size() << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	return forReturn;
}

std::vector<oneReadPair> readSimulator::simulate_paired_reads_from_edgePath(std::vector<Edge*> edgePath, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd, bool perfectly, bool is2nd)
{
	std::vector<oneReadPair> forReturn;

	Graph* g = edgePath.front()->From->g;
	std::string edgePath_string;
	std::vector<unsigned int> edgePath_string_originLevel;
	for(unsigned int eI = 0; eI < edgePath.size(); eI++)
	{
		std::string edgeEmission = g->CODE.deCode(edgePath.at(eI)->locus_id, edgePath.at(eI)->emission);
		assert(edgeEmission.length() == 1);
		if(edgeEmission != "_")
		{
			edgePath_string.append(edgeEmission);
			edgePath_string_originLevel.push_back(eI);
		}
	}
	assert(edgePath_string.size() == edgePath_string_originLevel.size());

	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_paired_reads_from_edgePath(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;
	
	{
		boost::mt19937 rnd_gen;

		auto seed = boost::random::random_device()();
		rnd_gen.seed(seed);  
		
		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			if(is2nd)
			{
				returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies_2nd.at(position_in_read), rnd_gen);
			}
			else
			{
				returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);
			}


			assert(returnedQuality > 0);

			bool generateError;

			if(is2nd)
			{
				generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness_2nd.at(position_in_read).at(returnedQuality), rnd_gen);
			}
			else
			{
				generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);
			}
						
			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};
		
		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}
				
				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);
						
						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);
				
				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}
		

			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}
			
			assert((! success) || (coordinates_string.size() == read.size()));
			
			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));


				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);
				sampleRead(i + this->read_length + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					std::reverse(read2_qualities.begin(), read2_qualities.end());
					
					std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i);
					std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize);

					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					std::vector<int> read1_coordinates_edgePath;
					std::vector<int> read2_coordinates_edgePath;
					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
					{
						int c = read1_coordinates_string.at(cI);
						if(c == -1)
						{
							read1_coordinates_edgePath.push_back(c);
						}
						else
						{
							read1_coordinates_edgePath.push_back(edgePath_string_originLevel.at(c));
						}
					}
					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
					{
						int c = read2_coordinates_string.at(cI);
						if(c == -1)
						{
							read2_coordinates_edgePath.push_back(c);
						}
						else
						{
							read2_coordinates_edgePath.push_back(edgePath_string_originLevel.at(c));
						}
					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;
					
					std::reverse(read2_coordinates_string.begin(), read2_coordinates_string.end());
					std::reverse(read2_coordinates_edgePath.begin(), read2_coordinates_edgePath.end());

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());
					assert(r2.coordinates_string.size() == r2.sequence.size());
					assert(r1.coordinates_string.size() == r1.coordinates_edgePath.size());
					assert(r2.coordinates_string.size() == r2.coordinates_edgePath.size());
					  
					oneReadPair rP(r1, r2, jumpSize);

					if(Utilities::oneBernoulliTrial(0.5, rnd_gen))
					{
						rP.invert();
					}

					forReturn.push_back(rP);
				}
			}
		}
		  
		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_edgePath(..): Simulated " << forReturn.size() << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	return forReturn;
}

// see header
size_t readSimulator::simulate_paired_reads_from_string(std::string readNamePrefix, std::string& s, double expected_haploid_coverage, std::vector<std::pair<std::ofstream*, std::ofstream*>>& output_FHs_perThread, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd)
{
	std::vector<oneReadPair> forReturn;

	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2
//	std::cout << "Poisson start rate: " << poissonStartRate << "\n";

	long long firstPosition = 0;
	long long lastPosition = s.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_reads_from_string(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	omp_set_num_threads(threads);

	std::vector< std::vector<oneReadPair> > forReturn_perThread;
	forReturn_perThread.resize(threads);

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;

	auto print_one_readPair = [&output_FHs_perThread] (oneReadPair& rP, unsigned int threadI) -> void
	{
		*(output_FHs_perThread.at(threadI).first) << "@" << rP.reads.first.name << "\n";
		*(output_FHs_perThread.at(threadI).first) << rP.reads.first.sequence << "\n";
		*(output_FHs_perThread.at(threadI).first) << "+" << "\n";
		*(output_FHs_perThread.at(threadI).first) << rP.reads.first.quality << "\n";

		*(output_FHs_perThread.at(threadI).second) << "@" << rP.reads.second.name << "\n";
		*(output_FHs_perThread.at(threadI).second) << rP.reads.second.sequence << "\n";
		*(output_FHs_perThread.at(threadI).second) << "+" << "\n";
		*(output_FHs_perThread.at(threadI).second) << rP.reads.second.quality << "\n";
	};

	#pragma omp parallel
	{
		assert(omp_get_num_threads() == (int)threads);
		int thisThread = omp_get_thread_num();

		boost::mt19937 rnd_gen;

		auto seed = boost::random::random_device()();
		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);
			if(generateError)
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;

		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= s.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < s.size());

				sampleOneBase(base, s.at(index_into_baseString), base_for_read, quality_for_read);
				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;

				index_into_baseString++;
			}

			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		#pragma omp for
		for(long long i = 0; i < lastPosition; i++)
		{

			assert(omp_get_num_threads() == (int)threads);
			assert(thisThread == omp_get_thread_num());

			if((i % 10000) == 0)
			{
				#pragma omp critical
				{
					if(omp_get_thread_num() == 0)
					{
						double approx_all = thread_generated_bases * omp_get_num_threads();
						std::cout << "\r" << "Generated bases this thread: " <<  thread_generated_bases << " x total: " << approx_all << std::flush;
					}
				}
			}
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));
				

				std::string read1; std::string read1_qualities; bool read1_success;
				std::string read2; std::string read2_qualities; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_success);
				sampleRead(i + jumpSize, read2, read2_qualities, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					
					
					std::string read1_name = readNamePrefix + readName_field_separator + "p1" + readName_field_separator + Utilities::ItoStr(i);
					std::string read2_name = readNamePrefix + readName_field_separator + "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize);

					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					oneReadPair rP(r1, r2, jumpSize);

					print_one_readPair(rP, thisThread);

					// forReturn_perThread.at(thisThread).push_back(rP);
				}
			}
		}

		#pragma omp critical
		{
//			std::cout << "Thread " << thisThread << "; read pairs " << thread_read_pairs << "; thread bases: " << thread_generated_bases << "; thread errors: " << thread_generated_errors << "; thread INDELs: " << thread_indel_events << "\n" << std::flush;

			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}
	
	std::cout << "\n";


//	std::cout << "readSimulator::simulate_paired_reads_from_string() summary\n";
//	std::cout << "\tUnderlying string length: " << s.length() << "\n";
//	std::cout << "\tGenerated bases: " << global_generated_bases << "\n";
//	std::cout << "\tGenerated errors: " << global_generated_errors << "\n";
//	std::cout << "\tGenerated INDEL events: " << global_indel_events << "\n";
//	std::cout << "\tAchieved approximate coverage: " << global_generated_bases/double(s.length()) << "\n\n";

	return global_generated_bases;

//	for(unsigned int tI = 0; tI < threads; tI++)
//	{
//		forReturn.insert(forReturn.end(), forReturn_perThread.at(tI).begin(), forReturn_perThread.at(tI).end());
//	}
//
//	for(unsigned int rI = 0; rI < forReturn.size(); rI++)
//	{
//		oneReadPair& p = forReturn.at(rI);
//		std::cout << "Read pair #" << rI << "\n======================================\n";
//		std::cout << p.reads.first.name << "\n";
//		std::cout << p.reads.first.sequence << "\n";
//		std::cout << p.reads.first.quality << "\n";
//		std::cout << p.diff_starting_coordinates << "\n";
//		std::cout << p.reads.second.name << "\n";
//		std::cout << p.reads.second.sequence << "\n";
//		std::cout << p.reads.second.quality << "\n\n";
//	}

//	return forReturn;
}
