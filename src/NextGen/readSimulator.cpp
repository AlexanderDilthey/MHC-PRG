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

#include "readSimulator.h"
#include "../Utilities.h"

#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>

std::string readName_field_separator = "|||";

readSimulator::readSimulator(std::string qualityMatrixFile, unsigned int readLength) {

	std::ifstream matrixStream;
	matrixStream.open (qualityMatrixFile.c_str(), std::ios::in);

	read_quality_frequencies.resize(readLength);
	read_quality_correctness.resize(readLength);

	read_INDEL_freq.resize(readLength, 0);
	std::vector<double> read_nonINDEL_freq;
	read_nonINDEL_freq.resize(readLength, 0);

	if(matrixStream.is_open())
	{
		std::string line;
		std::string currentIdentifier;
		std::string currentSequence;

		size_t lineCounter = 0;
		bool found_read_length = false;
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
				if(Utilities::StrtoI(line_fields.at(0)) == (int)readLength)
				{
					found_read_length = true;
				}
				else
				{
					continue;
				}

				char quality = line_fields.at(1).at(0);
				unsigned int positionInRead = Utilities::StrtoI(line_fields.at(2));
				int datapoints = Utilities::StrtoI(line_fields.at(3));

				assert(positionInRead < readLength);


				if(quality == 0)
				{
					read_INDEL_freq.at(positionInRead) += datapoints;
				}
				else
				{
					read_nonINDEL_freq.at(positionInRead) += datapoints;

					if(read_quality_frequencies.at(positionInRead).count(quality) == 0)
					{
						read_quality_frequencies.at(positionInRead)[quality] = 0;
					}
					read_quality_frequencies.at(positionInRead)[quality] += datapoints;

					if(!(read_quality_correctness.at(positionInRead).count(quality) == 0))
					{
						throw std::runtime_error("readSimulator::readSimulator(): file "+qualityMatrixFile+", double-defined quality. Line "+Utilities::ItoStr(lineCounter)+".\n"+line+"\n");
					}
					assert(read_quality_correctness.at(positionInRead).count(quality) == 0);
					read_quality_correctness.at(positionInRead)[quality] += Utilities::StrtoD(line_fields.at(5));
				}
			}
		}
		matrixStream.close();

		if(! found_read_length)
		{
			throw std::runtime_error("Cannot find the right read length "+Utilities::ItoStr(readLength)+" in file "+qualityMatrixFile);
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
				read_INDEL_freq.at(i) = 1e-8;
			}
		}
	}
	else
	{
		throw std::runtime_error("readSimulator::readSimulator(): Cannot open file "+qualityMatrixFile);
	}

	read_length = readLength;
	threads = 40;
	paranoid = true;
}


std::vector<oneReadPair> readSimulator::simulate_paired_reads_from_edgePath(std::vector<Edge*> edgePath, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd)
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

						coordinates_string.push_back(-1);
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
				coordinates_string.push_back(index_into_baseString);

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

		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));


				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);
				sampleRead(i + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);

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

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;

					oneReadPair rP(r1, r2, jumpSize);

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
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n";
	std::cout << "\n" << std::flush;

	return forReturn;
}

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
