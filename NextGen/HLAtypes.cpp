/*
 * HLAtypes.cpp
 *
 *  Created on: 20.06.2014
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
#include <limits>
#include <unordered_map>

#include "HLAtypes.h"


#include "readSimulator.h"
#include "../Utilities.h"
#include <boost/math/distributions/poisson.hpp>

#include "../hash/deBruijn/DeBruijnGraph.h"
#include "../GraphAligner/GraphAlignerAffine.h"
#include "../GraphAlignerUnique/GraphAlignerUnique.h"

#include <stdio.h>
#include "dirent.h"
#include <boost/math/distributions/chi_squared.hpp>


// constants


void fillAS();

double insertionP = 0.001;  
double deletionP = 0.001;
double log_likelihood_insertion = log(insertionP);
double log_likelihood_insertion_actualAllele = log_likelihood_insertion + log(1.0/4.0);

double log_likelihood_deletion = log(deletionP);
double log_likelihood_nonInsertion = log(1 - insertionP);
double log_likelihood_nonDeletion = log(1 - deletionP);
double log_likelihood_match_mismatch = log(1 - insertionP - deletionP);

// int max_mismatches_perRead = 2;
double min_alignmentFraction_OK = 0.96; // measures all alignment positions but graph AND sequence gaps, separately for both reads
double min_oneRead_weightedCharactersOK = 0.995; // one read, mismatches downweighted by quality
// double min_bothReads_weightedCharactersOK = 0.985; // both reads, mismatches downweighted by quality
// double min_bothReads_weightedCharactersOK = 0.95; // todo reinstate
double min_bothReads_weightedCharactersOK = 0.0;

double minimumMappingQuality = 0.0;
double minimumPerPositionMappingQuality = 0.7;

bool combineReadAndBaseLikelihoods = false;


using namespace boost::math::policies;
using namespace boost::math;

typedef boost::math::poisson_distribution< double, policy < discrete_quantile < integer_round_inwards> > > poisson_up;
poisson_up poisson1(1);

bool veryConservativeReadLikelihoods = true;

std::map<std::string, std::string> alleles_to_G;
std::set<std::string> G_loci;

std::map<std::string, std::string> codon2AS;
std::map<std::string, std::vector<std::string> > AS2codon;

/*
double logAvg(double a, double b)
{
	if(a > b)
	{
		return(log(0.5) + (log(1 + exp(b - a)) + a));
		
	} else
	{
		return(log(0.5) + (log(1 + exp(a - b)) + b));
	}
}
*/

class oneExonPosition {
public:
	unsigned int positionInExon;
	int graphLevel;
	std::string genotype;
	std::string alignment_edgelabels;
	std::string qualities;

	std::string thisRead_ID;
	double thisRead_fractionOK;
	double thisRead_WeightedCharactersOK;

	std::string pairedRead_ID;
	double pairedRead_fractionOK;
	double pairedRead_WeightedCharactersOK;

	std::string read1_ID;

	bool pairs_strands_OK;
	double pairs_strands_distance;

	double mapQ;
	double mapQ_genomic;
	double mapQ_position;
};


void oneReadAlignment_2_exonPositions_paired(seedAndExtend_return_local& alignment, oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, seedAndExtend_return_local& paired_alignment, oneRead& paired_read, int read_1_or_2, const std::vector<int>& combined_exon_sequences_graphLevels, const std::map<int, unsigned int>& graphLevel_2_exonPosition);
void oneReadAlignment_2_exonPositions_unpaired(seedAndExtend_return_local& alignment, oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, const std::vector<int>& combined_exon_sequences_graphLevels, const std::map<int, unsigned int>& graphLevel_2_exonPosition);

auto countMismatchesInExon(std::vector<oneExonPosition>& exonPositions) -> int {
	int forReturn = 0;
	for(unsigned int i = 0; i < exonPositions.size(); i++)
	{
		if(exonPositions.at(i).genotype != exonPositions.at(i).alignment_edgelabels)
		{
			forReturn++;
		}
	}
	return forReturn;
};


double alignmentFractionOK (seedAndExtend_return_local& r)
{
	int positions_OK = 0;
	int positions_checked = 0;

	for(unsigned int i = 0; i < r.graph_aligned.size(); i++)
	{
		if((r.graph_aligned.at(i) == '_') && (r.sequence_aligned.at(i) == '_'))
		{
			continue;
		}

		positions_checked++;
		if(r.graph_aligned.at(i) == r.sequence_aligned.at(i))
		{
			positions_OK++;
		}
	}
	return double(positions_OK)/double(positions_checked);
};

bool debugOpenMP = false;
auto convertLogVectorToP = [](std::vector<double>& v) -> void {
	if((v.size() > 10000) && (omp_get_max_threads() > 1))
	{
		std::vector<double> v2;
	
		if(debugOpenMP)
		{
			v2 = v;
			double v_max = 0;
			for(unsigned int i = 0; i < v2.size(); i++)
			{
				if((i == 0) || (v_max < v2.at(i)))
				{
					v_max = v2.at(i);
				}
			}
			for(unsigned int i = 0; i < v2.size(); i++)
			{
				v2.at(i) = exp(v2.at(i) - v_max);
				assert(v2.at(i) >= 0);
				assert(v2.at(i) <= 1);
			}		
		}	
		
		
		{
			int T = omp_get_max_threads();
			int chunk_size = v.size() / T;
		
			std::vector<double> v_max_pT; 
			v_max_pT.resize(omp_get_max_threads(), 0);
					
			#pragma omp parallel
			{
				assert(omp_get_num_threads() == T);
			
				int thisThread = omp_get_thread_num();
				int firstState = thisThread * chunk_size;
				int lastState = (thisThread+1) * chunk_size - 1;
				if(thisThread == (omp_get_num_threads()-1))
				{
					lastState = (v.size()-1);
				}
				
				for(int i = firstState; i <= lastState; i++)
				{
					if((i == firstState) || (v_max_pT.at(thisThread) < v.at(i)))
					{
						v_max_pT.at(thisThread) = v.at(i);
					}			
				}
			}
			
			double v_max = 0;
			for(unsigned int i = 0; i < v_max_pT.size(); i++)
			{
				if((i == 0) || (v_max < v_max_pT.at(i)))
				{
					v_max = v_max_pT.at(i);
				}
			}

			#pragma omp parallel for
			for(unsigned int i = 0; i < v.size(); i++)
			{
				v.at(i) = exp(v.at(i) - v_max);
				assert(v.at(i) >= 0);
				assert(v.at(i) <= 1);
			}
			
			if(debugOpenMP)
			{
				for(unsigned int i = 0; i < v.size(); i++)
				{
					if(!(abs(v.at(i) - v2.at(i)) < 1e-5))
					{
						#pragma omp critical
						{					
							std::cerr << "Problem 1!\n";
							std::cerr << "i: " << i << "\n";
							std::cerr << "v.at(i): " << v.at(i)  << "\n";
							std::cerr << "v2.at(i): " << v2.at(i)  << "\n";
							std::cerr << "\n" << std::flush;
						}
					}
					
					assert(abs(v.at(i) - v2.at(i)) < 1e-5);
				}			
				// std::cout << "OK1" << std::flush;
			}
		}

	}
	else
	{
		double v_max = 0;
		for(unsigned int i = 0; i < v.size(); i++)
		{
			if((i == 0) || (v_max < v.at(i)))
			{
				v_max = v.at(i);
			}
		}
		for(unsigned int i = 0; i < v.size(); i++)
		{
			v.at(i) = exp(v.at(i) - v_max);
			assert(v.at(i) >= 0);
			assert(v.at(i) <= 1);
		}
	}
};

auto normalizeVector = [](std::vector<double>& v) -> void {
	if((v.size() > 10000) && (omp_get_max_threads() > 1))
	{
		std::vector<double> v2;
		
		if(debugOpenMP)
		{
			v2 = v;
			double v_sum = 0;
			for(unsigned int i = 0; i < v2.size(); i++)
			{
				v_sum += v2.at(i);
			}
			assert(v_sum > 0);
			for(unsigned int i = 0; i < v2.size(); i++)
			{
				v2.at(i) = v2.at(i) / v_sum;
			}			
		}
		
		{
			int T = omp_get_max_threads();
			std::vector<double> sums; 
			sums.resize(T, 0);
			
			size_t M = v.size();
			#pragma omp parallel for
			for(unsigned int i = 0; i < M; i++)
			{
				assert(omp_get_num_threads() == T);
				sums.at(omp_get_thread_num()) += v.at(i);
			}
			
			double v_sum = 0;
			for(int i = 0; i < (int)sums.size(); i++)
			{
				v_sum += sums.at(i);
			}

			assert(v_sum > 0);
			#pragma omp parallel for		
			for(unsigned int i = 0; i < v.size(); i++)
			{
				v.at(i) = v.at(i) / v_sum;				
			}
			
			if(debugOpenMP)
			{
				for(unsigned int i = 0; i < v.size(); i++)
				{
					if(!(abs(v.at(i) - v2.at(i)) < 1e-5))
					{
						#pragma omp critical
						{
							std::cerr << "Problem 2!\n";
							std::cerr << "i: " << i << "\n";
							std::cerr << "v.at(i): " << v.at(i)  << "\n";
							std::cerr << "v2.at(i): " << v2.at(i)  << "\n";
							std::cerr << "\n" << std::flush;
						}
					}
					assert(abs(v.at(i) - v2.at(i)) < 1e-5);
				}			
				// std::cout << "OK2" << std::flush;
			}			
		}
	}
	else
	{
		double v_sum = 0;
		for(unsigned int i = 0; i < v.size(); i++)
		{
			v_sum += v.at(i);
		}
		assert(v_sum > 0);
		for(unsigned int i = 0; i < v.size(); i++)
		{
			v.at(i) = v.at(i) / v_sum;
		}
	}
};

std::vector<oneExonPosition> removeDoublePositionsFromRead(const std::vector<oneExonPosition>& positions);

std::string printPerc (double v1, double v2)
{
	double perc = (v1/v2) * 100;
	return Utilities::DtoStr(perc);
};

std::pair<double, double> meanMedian (std::vector<double> L)
{
	std::sort(L.begin(), L.end());
	double S = 0;
	for(unsigned int i = 0; i < L.size(); i++)
	{
		S += L.at(i);
		// std::cout << L.at(i) << " ";
	}
	double mean = 0;
	double median = 0;
	if(L.size() > 0)
	{
		mean = S / (double)L.size();
		unsigned int medium_index = L.size() / 2;
		median = L.at(medium_index);
	}
	return make_pair(mean, median);
};



void fill_loci_2_exons(std::map<std::string, std::vector<std::string> >& loci_2_exons)
{
	std::vector<std::string> exons_A = {"exon_2", "exon_3"};
	loci_2_exons["A"] = exons_A;

	std::vector<std::string> exons_B = {"exon_2", "exon_3"};
	loci_2_exons["B"] = exons_B;

	std::vector<std::string> exons_C = {"exon_2", "exon_3"};
	loci_2_exons["C"] = exons_C;

	std::vector<std::string> exons_DQA1 = {"exon_2"};
	loci_2_exons["DQA1"] = exons_DQA1;

	std::vector<std::string> exons_DQB1 = {"exon_2"};
	loci_2_exons["DQB1"] = exons_DQB1;

	std::vector<std::string> exons_DRB1 = {"exon_2"};
	loci_2_exons["DRB1"] = exons_DRB1;

	loci_2_exons["DPA1"] = {"exon_2"};
	loci_2_exons["DPB1"] = {"exon_2"};
	loci_2_exons["DRA"] = {"exon_2"};
	loci_2_exons["DRB3"] = {"exon_2"};
	loci_2_exons["DRB4"] = {"exon_2"};

	loci_2_exons["E"] = {"exon_2", "exon_3"};
	loci_2_exons["F"] = {"exon_2", "exon_3"};
	loci_2_exons["G"] = {"exon_2", "exon_3"};
	loci_2_exons["H"] = {"exon_2", "exon_3"};
	loci_2_exons["J"] = {"exon_2", "exon_3"};
	loci_2_exons["K"] = {"exon_2", "exon_3"};
	loci_2_exons["L"] = {"exon_2", "exon_3"};
	loci_2_exons["V"] = {"exon_2", "exon_3"};


}

double alignmentWeightedOKFraction(oneRead& underlyingRead, seedAndExtend_return_local& alignment)
{
	int indexIntoOriginalReadData = -1;

	int totalMismatches = 0;
	double weightedMismatches = 0;

	for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
	{
		std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
		std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
		int graphLevel = alignment.graph_aligned_levels.at(cI);

		if(sequenceCharacter != "_")
		{
			indexIntoOriginalReadData++;
			int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
			if(alignment.reverse)
			{
				indexIntoOriginalReadData_correctlyAligned = underlyingRead.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
			}
			assert(indexIntoOriginalReadData_correctlyAligned >= 0);
			assert(indexIntoOriginalReadData_correctlyAligned < (int)underlyingRead.sequence.length());;

			std::string underlyingReadCharacter = underlyingRead.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
			if(alignment.reverse)
			{
				underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
			}
			assert(underlyingReadCharacter == sequenceCharacter);

			if(graphCharacter == "_")
			{
				totalMismatches++;
				weightedMismatches++;
			}
			else
			{
				// two well-defined characters
				char qualityCharacter = underlyingRead.quality.at(indexIntoOriginalReadData_correctlyAligned);
				double pCorrect = Utilities::PhredToPCorrect(qualityCharacter);
				assert((pCorrect >= 0) && (pCorrect <= 1));
				if(sequenceCharacter == graphCharacter)
				{
					// match!
				}
				else
				{
					weightedMismatches += pCorrect;
					totalMismatches++;
				}
			}

		}
		else
		{
			assert(sequenceCharacter == "_");
			if(graphCharacter == "_")
			{
				// sequence gap, graph gap - likelihood 1
				assert(graphLevel != -1);
			}
			else
			{
				// sequence gap, graph non gap - deletion in sequence
				totalMismatches++;
				weightedMismatches++;
			}
		}
	}

	double readLength = underlyingRead.sequence.length();
	assert(totalMismatches >= weightedMismatches);

	return  (1.0 - (weightedMismatches / readLength));
};


double read_likelihood_per_position(const std::string& exonGenotypeR, const std::string& readGenotypeR, const std::string& readQualities, const int& alignmentGraphLevel)
{
	// not exact, only approximation, but seems to work well empirically

	double log_likelihood = 0;
	bool verbose = 0;

	// if(!(exonGenotype.length() == 1))
	// {
		// std::cerr << "! (exonGenotype.length() == 1)" << "\n";
		// std::cerr << exonGenotype.length()  << "\n";
		// std::cerr << exonGenotype << "\n" << std::flush;
	// }
	// 
	
	assert(exonGenotypeR.length() >= 1);
	assert(readGenotypeR.length() >= 1);
	
	for(unsigned int i = 0; i < exonGenotypeR.length(); i++)
	{
		std::string exonGenotype = exonGenotypeR.substr(i, 1);
		
		std::string readGenotype;
		if(i == (exonGenotypeR.length() - 1))
		{
			if(i < readGenotypeR.length())
			{
				int extractLength = (int)readGenotypeR.length() - (int)i;
				assert(extractLength > 0);
				readGenotype = readGenotypeR.substr(i, extractLength);
				assert((int)readGenotype.length() == extractLength);
			}
			else
			{
				readGenotype = "_";
			}
		}
		else
		{
			if(i < readGenotypeR.length())
			{
				readGenotype = readGenotypeR.substr(i, 1);
			}
			else
			{
				readGenotype = "_";			
			}
		}
		
		assert(exonGenotype.length() == 1);		
		assert(readGenotype.length() >= 1);

		unsigned int l_diff = readGenotype.length() - exonGenotype.length();
		if(exonGenotype == "_")
		{
			// assert(l_diff == 0);
			if(readGenotype == "_")
			{
				assert(alignmentGraphLevel != -1);
				// likelihood 1 - intrinsic graph gap

				if(verbose)
				{
					std::cout << "\t\t" << "Intrinsic graph gap" << "\n";
				}
				
				// log_likelihood += log_likelihood_nonInsertion;

			}
			else
			{
				if(verbose)
				{
					std::cout << "\t\t" << "Insertion " << (1 + l_diff) << "\n";
				}
				
				assert(l_diff >= 0);	
				double p_pdf = pdf(poisson1, l_diff);
				double log_from_pdf;
				if(p_pdf == 0)
				{
					log_from_pdf = log(1e-30);
				}
				else
				{
					log_from_pdf = log(p_pdf);
				}
				
				log_likelihood += (log_likelihood_insertion_actualAllele + log_from_pdf);
			}
		}
		else
		{

			// score from first position match
			if(readGenotype.substr(0, 1) == "_")
			{
				log_likelihood += log_likelihood_deletion;

				if(verbose)
				{
					std::cout << "\t\t" << "Deletion" << "\n";
				}
			}
			else
			{
				if(readGenotype.length() > 1)
				{
					std::string readGenotype_after1 = readGenotype.substr(1);
					assert(readGenotype_after1.find("_") == std::string::npos);
				}
				
				log_likelihood += log_likelihood_match_mismatch;
				
				assert(readQualities.length());
				double pCorrect = Utilities::PhredToPCorrect(readQualities.at(0));
				if(veryConservativeReadLikelihoods)
				{
					if(pCorrect > 0.999)
						pCorrect = 0.999;
				}
				assert((pCorrect >= 0) && (pCorrect <= 1));

				if(pCorrect == 0)
				{
					pCorrect = 0.001;
				}
				
				// if(!(pCorrect >= 0.25))
				// {
					// std::cerr << "pCorrect = " << pCorrect << "\n" << std::flush;
				// }
				// assert(pCorrect >= 0.25);

				if(exonGenotype == readGenotype.substr(0, 1))
				{
					if(verbose)
					{
						std::cout << "\t\t" << "Match " << pCorrect << "\n";
					}

					log_likelihood += log(pCorrect);
				}
				else
				{
					double pIncorrect = (1 - pCorrect)*(1.0/3.0);
					assert(pIncorrect <= 0.75);
					assert((pIncorrect > 0) && (pIncorrect < 1));
					log_likelihood += log(pIncorrect);

					if(verbose)
					{
						std::cout << "\t\t" << "Mismatch " << pIncorrect << "\n";
					}
				}
			}
			
			// if read allele is longer
			if(l_diff == 0)
			{
				// log_likelihood += log_likelihood_nonInsertion;
			}
			else
			{
				assert(l_diff >= 1);
				double p_pdf = pdf(poisson1, l_diff - 1);
				double log_from_pdf;
				if(p_pdf == 0)
				{
					log_from_pdf = log(1e-30);
				}
				else
				{
					log_from_pdf = log(p_pdf);
				}
				log_likelihood += (log_likelihood_insertion_actualAllele + log_from_pdf);				
			}

			if(l_diff > 0)
			{
				if(verbose)
				{
					std::cout << "\t\t" << "Insertion " << l_diff << "\n";
				}
			}
		}		
	}


	if(verbose)
	{
		std::cout << "\t\t" << "Running log likelihood: " << log_likelihood << "\n";
	}

	double l = exp(log_likelihood);
	
	if(!((l >= 0) && (l <= 1)))  
	{
		std::cerr << "l" << ": " << l << "\n";
		std::cerr << "exonGenotypeR" << ": " << exonGenotypeR << "\n";
		std::cerr << "readGenotypeR" << ": " << readGenotypeR << "\n";
		std::cerr << "readQualities" << ": " << readQualities << "\n";
		std::cerr << "alignmentGraphLevel" << ": " << alignmentGraphLevel << "\n";
		
		std::cerr << std::flush;
	}
	
	assert(l >= 0);
	assert(l <= 1);
	
	return log_likelihood;
}

void simulateHLAreads_perturbHaplotype(std::vector<std::string>& haplotype, std::vector<int>& perturbedPositions)
{
	// std::cout << "simulateHLAreads_perturbHaplotype:\n";
	// std::cout << "\t" << "haplotype.size()" << ": " << haplotype.size() << "\n";
	// std::cout << "\t" << "perturbedPositions.size()" << ": " << perturbedPositions.size() << "\n";
	// std::cout << std::flush;
	
	perturbedPositions.clear();
	
	double rate_perturbation = 0.01;
	for(unsigned int pI = 0; pI < haplotype.size(); pI++)
	{
		// std::cout << pI << " / " << haplotype.size() << "\n";
		// std::cout << "\t" << haplotype.at(pI).size() << "\n";
		// std::cout << std::flush;
		
		bool madePerturbation = false;
		if(Utilities::randomDouble() >= (1 - rate_perturbation))
		{
			int type_of_event = Utilities::randomNumber(100);

			if(type_of_event <= 90)
			{
				// SNV
				std::string newAllele;
				newAllele.push_back(Utilities::randomNucleotide());
				if(haplotype.at(pI) != newAllele)
				{
					madePerturbation = true;
				}
				haplotype.at(pI) = newAllele;
			}
			else if((type_of_event > 90) && (type_of_event <= 95))
			{
				// deletion
				if(haplotype.at(pI) != "_")
				{
					madePerturbation = true;
				}				
				haplotype.at(pI) = "_";
			}
			else
			{
				// insertion
				std::string newAllele = haplotype.at(pI);
				int length = Utilities::randomNumber(1)+1;
				// std::cout << "Length: " << length << "\n";
				std::cout << std::flush;
				for(int l = 0; l <= length; l++)
				{
					newAllele.push_back(Utilities::randomNucleotide());
				}
				haplotype.at(pI) = newAllele;
				madePerturbation = true;
			}
		}
		
		if(madePerturbation)
		{
			perturbedPositions.push_back(pI);
		}
	}
}

void simulateHLAreads(std::string graphDir, int nIndividuals, bool exon23, bool perturbHaplotypes, bool readError, std::string outputDirectory, std::string qualityMatrixFile, int readLength, double insertSize_mean, double insertSize_sd, double haploidCoverage)
{
	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	if(! Utilities::directoryExists(outputDirectory))
	{
		Utilities::makeDir(outputDirectory);
	}

	std::set<std::string> loci;
	std::map<std::string, std::vector<std::string> > files_per_locus;
	std::map<std::string, std::vector<std::string> > files_per_locus_type;
	std::map<std::string, std::vector<std::string> > files_per_locus_number;
	
	// loci_2_exons
	std::map<std::string, std::vector<std::string> > loci_2_exons;
	fill_loci_2_exons(loci_2_exons);
	
	// find files and loci

	std::vector<std::string> files_in_order;
	std::vector<std::string> files_in_order_type;
	std::vector<int> files_in_order_number;

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			std::string file_locus = split_by_underscore.at(0);
			if(file_locus == "before")
			{
				continue;
			}
			if(file_locus == "after")
			{
				continue;
			}
			loci.insert(file_locus);

			std::string filePath = graphDir + "/" + line;

			files_per_locus[file_locus].push_back(filePath);
			files_per_locus_type[file_locus].push_back(split_by_underscore.at(2));
			
			if(split_by_underscore.size() >= 4)
			{
				std::vector<std::string> split_by_colon = Utilities::split(split_by_underscore.at(3), ".");
				assert(split_by_colon.size() == 2);				
				files_per_locus_number[file_locus].push_back(split_by_colon.at(0));
			}
			else
			{
				files_per_locus_number[file_locus].push_back("");
			}
			
		}
	}
	assert(loci.size() > 0);


	std::cout << "simulateHLAreads(..): Found " << loci.size() << " loci.\n";

	// find available types

	std::map<std::string, std::vector<std::string> > loci_availableTypes;
	for(std::set<std::string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
	{
		std::string locus = *locusIt;

		std::string arbitraryIntronFile;
		if(exon23 && loci_2_exons.count(locus))
		{				
			for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
			{
				if((files_per_locus_type.at(locus).at(fI) == "exon") && (files_per_locus_number.at(locus).at(fI) == "2"))
				{
					arbitraryIntronFile = files_per_locus.at(locus).at(fI);
					break;
				}
			}		
		}
		else
		{
			for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
			{
				if(files_per_locus_type.at(locus).at(fI) == "intron")
				{
					arbitraryIntronFile = files_per_locus.at(locus).at(fI);
					break;
				}
			}
			
			if(arbitraryIntronFile.length() == 0)
			{
				for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
				{
					if(files_per_locus_type.at(locus).at(fI) == "exon")
					{
						arbitraryIntronFile = files_per_locus.at(locus).at(fI);
						break;
					}
				}		
			}		
		}
		
		if(!(arbitraryIntronFile.length()))
		{
			std::cerr << "Cannot find an arbitrary intron/exon file for locus " << locus << "\n" << std::flush;
		}
		assert(arbitraryIntronFile.length());
		
		std::vector<std::string> availableTypes;

		std::ifstream fileInputStream;
		fileInputStream.open(arbitraryIntronFile.c_str());
		if(!fileInputStream.is_open())
		{
			std::cerr << "Cannot open file " << arbitraryIntronFile << "\n" << std::flush;
		}
		assert(fileInputStream.is_open());
		unsigned int lI = 0;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			std::vector<std::string> line_fields = Utilities::split(line, " ");
			if(lI == 0)
			{
				if(!(line_fields.size() > 0))
				{
					std::cerr << "First line of file weird: " << arbitraryIntronFile << "\n" << std::flush;
				}
				assert(line_fields.size() > 0);
				assert(line_fields.at(0) == "IndividualID");
			}
			else
			{
				if(line.length())
				{
					std::string type = line_fields.at(0);
					availableTypes.push_back(type);
				}
			}
			lI++;
		}
		fileInputStream.close();

		loci_availableTypes[locus] = availableTypes;
	}

	// find type haplotypes

	std::map<std::string, std::set<std::string> > removeTypes;

	std::map<std::string, std::map<std::string, std::vector<std::string>> > loci_types_haplotypes;
	std::map<std::string, std::vector<std::string> > loci_graphLevelIDs;

	for(std::set<std::string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
	{
		std::string locus = *locusIt;
		std::set<std::string> availableTypes_set(loci_availableTypes.at(locus).begin(), loci_availableTypes.at(locus).end());

		// std::cerr << "Locus: " << locus << "\n";
		
		for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
		{
			std::ifstream fileInputStream;
			fileInputStream.open(files_per_locus.at(locus).at(fI).c_str());
			assert(fileInputStream.is_open());
			std::vector<std::string> file_lines;
			while(fileInputStream.good())
			{
				std::string line;
				std::getline(fileInputStream, line);
				Utilities::eraseNL(line);
				file_lines.push_back(line);
			}
			fileInputStream.close();

			std::string firstLine = file_lines.at(0);
			std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
			assert(firstLine_fields.at(0) == "IndividualID");

			std::vector<std::string> graph_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
			loci_graphLevelIDs[locus].insert(loci_graphLevelIDs[locus].end(), graph_level_names.begin(), graph_level_names.end());

			// std::cerr << "\tfI = " << fI << "\n";
			if(exon23 && loci_2_exons.count(locus))
			{
				
				std::cout << "Locus " << locus << " in exon23 mode!\n";
				
				bool relevantExon = false;

				if(files_per_locus_type.at(locus).at(fI) == "exon")
				{
					
					if(loci_2_exons.at(locus).size() == 1)
					{
						assert(loci_2_exons.at(locus).at(0) == "exon_2");
						if(files_per_locus_number.at(locus).at(fI) == "2")
						{
							relevantExon = true;
						}
					}
					
					if(loci_2_exons.at(locus).size() == 2)
					{
						assert(loci_2_exons.at(locus).at(0) == "exon_2");
						assert(loci_2_exons.at(locus).at(1) == "exon_3");
						
						if(files_per_locus_number.at(locus).at(fI) == "2")
						{
							relevantExon = true;
						}
						
						if(files_per_locus_number.at(locus).at(fI) == "3")
						{
							relevantExon = true;
						}										
					}
				}
					
							
				std::cout << "\tFile " << files_per_locus.at(locus).at(fI) << " relevant: " << relevantExon << "\n";
							
				for(unsigned int lI = 1; lI < file_lines.size(); lI++)
				{
					if(file_lines.at(lI).length())
					{
						std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
						assert(line_fields.size() == firstLine_fields.size());
						std::string HLA_type = line_fields.at(0);
						std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

						if(relevantExon)
						{
							if(availableTypes_set.count(HLA_type))
							{
								if(loci_types_haplotypes[locus].count(HLA_type) == 0)
								{
									loci_types_haplotypes[locus][HLA_type].resize(0);
								}

								loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());
								
								// std::cout << "\t\t individual exon for type: " << HLA_type << "\n";

												
								// std::cerr << "\t\t" << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " / " << loci_graphLevelIDs[locus].size() << "\n";			
							}
						}
						else
						{
							bool noStars = true;
							for(unsigned int aI = 0; aI < line_alleles.size(); aI++)
							{
								std::string& S = line_alleles.at(aI);
								if(S.find("*") != std::string::npos)
								{
									noStars = false;
									break;
								}
							}
							
							if(lI == (file_lines.size() - 1))
							{
								noStars = true;								
								std::cerr << "Locus " << locus << " file " << files_per_locus.at(locus).at(fI) << ": No line without stars -- force " << lI << " / " << file_lines.size() << "\n" << std::flush;
							}
							
							if(noStars)
							{
								for(std::set<std::string>::iterator typeIt = availableTypes_set.begin(); typeIt != availableTypes_set.end(); typeIt++)
								{
									std::string HLA_type = *typeIt;
									if(loci_types_haplotypes[locus].count(HLA_type) == 0)
									{
										loci_types_haplotypes[locus][HLA_type].resize(0);
									}

									loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());																	
									
									// std::cout << "\t\t identical sequence for type: " << HLA_type << " (line " << lI << ")\n";
								}	

								break;
							}
						}
					}
				}
			}
			else
			{								
				for(unsigned int lI = 1; lI < file_lines.size(); lI++)
				{
					if(file_lines.at(lI).length())
					{
						std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
						assert(line_fields.size() == firstLine_fields.size());
						std::string HLA_type = line_fields.at(0);
						std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

						if((files_per_locus_type.at(locus).at(fI) == "intron") || (files_per_locus_type.at(locus).at(fI) == "exon"))
						{
							if(availableTypes_set.count(HLA_type))
							{
								if(loci_types_haplotypes[locus].count(HLA_type) == 0)
								{
									loci_types_haplotypes[locus][HLA_type].resize(0);
								}

								loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());
									
								// std::cerr << "\t\t" << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " / " << loci_graphLevelIDs[locus].size() << "\n";			
							}
						}
						else
						{
							// this must be a pre- or post-padding sequence
							assert((files_per_locus_type.at(locus).at(fI) == "paddingLeft.txt") || (files_per_locus_type.at(locus).at(fI) == "paddingRight.txt"));
							
							bool noStars = true;
							for(unsigned int aI = 0; aI < line_alleles.size(); aI++)
							{
								std::string& S = line_alleles.at(aI);
								if(S.find("*") != std::string::npos)
								{
									noStars = false;
									break;
								}
							}
							
							if(lI == (file_lines.size() - 1))
							{
								noStars = true;								
								std::cerr << "Locus " << locus << " file " << files_per_locus.at(locus).at(fI) << ": No line without stars -- force " << lI << " / " << file_lines.size() << "\n" << std::flush;
							}
							
							// std::cout << "TYPE: " << files_per_locus_type.at(locus).at(fI) << "\n";
							if(noStars)
							{
								for(std::set<std::string>::iterator typeIt = availableTypes_set.begin(); typeIt != availableTypes_set.end(); typeIt++)
								{
									std::string HLA_type = *typeIt;
									if(loci_types_haplotypes[locus].count(HLA_type) == 0)
									{
										loci_types_haplotypes[locus][HLA_type].resize(0);
									}

									loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());
									
									// std::cerr << "\t\t" << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " / " << loci_graphLevelIDs[locus].size() << "\n";								
								}
								break;
							}
						}
					}
				}
			}
		}
	}

	// check that all haplotypes make sense

	for(std::set<std::string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
	{
		std::string locus = *locusIt;
		
		std::cout << "Locus " << locus << "\n";
		
		assert(loci_availableTypes.at(locus).size() > 0);
		for(unsigned int tI = 0; tI < loci_availableTypes.at(locus).size(); tI++)
		{
			std::string HLA_type = loci_availableTypes.at(locus).at(tI);
			
			for(unsigned int gI = 0; gI < loci_types_haplotypes.at(locus).at(HLA_type).size(); gI++)
			{
				std::string& S = loci_types_haplotypes.at(locus).at(HLA_type).at(gI);
				if(S.find("*") != std::string::npos)
				{
					if(removeTypes[locus].count(HLA_type) == 0)
					{
						removeTypes[locus].insert(HLA_type);					
						std::cerr << "Locus " << locus << " allele " << HLA_type << " out because of stars.\n\n" << std::flush;
						std::cerr << gI << "/" << loci_types_haplotypes.at(locus).at(HLA_type).size() << ": " <<  S << "\n\n" << std::flush;
					}
				}
			}
			
			if(!(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size()))
			{
				std::cerr << "Length problem for locus " << locus << " allele " << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " vs " << loci_graphLevelIDs.at(locus).size() << "\n" << std::flush;
				removeTypes[locus].insert(HLA_type);
				continue;
			}	
			assert(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size());
		}
		
		if(removeTypes.count(locus) && (removeTypes.at(locus).size() > 0))
		{
			std::cerr << "\t\t" << locus << " before removal: " << loci_availableTypes.at(locus).size() << " types.\n";
			
			std::vector<std::string> availableTypes_new;
			for(unsigned int tI = 0; tI < loci_availableTypes.at(locus).size(); tI++)
			{
				std::string HLA_type = loci_availableTypes.at(locus).at(tI);
				if(removeTypes.at(locus).count(HLA_type) == 0)
				{
					availableTypes_new.push_back(HLA_type);
				}
			}
			loci_availableTypes.at(locus) = availableTypes_new;
			
			std::cerr << "\t\t\t" << locus << " after removal: " << loci_availableTypes.at(locus).size() << " types.\n" << std::flush;
			
			assert(loci_availableTypes.at(locus).size() > 0);
		}
	}

	// start simulations
	std::string outputFile_trueHLA = outputDirectory + "/trueHLA.txt";
	std::string outputFile_trueHaplotypes = outputDirectory + "/trueHaplotypes.txt";
	std::string outputFile_trueHaplotypes_perturbed = outputDirectory + "/trueHaplotypes.txt.perturbed";

	std::ofstream trueHLAstream;
	std::ofstream trueHaplotypesstream;
	std::ofstream trueHaplotypesPerturbedstream;

	trueHLAstream.open(outputFile_trueHLA.c_str());
	assert(trueHLAstream.is_open());
	trueHaplotypesstream.open(outputFile_trueHaplotypes.c_str());
	assert(trueHaplotypesstream.is_open());
	trueHaplotypesPerturbedstream.open(outputFile_trueHaplotypes_perturbed.c_str());
	assert(trueHaplotypesPerturbedstream.is_open());

	std::vector<std::string> truthFiles_headerFields;
	truthFiles_headerFields.push_back("IndividualID");
	truthFiles_headerFields.insert(truthFiles_headerFields.end(), loci.begin(), loci.end());
	
	trueHLAstream << Utilities::join(truthFiles_headerFields, " ") << "\n";
	trueHaplotypesstream << Utilities::join(truthFiles_headerFields, " ") << "\n";
	trueHaplotypesPerturbedstream << Utilities::join(truthFiles_headerFields, " ") << "\n";

	std::string outputFile_trueHaplotypes_fieldNames = outputDirectory + "/trueHaplotypes_fieldNames.txt";
	std::ofstream trueHaplotypes_fieldNames_stream;
	trueHaplotypes_fieldNames_stream.open(outputFile_trueHaplotypes_fieldNames.c_str());
	assert(trueHaplotypes_fieldNames_stream.is_open());
	for(std::set<std::string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
	{
		std::vector<std::string> lineFields;
		lineFields.push_back(*locusIt);
		lineFields.insert(lineFields.end(), loci_graphLevelIDs.at(*locusIt).begin(), loci_graphLevelIDs.at(*locusIt).end());
		trueHaplotypes_fieldNames_stream << Utilities::join(lineFields, " ") << "\n";
	}
	trueHaplotypes_fieldNames_stream.close();
 
	std::string outputFile_simulationDetails = outputDirectory + "/simulationDetails.txt";
	std::ofstream simulationDetailsStream;
	simulationDetailsStream.open(outputFile_simulationDetails.c_str());
	assert(simulationDetailsStream.is_open());
	simulationDetailsStream << "nIndividuals" << " " << nIndividuals << "\n";
	simulationDetailsStream << "exon23" << " " << exon23 << "\n";
	simulationDetailsStream << "perturbHaplotypes" << " " << perturbHaplotypes << "\n";
	simulationDetailsStream << "readError" << " " << readError << "\n";
	simulationDetailsStream << "haploidCoverage" << " " << haploidCoverage << "\n";
	simulationDetailsStream << "insertSize_mean" << " " << insertSize_mean << "\n";
	simulationDetailsStream << "insertSize_sd" << " " << insertSize_sd << "\n";
	simulationDetailsStream << "qualityMatrix" << " " << qualityMatrixFile << "\n";
	simulationDetailsStream << "readLength" << " " << readLength << "\n";

	simulationDetailsStream.close();

	assert(readLength > 0);
	readSimulator rS(qualityMatrixFile, (unsigned int)readLength);

	for(int sI = 0; sI < nIndividuals; sI++)
	{
		std::cout << "simulateHLAreads(..): Individual " << sI << " / " << nIndividuals << "\n" << std::flush;

		std::string outputFile_FASTQ_1 = outputDirectory + "/S_" + Utilities::ItoStr(sI) + "_1.fastq";
		std::string outputFile_FASTQ_2 = outputDirectory + "/S_" + Utilities::ItoStr(sI) + "_2.fastq";

		std::ofstream fastQStream_1;
		fastQStream_1.open(outputFile_FASTQ_1.c_str());
		assert(fastQStream_1.is_open());

		std::ofstream fastQStream_2;
		fastQStream_2.open(outputFile_FASTQ_2.c_str());
		assert(fastQStream_2.is_open());

		std::vector<std::string> outputFields_trueHLA;
		std::vector<std::string> outputFields_trueHaplotypes;
		std::vector<std::string> outputFields_trueHaplotypesPerturbed;

		outputFields_trueHLA.push_back("S_" + Utilities::ItoStr(sI));
		outputFields_trueHaplotypes.push_back("S_" + Utilities::ItoStr(sI));
		outputFields_trueHaplotypesPerturbed.push_back("S_" + Utilities::ItoStr(sI));
		
		for(std::set<std::string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
		{
			std::string locus = *locusIt;
			
			assert(loci_availableTypes.at(locus).size() > 0);
			
			std::string selectedType_1 = loci_availableTypes.at(locus).at(Utilities::randomNumber(loci_availableTypes.at(locus).size() - 1));
			std::string selectedType_2 = loci_availableTypes.at(locus).at(Utilities::randomNumber(loci_availableTypes.at(locus).size() - 1));

			std::cout << "\tLocus " << locus << " selected " << selectedType_1 << " / " << selectedType_2 << "\n" << std::flush;

			std::vector<std::string> haplotype_1 = loci_types_haplotypes.at(locus).at(selectedType_1);
			std::vector<std::string> haplotype_2 = loci_types_haplotypes.at(locus).at(selectedType_2);

			std::vector<int> perturbed_h1;
			std::vector<int> perturbed_h2;		

			if(perturbHaplotypes)
			{
				std::cout << "\t.. perturb...\n" << std::flush;			
				simulateHLAreads_perturbHaplotype(haplotype_1, perturbed_h1);
				simulateHLAreads_perturbHaplotype(haplotype_2, perturbed_h2); 
				std::cout << "\t\tdone " << haplotype_1.size() << " " << haplotype_2.size() << " " << perturbed_h1.size() << " " << perturbed_h2.size() << "...\n" << std::flush;			
			}

			std::string haplotype_1_noGaps = Utilities::join(haplotype_1, "");
			std::string haplotype_2_noGaps = Utilities::join(haplotype_2, "");

			haplotype_1_noGaps.erase(std::remove_if(haplotype_1_noGaps.begin(),haplotype_1_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), haplotype_1_noGaps.end());
			haplotype_2_noGaps.erase(std::remove_if(haplotype_2_noGaps.begin(),haplotype_2_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), haplotype_2_noGaps.end());

			std::vector<oneReadPair> simulatedReadPairs_h1 = rS.simulate_paired_reads_from_string(haplotype_1_noGaps, haploidCoverage, insertSize_mean, insertSize_sd, (! readError));
			std::vector<oneReadPair> simulatedReadPairs_h2 = rS.simulate_paired_reads_from_string(haplotype_2_noGaps, haploidCoverage, insertSize_mean, insertSize_sd, (! readError));

			auto print_one_readPair = [] (oneReadPair& rP, std::ofstream& output_1, std::ofstream& output_2) -> void
			{
				output_1 << "@" << rP.reads.first.name << "/1\n";
				output_1 << rP.reads.first.sequence << "\n";
				output_1 << "+" << "\n";
				output_1 << rP.reads.first.quality << "\n";

				output_2 << "@" << rP.reads.second.name << "/2\n";
				output_2 << rP.reads.second.sequence << "\n";
				output_2 << "+" << "\n";
				output_2 << rP.reads.second.quality << "\n";
			};

			for(unsigned int pI = 0; pI < simulatedReadPairs_h1.size(); pI++)
			{
				print_one_readPair(simulatedReadPairs_h1.at(pI), fastQStream_1, fastQStream_2);
			}

			for(unsigned int pI = 0; pI < simulatedReadPairs_h2.size(); pI++)
			{
				print_one_readPair(simulatedReadPairs_h2.at(pI), fastQStream_1, fastQStream_2);
			}


			auto removeStar = [](std::string t) -> std::string {
				std::vector<std::string> f = Utilities::split(t, "*");
				assert(f.size() == 2);
				return f.at(1);
			};
			
			outputFields_trueHLA.push_back(removeStar(selectedType_1) + "/" + removeStar(selectedType_2));
			outputFields_trueHaplotypes.push_back(Utilities::join(haplotype_1, ";") + "/" + Utilities::join(haplotype_2, ";"));
			
			std::vector<std::string> perturbed_h1_str;
			std::vector<std::string> perturbed_h2_str;
			for(unsigned int i = 0; i < perturbed_h1.size(); i++)
			{
				perturbed_h1_str.push_back(Utilities::ItoStr(perturbed_h1.at(i)));
			}
			for(unsigned int i = 0; i < perturbed_h2.size(); i++)
			{
				perturbed_h2_str.push_back(Utilities::ItoStr(perturbed_h2.at(i)));
			}			
			
			outputFields_trueHaplotypesPerturbed.push_back(Utilities::join(perturbed_h1_str, ";") + "/" + Utilities::join(perturbed_h2_str, ";"));
		}

		trueHLAstream << Utilities::join(outputFields_trueHLA, " ") << "\n" << std::flush;
		trueHaplotypesstream << Utilities::join(outputFields_trueHaplotypes, " ") << "\n" << std::flush;
		trueHaplotypesPerturbedstream << Utilities::join(outputFields_trueHaplotypesPerturbed, " ") << "\n" << std::flush;
		
		fastQStream_1.close();
		fastQStream_2.close();
	}

	trueHLAstream.close();
	trueHaplotypesstream.close();
	trueHaplotypesPerturbedstream.close();
}


std::set<std::string> getCompletelyDefinedHLAAlleles(std::string graphDir, std::string locus)
{
	std::vector<std::string> files_in_order;
	std::vector<std::string> files_in_order_type;
	std::vector<int> files_in_order_number;

	std::set<std::string> loci;
	std::map<std::string, std::vector<std::string> > files_per_locus;
	std::map<std::string, std::vector<std::string> > files_per_locus_type;
	std::map<std::string, std::vector<std::string> > files_per_locus_number;

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			std::string file_locus = split_by_underscore.at(0);
			if(file_locus == "before")
			{
				continue;
			}
			if(file_locus == "after")
			{
				continue;
			}
			loci.insert(file_locus);

			std::string filePath = graphDir + "/" + line;

			files_per_locus[file_locus].push_back(filePath);
			files_per_locus_type[file_locus].push_back(split_by_underscore.at(2));

			if(split_by_underscore.size() >= 4)
			{
				std::vector<std::string> split_by_colon = Utilities::split(split_by_underscore.at(3), ".");
				assert(split_by_colon.size() == 2);
				files_per_locus_number[file_locus].push_back(split_by_colon.at(0));
			}
			else
			{
				files_per_locus_number[file_locus].push_back("");
			}

		}
	}
	assert(loci.size() > 0);
	assert(loci.count(locus));


	std::string arbitraryIntronFile;
	for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
	{
		if(files_per_locus_type.at(locus).at(fI) == "intron")
		{
			arbitraryIntronFile = files_per_locus.at(locus).at(fI);
			break;
		}
	}

	if(arbitraryIntronFile.length() == 0)
	{
		for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
		{
			if(files_per_locus_type.at(locus).at(fI) == "exon")
			{
				arbitraryIntronFile = files_per_locus.at(locus).at(fI);
				break;
			}
		}
	}

	if(!(arbitraryIntronFile.length()))
	{
		std::cerr << "Cannot find an arbitrary intron/exon file for locus " << locus << "\n" << std::flush;
	}
	assert(arbitraryIntronFile.length());

	std::set<std::string> availableTypes;

	std::ifstream fileInputStream;
	fileInputStream.open(arbitraryIntronFile.c_str());
	if(!fileInputStream.is_open())
	{
		std::cerr << "Cannot open file " << arbitraryIntronFile << "\n" << std::flush;
	}
	assert(fileInputStream.is_open());
	unsigned int lI = 0;
	while(fileInputStream.good())
	{
		std::string line;
		std::getline(fileInputStream, line);
		Utilities::eraseNL(line);
		std::vector<std::string> line_fields = Utilities::split(line, " ");
		if(lI == 0)
		{
			if(!(line_fields.size() > 0))
			{
				std::cerr << "First line of file weird: " << arbitraryIntronFile << "\n" << std::flush;
			}
			assert(line_fields.size() > 0);
			assert(line_fields.at(0) == "IndividualID");
		}
		else
		{
			if(line.length())
			{
				std::string type = line_fields.at(0);
				availableTypes.insert(type);
			}
		}
		lI++;
	}
	fileInputStream.close();

	// find type haplotypes

	std::map<std::string, std::set<std::string> > removeTypes;

	std::map<std::string, std::map<std::string, std::vector<std::string>> > loci_types_haplotypes;
	std::map<std::string, std::vector<std::string> > loci_graphLevelIDs;

	for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
	{
		std::ifstream fileInputStream;
		fileInputStream.open(files_per_locus.at(locus).at(fI).c_str());
		assert(fileInputStream.is_open());
		std::vector<std::string> file_lines;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			file_lines.push_back(line);
		}
		fileInputStream.close();

		std::string firstLine = file_lines.at(0);
		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
		assert(firstLine_fields.at(0) == "IndividualID");

		std::vector<std::string> graph_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
		loci_graphLevelIDs[locus].insert(loci_graphLevelIDs[locus].end(), graph_level_names.begin(), graph_level_names.end());

		for(unsigned int lI = 1; lI < file_lines.size(); lI++)
		{
			if(file_lines.at(lI).length())
			{
				std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
				assert(line_fields.size() == firstLine_fields.size());
				std::string HLA_type = line_fields.at(0);
				std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

				if((files_per_locus_type.at(locus).at(fI) == "intron") || (files_per_locus_type.at(locus).at(fI) == "exon"))
				{
					if(availableTypes.count(HLA_type))
					{
						if(loci_types_haplotypes[locus].count(HLA_type) == 0)
						{
							loci_types_haplotypes[locus][HLA_type].resize(0);
						}

						loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());
					}
				}
				else
				{
					// this must be a pre- or post-padding sequence
					assert((files_per_locus_type.at(locus).at(fI) == "paddingLeft.txt") || (files_per_locus_type.at(locus).at(fI) == "paddingRight.txt"));

					bool noStars = true;
					for(unsigned int aI = 0; aI < line_alleles.size(); aI++)
					{
						std::string& S = line_alleles.at(aI);
						if(S.find("*") != std::string::npos)
						{
							noStars = false;
							break;
						}
					}

					if(lI == (file_lines.size() - 1))
					{
						noStars = true;
						std::cerr << "Locus " << locus << " file " << files_per_locus.at(locus).at(fI) << ": No line without stars -- force " << lI << " / " << file_lines.size() << "\n" << std::flush;
					}

					// std::cout << "TYPE: " << files_per_locus_type.at(locus).at(fI) << "\n";
					if(noStars)
					{
						for(std::set<std::string>::iterator typeIt = availableTypes.begin(); typeIt != availableTypes.end(); typeIt++)
						{
							std::string HLA_type = *typeIt;
							if(loci_types_haplotypes[locus].count(HLA_type) == 0)
							{
								loci_types_haplotypes[locus][HLA_type].resize(0);
							}

							loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());

							// std::cerr << "\t\t" << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " / " << loci_graphLevelIDs[locus].size() << "\n";
						}
						break;
					}
				}
			}
		}
	}

	// check that all haplotypes make sense
	for(std::set<std::string>::iterator HLA_typeIT = availableTypes.begin(); HLA_typeIT != availableTypes.end(); HLA_typeIT++)
	{
		std::string HLA_type = *HLA_typeIT;

		for(unsigned int gI = 0; gI < loci_types_haplotypes.at(locus).at(HLA_type).size(); gI++)
		{
			std::string& S = loci_types_haplotypes.at(locus).at(HLA_type).at(gI);
			if(S.find("*") != std::string::npos)
			{
				if(removeTypes[locus].count(HLA_type) == 0)
				{
					removeTypes[locus].insert(HLA_type);
					std::cerr << "Locus " << locus << " allele " << HLA_type << " out because of stars.\n\n" << std::flush;
					std::cerr << gI << "/" << loci_types_haplotypes.at(locus).at(HLA_type).size() << ": " <<  S << "\n\n" << std::flush;
				}
			}
		}

		if(!(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size()))
		{
			std::cerr << "Length problem for locus " << locus << " allele " << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " vs " << loci_graphLevelIDs.at(locus).size() << "\n" << std::flush;
			removeTypes[locus].insert(HLA_type);
			continue;
		}
		assert(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size());
	}

	if(removeTypes.count(locus) && (removeTypes.at(locus).size() > 0))
	{
		std::cerr << "\t\t" << locus << " before removal: " << availableTypes.size() << " types.\n";

		std::set<std::string> availableTypes_new;
		for(std::set<std::string>::iterator HLA_typeIT = availableTypes.begin(); HLA_typeIT != availableTypes.end(); HLA_typeIT++)
		{
			std::string HLA_type = *HLA_typeIT;
			if(removeTypes.at(locus).count(HLA_type) == 0)
			{
				availableTypes_new.insert(HLA_type);
			}
		}
		availableTypes = availableTypes_new;

		std::cerr << "\t\t\t" << locus << " after removal: " << availableTypes.size() << " types.\n" << std::flush;

		assert(availableTypes.size() > 0);
	}

	return availableTypes;
}


void HLAHaplotypeInference(std::string alignedReads_file, std::string graphDir, std::string outputDir, std::string sampleName, std::string loci_str, std::string starting_haplotypes_perLocus_1_str, std::string starting_haplotypes_perLocus_2_str, bool longUnpairedReads)
{
	
	std::cout << Utilities::timestamp() << "HLAHaplotypeInference(..): Start.\n";
	
	// std::cout << "omp_get_max_threads(): " << omp_get_max_threads() << "\n" << std::flush;
	
	std::cout << "\t\t" << "loci_str" << ": " << loci_str << "\n";
	std::cout << "\t\t" << "starting_haplotypes_perLocus_1_str" << ": " << starting_haplotypes_perLocus_1_str << "\n";
	std::cout << "\t\t" << "starting_haplotypes_perLocus_2_str" << ": " << starting_haplotypes_perLocus_2_str << "\n\n" << std::flush;

	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	std::vector<std::string> loci = Utilities::split(loci_str, ",");
	std::vector<std::string> starting_haplotypes_perLocus_1 = Utilities::split(starting_haplotypes_perLocus_1_str, ",");
	std::vector<std::string> starting_haplotypes_perLocus_2 = Utilities::split(starting_haplotypes_perLocus_2_str, ",");
	
	if(!((loci.size() == starting_haplotypes_perLocus_1.size()) && (loci.size() == starting_haplotypes_perLocus_2.size())))
	{
		std::cout << "Problem!\n";
		std::cout << "\t" << "loci.size()" << ": " << loci.size() << "\n";
		std::cout << "\t" << "starting_haplotypes_perLocus_1.size()" << ": " << starting_haplotypes_perLocus_1.size() << "\n";
		std::cout << "\t" << "starting_haplotypes_perLocus_2.size()" << ": " << starting_haplotypes_perLocus_2.size() << "\n";
		std::cout << "\t" << "loci_str" << ": " << loci_str << "\n";
		std::cout << "\t" << "starting_haplotypes_perLocus_1_str" << ": " << starting_haplotypes_perLocus_1_str << "\n";
		std::cout << "\t" << "starting_haplotypes_perLocus_2_str" << ": " << starting_haplotypes_perLocus_2_str << "\n";
		std::cout << "\n" << std::flush;
	}
	assert(loci.size() == starting_haplotypes_perLocus_1.size());
	assert(loci.size() == starting_haplotypes_perLocus_2.size());

	// translate location IDs to graph levels
	std::vector<std::string> graphLoci = readGraphLoci(graphDir);
	std::map<std::string, unsigned int> graphLocus_2_levels;
	for(unsigned int i = 0; i < graphLoci.size(); i++)
	{
		std::string locusID = graphLoci.at(i);
		assert(graphLocus_2_levels.count(locusID) == 0);
		graphLocus_2_levels[locusID] = i;
	}

	// load reads
	std::cout << Utilities::timestamp() << "HLAHaplotypeInference(..): Load reads.\n" << std::flush;
	std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> alignments_paired;
	std::vector<oneReadPair> alignments_originalReads_paired;

	std::vector<seedAndExtend_return_local> alignments_unpaired;
	std::vector<oneRead> alignments_originalReads_unpaired;

	double insertSize_mean;
	double insertSize_sd;

	if(longUnpairedReads)
	{
		read_longReadAlignments_fromFile(alignedReads_file, alignments_unpaired, alignments_originalReads_unpaired);
		assert(alignments_unpaired.size() == alignments_originalReads_unpaired.size());

		std::cout << Utilities::timestamp() << "HLAHaplotypeInference(..): Load reads -- done. Have " << alignments_unpaired.size() << " unpaired reads, long unpaired reads: " << longUnpairedReads << ".\n" << std::flush;

	}
	else
	{
		read_shortReadAlignments_fromFile(alignedReads_file, alignments_paired, alignments_originalReads_paired, insertSize_mean, insertSize_sd);
		assert(alignments_paired.size() == alignments_originalReads_paired.size());

		std::cout << Utilities::timestamp() << "HLAHaplotypeInference(..): Load reads -- done. Have " << alignments_paired.size() << " read pairs, long unpaired reads: " << longUnpairedReads << ", IS mean " << insertSize_mean << " / sd " << insertSize_sd << ".\n" << std::flush;
	}

	// read alignment statistics

	int alignmentStats_strandsValid = 0;
	int alignments_paired_perfect = 0;
	int alignments_paired_oneReadPerfect = 0;

	int alignments_unpaired_perfect = 0;

	int alignmentStats_strandsValid_and_distanceOK = 0;
	std::vector<double> alignmentStats_strandsValid_distances;
	double alignmentStats_paired_fractionOK_sum = 0;
	double alignmentStats_unpaired_fractionOK_sum = 0;

	for(unsigned int alignmentI = 0; alignmentI < alignments_paired.size(); alignmentI++)
	{

		std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignedReadPair = alignments_paired.at(alignmentI);

		if(alignedReadPair_strandsValid(alignedReadPair))
		{
			alignmentStats_strandsValid++;
			double pairsDistance = alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair);
			alignmentStats_strandsValid_distances.push_back(pairsDistance);
			if(abs(pairsDistance - insertSize_mean) <= (5 * insertSize_sd))
			{
				alignmentStats_strandsValid_and_distanceOK++;
			}
		}

		double fractionOK_1 = alignmentFractionOK(alignedReadPair.first);
		double fractionOK_2 = alignmentFractionOK(alignedReadPair.second);

		if(fractionOK_1 == 1)
		{
			alignments_paired_perfect++;
		}
		if(fractionOK_2 == 1)
		{
			alignments_paired_perfect++;
		}
		if((fractionOK_1 == 1) || (fractionOK_2 == 1))
		{
			alignments_paired_oneReadPerfect++;
		}

		alignmentStats_paired_fractionOK_sum += fractionOK_1;
		alignmentStats_paired_fractionOK_sum += fractionOK_2;
	}

	for(unsigned int alignmentI = 0; alignmentI < alignments_unpaired.size(); alignmentI++)
	{

		seedAndExtend_return_local& alignedRead = alignments_unpaired.at(alignmentI);

		double fractionOK = alignmentFractionOK(alignedRead);

		if(fractionOK == 1)
		{
			alignments_unpaired_perfect++;
		}

		alignmentStats_unpaired_fractionOK_sum += fractionOK;
	}

	// std::pair<double, double> alignmentStats_distance_meanMedian = meanMedian(alignmentStats_strandsValid_distances);
	// double alignmentStats_fractionOK_avg = (alignments.size() > 0) ? (alignmentStats_fractionOK_sum / (2.0* (double)alignments.size())) : 0;

	if(! Utilities::directoryExists("../tmp/hla"))
	{
		Utilities::makeDir("../tmp/hla");
	}

	std::string outputDirectory = outputDir + "/hla/" +sampleName;
	if(! Utilities::directoryExists(outputDirectory))
	{
		Utilities::makeDir(outputDirectory);
	}

	std::string outputFN_bestGuess_6dig = outputDirectory + "/R2_bestguess_6dig.txt";
	std::ofstream bestGuess_6dig_outputStream;
	bestGuess_6dig_outputStream.open(outputFN_bestGuess_6dig.c_str());
	assert(bestGuess_6dig_outputStream.is_open());
	bestGuess_6dig_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "AbsoluteDifferences" << "\t"<< "RelativeDifference" << "\n";

	std::string outputFN_bestGuess_8dig = outputDirectory + "/R2_bestguess_8dig.txt";
	std::ofstream bestGuess_8dig_outputStream;
	bestGuess_8dig_outputStream.open(outputFN_bestGuess_8dig.c_str());
	assert(bestGuess_8dig_outputStream.is_open());
	bestGuess_8dig_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "AbsoluteDifferences" << "\t" << "RelativeDifference" << "\n";


	for(unsigned int locusI = 0; locusI < loci.size(); locusI++)
	{
		std::string locus = loci.at(locusI);

		std::cout << Utilities::timestamp() << "HLAHaplotypeInference(..): Making inference for " << locus << "\n" << std::flush;

		std::vector<int> combined_sequences_graphLevels;
		std::vector<std::string> combined_sequences_graphLevels_individualType;
		std::vector<int> combined_sequences_graphLevels_individualTypeNumber;
		std::vector<int> combined_sequences_graphLevels_individualPosition;

		std::vector<std::string> combined_sequences_locusIDs;
		std::map<std::string, std::string> combined_incomplete_sequences_paddingSetToStars;
		std::vector<std::string> starting_haplotypes_combined_vec;

		read_HLA_alleles_for_haplotypeInference(
			graphDir,
			locus,
			combined_sequences_graphLevels,
			combined_sequences_graphLevels_individualType,
			combined_sequences_graphLevels_individualTypeNumber,
			combined_sequences_graphLevels_individualPosition,
			combined_sequences_locusIDs,
			combined_incomplete_sequences_paddingSetToStars,
			graphLocus_2_levels,
			false,
			false,
			true
		);

		std::map<std::string, std::string> combined_sequences;

		read_HLA_alleles_for_haplotypeInference(
			graphDir,
			locus,
			combined_sequences_graphLevels,
			combined_sequences_graphLevels_individualType,
			combined_sequences_graphLevels_individualTypeNumber,
			combined_sequences_graphLevels_individualPosition,
			combined_sequences_locusIDs,
			combined_sequences,
			graphLocus_2_levels,
			true,
			false,
			false
		);

		assert(combined_sequences.size());

		std::string starting_haplotypes_1_str = starting_haplotypes_perLocus_1.at(locusI);
		std::string starting_haplotypes_2_str = starting_haplotypes_perLocus_2.at(locusI);
		std::vector<std::string> starting_haplotypes_1_vec_preTranslation = Utilities::split(starting_haplotypes_1_str, ";");
		std::vector<std::string> starting_haplotypes_2_vec_preTranslation = Utilities::split(starting_haplotypes_2_str, ";");

		// translate predefined types into full-length types
		auto translatePartialTypesToFull = [&](std::vector<std::string> partialTypes) -> std::vector<std::string>
		{
			std::vector<std::string> forReturn;
			for(unsigned int i = 0; i < partialTypes.size(); i++)
			{
				std::string partialType = partialTypes.at(i);
				if(combined_sequences.count(partialType))
				{
					forReturn.push_back(partialType);
				}
			}

			if(forReturn.size() > 0)
			{
				std::cout << Utilities::timestamp() << "\t\ttranslatePartialTypesToFull(..): Input set of size " << partialTypes.size() << " immediately mapped onto output set of size " << forReturn.size() << ".\n";

				return forReturn;
			}
			else
			{
				std::map<std::string, int> completeTypes_editDistances;
				for(std::map<std::string,std::string>::iterator completeTypeIt = combined_sequences.begin(); completeTypeIt != combined_sequences.end(); completeTypeIt++)
				{
					std::string completeType = completeTypeIt->first;

					for(unsigned int i = 0; i < partialTypes.size(); i++)
					{
						std::string partialType = partialTypes.at(i);

						if(! combined_incomplete_sequences_paddingSetToStars.count(partialType))
						{
							std::cerr << "Don't have type '" << partialType << "' at locus " << locus << "\n" << std::flush;
						}
						std::string completeTypeSequence = combined_sequences.at(completeType);
						std::string partialTypeSequence = combined_incomplete_sequences_paddingSetToStars.at(partialType);

						assert(completeTypeSequence.length() == partialTypeSequence.length());
						if(completeTypeSequence.length() == partialTypeSequence.length())
						{
							int editDistance = compute_Hamming_distance(completeTypeSequence, partialTypeSequence, true);
							if(completeTypes_editDistances.count(completeType) == 0)
							{
								completeTypes_editDistances[completeType] = editDistance;
							}
							else
							{
								if(editDistance < completeTypes_editDistances.at(completeType))
								{
									completeTypes_editDistances.at(completeType) = editDistance;
								}
							}
						}
					}
				}

				if(completeTypes_editDistances.size() == 0)
				{
					std::cout << Utilities::timestamp() << "\t\ttranslatePartialTypesToFull(..): Could not use a single complete type for edit distance check - use arbitrary complete type.";
					forReturn.push_back(combined_sequences.begin()->first);
				}
				else
				{
					std::vector<std::string> completeTypes_editDistances_keys;
					for(std::map<std::string, int>::iterator completeTypeIt = completeTypes_editDistances.begin(); completeTypeIt != completeTypes_editDistances.end(); completeTypeIt++)
					{
						completeTypes_editDistances_keys.push_back(completeTypeIt->first);
					}
					std::sort(completeTypes_editDistances_keys.begin(), completeTypes_editDistances_keys.end(), [&](std::string a, std::string b){return (completeTypes_editDistances.at(a) < completeTypes_editDistances.at(b));});

					if(completeTypes_editDistances_keys.size() > 1)
					{
						assert(completeTypes_editDistances.at(completeTypes_editDistances_keys.at(0)) <= completeTypes_editDistances.at(completeTypes_editDistances_keys.at(1)));
					}

					std::vector<std::string> equallyGoodEditDistances;
					int bestEditDistance = completeTypes_editDistances.at(completeTypes_editDistances_keys.at(0));
					for(unsigned int i = 0; i < completeTypes_editDistances_keys.size(); i++)
					{
						int competingEditDistance = completeTypes_editDistances.at(completeTypes_editDistances_keys.at(i));
						if(competingEditDistance == bestEditDistance)
						{
							equallyGoodEditDistances.push_back(completeTypes_editDistances_keys.at(i));
						}
						else
						{
							break;
						}
					}

					assert(equallyGoodEditDistances.size() >= 1);

					forReturn = equallyGoodEditDistances;

					std::cout << Utilities::timestamp() << "\t\ttranslatePartialTypesToFull(..): Input set of size " << partialTypes.size() << " edit-distance mapped onto output set of size " << forReturn.size() << ".\n";
				}

				return forReturn;
			}
		};

		std::vector<std::string> starting_haplotypes_1_vec = translatePartialTypesToFull(starting_haplotypes_1_vec_preTranslation);
		std::vector<std::string> starting_haplotypes_2_vec = translatePartialTypesToFull(starting_haplotypes_2_vec_preTranslation);

		std::set<std::string> starting_haplotypes_1_set(starting_haplotypes_1_vec.begin(), starting_haplotypes_1_vec.end());
		std::set<std::string> starting_haplotypes_2_set(starting_haplotypes_2_vec.begin(), starting_haplotypes_2_vec.end());
		std::set<std::string> starting_haplotypes_combined_set = starting_haplotypes_1_set;
		starting_haplotypes_combined_set.insert(starting_haplotypes_2_set.begin(), starting_haplotypes_2_set.end());
		
		starting_haplotypes_combined_vec.clear();
		starting_haplotypes_combined_vec.insert(starting_haplotypes_combined_vec.end(), starting_haplotypes_combined_set.begin(), starting_haplotypes_combined_set.end());

		assert(starting_haplotypes_1_set.size() > 0);
		assert(starting_haplotypes_2_set.size() > 0);

		std::cout << "\n\n" << Utilities::timestamp() << "\thave " << starting_haplotypes_combined_set.size() << " combined starting haplotypes\n" << std::flush;
		std::cout << "\t\t" << Utilities::join(starting_haplotypes_1_vec, ";") << "\n";
		std::cout << "\t\t" << Utilities::join(starting_haplotypes_2_vec, ";") << "\n\n" << std::flush;

		std::string file_AAmapping = graphDir + "/AAmapping/"+locus+".txt";
		// std::cerr << "File: " << file_AAmapping << "\n\n" << std::flush;
		// assert( 1 == 0 );

		std::map<std::string, std::string> graphLevel_2_AAid;
		if(Utilities::fileExists(file_AAmapping))
		{
			std::ifstream fileInputStream;
			fileInputStream.open(file_AAmapping.c_str());
			assert(fileInputStream.is_open());
			std::vector<std::string> file_lines;
			while(fileInputStream.good())
			{
				std::string line;
				std::getline(fileInputStream, line);
				Utilities::eraseNL(line);
				if(line.length())
				{
					std::vector<std::string> line_parts = Utilities::split(line, "\t");
					assert(line_parts.size() == 2);
					assert(graphLevel_2_AAid.count(line_parts.at(0)) == 0);
					graphLevel_2_AAid[line_parts.at(0)] = line_parts.at(1);
				}
			}
			fileInputStream.close();
		}
		
		std::set<std::string> non_full_length_types;
		for(std::set<std::string>::iterator haplotypeIt = starting_haplotypes_combined_set.begin(); haplotypeIt != starting_haplotypes_combined_set.end(); haplotypeIt++)
		{
			std::string haplotypeID = *haplotypeIt;
			assert(combined_sequences.count(haplotypeID) > 0);
			// std::cout << "Length comparison " << haplotypeID << ": "  << combined_sequences.at(haplotypeID).length() << " vs " << combined_sequences_locusIDs.size() << "\n" << std::flush;
			if(combined_sequences.at(haplotypeID).length() != combined_sequences_locusIDs.size())
			{
				non_full_length_types.insert(haplotypeID);
			}
		}
		assert(non_full_length_types.size() == 0);
//
//		for(std::set<std::string>::iterator haplotypeIt = non_full_length_types.begin(); haplotypeIt != non_full_length_types.end(); haplotypeIt++)
//		{
//			combined_sequences.erase(*haplotypeIt);
//		}
		
		for(std::set<std::string>::iterator haplotypeIt = starting_haplotypes_combined_set.begin(); haplotypeIt != starting_haplotypes_combined_set.end(); haplotypeIt++)
		{
			std::string haplotypeID = *haplotypeIt;
			assert(combined_sequences.count(haplotypeID) > 0);
			assert(combined_sequences.at(haplotypeID).length() == combined_sequences_locusIDs.size());
		}
		
		for(std::set<std::string>::iterator haplotypeIt = starting_haplotypes_1_set.begin(); haplotypeIt != starting_haplotypes_1_set.end(); haplotypeIt++)
		{
			assert(combined_sequences.count(*haplotypeIt));
		}		
		for(std::set<std::string>::iterator haplotypeIt = starting_haplotypes_2_set.begin(); haplotypeIt != starting_haplotypes_2_set.end(); haplotypeIt++)
		{
			assert(combined_sequences.count(*haplotypeIt));
		}		

		std::map<int, unsigned int> graphLevel_2_position;
		std::map<int, std::string> graphLevel_2_position_type;
		std::map<int, unsigned int> graphLevel_2_exonPosition_typeNumber;
		std::map<int, unsigned int> graphLevel_2_exonPosition_position;
		std::map<int, std::string> graphLevel_2_locusID;
		for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
		{
			int graphLevel = combined_sequences_graphLevels.at(pI);
			assert(graphLevel >= 0);
			graphLevel_2_position[graphLevel] = pI;
			graphLevel_2_position_type[graphLevel] = combined_sequences_graphLevels_individualType.at(pI);
			graphLevel_2_exonPosition_typeNumber[graphLevel] = combined_sequences_graphLevels_individualTypeNumber.at(pI);
			graphLevel_2_exonPosition_position[graphLevel] = combined_sequences_graphLevels_individualPosition.at(pI);
			graphLevel_2_locusID[graphLevel] = combined_sequences_locusIDs.at(pI);

		}

		std::cout << Utilities::timestamp() << "Have collected " << combined_sequences.size() << " sequences -- first level " << combined_sequences_graphLevels.front() << ", last level " << combined_sequences_graphLevels.back() << ".\n" << std::flush;

		// now transform reads into sequences specifying genotype value

		std::cout << Utilities::timestamp() << "Compute exon positions and specified genotypes from reads\n" << std::flush;

		std::vector< std::vector<oneExonPosition> > positions_fromReads;

		unsigned int readPairs_OK = 0;
		unsigned int readPairs_broken = 0;
		unsigned int reads_unpaired_OK = 0;

		for(unsigned int readPairI = 0; readPairI < alignments_paired.size(); readPairI++)
		{
			oneReadPair& originalReadPair = alignments_originalReads_paired.at(readPairI);
			std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignedReadPair = alignments_paired.at(readPairI);

			std::vector<oneExonPosition> read1_positions;
			std::vector<oneExonPosition> read2_positions;

			oneReadAlignment_2_exonPositions_paired(alignedReadPair.first, originalReadPair.reads.first, read1_positions, alignedReadPair.second, originalReadPair.reads.second, 1, combined_sequences_graphLevels, graphLevel_2_position);
			oneReadAlignment_2_exonPositions_paired(alignedReadPair.second, originalReadPair.reads.second, read2_positions, alignedReadPair.first, originalReadPair.reads.first, 2, combined_sequences_graphLevels, graphLevel_2_position);

			assert(alignedReadPair.first.mapQ_genomic != 2);
			double mapQ_thisAlignment = alignedReadPair.first.mapQ_genomic;

			if(
					alignedReadPair_strandsValid(alignedReadPair) &&
					(abs(alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)) &&
					(mapQ_thisAlignment >= minimumMappingQuality) &&
					((alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.first) >= min_bothReads_weightedCharactersOK) && (alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.second) >= min_bothReads_weightedCharactersOK))
			)
			{
				// good

				// std::cout << "\t\t" << "readPair " << readPairI << ", pairing OK.\n" << std::flush;

				std::vector<oneExonPosition> thisRead_positions = read1_positions;
				thisRead_positions.insert(thisRead_positions.end(), read2_positions.begin(), read2_positions.end());
				

								
				if(thisRead_positions.size() > 0)
				{
					thisRead_positions = removeDoublePositionsFromRead(thisRead_positions);
					positions_fromReads.push_back(thisRead_positions);
				}
				
				readPairs_OK++;
			}
		}

		for(unsigned int readI = 0; readI < alignments_unpaired.size(); readI++)
		{
			oneRead& originalRead = alignments_originalReads_unpaired.at(readI);
			seedAndExtend_return_local& alignedRead = alignments_unpaired.at(readI);

			std::vector<oneExonPosition> read_positions;

			oneReadAlignment_2_exonPositions_unpaired(alignedRead, originalRead, read_positions, combined_sequences_graphLevels, graphLevel_2_position);

			assert(alignedRead.mapQ_genomic != 2);
			double mapQ_thisAlignment = alignedRead.mapQ_genomic;

			if(
					(mapQ_thisAlignment >= minimumMappingQuality)
			)
			{
				// good

				// std::cout << "\t\t" << "readPair " << readPairI << ", pairing OK.\n" << std::flush;
				
				std::vector<oneExonPosition> thisRead_positions = read_positions;
				if(thisRead_positions.size() > 0)
					positions_fromReads.push_back(thisRead_positions);

				reads_unpaired_OK++;
			}
		}


		std::cout << Utilities::timestamp() << "Mapped reads to exons and introns. " << readPairs_OK << " pairs OK, " << readPairs_broken << " pairs broken, " << reads_unpaired_OK << " unpaired reads OK." << "\n" << std::flush;

		// Pileup of mapped reads

		std::vector<std::vector<oneExonPosition> > pileUpPerPosition;
		pileUpPerPosition.resize(combined_sequences_graphLevels.size());
		for(unsigned int positionSpecifierI = 0; positionSpecifierI < positions_fromReads.size(); positionSpecifierI++)
		{
			std::vector<oneExonPosition>& individualPositions = positions_fromReads.at(positionSpecifierI);
			for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
			{
				oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);

				// std::cout << "Read " << positionSpecifierI << " position " << positionI << " mapQ_position: " << onePositionSpecifier.mapQ_position << " (vs " << minimumPerPositionMappingQuality << ")\n" << std::flush;
				
				if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
				{
					continue;
				}

				unsigned int position = graphLevel_2_position.at(onePositionSpecifier.graphLevel);
				pileUpPerPosition.at(position).push_back(onePositionSpecifier);
			}
		}

		std::string fileName_pileUp = outputDirectory + "/R2_pileup_"+locus+".txt";
		std::ofstream pileUpStream;
		pileUpStream.open(fileName_pileUp.c_str());
		assert(pileUpStream.is_open());

		class haplotypeAlternative {
		protected:

			std::vector<std::string> h1;
			std::vector<std::string> h2;
			double LL1;
			double LL2;
			double LLPriors;
			double Pexternal;
			bool setPexternal;
			
			std::map<std::string, std::pair<double, double> > ll_per_read;
			std::map<std::string, double> averaged_ll_per_read;

			int ll_computed_until_position;

		public:
			haplotypeAlternative()
			{
				LL1 = 0;
				LL2 = 0;
				LLPriors = 0;
				ll_computed_until_position = -1;
				setPexternal = false;
			}

			const std::vector<std::string>& getH1()
			{
				return h1;
			}

			const std::vector<std::string>& getH2()
			{
				return h2;
			}

			const std::vector<std::string>& getH(unsigned int i)
			{
				assert((i == 1) || (i == 2));
				if(i == 1)
				{
					return h1;
				}
				else
				{
					return h2;
				}
			}


			std::pair<std::string, std::string> getHaplotypeAlleles(unsigned int pI)
			{
				return make_pair(h1.at(pI), h2.at(pI));
			}
			
			std::pair<std::string, std::string> getTrailingHaplotypeAlleles()
			{
				return make_pair(h1.back(), h2.back());
			}			
			
			std::pair<std::string, std::string> getTrailingHaplotypeAlleles(unsigned int pI)
			{	
				std::string h1Str;
				std::string h2Str;
				for(unsigned int i = pI; i < h1.size(); i++)
				{
					h1Str += h1.at(i);
					h2Str += h2.at(i);
				}
				
				return make_pair(h1Str, h2Str);
			}			
			

			void setRescaledP(double p)
			{
				setPexternal = true;
				Pexternal = p;
			}
			
			double getP()
			{
				assert(setPexternal);
				return Pexternal;
			}
			
			double getLL1()
			{
				return LL1 + LLPriors;
			}
			
			double getLL2()
			{
				return LL1 + LLPriors;
			}
			
			int getLL_computedUntil()
			{
				return ll_computed_until_position;
			}

			void updateLikelihood(const std::set<std::string>& alleles_from_h1, const std::set<std::string>& alleles_from_h2, const std::vector<oneExonPosition>& pileUpPerPosition)
			{
				const std::string& h1_underlying_position = h1.back();
				const std::string& h2_underlying_position = h2.back();
			
				bool verbose = (ll_computed_until_position > 502);
				verbose = false;
				
				double newAlleles_h1 = 4 - alleles_from_h1.size();
				double newAlleles_h2 = 4 - alleles_from_h2.size();
				if(newAlleles_h1 < 1)
					newAlleles_h1 = 1;
				if(newAlleles_h2 < 1)
					newAlleles_h2 = 1;					
				
				double existingAlleles_h1 = 4 - newAlleles_h1;
				double existingAlleles_h2 = 4 - newAlleles_h2;
				
				assert(existingAlleles_h1 >= 1);
				assert(existingAlleles_h1 <= 4);
				assert(existingAlleles_h2 >= 1);
				assert(existingAlleles_h2 <= 4);
				assert(newAlleles_h1 >= 1);
				assert(newAlleles_h1 <= 4);
				assert(newAlleles_h2 >= 1);
				assert(newAlleles_h2 <= 4);
				
				double p_new = 0.999;
				double log_existingAllele_h1 = log(p_new *  ( 1.0 / existingAlleles_h1) );
				double log_newAllele_h1 = log((1 - p_new) * (1.0 / newAlleles_h1) );
				double log_existingAllele_h2 = log(p_new *  ( 1.0 / existingAlleles_h2) );
				double log_newAllele_h2 = log((1 - p_new) * (1.0 / newAlleles_h2) );
							
				
				LLPriors += (alleles_from_h1.count(h1_underlying_position)) ? log_existingAllele_h1 : log_newAllele_h1;
				LLPriors += (alleles_from_h2.count(h2_underlying_position)) ? log_existingAllele_h2 : log_newAllele_h2;
				
				//first LL
				{
					assert(h1.size() == h2.size());
					if(!((int)h1.size() == (ll_computed_until_position + 2)))
					{
						std::cerr << "h1.size(): " << h1.size() << "\n";
						std::cerr << "ll_computed_until_position: " << ll_computed_until_position << "\n" << std::flush;
						
					}
					

						
					assert((int)h1.size() == (ll_computed_until_position + 2));

					if(verbose)
					{
						std::cout << "\t\\tAfter 'priors' update: " << LLPriors << "\n" << std::flush;
					}
					
					if(verbose)
					{
						std::cout << "\t\tBefore update LL1: " << LL1 << "\n" << std::flush;
					}
					
					std::set<std::string> read_IDs_modified;
					for(unsigned int rI = 0; rI < pileUpPerPosition.size(); rI++)
					{
						const oneExonPosition& r = pileUpPerPosition.at(rI);
						const std::string& readID = r.thisRead_ID;

						double position_likelihood_h1 = read_likelihood_per_position(h1_underlying_position, r.genotype, r.qualities, r.graphLevel);
						double position_likelihood_h2 = read_likelihood_per_position(h2_underlying_position, r.genotype, r.qualities, r.graphLevel);

						if(ll_per_read.count(readID) == 0)
						{
							ll_per_read[readID].first = 0;
							ll_per_read[readID].second = 0;
						}

						// double first_before = ll_per_read.at(readID).first;
						// double second_before = ll_per_read.at(readID).second;
						
						// std::cout << "\t\t\tRead " << readID << ": " << r.genotype << "\n";

						ll_per_read.at(readID).first += position_likelihood_h1;
						ll_per_read.at(readID).second += position_likelihood_h2;

						// std::cout << "\t\t\t\tOn H1 [" << h1_underlying_position << "]: " << first_before << " -> " << ll_per_read.at(readID).first << "\n";
						// std::cout << "\t\t\t\tOn H2 [" << h2_underlying_position << "]: " << second_before << " -> " << ll_per_read.at(readID).second << "\n";
						
						read_IDs_modified.insert(readID);
						
						
					}

					for(std::set<std::string>::iterator readIDit = read_IDs_modified.begin(); readIDit != read_IDs_modified.end(); readIDit++)
					{
						if(averaged_ll_per_read.count(*readIDit))
						{
							LL1 -= averaged_ll_per_read.at(*readIDit);
						}
						else
						{
							averaged_ll_per_read[*readIDit] = 0;
						}

						// std::cout << *readIDit << "\n";
						// std::cout << "\t" << ll_per_read.at(*readIDit).first << "\n";
						// std::cout << "\t" << ll_per_read.at(*readIDit).second << "\n" << std::flush;
						
						// double new_average_p = 0.5 * exp(ll_per_read.at(*readIDit).first) + 0.5 * exp(ll_per_read.at(*readIDit).second);
						
						// std::cout << "\t" << new_average_p << "\n" << std::flush;
						
						// assert(new_average_p >= 0);
						// assert(new_average_p <= 1);

						double new_average_logP = Utilities::logAvg(ll_per_read.at(*readIDit).first, ll_per_read.at(*readIDit).second);
						// if(new_average_p == 0)
						// {
							// new_average_logP = -1e100;
						// }
						// std::cout << new_average_logP << " " << new_average_p << "\n" << std::flush;
						LL1 += new_average_logP;
						averaged_ll_per_read.at(*readIDit) = new_average_logP;
					}

					if(verbose)
					{
						std::cout << "\t\tLL1 After complete update: " << LL1 << "\n" << std::flush;
					}
				}
				
				//LL2
				{
				
					if(verbose)
					{
						std::cout << "\t\tBefore update LL2: " << LL2 << "\n" << std::flush;
					}	
					
					double thisPosition_LL2 = 0;
					for(unsigned int rI = 0; rI < pileUpPerPosition.size(); rI++)
					{
						const oneExonPosition& r = pileUpPerPosition.at(rI);
						// const std::string& readID = r.thisRead_ID;
					
						double position_likelihood_h1 = read_likelihood_per_position(h1_underlying_position, r.genotype, r.qualities, r.graphLevel);
						double position_likelihood_h2 = read_likelihood_per_position(h2_underlying_position, r.genotype, r.qualities, r.graphLevel);
						
						double combined_LL =  Utilities::logAvg(position_likelihood_h1, position_likelihood_h2);
						double combined_L = exp(combined_LL);
						if(!((combined_L >= 0) && (combined_L <= 1)))
						{
							std::cerr << "h1_underlying_position" << ": " << h1_underlying_position << "\n";
							std::cerr << "h2_underlying_position" << ": " << h2_underlying_position << "\n";
							std::cerr << "r.genotype" << ": " << r.genotype << "\n";
							std::cerr << "r.qualities" << ": " << r.qualities << "\n";
							std::cerr << "r.graphLevel" << ": " << r.graphLevel << "\n\n";
						
							std::cerr << "position_likelihood_h1" << ": " << position_likelihood_h1 << "\n";
							std::cerr << "position_likelihood_h2" << ": " << position_likelihood_h2 << "\n";
							std::cerr << "combined_L" << ": " << combined_L << "\n";
							std::cerr << std::flush;
						}
						assert(combined_L >= 0);
						assert(combined_L <= 1);
						
						thisPosition_LL2 += combined_LL;
					}
					
					LL2 += thisPosition_LL2;
					
					/*
					
					auto compute_p_isUnderlyingAllele = [](std::string rGenotype, std::string rQualities, int rGraphLevel, std::string underlyingAllele) -> double {
					
						return exp(read_likelihood_per_position(underlyingAllele, rGenotype, rQualities, rGraphLevel));						
					};
					
					auto compute_p_isOneOfUnderlyingAllele = [](std::string rGenotype, std::string rQualities, int rGraphLevel, std::string underlyingAllele1, std::string underlyingAllele2) -> std::vector<double> {
					
						double a1 = exp(read_likelihood_per_position(underlyingAllele1, rGenotype, rQualities, rGraphLevel));
						double a2 = exp(read_likelihood_per_position(underlyingAllele2, rGenotype, rQualities, rGraphLevel));
						assert((a1 + a2) <= 1);
						double a3 = 1 - a1 - a2;
						
						std::vector<double> forReturn = {a1, a2, a3};
						return forReturn;
					};					
					
					double genericError = 0.01;
					double thisPosition_LL2 = 0;
					for(unsigned int rI = 0; rI < pileUpPerPosition.size(); rI++)
					{
						const oneExonPosition& r = pileUpPerPosition.at(rI);
						const std::string& readID = r.thisRead_ID;
					
						double thisRead_L = 0;
						if(h1_underlying_position == h2_underlying_position)
						{
							double p_is_underlyingAllele = compute_p_isUnderlyingAllele(r.genotype, r.qualities, r.graphLevel, h1_underlying_position);
							assert(p_is_underlyingAllele >= 0);
							assert(p_is_underlyingAllele <= 1);
							
							thisRead_L = genericError * (1 - p_is_underlyingAllele) + (1 - genericError) * (p_is_underlyingAllele);
							assert(thisRead_L >= 0);
							assert(thisRead_L <= 1);
						}
						else
						{
							std::vector<double> p_underlyingAlleles = compute_p_isOneOfUnderlyingAllele(r.genotype, r.qualities, r.graphLevel, h1_underlying_position, h2_underlying_position);
							
							assert(p_underlyingAlleles.at(0) >= 0);
							assert(p_underlyingAlleles.at(0) <= 1);							
							assert(p_underlyingAlleles.at(1) >= 0);
							assert(p_underlyingAlleles.at(1) <= 1);							
							assert(p_underlyingAlleles.at(2) >= 0);
							assert(p_underlyingAlleles.at(2) <= 1);		
							assert(abs(p_underlyingAlleles.at(0) + p_underlyingAlleles.at(1) + p_underlyingAlleles.at(2) - 1) < 1e-5);
							
							thisRead_L = genericError * p_underlyingAlleles.at(2) + (1 - genericError)*0.5 * p_underlyingAlleles.at(0) + (1 - genericError)*0.5 * p_underlyingAlleles.at(1);
							assert(thisRead_L >= 0);
							assert(thisRead_L <= 1);							
						}
						
						assert(thisRead_L > 0);
						assert(thisRead_L <= 1);
						double thisRead_LL = log(thisRead_L);
						
						thisPosition_LL2 += thisRead_LL;
					}
					
					
					
					LL2 += thisPosition_LL2;
					*/
					
					if(verbose)
					{
						std::cout << "\t\tAfter update LL2: " << LL2 << "\n" << std::flush;
					}				
										
				}
				
				ll_computed_until_position++;
				setPexternal = false;
			}

			void extendHaplotypes(const std::pair<std::string, std::string>& extension)
			{
				h1.push_back(extension.first);
				h2.push_back(extension.second);
				setPexternal = false;
			}
		};

		class haplotypeAlternatives {
		public:
			std::list<haplotypeAlternative> runningAlternatives;
			
			int completedLevel;

			haplotypeAlternatives()
			{
				completedLevel = -1;
			}

			void printAlternativesAndP(unsigned int pI)
			{	
				unsigned int aI = 0;

				std::cout << "Haplotype alternative print starting at position " << pI << "\n";
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					double P = alternativeIt->getP();
					std::pair<std::string, std::string> haplos = alternativeIt->getTrailingHaplotypeAlleles(pI);
					
					std::cout << "\t" << "Alternative " << aI << " relative likelihood " << P << "\n";
					std::cout << "\t\t" << haplos.first << "\n";
					std::cout << "\t\t" << haplos.second << "\n";
					
					std::cout << std::flush;
					aI++;
				}					
			}
			
			haplotypeAlternative getBestAlternative()
			{
				unsigned int l_max = 0;
				double ll_max = 0;

				unsigned int aI = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					if((alternativeIt == runningAlternatives.begin()) || (alternativeIt->getP() > ll_max))
					{
						l_max = aI;
						ll_max = alternativeIt->getP();
					}
					aI++;
				}

				haplotypeAlternative forReturn;
				aI = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					if(aI == l_max)
					{
						forReturn = *alternativeIt;
						break;
					}
					aI++;
				}

				return forReturn;
			}

			void processNewAlternatives(const std::vector<std::pair<std::string, std::string>>& alternatives)
			{
				if(completedLevel == -1)
				{
					for(unsigned int aI = 0; aI < alternatives.size(); aI++)
					{
						haplotypeAlternative A;
						A.extendHaplotypes(alternatives.at(aI));
						runningAlternatives.push_back(A);
					}
				}
				else
				{
					if(alternatives.size() == 1)
					{
						for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
						{
							alternativeIt->extendHaplotypes(alternatives.front());
						}
					}
					else
					{
						std::list<haplotypeAlternative> runningAlternatives_copy = runningAlternatives;
						for(unsigned int aI = 0; aI < alternatives.size(); aI++)
						{
							if(aI == 0)
							{
								for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
								{
									alternativeIt->extendHaplotypes(alternatives.front());
								}
							}
							else
							{
								for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives_copy.begin(); alternativeIt != runningAlternatives_copy.end(); alternativeIt++)
								{
									haplotypeAlternative newAlternative = *alternativeIt;
									newAlternative.extendHaplotypes(alternatives.at(aI));
									runningAlternatives.push_back(newAlternative);
								}
							}
						}
					}
				}

				completedLevel++;
			}

			void updateLikelihood(std::set<std::string>& alleles_from_h1, std::set<std::string>& alleles_from_h2, std::vector<oneExonPosition>& pileUpPerPosition)
			{
				unsigned int aI = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					if(completedLevel > 502)
					{
						// std::cout << "Update likelihood for alternative " << aI << "\n" << std::flush;
					}
					alternativeIt->updateLikelihood(alleles_from_h1, alleles_from_h2, pileUpPerPosition);
					
					aI++;
				}
				
				setNormalizedLikelihoods();
			}
			
			void setNormalizedLikelihoods()
			{
				bool verbose = ((completedLevel >= 515) && (completedLevel <= 520));
				verbose = false; // todo activate
				
				if(verbose)
					std::cout << "Enter setNormalizedLikelihoods()...\n" << std::flush;
					
				std::vector<double> scaled_combined_likelihoods_per_index_LL1;
				std::vector<double> scaled_combined_likelihoods_per_index_LL2;
				
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					scaled_combined_likelihoods_per_index_LL1.push_back(alternativeIt->getLL1());
					scaled_combined_likelihoods_per_index_LL2.push_back(alternativeIt->getLL2());					
				}
				
				convertLogVectorToP(scaled_combined_likelihoods_per_index_LL1);
				convertLogVectorToP(scaled_combined_likelihoods_per_index_LL2);
				
				normalizeVector(scaled_combined_likelihoods_per_index_LL1);
				normalizeVector(scaled_combined_likelihoods_per_index_LL2);
			
				std::vector<double> scaled_combined_likelihoods_per_index;
				for(unsigned int i = 0; i < scaled_combined_likelihoods_per_index_LL1.size(); i++)
				{
					assert(scaled_combined_likelihoods_per_index_LL1.at(i) >= 0);
					assert(scaled_combined_likelihoods_per_index_LL1.at(i) <= 1);
					assert(scaled_combined_likelihoods_per_index_LL2.at(i) >= 0);
					assert(scaled_combined_likelihoods_per_index_LL2.at(i) <= 1);
					
					double l_combined = 0.5 * scaled_combined_likelihoods_per_index_LL1.at(i) + 0.5 * scaled_combined_likelihoods_per_index_LL2.at(i);
					
					if(verbose)
						std::cout << "Alternative " << i << ": " <<  scaled_combined_likelihoods_per_index_LL1.at(i) << " x " << scaled_combined_likelihoods_per_index_LL2.at(i) << "\n" << std::flush;
					
					
					assert(l_combined >= 0);
					assert(l_combined <= 1);
					
					scaled_combined_likelihoods_per_index.push_back(l_combined);
				}	
				
				normalizeVector(scaled_combined_likelihoods_per_index);			

				unsigned int aI = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					alternativeIt->setRescaledP(scaled_combined_likelihoods_per_index.at(aI));
					aI++;
				}		
				assert(aI == scaled_combined_likelihoods_per_index.size());
			}

			std::map<std::string, double> getBackwardConfidence(unsigned int pI_start)
			{
				double max_LL = -1;
				int pI_stop = -1;

				assert(runningAlternatives.size() > 0);
				
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					int LL_available_to = alternativeIt->getLL_computedUntil();
					assert((int)LL_available_to >= (int)pI_start);
					if(alternativeIt == runningAlternatives.begin())
					{
						pI_stop = LL_available_to;
					}
					else
					{
						assert(pI_stop == LL_available_to);
					}

					if((alternativeIt == runningAlternatives.begin()) || (alternativeIt->getP() > max_LL))
					{
						max_LL = alternativeIt->getP();
					}
				}

				std::vector<double> relative_P;
				relative_P.resize(runningAlternatives.size());
				unsigned int aI = 0;
				double relative_P_sum = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					relative_P.at(aI) =  alternativeIt->getP(); //exp(alternativeIt->getLL1() - max_LL);
					relative_P_sum += relative_P.at(aI);
					aI++;
				}

				assert(relative_P_sum != 0);
				for(unsigned int i = 0; i < relative_P.size(); i++)
				{
					relative_P.at(i) = relative_P.at(i) / relative_P_sum;
				}

				std::map<std::string, double> forReturn;
				aI = 0;
				double check_sum = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					std::vector<std::string> h1;
					std::vector<std::string> h2;

					assert((int)pI_start <= pI_stop);
					for(int pI = pI_start; pI <= pI_stop; pI++)
					{
						std::pair<std::string, std::string> hA = alternativeIt->getHaplotypeAlleles(pI);
						h1.push_back(hA.first);
						h2.push_back(hA.second);
					}

					std::string h1_str = Utilities::join(h1, ";");
					std::string h2_str = Utilities::join(h2, ";");

					std::string h_comd = h1_str + "|" + h2_str;

					if(! forReturn.count(h_comd))
					{
						forReturn[h_comd] = 0;
					}

					forReturn.at(h_comd) += relative_P.at(aI);
					
					check_sum += relative_P.at(aI);
					assert(forReturn.at(h_comd) >= 0);
					if(!((forReturn.at(h_comd) - 1) <= 1e-4))
					{
						std::cerr << "! forReturn.at(h_comd) <= 1" << ": " << forReturn.at(h_comd) << "\n" << std::flush;
					}
					assert((forReturn.at(h_comd) - 1) <= 1e-4);
					
					aI++;
				}
				
				if(!((abs(1 - check_sum) < 1e-9)))
				{
					std::cerr << "! " << "(abs(1 - check_sum) < 1e-9)" << "\n";
					std::cerr << check_sum << "\n";
					
					std::cerr << "\n\n" << std::flush;
				}
				
				assert(abs(1 - check_sum) < 1e-9);

				return forReturn;
			}

			void pruneAlternatives()
			{
			
				bool verbose = ((completedLevel >= 589) && (completedLevel <= 596));
				verbose = false;
				
				int maximumPostpruning = 32;				

				std::vector<double> P_per_index;
				std::vector<std::pair<std::string, std::string>> trailingAlleles_per_index;
				double P_max = 0;
				unsigned int aI = 0;
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					double thisAlternative_P = alternativeIt->getP();
					// assert(thisAlternative_P >= 0);
					// assert(thisAlternative_P <= 1);
					
					P_per_index.push_back(thisAlternative_P);
					
					if((aI == 0) || (thisAlternative_P > P_max))
					{
						P_max = thisAlternative_P;
					}
					
					if(verbose)
						trailingAlleles_per_index.push_back(alternativeIt->getTrailingHaplotypeAlleles());
					
					aI++;
				}
				assert(P_max > 0);
				
				std::vector<unsigned int> sorted_indices;
				for(unsigned int i = 0; i < P_per_index.size(); i++)
				{
					sorted_indices.push_back(i);
				}
				
				std::sort(sorted_indices.begin(), sorted_indices.end(),
					[&P_per_index](unsigned int a, unsigned int b){return (P_per_index.at(a) < P_per_index.at(b));}
				);
				std::reverse(sorted_indices.begin(), sorted_indices.end());
				
				assert(P_per_index.at(sorted_indices.at(0)) == P_max);
				if(sorted_indices.size() > 1)
				{
					if(!(P_per_index.at(sorted_indices.at(0)) >= P_per_index.at(sorted_indices.at(1))))
					{
						std::cerr << "P_per_index.at(sorted_indices.at(0): " << P_per_index.at(sorted_indices.at(0)) << "\n";
						std::cerr << "P_per_index.at(sorted_indices.at(1)): " << P_per_index.at(sorted_indices.at(1)) << "\n";
						std::cerr << std::flush;
					}
					assert(P_per_index.at(sorted_indices.at(0)) >= P_per_index.at(sorted_indices.at(1)));
				}
				
				std::vector<int> P_per_index_sortIndex;
				P_per_index_sortIndex.resize(P_per_index.size(), -1);
				for(unsigned int i = 0; i < sorted_indices.size(); i++)
				{
					unsigned int originalIndex = sorted_indices.at(i);
					P_per_index_sortIndex.at(originalIndex) = i;
				}

				std::vector<bool> keep_per_index;
				for(unsigned int i = 0; i < P_per_index.size(); i++)
				{
					int sortedIndex = P_per_index_sortIndex.at(i);
					assert(sortedIndex >= 0);
					
					
					bool keep = true;
					if(sortedIndex >= maximumPostpruning)
					{
						keep = false;
					}
					
					double thisAlternative_P = P_per_index.at(i);
					if(P_max != 0)
					{
						if(thisAlternative_P != 0)
						{
							keep = (keep && ((log(P_max) - log(thisAlternative_P)) < 10));
						}
						else
						{
							keep = false;
						}
					}

					keep_per_index.push_back(keep);
					
					// std::cout << i << " / " << P_per_index.size() << ": " << keep << "\n" << std::flush ;
					
					if(verbose)
					{
						std::string trailing = "Ending in " + trailingAlleles_per_index.at(i).first + " / " + trailingAlleles_per_index.at(i).second;
						std::cout << "\t" << i << " P: " << thisAlternative_P << " (sorted " << sortedIndex << "/" << P_per_index_sortIndex.size() << ") - keep: " << keep << "\t\t" << trailing << "\n" << std::flush;
					}
				}
				
				aI = 0;
				std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin();
				while(alternativeIt != runningAlternatives.end())
				{
					std::list<haplotypeAlternative>::iterator thisIt = alternativeIt;
					alternativeIt++;

					if(keep_per_index.at(aI) == false)
					{
						runningAlternatives.erase(thisIt);
					}
					aI++;
				}

	
				assert(runningAlternatives.size() > 0);
				if(!((int)runningAlternatives.size() <= maximumPostpruning))
				{
					std::cerr << "! ((int)runningAlternatives.size() <= maximumPostpruning) " << "\n";
					std::cerr << "runningAlternatives.size(): " << runningAlternatives.size() << "\n";
					std::cerr << "maximumPostpruning: " << maximumPostpruning << "\n" << std::flush;
					
				}
				assert((int)runningAlternatives.size() <= maximumPostpruning);
				
				setNormalizedLikelihoods();
								
	
			/*	
				
				
				std::vector<std::pair<unsigned int, double> > likelihood_per_index;
				
				if(verbose)
					std::cout << "PRUNING:\n";
				for(std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin(); alternativeIt != runningAlternatives.end(); alternativeIt++)
				{
					likelihood_per_index.push_back(make_pair(aI, alternativeIt->getP()));
					// std::cout << "\t" << aI << "\t" << alternativeIt->getP() << "\n" << std::flush;
					aI++;
				}

				std::sort(likelihood_per_index.begin(), likelihood_per_index.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b){return (a.second < b.second);});
				
				if(likelihood_per_index.size() > 1)
				{
					if(!((likelihood_per_index.at(likelihood_per_index.size() - 1).second >= likelihood_per_index.at(likelihood_per_index.size() - 2).second)))
					{
						std::cout << "likelihood_per_index.at(likelihood_per_index.size() - 1).second" << ": " << likelihood_per_index.at(likelihood_per_index.size() - 1).second << "\n";
						std::cout << "likelihood_per_index.at(likelihood_per_index.size() - 2).second" << ": " << likelihood_per_index.at(likelihood_per_index.size() - 2).second << "\n" << std::flush;
					}
					assert(likelihood_per_index.at(likelihood_per_index.size() - 1).second >= likelihood_per_index.at(likelihood_per_index.size() - 2).second);
				}

				double ll_max = 10;
				double ll_min;
				double ll_cutoff = 1;

				std::vector<unsigned int> alternatives_sortedPosition;
				alternatives_sortedPosition.resize(runningAlternatives.size());
				for(unsigned int i = 0; i < likelihood_per_index.size(); i++)
				{
					int invertedI = likelihood_per_index.size() - i - 1;
					alternatives_sortedPosition.at(likelihood_per_index.at(i).first) = invertedI;

					if(i == 0)
					{
						ll_min = likelihood_per_index.at(i).second;
					}
					if(i == (likelihood_per_index.size() - 1))
					{
						ll_max = likelihood_per_index.at(i).second;
					}
				}
				assert(ll_max != 10);
				
				int post10Removal = 0;
				for(unsigned int i = 0; i < likelihood_per_index.size(); i++)
				{
					double LL = likelihood_per_index.at(i).second;
					assert(LL >= 0);
					assert(LL <= 1);
					if(log(LL) >= (log(ll_max) - 10))
					{
						post10Removal++;
					}
				}				
				
				bool activateOneThirdRemoval = (post10Removal > maximumPostpruning);
				int pruningFactor = 3;
				if(activateOneThirdRemoval)
				{
					pruningFactor = (runningAlternatives.size() / maximumPostpruning) + 1;
					if(pruningFactor < 3)
					{
						pruningFactor = 3;
					}
					for(unsigned int i = 0; i < likelihood_per_index.size(); i++)
					{
						// std::cout << "i: " << likelihood_per_index.at(i).second << "\n" << std::flush;
						int invertedI = likelihood_per_index.size() - i - 1;
						alternatives_sortedPosition.at(likelihood_per_index.at(i).first) = invertedI;

						if(i == 0)
						{
							ll_min = likelihood_per_index.at(i).second;
						}
						if(i == (likelihood_per_index.size() - 1))
						{
							ll_max = likelihood_per_index.at(i).second;
						}
						if((ll_cutoff == 1) && (invertedI > int((double)likelihood_per_index.size() * (1.0/(double)pruningFactor))))
						{
							ll_cutoff = likelihood_per_index.at(i).second;
						}
					}				
				}
				
				assert((ll_cutoff == 1) || (ll_max >= ll_cutoff));
				assert(ll_max >= ll_min);
				assert((ll_cutoff == 1) || (ll_cutoff >= ll_min));

				assert(ll_max >= 0);
				assert(ll_max <= 1);
				
				// std::cout << "Before pruning: " << runningAlternatives.size() << "\n\n";
				if(activateOneThirdRemoval)
				{
					// std::cout << "Likelihood pruning: remove 1 - 1/" << pruningFactor << " between min = " << ll_min << " and " << ll_max << " with cutoff " << ll_cutoff << "\n" << std::flush;
				}
				else
				{
					// std::cout << "Likelihood pruning: remove ll -10 between min = " << ll_min << " and " << ll_max << " without one-third-removal\n" << std::flush;				
				}
				
				aI = 0;
				std::list<haplotypeAlternative>::iterator alternativeIt = runningAlternatives.begin();
				while(alternativeIt != runningAlternatives.end())
				{
					std::list<haplotypeAlternative>::iterator thisIt = alternativeIt;
					alternativeIt++;
					
					// if( 0 || ((completedLevel > 500) && (completedLevel < 510)))
					// { 
						// std::pair<std::string, std::string> trailingAlleles = thisIt->getTrailingHaplotypeAlleles();
						// std::cout << "CompletedLevel " << completedLevel << ": Pruning " << aI << " with relative P " << thisIt->getP() << " vs max " << ll_max << " (log diff " << (log(ll_max) - log(alternativeIt->getP())) << "): Ending in " << trailingAlleles.first << " / " << trailingAlleles.second << "\n" << std::flush;						
					// }
					// if(completedLevel == 510)
					// {
						// assert( 1 == 0 );
					// }
					
					
					if((log(thisIt->getP()) < (log(ll_max) - 10)) || (activateOneThirdRemoval && ((int)alternatives_sortedPosition.at(aI) > int((double)likelihood_per_index.size() * (1.0/double(pruningFactor))))))
					{
						runningAlternatives.erase(thisIt);
					}
					aI++;
				}
				
				// std::cout << "After pruning: " << runningAlternatives.size() << "\n";
				
				*/						
			}
		};

		haplotypeAlternatives runningHaplotypes;

		std::vector<unsigned int> combined_sequences_confidenceIndex;
		std::vector<std::pair<unsigned int, unsigned int> > confidenceIntervals_boundaries;
		std::vector<std::map<std::string, double> > pruning_confidences;

		unsigned int currentConfidenceInterval_start = 0;
		unsigned int currentConfidenceInterval_index = 0;

		std::vector<std::vector<std::string> > considered_pairs_perPosition;
		std::vector<std::vector<std::string> > considered_h1_perPosition;
		std::vector<std::vector<std::string> > considered_h2_perPosition;
		
		for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
		{
			if((pI % 100) == 0)
			{
				// std::cout << "\r" << "pI = " << pI << std::flush;
			}
			
			bool verbose = ((pI >= 948) && (pI <= 952));
			verbose = false;
			
			std::set<std::string> alleles_from_h1;
			std::set<std::string> alleles_from_h2;
			std::set<std::string> alleles_from_reads;

			// compute sets of plausible alleles

			for(std::set<std::string>::iterator haplotypeIDit = starting_haplotypes_1_set.begin(); haplotypeIDit != starting_haplotypes_1_set.end(); haplotypeIDit++)
			{
				std::string haplotypeID = *haplotypeIDit;
				std::string allele = combined_sequences.at(haplotypeID).substr(pI, 1);
				alleles_from_h1.insert(allele);
			}
			for(std::set<std::string>::iterator haplotypeIDit = starting_haplotypes_2_set.begin(); haplotypeIDit != starting_haplotypes_2_set.end(); haplotypeIDit++)
			{
				std::string haplotypeID = *haplotypeIDit;
				std::string allele = combined_sequences.at(haplotypeID).substr(pI, 1);
				alleles_from_h2.insert(allele);
			}
			for(unsigned int pileUpI = 0; pileUpI < pileUpPerPosition.at(pI).size(); pileUpI++)
			{
				oneExonPosition& onePositionSpecifier = pileUpPerPosition.at(pI).at(pileUpI);
				if(onePositionSpecifier.genotype.find("N") == std::string::npos)
				{
					alleles_from_reads.insert(onePositionSpecifier.genotype);
				}				
				// alleles_from_reads.insert(onePositionSpecifier.genotype);
			}
			
			for(std::set<std::string>::iterator alleleIt = alleles_from_h1.begin(); alleleIt != alleles_from_h1.end(); alleleIt++)
			{
				std::string allele = *alleleIt;
				alleles_from_reads.insert(allele);
			}
			
			for(std::set<std::string>::iterator alleleIt = alleles_from_h2.begin(); alleleIt != alleles_from_h2.end(); alleleIt++)
			{
				std::string allele = *alleleIt;
				alleles_from_reads.insert(allele);
			}
						
			
			std::vector<std::string> alleles_from_reads_vec(alleles_from_reads.begin(), alleles_from_reads.end());

			std::vector<std::pair<std::string, std::string>> possible_haplotype_extensions;

			std::vector<std::string> considered_pairs_str;
			for(unsigned int a1 = 0; a1 < alleles_from_reads.size(); a1++)
			{
				for(unsigned int a2 = 0; a2 < alleles_from_reads.size(); a2++)
				{
					std::pair<std::string, std::string> thisCombination = make_pair(alleles_from_reads_vec.at(a1), alleles_from_reads_vec.at(a2));
					possible_haplotype_extensions.push_back(thisCombination);
					considered_pairs_str.push_back(thisCombination.first + "/" + thisCombination.second);
				}
			}
			
			considered_pairs_perPosition.push_back(considered_pairs_str);
			considered_h1_perPosition.push_back(std::vector<std::string>(alleles_from_h1.begin(), alleles_from_h1.end()));
			considered_h2_perPosition.push_back(std::vector<std::string>(alleles_from_h2.begin(), alleles_from_h2.end()));
			
			if(verbose)
			{
			
			std::cout << "\t\tPosition " << pI << " / " << combined_sequences_graphLevels.size() << ": " << alleles_from_reads.size() << " alleles from " << pileUpPerPosition.at(pI).size() << " reads [" << Utilities::join(alleles_from_reads_vec, ", ") << "] and " << possible_haplotype_extensions.size() << " possible haplotype extensions.\n" << std::flush;
			
			
			
			
			std::cout << "\t\t\tFrom reads:\n";
			for(std::set<std::string>::iterator alleleIt = alleles_from_reads.begin(); alleleIt != alleles_from_reads.end(); alleleIt++)
			{
				std::cout << "\t\t\t\t" << *alleleIt << "\n";
			}
			std::cout << "\n\t\t\tH1:\n";
			
			for(std::set<std::string>::iterator alleleIt = alleles_from_h1.begin(); alleleIt != alleles_from_h1.end(); alleleIt++)
			{
				std::cout << "\t\t\t\t" << *alleleIt << "\n";
			}
			std::cout << "\n\t\t\tH2:\n";
			for(std::set<std::string>::iterator alleleIt = alleles_from_h2.begin(); alleleIt != alleles_from_h2.end(); alleleIt++)
			{
				std::cout << "\t\t\t\t" << *alleleIt << "\n";
			}						
			}
			// if(possible_haplotype_extensions.size() == 0)
			// {	
				// std::set<std::string> alleles_from_haplotypes = alleles_from_h1;
				// alleles_from_haplotypes.insert(alleles_from_h2.begin(), alleles_from_h2.end());
				// assert(alleles_from_haplotypes.size() > 0);
				
				// std::vector<std::string> all_alleles(alleles_from_haplotypes.begin(), alleles_from_haplotypes.end());
			
				// std::cout << "\t\tPosition " << pI << " / " << combined_sequences_graphLevels.size() << ": No reads, use haplotype alleles as possible!\n" << std::flush;
			
				// for(unsigned int a1 = 0; a1 < all_alleles.size(); a1++)
				// {
					// for(unsigned int a2 = 0; a2 < all_alleles.size(); a2++)
					// {
						// std::pair<std::string, std::string> thisCombination = make_pair(all_alleles.at(a1), all_alleles.at(a2));
						// possible_haplotype_extensions.push_back(thisCombination);
					// }
				// }					
			// }	
			
			assert(possible_haplotype_extensions.size() > 0);
			
			// std::cout << "\t\tPosition " << pI << " / " << combined_sequences_graphLevels.size() << ": " << alleles_from_reads.size() << " alleles from " << pileUpPerPosition.at(pI).size() << " reads [" << Utilities::join(alleles_from_reads_vec, ", ") << "] and " << possible_haplotype_extensions.size() << " possible haplotype extensions.\n" << std::flush;
			
			// std::cout << "\t\t\tH1:\n";
			// for(std::set<std::string>::iterator alleleIt = alleles_from_h1.begin(); alleleIt != alleles_from_h1.end(); alleleIt++)
			// {
				// std::cout << "\t\t\t\t" << *alleleIt << "\n";
			// }
			// std::cout << "\n\t\t\tH2:\n";
			// for(std::set<std::string>::iterator alleleIt = alleles_from_h2.begin(); alleleIt != alleles_from_h2.end(); alleleIt++)
			// {
				// std::cout << "\t\t\t\t" << *alleleIt << "\n";
			// }			
			
			std::cout << std::flush;
			
			// if(pI == 507)
			// {
				// assert( 1 == 0 );
			// }
			
			// if(pI > 502)
			// {
				// std::cout << "Position " << pI << " BEFORE amending " << "\n";
				// runningHaplotypes.printAlternativesAndP(500);
			// }
				

			runningHaplotypes.processNewAlternatives(possible_haplotype_extensions);
			
			// if(pI > 500)
			// {
				// std::cout << "Position " << pI << " AFTER amending, BEFORE LL update " << "\n";
				// runningHaplotypes.printAlternativesAndP(500);
			// }
			
			runningHaplotypes.updateLikelihood(alleles_from_h1, alleles_from_h2, pileUpPerPosition.at(pI));

			// if(pI > 500)
			// {
				// std::cout << "Position " << pI << " AFTER LL update " << "\n";
				// runningHaplotypes.printAlternativesAndP(500);
			// }
			
			// std::cout << "\t\tPosition " << pI << " / " << combined_sequences_graphLevels.size() << ": " << runningHaplotypes.runningAlternatives.size() << " alternatives.\n" << std::flush;

			
			combined_sequences_confidenceIndex.push_back(currentConfidenceInterval_index);

			bool pruneHere = true;


			std::string graphLevelID = combined_sequences_locusIDs.at(pI);
			// std::cerr << "graphLevelID: " << graphLevelID << "\n" << std::flush;
			// assert( 1 == 0);

			std::string next_graphLevelID;
			if(pI < (combined_sequences_graphLevels.size() - 1))
			{
				next_graphLevelID = combined_sequences_locusIDs.at(pI+1);
								
				if(	graphLevel_2_AAid.count(graphLevelID) &&
					graphLevel_2_AAid.count(next_graphLevelID) &&
					(graphLevel_2_AAid.at(graphLevelID) == graphLevel_2_AAid.at(next_graphLevelID)) )
				{
					pruneHere = false;
				}
				
				// std::cerr << graphLevelID << " " << next_graphLevelID << " " << graphLevel_2_AAid[graphLevelID] << " " << graphLevel_2_AAid[next_graphLevelID] << "\n";				
			}			

			// std::cerr << "Pruning " << pI << " " << graphLevelID << ": " << pruneHere << " index " << combined_sequences_confidenceIndex.back() << "   # running: " << runningHaplotypes.runningAlternatives.size() <<  "\n" << std::flush; // 
			
			if(pruneHere)
			{
				unsigned int currentConfidenceInterval_stop = pI;
				confidenceIntervals_boundaries.push_back(make_pair(currentConfidenceInterval_start, currentConfidenceInterval_stop));
	
				// std::cerr << "\t\t" << currentConfidenceInterval_start << " " << currentConfidenceInterval_stop << "\n" << std::flush;

			
				std::map<std::string, double> localConfidences = runningHaplotypes.getBackwardConfidence(currentConfidenceInterval_start);
				pruning_confidences.push_back(localConfidences);


				runningHaplotypes.pruneAlternatives();

				// std::cout << "\t\t\t post-pruning: " << runningHaplotypes.runningAlternatives.size() << "\n" << std::flush;
				
				assert(runningHaplotypes.runningAlternatives.size() > 0);
				
				currentConfidenceInterval_index++;
				currentConfidenceInterval_start = currentConfidenceInterval_stop + 1;
					
				// if(pI > 500)
				// {
					// std::cout << "Position " << pI << " AFTER pruning " << "\n";
					// runningHaplotypes.printAlternativesAndP(500);
				// }
			}
			else
			{
				if(runningHaplotypes.runningAlternatives.size() > 9000)
				{
					// prune without telling anyone
					runningHaplotypes.pruneAlternatives();					
				}
			}
			
				if(pileUpPerPosition.at(pI).size())
				{
					// int exonPos = exonPosIt->first;
					std::vector<oneExonPosition> pileUp = pileUpPerPosition.at(pI);

					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(pI));
					fieldsPerLine.push_back(Utilities::ItoStr(pileUp.size()));

					std::vector<std::string> piledUpGenotypes;

					for(unsigned int pileUpI = 0; pileUpI < pileUp.size(); pileUpI++)
					{
						oneExonPosition piledPosition = pileUp.at(pileUpI);

						std::vector<std::string> qualities_as_strings;
						for(unsigned int qI = 0; qI < piledPosition.qualities.size(); qI++)
						{
							char qC = piledPosition.qualities.at(qI);
							int qC_i = qC;
							qualities_as_strings.push_back(Utilities::ItoStr(qC_i));
						}

						std::string pileUpString = piledPosition.genotype
							+ " (" + Utilities::join(qualities_as_strings, ", ") + ")"
							+ " ["
							// + Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "
							// + Utilities::DtoStr(piledPosition.thisRead_fractionOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_fractionOK) + " | "
							// + Utilities::ItoStr(piledPosition.pairs_strands_OK) + " "
							+ Utilities::DtoStr(piledPosition.pairs_strands_distance) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ_position) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ) + " "
							+ Utilities::DtoStr(piledPosition.mapQ_genomic) + " | "
							+ Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							+ Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "						
							+ piledPosition.thisRead_ID + " "
							+ piledPosition.pairedRead_ID
							+ "]";

						piledUpGenotypes.push_back(pileUpString);
					}

					fieldsPerLine.push_back(Utilities::join(piledUpGenotypes, ", "));


					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\n";
				}
				else
				{
					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(pI));
					fieldsPerLine.push_back(Utilities::ItoStr(0));
					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\n";				
				}				
		}
		
		pileUpStream.close();
		
		// std::cout << "\n";

		if(currentConfidenceInterval_start != combined_sequences_graphLevels.size())
		{
			unsigned int currentConfidenceInterval_stop = (combined_sequences_graphLevels.size() - 1);
			confidenceIntervals_boundaries.push_back(make_pair(currentConfidenceInterval_start, currentConfidenceInterval_stop));

			std::map<std::string, double> localConfidences = runningHaplotypes.getBackwardConfidence(currentConfidenceInterval_start);
			pruning_confidences.push_back(localConfidences);

			currentConfidenceInterval_index++;
		}

		assert(currentConfidenceInterval_index >= 1);
		assert(combined_sequences_confidenceIndex.size() == combined_sequences_graphLevels.size());
		assert(confidenceIntervals_boundaries.size() == (combined_sequences_confidenceIndex.back() + 1));
		assert(confidenceIntervals_boundaries.size() == pruning_confidences.size());


		std::string outputFN_bestHaplotypeGuess = outputDirectory + "/R2_haplotypes_bestguess_"+locus+".txt";
		std::ofstream haplotypeBestGuess_outputStream;
		haplotypeBestGuess_outputStream.open(outputFN_bestHaplotypeGuess.c_str());
		assert(haplotypeBestGuess_outputStream.is_open());

		haplotypeAlternative bestHaplotype = runningHaplotypes.getBestAlternative();

		std::vector< std::string > h1_forConfidence_perPosition;
		std::vector< std::string > h2_forConfidence_perPosition;
		std::vector< double > h1_confidence_perPosition;
		std::vector< double > h2_confidence_perPosition;

		std::vector< double > h1h2_gt_confidence_perPosition;
		
		for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
		{
			std::vector<std::string> h1_forConfidence;
			std::vector<std::string> h2_forConfidence;

			unsigned int confidenceIndex = combined_sequences_confidenceIndex.at(pI);

			for(unsigned int pI = confidenceIntervals_boundaries.at(confidenceIndex).first; pI <= confidenceIntervals_boundaries.at(confidenceIndex).second; pI++)
			{
				std::pair<std::string, std::string> hA = bestHaplotype.getHaplotypeAlleles(pI);
				h1_forConfidence.push_back(hA.first);
				h2_forConfidence.push_back(hA.second);
			}

			std::string h1_forConfidence_str = Utilities::join(h1_forConfidence, ";");
			std::string h2_forConfidence_str = Utilities::join(h2_forConfidence, ";");

			std::string h_forConfidence_combd = h1_forConfidence_str + "|" + h2_forConfidence_str;

			double h1_confidence = 0;
			double h2_confidence = 0;
			double gt_confidence = 0;
			
			for(std::map<std::string, double>::iterator confidenceIt = pruning_confidences.at(confidenceIndex).begin(); confidenceIt != pruning_confidences.at(confidenceIndex).end(); confidenceIt++)
			{
				std::string alleles_str = confidenceIt->first;
				std::vector<std::string> alleles = Utilities::split(alleles_str, "|");
				assert(alleles.size() == 2);
				
				if(alleles.at(0) == h1_forConfidence_str)
				{	
					h1_confidence += confidenceIt->second;
				}
				if(alleles.at(1) == h2_forConfidence_str)
				{	
					h2_confidence += confidenceIt->second;
				}
				
				if(((alleles.at(0) == h1_forConfidence_str) && (alleles.at(1) == h2_forConfidence_str)) ||
				   ((alleles.at(1) == h1_forConfidence_str) && (alleles.at(0) == h2_forConfidence_str)))
				{
					gt_confidence += confidenceIt->second;
				}
				
				
			}
			
			h1_confidence_perPosition.push_back(h1_confidence);
			h2_confidence_perPosition.push_back(h2_confidence);
			h1h2_gt_confidence_perPosition.push_back(gt_confidence);
			
			h1_forConfidence_perPosition.push_back(h1_forConfidence_str);
			h2_forConfidence_perPosition.push_back(h2_forConfidence_str);
		}
		
		assert(h1_confidence_perPosition.size() == combined_sequences_graphLevels.size());
		assert(h2_confidence_perPosition.size() == combined_sequences_graphLevels.size());
		assert(h1_forConfidence_perPosition.size() == combined_sequences_graphLevels.size());
		assert(h2_forConfidence_perPosition.size() == combined_sequences_graphLevels.size());
		
		std::vector< std::vector<std::string> > AAs_h1;
		std::vector< std::vector<double> > AA_confidences_h1;
		std::vector< std::vector<std::string> > AAs_h2;
		std::vector< std::vector<double>> AA_confidences_h2;

		AAs_h1.resize(combined_sequences_graphLevels.size());
		AA_confidences_h1.resize(combined_sequences_graphLevels.size());
		AAs_h2.resize(combined_sequences_graphLevels.size());
		AA_confidences_h2.resize(combined_sequences_graphLevels.size());
		
		auto processStretch = [&](int start, int stop) -> void {
			fillAS();
			
			int L = stop - start + 1;
			assert(L >= 3);
			
			std::string h1;
			std::string h2;
			
			std::vector<int> h1_origins;
			std::vector<int> h2_origins;

			std::vector<double> h1_confidences;
			std::vector<double> h2_confidences;
			
			std::vector<int> h1_confidences_indices;
			std::vector<int> h2_confidences_indices;
			
			for(int i = start; i <= stop; i++)
			{
				std::pair<std::string, std::string> alleles = bestHaplotype.getHaplotypeAlleles(i);
				h1 += alleles.first;
				h2 += alleles.second;
				
				for(unsigned int j = 0; j < alleles.first.length(); j++)
				{
					h1_origins.push_back(i);
					h1_confidences.push_back(h1_confidence_perPosition.at(i));
					h1_confidences_indices.push_back(combined_sequences_confidenceIndex.at(i));
					
					
				}
				
				for(unsigned int j = 0; j < alleles.second.length(); j++)
				{
					h2_origins.push_back(i);
					h2_confidences.push_back(h2_confidence_perPosition.at(i));					
					h2_confidences_indices.push_back(combined_sequences_confidenceIndex.at(i));							
				}
			}
			
			assert(h1.length() >= 3);
			assert(h2.length() >= 3);
			assert(h1_origins.size() == h1.length());
			assert(h2_origins.size() == h2.length());
			
			std::string AA_string;
			int AA_start = 0;
			double confidence = 1;
			
			std::set<int> includedConfidences;
			
			for(int i = 0; i < (int)h1.length(); i++)
			{
				std::string c = h1.substr(i, 1);
				if(c != "_")
				{
					AA_string += c;
				}
				
				if(includedConfidences.count(h1_confidences_indices.at(i)) == 0)
				{
					includedConfidences.insert(h1_confidences_indices.at(i));
					confidence *= h1_confidences.at(i);
				}
				
				if(AA_string.length() == 3)
				{
					int pI_pos = h1_origins.at(AA_start);
					std::string correspondingAA = "X";
					if(codon2AS.count(AA_string))
					{
						correspondingAA = codon2AS.at(AA_string);
					}
					correspondingAA += "[" + AA_string + "]";
					
					AAs_h1.at(pI_pos).push_back(correspondingAA);
					AA_confidences_h1.at(pI_pos).push_back(confidence);
					
					AA_start = i + 1;
					AA_string = "";
					confidence = 1;
					includedConfidences.clear();
				}
			}
			
			AA_string = "";			
			AA_start = 0;
			confidence = 1;
			includedConfidences.clear();
			
			for(int i = 0; i < (int)h2.length(); i++)
			{
				std::string c = h2.substr(i, 1);
				if(c != "_")
				{
					AA_string += c;
				}
				
				if(includedConfidences.count(h2_confidences_indices.at(i)) == 0)
				{
					includedConfidences.insert(h2_confidences_indices.at(i));
					confidence *= h2_confidences.at(i);
				}
				
				if(AA_string.length() == 3)
				{
					int pI_pos = h2_origins.at(AA_start);
					std::string correspondingAA = "X";
					if(codon2AS.count(AA_string))
					{
						correspondingAA = codon2AS.at(AA_string);
					}
					correspondingAA += "[" + AA_string + "]";
					
					AAs_h2.at(pI_pos).push_back(correspondingAA);
					AA_confidences_h2.at(pI_pos).push_back(confidence);
					
					AA_start = i + 1;
					AA_string = "";
					confidence = 1;
					includedConfidences.clear();					
				}
			}			
		};
		
		bool inAA = false;
		int AAstretch_pI_start = -1;
		int AAstretch_pI_stop = -1;
		for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
		{
			int pI_stop = pI;
			while((
				(pI_stop + 1) < (int)combined_sequences_graphLevels.size()) &&
				(combined_sequences_confidenceIndex.at(pI_stop + 1) == combined_sequences_confidenceIndex.at(pI)))
			{
				pI_stop++;
			}
			
			int L = pI_stop - pI + 1;
			
			// std::cerr << pI << " " << pI_stop << " " << combined_sequences_locusIDs.at(pI) << "\n" << std::flush; 
			
			if(graphLevel_2_AAid.count(combined_sequences_locusIDs.at(pI)))
			{
				// std::cerr << "\t" << graphLevel_2_AAid.at(combined_sequences_locusIDs.at(pI)) << "\n" << std::flush;
			}
			
			if(L > 1)
			{
				if(!((L >= 3)))
				{	
					std::cerr << "Error\n";
					std::cerr << "L" << ": " << L << "\n";
					std::cerr << "inAA" << ": " << inAA << "\n";	
					std::cerr << "AAstretch_pI_start" << ": " << AAstretch_pI_start << "\n" << std::flush;			

					std::string graphLevelID = combined_sequences_locusIDs.at(pI);
					
					std::cerr << "graphLevelID" << ": " << graphLevelID << "\n" << std::flush;
					
					std::cerr << "graphLevel_2_AAid" << ": " << graphLevel_2_AAid.at(graphLevelID) << "\n" << std::flush;
					
					
					std::cerr << std::flush;
				}
				assert(L >= 3);
				if(! inAA)
				{
					inAA = true;
					AAstretch_pI_start = pI;
				}
				AAstretch_pI_stop = pI_stop;				
			}
			else
			{
				if(inAA)
				{
					processStretch(AAstretch_pI_start, AAstretch_pI_stop);
				}
				inAA = false;
			}
			
			pI += (L - 1);
		}
		if(inAA)
		{
			AAstretch_pI_stop = combined_sequences_graphLevels.size() - 1;
			processStretch(AAstretch_pI_start, AAstretch_pI_stop);
		}
		
		for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
		{
			std::pair<std::string, std::string> alleles = bestHaplotype.getHaplotypeAlleles(pI);

			std::vector<std::string> fieldsForLine;

			fieldsForLine.push_back(Utilities::ItoStr(pI));
			fieldsForLine.push_back(Utilities::ItoStr(combined_sequences_graphLevels.at(pI)));
			fieldsForLine.push_back(combined_sequences_graphLevels_individualType.at(pI));
			fieldsForLine.push_back(Utilities::ItoStr(combined_sequences_graphLevels_individualTypeNumber.at(pI)));
			fieldsForLine.push_back(Utilities::ItoStr(combined_sequences_graphLevels_individualPosition.at(pI)));

			fieldsForLine.push_back(alleles.first);
			fieldsForLine.push_back(alleles.second);

			fieldsForLine.push_back(Utilities::DtoStr(h1h2_gt_confidence_perPosition.at(pI)));
			
			unsigned int confidenceIndex = combined_sequences_confidenceIndex.at(pI);
			fieldsForLine.push_back(Utilities::ItoStr(confidenceIndex));

			fieldsForLine.push_back(h1_forConfidence_perPosition.at(pI));
			
			fieldsForLine.push_back(Utilities::DtoStr(h1_confidence_perPosition.at(pI)));
			
			fieldsForLine.push_back(h2_forConfidence_perPosition.at(pI));

			fieldsForLine.push_back(Utilities::DtoStr(h2_confidence_perPosition.at(pI)));
			
			std::vector<std::string> h1_AAs = AAs_h1.at(pI);
			std::vector<std::string> h2_AAs = AAs_h2.at(pI);
			
			std::vector<double> h1_AAs_confidences = AA_confidences_h1.at(pI);
			std::vector<double> h2_AAs_confidences = AA_confidences_h2.at(pI);
			std::vector<std::string> h1_AAs_confidences_str;
			std::vector<std::string> h2_AAs_confidences_str;

			for(unsigned int i = 0; i < h1_AAs_confidences.size(); i++)
			{
				h1_AAs_confidences_str.push_back(Utilities::DtoStr(h1_AAs_confidences.at(i)));
			}

			for(unsigned int i = 0; i < h2_AAs_confidences.size(); i++)
			{
				h2_AAs_confidences_str.push_back(Utilities::DtoStr(h2_AAs_confidences.at(i)));
			}
			
			assert(h1_AAs_confidences_str.size() == h1_AAs.size());
			assert(h2_AAs_confidences_str.size() == h2_AAs.size());

			
			std::vector<std::string> pileUpVec;
			for(unsigned int pileUpI = 0; pileUpI < pileUpPerPosition.at(pI).size(); pileUpI++)
			{
				//pileUpVec.push_back(pileUpPerPosition.at(pI).at(pileUpI).genotype+"["+pileUpPerPosition.at(pI).at(pileUpI).thisRead_ID+"]"); // .qualities
				pileUpVec.push_back(pileUpPerPosition.at(pI).at(pileUpI).genotype);
			}
			std::string pileUpString = Utilities::join(pileUpVec, "-");
			
			fieldsForLine.push_back(Utilities::ItoStr(pileUpPerPosition.at(pI).size()));
			fieldsForLine.push_back(pileUpString);
			
			fieldsForLine.push_back(Utilities::join(considered_pairs_perPosition.at(pI), ";"));

			fieldsForLine.push_back(Utilities::join(considered_h1_perPosition.at(pI), ";"));
			
			fieldsForLine.push_back(Utilities::join(considered_h2_perPosition.at(pI), ";"));			
						
			fieldsForLine.push_back(Utilities::join(h1_AAs, ";"));
			fieldsForLine.push_back(Utilities::join(h1_AAs_confidences_str, ";"));
			fieldsForLine.push_back(Utilities::join(h2_AAs, ";"));
			fieldsForLine.push_back(Utilities::join(h2_AAs_confidences_str, ";"));
			
			haplotypeBestGuess_outputStream << Utilities::join(fieldsForLine, "\t") << "\n";
		}

		haplotypeBestGuess_outputStream.close();

		// find closest 6- and 8-digit matches

		assert(bestHaplotype.getH1().size() == combined_sequences_graphLevels.size());
		assert(bestHaplotype.getH2().size() == combined_sequences_graphLevels.size());

		std::map<std::string, std::vector<std::string>> combined_sequences_6_8_digit_inference;
		read_HLA_alleles_for_6_8_digits(graphDir, locus, combined_sequences_graphLevels, graphLocus_2_levels, combined_sequences_6_8_digit_inference);

		std::string outputFN_similarReferenceHaplotypes = outputDirectory + "/R2_haplotypes_similarReferenceHaplotypes_"+locus+".txt";
		std::ofstream similarReferenceHaplotypesStream;
		similarReferenceHaplotypesStream.open(outputFN_similarReferenceHaplotypes.c_str());
		assert(similarReferenceHaplotypesStream.is_open());

		similarReferenceHaplotypesStream 	<< "Locus" << "\t"
											<< "6DigitOnly" << "\t"
											<< "Haplotype" << "\t"
											<< "Alternative" << "\t"
											<< "ReferenceAllele" << "\t"
											<< "AbsoluteDifferences" << "\t"
											<< "RelativeDifference" << "\n";

		for(unsigned int haplotypeI = 1; haplotypeI <= 2; haplotypeI++)
		{
			for(unsigned int dig_6 = 0; dig_6 <= 1; dig_6++)
			{
				std::vector<std::string> h = bestHaplotype.getH(haplotypeI);

				int comparisons = 0;
				int differences = 0;

				if(dig_6 == 1)
				{
					assert(h.size() == combined_sequences_graphLevels.size());
					int masked_characters = 0;
					for(unsigned int pI = 0; pI < combined_sequences_graphLevels.size(); pI++)
					{
						if(combined_sequences_graphLevels_individualType.at(pI) != "exon")
						{
							h.at(pI) = "*";
							masked_characters++;
						}
					}
					assert(masked_characters > 0); // deactivate later for loci without exons
				}
				std::map<std::string, double> referenceHaplotypes_absoluteDifferences;
				std::map<std::string, double> referenceHaplotypes_differenceRatio;

				for(std::map<std::string, std::vector<std::string>>::iterator referenceHaplotypeIt = combined_sequences_6_8_digit_inference.begin(); referenceHaplotypeIt != combined_sequences_6_8_digit_inference.end(); referenceHaplotypeIt++)
				{
					compute_weird_Edit_distance(h, referenceHaplotypeIt->second, comparisons, differences);

					double diff_ratio = 0;
					if(comparisons != 0)
					{
						diff_ratio = (double)differences/(double)comparisons;
					}

					referenceHaplotypes_absoluteDifferences[referenceHaplotypeIt->first] = differences;
					referenceHaplotypes_differenceRatio[referenceHaplotypeIt->first] = diff_ratio;
				}

				std::vector<std::string> referenceHaplotypes_keys;
				for(std::map<std::string, double>::iterator completeTypeIt = referenceHaplotypes_absoluteDifferences.begin(); completeTypeIt != referenceHaplotypes_absoluteDifferences.end(); completeTypeIt++)
				{
					referenceHaplotypes_keys.push_back(completeTypeIt->first);
				}
				std::sort(referenceHaplotypes_keys.begin(), referenceHaplotypes_keys.end(), [&](std::string a, std::string b){return (referenceHaplotypes_absoluteDifferences.at(a) < referenceHaplotypes_absoluteDifferences.at(b));});

				if(referenceHaplotypes_keys.size() > 1)
				{
					assert(referenceHaplotypes_absoluteDifferences.at(referenceHaplotypes_keys.at(0)) <= referenceHaplotypes_absoluteDifferences.at(referenceHaplotypes_keys.at(1)));
				}

				int print_index = (referenceHaplotypes_keys.size() > 5) ? 5 : referenceHaplotypes_keys.size();


				for(int printI = 1; printI <= print_index; printI++)
				{
					std::string printAllele = referenceHaplotypes_keys.at(printI-1);

					similarReferenceHaplotypesStream 	<< locus << "\t"
														<< dig_6 << "\t"
														<< haplotypeI << "\t"
														<< printI << "\t"
														<< printAllele << "\t"
														<< referenceHaplotypes_absoluteDifferences.at(printAllele) << "\t"
														<< referenceHaplotypes_differenceRatio.at(printAllele) << "\n";

					if(printI == 1)
					{
						((dig_6 == 0) ? bestGuess_8dig_outputStream : bestGuess_6dig_outputStream )
														<< locus << "\t"
														<< haplotypeI << "\t"
														<< printAllele << "\t"
														<< referenceHaplotypes_absoluteDifferences.at(printAllele) << "\t"
														<< referenceHaplotypes_differenceRatio.at(printAllele) << "\n";
					}
				}
			}
		}

		similarReferenceHaplotypesStream.close();

	}


	std::string outputFN_parameters = outputDirectory + "/R2_parameters.txt";
	std::ofstream outputFN_parameters_outputStream;
	outputFN_parameters_outputStream.open(outputFN_parameters.c_str());
	assert(outputFN_parameters_outputStream.is_open());
	outputFN_parameters_outputStream << "loci_str" << " = " << loci_str << "\n";
	outputFN_parameters_outputStream << "starting_haplotypes_perLocus_1_str" << " = " << starting_haplotypes_perLocus_1_str << "\n";
	outputFN_parameters_outputStream << "starting_haplotypes_perLocus_2_str" << " = " << starting_haplotypes_perLocus_2_str << "\n";
	outputFN_parameters_outputStream << "veryConservativeReadLikelihoods" << " = " << veryConservativeReadLikelihoods << "\n";
	outputFN_parameters_outputStream.close();
}

void HLATypeInference(std::string alignedReads_file, std::string graphDir, std::string outputDir, std::string sampleName, bool restrictToFullHaplotypes, std::string& forReturn_lociString, std::string& forReturn_starting_haplotype_1, std::string& forReturn_starting_haplotype_2, bool longUnpairedReads, bool MiSeq250bp)
{
	std::string graph = graphDir + "/graph.txt";
	assert(Utilities::fileReadable(graph));

	// define loci
	std::vector<std::string> loci;
	//std::vector<std::string> loci = {"A"}; // todo activate later

	bool allLoci = true;
	if(allLoci)
	{
		loci = {"A", "B", "C", "DQA1", "DQB1", "DRB1", "DPA1", "DPB1", "DRA", "DRB3", "DRB4", "F", "G", "H", "J", "K", "L", "V"};
		// loci = {"E", "F", "G", "H", "J", "K", "L", "V"};
		// loci = {"F", "G", "H", "J", "K", "L", "V"};
	}
	else
	{
		loci = {"A", "B", "C", "DQA1", "DQB1", "DRB1"};
	}
	
	int threshold_reportColumn_forPresenceOfUnaccountedAlleles_minCoverage = 30;
	double threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction = 0.2;
	bool highCoverage_filter_alleles = true;
	int highCoverage_minCoverage = 100;
	double highCoverage_minAlleleFreq = 0.2;

	forReturn_lociString = Utilities::join(loci, ",");

	// define locus -> exon
	std::map<std::string, std::vector<std::string> > loci_2_exons;
	fill_loci_2_exons(loci_2_exons);

	int HLATypeInference_totalBases_used = 0;
	int HLAtypeInference_totalColumns = 0;
	
	// function to find right exon file
	std::vector<std::string> files_in_graphDir = filesInDirectory(graphDir);
	auto find_file_for_exon = [&](std::string locus, std::string exon) -> std::string
	{
		std::string forReturn;
		for(unsigned int fI = 0; fI < files_in_graphDir.size(); fI++)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(files_in_graphDir.at(fI), "_");
			// if(!(split_by_underscore.size() >= 3))
			// {
				// std::cerr << "Can't split according to underscores: " << files_in_graphDir.at(fI) << "\n" << std::flush;
			// }
			// assert(split_by_underscore.size() >= 3);

			if(split_by_underscore.size() >= 4)
			{
				if(split_by_underscore.at(split_by_underscore.size()-4).substr(split_by_underscore.at(split_by_underscore.size()-4).length() - (locus.length()+1)) == ("/"+locus))
				{
					if((split_by_underscore.at(split_by_underscore.size()-2)+"_"+split_by_underscore.at(split_by_underscore.size()-1)) == (exon + ".txt"))
					{
						forReturn = files_in_graphDir.at(fI);
					}
				}
			}
		}
		if(! forReturn.length())
		{
			std::cerr << "find_file_for_exon -- problem -- " << locus << " -- " << exon << "\n";
			std::cerr << Utilities::join(files_in_graphDir, " ") << "\n" << std::flush;
			assert(forReturn.length());
		}
		return forReturn;
	};

	// translate location IDs to graph levels
	std::vector<std::string> graphLoci = readGraphLoci(graphDir);
	std::map<std::string, unsigned int> graphLocus_2_levels;
	for(unsigned int i = 0; i < graphLoci.size(); i++)
	{
		std::string locusID = graphLoci.at(i);
		assert(graphLocus_2_levels.count(locusID) == 0);
		graphLocus_2_levels[locusID] = i;
	}

	// load reads
	std::cout << Utilities::timestamp() << "HLATypeInference(..): Load reads.\n" << std::flush;
	std::vector<std::pair<seedAndExtend_return_local, seedAndExtend_return_local>> alignments_paired;
	std::vector<seedAndExtend_return_local> alignments_unpaired;

	std::vector<oneReadPair> alignments_originalReads_paired;
	std::vector<oneRead> alignments_originalReads_unpaired;

	double insertSize_mean;
	double insertSize_sd;

	if(longUnpairedReads)
	{
		read_longReadAlignments_fromFile(alignedReads_file, alignments_unpaired, alignments_originalReads_unpaired);
		assert(alignments_unpaired.size() == alignments_originalReads_unpaired.size());

		std::cout << Utilities::timestamp() << "HLATypeInference(..): Load reads -- done. Have " << alignments_unpaired.size() << " unpaired reads, long unpaired reads: " << longUnpairedReads << "\n" << std::flush;
	}
	else
	{
		read_shortReadAlignments_fromFile(alignedReads_file, alignments_paired, alignments_originalReads_paired, insertSize_mean, insertSize_sd);
		assert(alignments_paired.size() == alignments_originalReads_paired.size());

		std::cout << Utilities::timestamp() << "HLATypeInference(..): Load reads -- done. Have " << alignments_paired.size() << " read pairs, long unpaired reads: " << longUnpairedReads << ", IS mean " << insertSize_mean << " / sd " << insertSize_sd << ".\n" << std::flush;
	}

	int k_for_kMer_index = 31;
	std::unordered_map<std::string, int> kMer_counts;
	// build kMer index of reads
	{
		auto add_sequence_to_index = [&kMer_counts, &k_for_kMer_index](const std::string& S) -> void {
			std::vector<std::string> kMers = partitionStringIntokMers(S, k_for_kMer_index);
			for(const auto& kMer : kMers)
			{
				std::string kMerKey = kMer_canonical_representation(kMer);
				if(kMer_counts.count(kMerKey) == 0)
				{
					kMer_counts[kMerKey] = 0;
				}
				kMer_counts.at(kMerKey)++;
			}
		};

		for(auto rP : alignments_originalReads_paired)
		{
			add_sequence_to_index(rP.reads.first.sequence);
			add_sequence_to_index(rP.reads.second.sequence);
		}

		for(auto uP : alignments_originalReads_unpaired)
		{
			add_sequence_to_index(uP.sequence);
		}
		
		std::cout << "Have " << kMer_counts.size() << " " << k_for_kMer_index << "-mers in index.\n" << std::flush;
	}

	// read alignment statistics

	int alignmentStats_strandsValid = 0;
	int alignments_paired_perfect = 0;
	int alignments_paired_oneReadPerfect = 0;
	int alignments_unpaired_perfect = 0;
	int alignmentStats_strandsValid_and_distanceOK = 0;
	std::vector<double> alignmentStats_strandsValid_distances;
	double alignmentStats_paired_fractionOK_sum = 0;
	double alignmentStats_unpaired_fractionOK_sum = 0;

	for(unsigned int alignmentI = 0; alignmentI < alignments_paired.size(); alignmentI++)
	{

		std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignedReadPair = alignments_paired.at(alignmentI);
		
		if(alignedReadPair_strandsValid(alignedReadPair))
		{
			alignmentStats_strandsValid++;
			double pairsDistance = alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair);
			alignmentStats_strandsValid_distances.push_back(pairsDistance);
			if(abs(pairsDistance - insertSize_mean) <= (5 * insertSize_sd))
			{
				alignmentStats_strandsValid_and_distanceOK++;
			}
		}

		double fractionOK_1 = alignmentFractionOK(alignedReadPair.first);
		double fractionOK_2 = alignmentFractionOK(alignedReadPair.second);

		if(fractionOK_1 == 1)
		{
			alignments_paired_perfect++;
		}
		if(fractionOK_2 == 1)
		{
			alignments_paired_perfect++;
		}
		if((fractionOK_1 == 1) || (fractionOK_2 == 1))
		{
			alignments_paired_oneReadPerfect++;
		}

		alignmentStats_paired_fractionOK_sum += fractionOK_1;
		alignmentStats_paired_fractionOK_sum += fractionOK_2;
	}
	

	for(unsigned int alignmentI = 0; alignmentI < alignments_unpaired.size(); alignmentI++)
	{

		seedAndExtend_return_local& alignedRead = alignments_unpaired.at(alignmentI);


		double fractionOK_1 = alignmentFractionOK(alignedRead);

		if(fractionOK_1 == 1)
		{
			alignments_unpaired_perfect++;
		}

		alignmentStats_unpaired_fractionOK_sum += fractionOK_1;
	}

	std::pair<double, double> alignmentStats_distance_meanMedian = meanMedian(alignmentStats_strandsValid_distances);
	double alignmentStats_paired_fractionOK_avg = (alignments_paired.size() > 0) ? (alignmentStats_paired_fractionOK_sum / (2.0* (double)alignments_paired.size())) : 0;
	double alignmentStats_unpaired_fractionOK_avg = (alignments_unpaired.size() > 0) ? (alignmentStats_unpaired_fractionOK_sum / (2.0* (double)alignments_unpaired.size())) : 0;


	if(! Utilities::directoryExists(outputDir + "/hla"))
	{
		Utilities::makeDir( outputDir + "/hla");
	}

	std::string outputDirectory = outputDir + "/hla/" + sampleName; 

	if(! Utilities::directoryExists(outputDirectory))
	{
		Utilities::makeDir(outputDirectory);
	}



	std::ofstream summaryStatisticsStream;
	std::string summaryStatisticsFilename = outputDirectory + "/" + "summaryStatistics.txt";
	summaryStatisticsStream.open(summaryStatisticsFilename.c_str());
	assert(summaryStatisticsStream.is_open());
	summaryStatisticsStream << "\nRead alignment statistics:\n";
	summaryStatisticsStream << "\t - Total number (paired) alignments:                 " << alignments_paired.size() << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK:                  " << alignmentStats_strandsValid << " (" << printPerc(alignmentStats_strandsValid, alignments_paired.size()) << "%)\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK && distance OK:   " << alignmentStats_strandsValid_and_distanceOK << " (" << printPerc(alignmentStats_strandsValid_and_distanceOK, alignments_paired.size()) << "%)\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK, mean distance:   " << alignmentStats_distance_meanMedian.first << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK, median distance: " << alignmentStats_distance_meanMedian.second << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, average fraction alignment OK:   " << alignmentStats_paired_fractionOK_avg << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, at least one alignment perfect:   " << alignments_paired_oneReadPerfect << "\n";
	summaryStatisticsStream << "\t\t - Single alignments, perfect (total):   " << alignments_paired_perfect << " (" << alignments_paired.size()*2 << ")\n";
	summaryStatisticsStream << "\t - Total number (unpaired) alignments:                 " << alignments_unpaired.size() << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, average fraction alignment OK:   " << alignmentStats_unpaired_fractionOK_avg << "\n";
	summaryStatisticsStream << "\t\t - Single alignments, perfect (total):   " << alignments_unpaired_perfect << " (" << alignments_unpaired.size()*2 << ")\n";
	summaryStatisticsStream.close();

	std::string outputFN_bestGuess = outputDirectory + "/R1_bestguess.txt";
	std::ofstream bestGuess_outputStream;
	bestGuess_outputStream.open(outputFN_bestGuess.c_str());
	assert(bestGuess_outputStream.is_open());
	std::string fieldName_unaccountedAlleles = "NColumns_UnaccountedAllele_fGT" + Utilities::DtoStr(threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction);

	// bestGuess_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << "LocusMinimumColumnErrorP" << "\n";
	bestGuess_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << fieldName_unaccountedAlleles << "\n";

	std::string outputFN_bestGuess_G = outputDirectory + "/R1_bestguess_G.txt";
	std::ofstream bestGuess_G_outputStream;
	bestGuess_G_outputStream.open(outputFN_bestGuess_G.c_str());
	assert(bestGuess_G_outputStream.is_open());
	// bestGuess_G_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << "LocusMinimumColumnErrorP" << "\t" << "perfectG" << "\n";
	bestGuess_G_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << fieldName_unaccountedAlleles << "\t" << "perfectG" << "\n";

	std::vector<std::string> forReturn_starting_haplotype_1_vec;
	std::vector<std::string> forReturn_starting_haplotype_2_vec;

	for(unsigned int locusI = 0; locusI < loci.size(); locusI++)
	{
		std::string locus = loci.at(locusI);
		std::set<std::string> utilized_reads;
		
		int HLATypeInference_thisLocus_bases_used = 0;

		std::string outputFN_allPairs = outputDirectory + "/R1_PP_"+locus+"_pairs.txt";
		std::string outputFN_columnError = outputDirectory + "/R1_columnIncompatibilities_"+locus+".txt";

		std::cout << Utilities::timestamp() << "HLATypeInference(..): Making inference for " << locus << "\n" << std::flush;

		std::vector<int> combined_exon_sequences_graphLevels;
		std::vector<int> combined_exon_sequences_graphLevels_individualExon;
		std::vector<int> combined_exon_sequences_graphLevels_individualExonPosition;

		std::vector<std::string> combined_exon_sequences_locusIDs;
		std::map<std::string, std::string> combined_exon_sequences;

		std::set<std::string> completeDefinedTypes;
		if(restrictToFullHaplotypes)
		{
			completeDefinedTypes = getCompletelyDefinedHLAAlleles(graphDir, locus);
			std::cout << Utilities::timestamp() << "HLATypeInference(..): restrictToFullHaplotypes in force, restrict to " << completeDefinedTypes.size() << " types.\n" << std::flush;
		}

		int thisLocus_totalColumns = 0;

		std::map<int, int> exon_lengths;
		
		for(unsigned int exonI = 0; exonI < loci_2_exons.at(locus).size(); exonI++)
		{
			std::string exonID = loci_2_exons.at(locus).at(exonI);
			std::cout << Utilities::timestamp() << "\tLocus" << locus << ", exon " << exonID << "\n" << std::flush;

			std::string exonFile = find_file_for_exon(locus, exonID);
			if(! Utilities::fileReadable(exonFile))
			{
				std::cerr << "HLATypeInference(..): Locus " << locus << ", exon " << exonID << ": Can't read file " << exonFile << "\n";
			}
			assert(Utilities::fileReadable(exonFile));

			std::ifstream exonInputStream;
			exonInputStream.open(exonFile.c_str());
			assert(exonInputStream.is_open());
			std::vector<std::string> exon_lines;
			while(exonInputStream.good())
			{
				std::string line;
				std::getline(exonInputStream, line);
				Utilities::eraseNL(line);
				exon_lines.push_back(line);
			}
			exonInputStream.close();

			std::string firstLine = exon_lines.at(0);
			std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
			assert(firstLine_fields.at(0) == "IndividualID");

			std::vector<std::string> exon_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
			std::string first_graph_locusID = exon_level_names.front();
			std::string last_graph_locusID = exon_level_names.back();

			assert(graphLocus_2_levels.count(first_graph_locusID));
			assert(graphLocus_2_levels.count(last_graph_locusID));

			unsigned int first_graph_level = graphLocus_2_levels.at(first_graph_locusID);
			unsigned int last_graph_level = graphLocus_2_levels.at(last_graph_locusID);

			std::cout << Utilities::timestamp() << "\tLocus" << locus << ", exon " << exonID << ": from " << first_graph_locusID << " (" << first_graph_level << ") to " << last_graph_locusID << " (" << last_graph_level << ").\n" << std::flush;

			assert(last_graph_level > first_graph_level);
			unsigned int expected_allele_length = last_graph_level - first_graph_level + 1;
			if(!(exon_level_names.size() == expected_allele_length))
			{
				std::cerr << "For locus " << locus << " exon " << exonID << " (" << exonFile << "), we have a problem with expected graph length.\n";
				std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
				std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
				std::cerr << std::flush;
			}
			assert(exon_level_names.size() == expected_allele_length);

			thisLocus_totalColumns += exon_level_names.size();
			
			combined_exon_sequences_locusIDs.insert(combined_exon_sequences_locusIDs.end(), exon_level_names.begin(), exon_level_names.end());
			for(unsigned int lI = 0; lI < expected_allele_length; lI++)
			{
				unsigned int graphLevel = first_graph_level + lI;
				assert(graphLocus_2_levels.at(exon_level_names.at(lI)) == graphLevel);
				combined_exon_sequences_graphLevels.push_back(graphLevel);
				combined_exon_sequences_graphLevels_individualExon.push_back(exonI);
				combined_exon_sequences_graphLevels_individualExonPosition.push_back(lI);
			}
			
			exon_lengths[exonI] = expected_allele_length;


			for(unsigned int lI = 1; lI < exon_lines.size(); lI++)
			{
				if(exon_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(exon_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					if(restrictToFullHaplotypes)
					{
						if(completeDefinedTypes.count(HLA_type) == 0)
						{
							continue;
						}
					}
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					std::string HLA_type_sequence = Utilities::join(line_alleles, "");

					if(can_translateToG_locus(locus))
					{
						// this is an initial test!
						bool _ignore;
						translate_allele_list_to_G_allele({HLA_type}, _ignore);
					}
					if(exonI == 0)
					{
						assert(combined_exon_sequences.count(HLA_type) == 0);
						combined_exon_sequences[HLA_type] = HLA_type_sequence;
					}
					else
					{
						assert(combined_exon_sequences.count(HLA_type));
						combined_exon_sequences.at(HLA_type) += HLA_type_sequence;
					}
				}
			}
			assert(combined_exon_sequences.size() > 0);
		}

		std::map<int, unsigned int> graphLevel_2_exonPosition;
		std::map<int, unsigned int> graphLevel_2_exonPosition_individualExon;
		std::map<int, unsigned int> graphLevel_2_exonPosition_individualExonPosition;
		for(unsigned int pI = 0; pI < combined_exon_sequences_graphLevels.size(); pI++)
		{
			int graphLevel = combined_exon_sequences_graphLevels.at(pI);
			assert(graphLevel >= 0);
			graphLevel_2_exonPosition[graphLevel] = pI;
			graphLevel_2_exonPosition_individualExon[graphLevel] = combined_exon_sequences_graphLevels_individualExon.at(pI);
			graphLevel_2_exonPosition_individualExonPosition[graphLevel] = combined_exon_sequences_graphLevels_individualExonPosition.at(pI);

		}

		std::cout << Utilities::timestamp() << "Have collected " << combined_exon_sequences.size() << " sequences -- first level " << combined_exon_sequences_graphLevels.front() << ", last level " << combined_exon_sequences_graphLevels.back() << ".\n" << std::flush;
		
		HLAtypeInference_totalColumns += thisLocus_totalColumns;
		
		std::map<std::string, unsigned int> HLAtype_2_clusterID;
		std::vector<std::set<std::string>> HLAtype_clusters;
		std::map<std::string, unsigned int> sequence_2_cluster;
		std::vector<std::string> cluster_2_sequence;

		for(std::map<std::string, std::string>::iterator HLAtypeIt = combined_exon_sequences.begin(); HLAtypeIt != combined_exon_sequences.end(); HLAtypeIt++)
		{
			std::string HLAtypeID = HLAtypeIt->first;
			std::string sequence = HLAtypeIt->second;
			if(sequence_2_cluster.count(sequence))
			{
				unsigned int cluster = sequence_2_cluster.at(sequence);
				HLAtype_clusters.at(cluster).insert(HLAtypeID);
				HLAtype_2_clusterID[HLAtypeID] = cluster;
			}
			else
			{
				std::set<std::string> newCluster;
				newCluster.insert(HLAtypeID);
				HLAtype_clusters.push_back(newCluster);
				unsigned int newClusterID = HLAtype_clusters.size() - 1;

				assert(sequence_2_cluster.count(sequence) == 0);
				sequence_2_cluster[sequence] = newClusterID;

				HLAtype_2_clusterID[HLAtypeID] = newClusterID;
			}
		}

		std::cout << Utilities::timestamp() << "Clustered into " << HLAtype_clusters.size() << " identical (over exons considered) clusters." << "\n" << std::flush;


		for(unsigned int clusterI = 0; clusterI < HLAtype_clusters.size(); clusterI++)
		{
			// std::cout << "\t\t\tcluster " << clusterI << "\n";
			for(std::set<std::string>::iterator typeIt = HLAtype_clusters.at(clusterI).begin(); typeIt != HLAtype_clusters.at(clusterI).end(); typeIt++)
			{
				std::string HLAtype = *typeIt;
				// std::cout << "\t\t\t\t" << HLAtype << "\n";
				assert(HLAtype_2_clusterID.at(HLAtype) == clusterI);

				if(typeIt == HLAtype_clusters.at(clusterI).begin())
				{
					std::string sequence = combined_exon_sequences.at(HLAtype);
					cluster_2_sequence.push_back(sequence);
				}
			}
		}

		// now transform reads into sequences specifying exon genotype value


		std::cout << Utilities::timestamp() << "Compute exon positions and specified genotypes from reads\n" << std::flush;

		std::vector< std::vector<oneExonPosition> > exonPositions_fromReads;



		unsigned int readPairs_OK = 0;
		unsigned int readPairs_broken = 0;

		for(unsigned int readPairI = 0; readPairI < alignments_paired.size(); readPairI++)
		{
			oneReadPair& originalReadPair = alignments_originalReads_paired.at(readPairI);
			std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& alignedReadPair = alignments_paired.at(readPairI);

			std::vector<oneExonPosition> read1_exonPositions;
			std::vector<oneExonPosition> read2_exonPositions;

			oneReadAlignment_2_exonPositions_paired(alignedReadPair.first, originalReadPair.reads.first, read1_exonPositions, alignedReadPair.second, originalReadPair.reads.second, 1, combined_exon_sequences_graphLevels, graphLevel_2_exonPosition);
			oneReadAlignment_2_exonPositions_paired(alignedReadPair.second, originalReadPair.reads.second, read2_exonPositions, alignedReadPair.first, originalReadPair.reads.first, 2, combined_exon_sequences_graphLevels, graphLevel_2_exonPosition);

			// if(originalReadPair.reads.first.name == "@@A819GMABXX:8:2204:2901:85228#GATCAGAT/1")
			// {
				// std::cerr << "!!!!!!!!!!!!!!!!!" << "\n";
				// std::cerr << "Locus: " << locus << "\n";
				// std::cerr << "Read name: " << originalReadPair.reads.first.name << "\n";
				// std::cerr << "countMismatchesInExon(read1_exonPositions): " << countMismatchesInExon(read1_exonPositions) << "\n" << std::flush;
				// std::cerr << "Alignment:\n";
				// std::cerr << "\t" << alignedReadPair.first.graph_aligned << "\n";
				// std::cerr << "\t" << alignedReadPair.first.sequence_aligned << "\n";
				// std::cerr << "!!!!!!!!!!!!!!!!!" << "\n" << std::flush;
			// }

			/*
			if(
					alignedReadPair_strandsValid(alignedReadPair) &&
					(abs(alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)) &&
					// (countMismatchesInExon(read1_exonPositions) < max_mismatches_perRead) &&
					// (countMismatchesInExon(read2_exonPositions) < max_mismatches_perRead)
					(alignmentFractionOK(alignedReadPair.first) >= min_alignmentFraction_OK) &&
					(alignmentFractionOK(alignedReadPair.second) >= min_alignmentFraction_OK) &&
					((alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.first) >= min_oneRead_weightedCharactersOK) || (alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.second) >= min_oneRead_weightedCharactersOK)) &&
					(alignedReadPair.first.mapQ >= minimumMappingQuality)
			)
			*/
			
			assert(alignedReadPair.first.mapQ_genomic != 2);
			double mapQ_thisAlignment = alignedReadPair.first.mapQ_genomic;
			if(
					alignedReadPair_strandsValid(alignedReadPair) &&
					(abs(alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)) &&
					(mapQ_thisAlignment >= minimumMappingQuality) &&
					((alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.first) >= min_bothReads_weightedCharactersOK) && (alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.second) >= min_bothReads_weightedCharactersOK))
					)  			
			{
				// good

				// std::cout << "\t\t" << "readPair " << readPairI << ", pairing OK.\n" << std::flush;

				std::vector<oneExonPosition> thisRead_exonPositions = read1_exonPositions;
				thisRead_exonPositions.insert(thisRead_exonPositions.end(), read2_exonPositions.begin(), read2_exonPositions.end());

				if(thisRead_exonPositions.size() > 0)
				{
					thisRead_exonPositions = removeDoublePositionsFromRead(thisRead_exonPositions);
					exonPositions_fromReads.push_back(thisRead_exonPositions);
				}

				readPairs_OK++;
			}
			else
			{
				// bad

				// std::cout << "\t\t" << "readPair " << readPairI << "/" < < alignments.size() << ", pairing FAILED.\n" << std::flush;

				if(alignedReadPair_strandsValid(alignedReadPair) &&
				(abs(alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)))
				{
					std::cout << "REJECTED MAPQ " << alignedReadPair.first.mapQ << " GENOMIC " << alignedReadPair.first.mapQ_genomic << "\n" << std::flush;
				}

				// (countMismatchesInExon(read1_exonPositions) < max_mismatches_perRead) &&
				// if((read1_exonPositions.size() > 0) && (countMismatchesInExon(read1_exonPositions) < max_mismatches_perRead))
				if(0 && (read1_exonPositions.size() > 0) && (alignmentFractionOK(alignedReadPair.first) >= min_alignmentFraction_OK) && (alignmentFractionOK(alignedReadPair.second) >= min_alignmentFraction_OK))
					exonPositions_fromReads.push_back(read1_exonPositions);

				// if((read2_exonPositions.size() > 0) && (countMismatchesInExon(read2_exonPositions) < max_mismatches_perRead))
				if(0 && (read2_exonPositions.size() > 0) && (alignmentFractionOK(alignedReadPair.first) >= min_alignmentFraction_OK) && (alignmentFractionOK(alignedReadPair.second) >= min_alignmentFraction_OK))
					exonPositions_fromReads.push_back(read2_exonPositions);

				readPairs_broken++;
			}
		}

		for(unsigned int readI = 0; readI < alignments_unpaired.size(); readI++)
		{
			oneRead& originalRead = alignments_originalReads_unpaired.at(readI);
			seedAndExtend_return_local& alignedRead = alignments_unpaired.at(readI);

			std::vector<oneExonPosition> read_exonPositions;
			oneReadAlignment_2_exonPositions_unpaired(alignedRead, originalRead, read_exonPositions, combined_exon_sequences_graphLevels, graphLevel_2_exonPosition);


			double mapQ_thisAlignment = alignedRead.mapQ_genomic;
			if(mapQ_thisAlignment >= minimumMappingQuality)
			{
				// good
				
				std::vector<oneExonPosition> thisRead_exonPositions = read_exonPositions;
				
				if(thisRead_exonPositions.size() > 0)
					exonPositions_fromReads.push_back(read_exonPositions);

				readPairs_OK++;
			}
			else
			{
				// bad

				std::cout << "REJECTED MAPQ " << alignedRead.mapQ << " GENOMIC " << alignedRead.mapQ_genomic << "\n" << std::flush;

				readPairs_broken++;
			}
		}

		std::cout << Utilities::timestamp() << "Mapped reads to exons. " << readPairs_OK << " pairs OK, " << readPairs_broken << " pairs broken -- " << exonPositions_fromReads.size() << " reads relevant to this locus.\n" << std::flush;

		// Pileup of mapped reads
		std::set<size_t> ignore_exonPositions_fromReads;

		std::map<unsigned int, std::map<std::string, int>> perPosition_allele_counts_postFiltering;
		std::map<unsigned int, std::set<std::string>> perPosition_ignore_alleles;
		{
			unsigned int n_relevant_positions = graphLevel_2_exonPosition_individualExon.size();
			size_t coverage_over_relevant_positions = 0;
			std::set<size_t> read_with_good_positions;

			std::map<unsigned int, std::map<std::string, int>> perPosition_allele_counts;
			for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
			{
				std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
				for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
				{
					oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);
					assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
					if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
					{
						continue;
					}

					int position = onePositionSpecifier.positionInExon;
					std::string allele = onePositionSpecifier.genotype;

					if(perPosition_allele_counts[position].count(allele) == 0)
					{
						perPosition_allele_counts[position][allele] = 0;
					}

					perPosition_allele_counts.at(position).at(allele)++;
				}
			}

			int positions_sufficientCoverage = 0;
			int positions_with_kickedOutAlleles = 0;
			int kickedOutAlleles = 0;
			if(highCoverage_filter_alleles)
				std::cout << "Locus " << locus << ", highCoverage_filter_alleles active, evaluate " << perPosition_allele_counts.size() << " positions.\n" << std::flush;
			for(auto position : perPosition_allele_counts)
			{
				int count_position = 0;
				for(auto allele : perPosition_allele_counts.at(position.first))
				{
					count_position += allele.second;
				}
				if(count_position >= highCoverage_minCoverage)
				{
					positions_sufficientCoverage++;
					bool kickedOutAtLeastOneAllele = false;
					for(auto allele : perPosition_allele_counts.at(position.first))
					{
						double aF = (double)allele.second / (double)count_position;
						if((aF < highCoverage_minAlleleFreq) && highCoverage_filter_alleles)
						{
							perPosition_ignore_alleles[position.first].insert(allele.first);
							kickedOutAlleles += allele.second;
							kickedOutAtLeastOneAllele = true;
						}
						else
						{
							perPosition_allele_counts_postFiltering[position.first][allele.first] = allele.second;
						}
					}
					if(kickedOutAtLeastOneAllele)
					{
						positions_with_kickedOutAlleles++;
					}
				}
			}

			if(highCoverage_filter_alleles)
			{
				std::cout << "\tPositions with coverage high enough for filtering: " << positions_sufficientCoverage << "\n";
				std::cout << "\tPositions with at least one allele removed with aF <= " << highCoverage_minAlleleFreq << ": " << positions_with_kickedOutAlleles << "\n";
				std::cout << "\tTotal number of (individual-read) alleles removed: " << kickedOutAlleles << "\n";
				std::cout << "\n" << std::flush;
			}
		}

		std::map<int, std::map<int, std::vector<oneExonPosition> > > pileUpPerPosition;
		for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
		{
			if(ignore_exonPositions_fromReads.count(readI))
			{
				continue;
			}

			std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
			for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
			{
				oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);
				
				assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
				if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
				{
					continue;
				}

				if(perPosition_ignore_alleles.count(onePositionSpecifier.positionInExon) && (perPosition_ignore_alleles.at(onePositionSpecifier.positionInExon).count(onePositionSpecifier.genotype)))
				{
					continue;
				}
					
				int individualExon = graphLevel_2_exonPosition_individualExon.at(onePositionSpecifier.graphLevel);
				int individualExonPosition = graphLevel_2_exonPosition_individualExonPosition.at(onePositionSpecifier.graphLevel);

				pileUpPerPosition[individualExon][individualExonPosition].push_back(onePositionSpecifier);
			}
		}

		std::string fileName_pileUp = outputDirectory + "/R1_pileup_"+locus+".txt";
		std::ofstream pileUpStream;
		pileUpStream.open(fileName_pileUp.c_str());
		assert(pileUpStream.is_open());
		

		for(std::map<int, std::map<int, std::vector<oneExonPosition> > >::iterator exonIt = pileUpPerPosition.begin(); exonIt != pileUpPerPosition.end(); exonIt++)
		{
			int exon = exonIt->first;
			int exonL = exon_lengths.at(exon);
			for(int exonPos = 0; exonPos < exonL; exonPos++)
			{
				if(pileUpPerPosition.at(exon).count(exonPos))
				{
					// int exonPos = exonPosIt->first;
					std::vector<oneExonPosition> pileUp = pileUpPerPosition.at(exon).at(exonPos);

					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(exon));
					fieldsPerLine.push_back(Utilities::ItoStr(exonPos));
					fieldsPerLine.push_back(Utilities::ItoStr(pileUp.size()));

					std::vector<std::string> piledUpGenotypes;

					for(unsigned int pI = 0; pI < pileUp.size(); pI++)
					{
						oneExonPosition piledPosition = pileUp.at(pI);

						std::vector<std::string> qualities_as_strings;
						for(unsigned int qI = 0; qI < piledPosition.qualities.size(); qI++)
						{
							char qC = piledPosition.qualities.at(qI);
							int qC_i = qC;
							qualities_as_strings.push_back(Utilities::ItoStr(qC_i));
						}

						std::string pileUpString = piledPosition.genotype
							+ " (" + Utilities::join(qualities_as_strings, ", ") + ")"
							+ " ["
							// + Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "
							// + Utilities::DtoStr(piledPosition.thisRead_fractionOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_fractionOK) + " | "
							// + Utilities::ItoStr(piledPosition.pairs_strands_OK) + " "
							+ Utilities::DtoStr(piledPosition.pairs_strands_distance) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ_position) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ) + " "
							+ Utilities::DtoStr(piledPosition.mapQ_genomic) + " | "
							+ Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							+ Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "						
							+ piledPosition.thisRead_ID + " "
							+ piledPosition.pairedRead_ID
							+ "]";

						utilized_reads.insert(piledPosition.thisRead_ID);
							
						piledUpGenotypes.push_back(pileUpString);
					}

					fieldsPerLine.push_back(Utilities::join(piledUpGenotypes, ", "));


					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\n";
				}
				else
				{
					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(exon));
					fieldsPerLine.push_back(Utilities::ItoStr(exonPos));
					fieldsPerLine.push_back(Utilities::ItoStr(0));
					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\n";				
				}				
			}
		}
		pileUpStream.close();

		std::string fileName_readIDs = outputDirectory + "/R1_readIDs_"+locus+".txt";
		std::ofstream readIDsStream;
		readIDsStream.open(fileName_readIDs.c_str());
		assert(readIDsStream.is_open());
		for(std::set<std::string>::iterator readIDit = utilized_reads.begin(); readIDit != utilized_reads.end(); readIDit++)
		{
			readIDsStream << *readIDit << "\n";
		}
		readIDsStream.close();
		
		// likelihoods for reads

		std::cout << Utilities::timestamp() << "Compute likelihoods for all exon-overlapping reads (" << exonPositions_fromReads.size() << "), conditional on underlying exons.\n" << std::flush;
		std::vector<std::vector<double> > likelihoods_perCluster_perRead;
		std::vector<std::vector<double> > likelihoods_perCluster_perObservedBase;
		std::vector<std::vector<int> > mismatches_perCluster_perRead;

		likelihoods_perCluster_perRead.resize(HLAtype_clusters.size());
		mismatches_perCluster_perRead.resize(HLAtype_clusters.size());
		likelihoods_perCluster_perObservedBase.resize(HLAtype_clusters.size());


		// std::cout << "\n\n" << HLAtype_2_clusterID.at("A*30:73N") << " " << HLAtype_2_clusterID.at("A*29:02:08") << "\n" << std::flush;
		// assert( 1 == 0 );

		std::set<int> printClusters;
		// printClusters.insert(1789);
		// printClusters.insert(1646);


		// can be removed later, will lead to problems during final normalization
//		assert(exonPositions_fromReads.size() > 0);
		for(unsigned int clusterI = 0; clusterI < HLAtype_clusters.size(); clusterI++)
		{
			std::vector<std::string> typesInCluster(HLAtype_clusters.at(clusterI).begin(), HLAtype_clusters.at(clusterI).end());
			std::string clusterName = Utilities::join(typesInCluster, "|");

			bool verbose = printClusters.count(clusterI);
			// verbose = (clusterI == 0);
			
			if(verbose)
			{
				std::cout << "CLUSTER " << clusterI << " " << clusterName << "\n";
			}

			std::string& clusterSequence = cluster_2_sequence.at(clusterI);

			likelihoods_perCluster_perRead.at(clusterI).resize(exonPositions_fromReads.size(), log(10));
			mismatches_perCluster_perRead.at(clusterI).resize(exonPositions_fromReads.size(), 0);

			// this is really readI
			for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
			{			
				if(ignore_exonPositions_fromReads.count(readI))
				{
					continue;
				}

				std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
				double log_likelihood_read = 0;
				int mismatches = 0;

				std::string readID;
				std::string read1_ID;
				if(individualPositions.size() > 0)
				{
					readID = individualPositions.at(0).thisRead_ID;
					read1_ID = individualPositions.at(0).read1_ID;

//					double mapQ_thisAlignment = (individualPositions.at(0).mapQ_genomic != 2) ? individualPositions.at(0).mapQ_genomic : individualPositions.at(0).mapQ;
//					assert((mapQ_thisAlignment >= 0) && (mapQ_thisAlignment <= 1));

					// log_likelihood_read += log(mapQ_thisAlignment);
					assert(! combineReadAndBaseLikelihoods); // this would not make sense, or at least one would have to think more about it
	
				}

				// bool verbose = ((clusterI == 1127) && (readID == "@@B81EP5ABXX:8:2208:11879:23374#GATCAGAT/2") && 0);

				if(verbose)
					std::cout << "Likelihood calculation for read " << readID << " / cluster " << clusterI << "\n";


				for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
				{
					oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);

					assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
					if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
					{
						if(verbose)
						{
							std::cout << "\t" << positionI << " mapQ_position " << onePositionSpecifier.mapQ_position << " too low (below " << minimumPerPositionMappingQuality << ").\n" << std::flush;
						}					
						continue;
					}
					
					if(perPosition_ignore_alleles.count(onePositionSpecifier.positionInExon) && (perPosition_ignore_alleles.at(onePositionSpecifier.positionInExon).count(onePositionSpecifier.genotype)))
					{
						continue;
					}

					if(clusterI == 0)
					{
						HLATypeInference_thisLocus_bases_used++;
						HLATypeInference_totalBases_used++;
					}
				
					double log_likelihood_position = 0;

					std::string exonGenotype = clusterSequence.substr(onePositionSpecifier.positionInExon, 1);
					std::string readGenotype = onePositionSpecifier.genotype;
					std::string readQualities = onePositionSpecifier.qualities;

					if(verbose)
					{
						std::cout << "\t" << positionI << " exon pos " << onePositionSpecifier.positionInExon << ": " << exonGenotype << " " << readGenotype << "\n" << std::flush;
					}

					assert(exonGenotype.length() == 1);
					assert(readGenotype.length() >= 1);
					unsigned int l_diff = readGenotype.length() - exonGenotype.length();
					if(exonGenotype == "_")
					{
						// assert(l_diff == 0);
						if(readGenotype == "_")
						{
							assert(onePositionSpecifier.graphLevel != -1);
							// likelihood 1 - intrinsic graph gap

							if(verbose)
							{
								std::cout << "\t\t" << "Intrinsic graph gap" << "\n";
							}

						}
						else
						{
							if(verbose)
							{
								std::cout << "\t\t" << "Insertion " << (1 + l_diff) << "\n";
							}

							assert(readGenotype.find("_") == std::string::npos);
							log_likelihood_position += (log_likelihood_insertion_actualAllele * (1 + l_diff));
						}
					}
					else
					{

						if(readGenotype.length() > 1)
						{
							std::string readGenotype_after1 = readGenotype.substr(1);
							assert(readGenotype_after1.find("_") == std::string::npos);
						}
						// score from first position match
						if(readGenotype.substr(0, 1) == "_")
						{
							log_likelihood_position += log_likelihood_deletion;

							if(verbose)
							{
								std::cout << "\t\t" << "Deletion" << "\n";
							}
						}
						else
						{
							log_likelihood_position += log_likelihood_match_mismatch;
							
							assert(readQualities.length());
							double pCorrect = Utilities::PhredToPCorrect(readQualities.at(0));
							if(veryConservativeReadLikelihoods)
							{
								if(pCorrect > 0.999)
									pCorrect = 0.999;
							}
							assert((pCorrect >= 0) && (pCorrect <= 1));

							if(pCorrect == 0)
							{
								pCorrect = 0.001;
							}
							
							// if(!(pCorrect >= 0.25))
							// {
								// std::cerr << "pCorrect = " << pCorrect << "\n" << std::flush;
							// }
							// assert(pCorrect >= 0.25);
					
							if(exonGenotype == readGenotype.substr(0, 1))
							{
								if(verbose)
								{
									std::cout << "\t\t" << "Match " << pCorrect << "\n";
								}

								log_likelihood_position += log(pCorrect);
							}
							else
							{
								double pIncorrect = (1 - pCorrect)*(1.0/3.0);
								assert(pIncorrect <= 0.75);
								assert((pIncorrect > 0) && (pIncorrect < 1));
								log_likelihood_position += log(pIncorrect);

								if(verbose)
								{
									std::cout << "\t\t" << "Mismatch " << pIncorrect << "\n";
								}
							}
						}
						// if read allele is longer
						log_likelihood_position += (log_likelihood_insertion_actualAllele * l_diff);

						if(l_diff > 0)
						{
							if(verbose)
							{
								std::cout << "\t\t" << "Insertion " << l_diff << "\n";
							}
						}
					}

					if(readGenotype != "_")
					{
						if(readGenotype != exonGenotype)
						{
							mismatches++;
						}
					}


					log_likelihood_read += log_likelihood_position;

					if(verbose)
					{
						std::cout << "\t\t" << "Running log likelihood: " << log_likelihood_read << "\n";
					}


					likelihoods_perCluster_perObservedBase.at(clusterI).push_back(log_likelihood_position);
				}

				if(individualPositions.size() > 0)
				{
					if((clusterI == 1160) || (clusterI == 1127))
					{
						// std::cout << locus << " " << readID << " cluster " << clusterI << ": " << log_likelihood << "\n" << std::flush;
					}
				}

				// std::cout << "cluster " << clusterI << ", position sequence " << positionSpecifierI << ": " << log_likelihood << "\n";

				assert(exp(log_likelihood_read) >= 0);
				assert(exp(log_likelihood_read) <= 1);
				
				likelihoods_perCluster_perRead.at(clusterI).at(readI) = log_likelihood_read;
				mismatches_perCluster_perRead.at(clusterI).at(readI) = mismatches;
			}

			assert(likelihoods_perCluster_perRead.at(clusterI).size() == exonPositions_fromReads.size());
			assert(mismatches_perCluster_perRead.at(clusterI).size() == exonPositions_fromReads.size());
		}

		// std::cout << Utilities::timestamp() << "Compute normalized likelihoods (over multiple alignments, if there are any)." << std::flush;
		//  std::cout << "\t" << "exonPositions_fromReads.size(): " << exonPositions_fromReads.size() << "\n" << std::flush;
		
		/*
		std::map<std::string, std::set<unsigned int> > readID_noAlignment_2_readIndex;
		std::map<std::string, unsigned int > readID_noAlignment_2_idx;
		for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
		{
			std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
			if(individualPositions.size() > 0)
			{
				std::string read1_ID = individualPositions.at(0).read1_ID;
				std::vector<std::string> readID_components = Utilities::split(read1_ID, ":");
				assert(readID_components.back().substr(0, 1) == "A");
				std::vector<std::string> readID_components_noA(readID_components.begin(), readID_components.end() - 1);
				assert(readID_components_noA.size() == (readID_components.size() - 1));

				std::string readID_noA = Utilities::join(readID_components_noA, ":");
				readID_noAlignment_2_readIndex[readID_noA].insert(readI);
			
				if(readID_noAlignment_2_idx.count(readID_noA) == 0)
				{
					unsigned int newIndex = readID_noAlignment_2_idx.size();
					readID_noAlignment_2_idx[readID_noA] = newIndex;
				}
			}
		}
		
		*/
		
		// std::cout << Utilities::timestamp() << "\tProcessed " << readID_noAlignment_2_idx.size() << " reads (no alignment)\n" << std::flush;



		std::vector<std::vector<double> > likelihoods_perCluster_perReads_allAlignments;
		
		/*


		unsigned maxNoAlignmentReadIdx = readID_noAlignment_2_idx.size() - 1;

		int reads_summed_more_than_1_alignment = 0;

		likelihoods_perCluster_perReads_allAlignments.resize(HLAtype_clusters.size());

		for(unsigned int clusterI = 0; clusterI < HLAtype_clusters.size(); clusterI++)
		{
			likelihoods_perCluster_perReads_allAlignments.at(clusterI).resize(maxNoAlignmentReadIdx + 1, 0);

			std::vector<std::string> typesInCluster(HLAtype_clusters.at(clusterI).begin(), HLAtype_clusters.at(clusterI).end());
			std::string clusterName = Utilities::join(typesInCluster, "|");

			for(std::map<std::string, std::set<unsigned int> >::iterator readIDit = readID_noAlignment_2_readIndex.begin(); readIDit != readID_noAlignment_2_readIndex.end(); readIDit++)
			{
				const std::set<unsigned int>& alignments_for_read = readIDit->second;
				double summed_mapQ = 0;
				std::vector<double> ll_across_alignments;
				ll_across_alignments.reserve(ll_across_alignments.size() + 1);
				// std::cout << "(ll_across_alignments.size() + 1)" << ": " << (ll_across_alignments.size() + 1) << "}
				for(std::set<unsigned int>::iterator readIndexit = alignments_for_read.begin(); readIndexit != alignments_for_read.end(); readIndexit++)
				{
					std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(*readIndexit);
					assert(individualPositions.at(0).mapQ_genomic != 2);

					double mapQ_genomic = individualPositions.at(0).mapQ_genomic;
					assert((mapQ_genomic > 0) && (mapQ_genomic <= 1));

					summed_mapQ += mapQ_genomic; 

					ll_across_alignments.push_back(log(mapQ_genomic) + likelihoods_perCluster_perRead.at(clusterI).at(*readIndexit));
				}

				if(!((summed_mapQ >= 0) && (summed_mapQ <= (1+1e-5))))
				{
					std::cerr << "! ((summed_mapQ >= 0) && (summed_mapQ <= (1+1e-5)))" << "\n";
					std::cerr << "summed_mapQ: " << summed_mapQ << "\n" << std::flush;
				}
				assert((summed_mapQ > 0) && (summed_mapQ <= (1+1e-5)));

				double component_outside = (1 - summed_mapQ) * 1; // likelihood is 1 for all non-captured regions

				if(component_outside > 0)
				{
					assert(component_outside > 0);
					assert(component_outside <= 1);
					
					ll_across_alignments.push_back(log(component_outside));
				}

				unsigned int idx_for_readID = readID_noAlignment_2_idx.at(readIDit->first);
				assert(likelihoods_perCluster_perReads_allAlignments.at(clusterI).at(idx_for_readID) == 0);
				likelihoods_perCluster_perReads_allAlignments.at(clusterI).at(idx_for_readID) = Utilities::LogSumLogPs(ll_across_alignments);

				double nonLog_P = exp(likelihoods_perCluster_perReads_allAlignments.at(clusterI).at(idx_for_readID));
				
				if(!((nonLog_P >= 0) && (nonLog_P <= (1+1e-5))))
				{
					std::cerr << "!((nonLog_P >= 0) && (nonLog_P <= 1))" << "\n";
					std::cerr << "nonLog_P" << ": " << nonLog_P << "\n";
					std::cerr << "likelihoods_perCluster_perReads_allAlignments.at(clusterI).at(idx_for_readID)" << ": " << likelihoods_perCluster_perReads_allAlignments.at(clusterI).at(idx_for_readID) << "\n";
					for(unsigned int i = 0; i < ll_across_alignments.size(); i++)
					{	
						std::cerr << "\t" << "ll_across_alignments" << " " << i << ": " << ll_across_alignments.at(i) << "\n";
					}		
					std::cerr << std::flush;
				}
				assert(nonLog_P >= 0);
				assert(nonLog_P <= (1+1e-5));
			}

			if(readID_noAlignment_2_readIndex.size() > 1)
			{
				reads_summed_more_than_1_alignment++;
			}
		}
		*/
		

		// std::cout << Utilities::timestamp() << "\tdone. " << reads_summed_more_than_1_alignment << " of " << exonPositions_fromReads.size() << " original reads summed over more than 1 alignment.\n" << std::flush;


		std::cout << Utilities::timestamp() << "Compute likelihoods for all exon cluster pairs (" << HLAtype_clusters.size() << "**2/2)\n" << std::flush;

		std::vector<double> LLs_completeReads;

		std::vector<double> LLs_observedBases;

		std::vector<double> Mismatches_avg;
		std::vector<double> Mismatches_min;

		std::vector<std::pair<unsigned int, unsigned int> > LLs_clusterIs;

		std::cout << "Threads: " << omp_get_num_threads() << "\n";
		
		
		size_t HLAtype_cluster_SIZE = HLAtype_clusters.size();
		
		#pragma omp parallel for schedule(dynamic)
		for(unsigned int clusterI1 = 0; clusterI1 < HLAtype_cluster_SIZE; clusterI1++)
		{	
			if(omp_get_thread_num() == 0)
			{
				std::cout << "\rClusterI1 = " << clusterI1 << std::flush;
			}
			
			std::vector<double> LLs_completeReads_perThread;

			std::vector<double> LLs_observedBases_perThread;

			std::vector<double> Mismatches_avg_perThread;
			std::vector<double> Mismatches_min_perThread;

			std::vector<std::pair<unsigned int, unsigned int> > LLs_clusterIs_perThread;
			// std::vector<unsigned int> LLs_completeReads_indices_perThread;
		
			for(unsigned int clusterI2 = clusterI1; clusterI2 < HLAtype_clusters.size(); clusterI2++)
			{

				double mismatches_sum_averages = 0;
				double mismatches_sum_min = 0;

				double pair_log_likelihood = 0;

				for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
				{
					if(ignore_exonPositions_fromReads.count(readI))
					{
						continue;
					}

					// now use allAlignments likelihoods instead!
//					std::string readID;
//					if(exonPositions_fromReads.at(readI).size() > 0)
//					{
//						readID = exonPositions_fromReads.at(readI).at(0).thisRead_ID;
//					}


					double LL_thisRead_cluster1 = likelihoods_perCluster_perRead.at(clusterI1).at(readI);
					double LL_thisRead_cluster2 = likelihoods_perCluster_perRead.at(clusterI2).at(readI);

					assert(exp(LL_thisRead_cluster1 <= 1));
					assert(exp(LL_thisRead_cluster2 <= 1));

					// double LL_average = log( (1 + (exp(LL_thisPositionSpecifier_cluster2+log(0.5)))/exp(LL_thisPositionSpecifier_cluster1+log(0.5)) )) + (LL_thisPositionSpecifier_cluster1+log(0.5));

					// double LL_average_2 = log(0.5 * exp(LL_thisRead_cluster1) + 0.5 * exp(LL_thisRead_cluster2));
					double LL_average_2_2 = Utilities::logAvg(LL_thisRead_cluster1, LL_thisRead_cluster2);
					
					// if((exp(LL_thisRead_cluster1) == 0) && (exp(LL_thisRead_cluster2) == 0))
					// {
						// LL_average_2_2 = -1 * std::numeric_limits<double>::infinity();
					// }
					
					// std::cerr << LL_average_2_2 << "\n";
					
					if(!((exp(LL_average_2_2) >= 0) && (exp(LL_average_2_2) <= 1)))
					{  
						#pragma omp critical
						{
							std::cerr << "Problem with readI " << readI << " conditional on clusters " << clusterI1 << "/" << clusterI2 << "\n";
							std::cerr << "LL_average_2_2" << ": " << LL_average_2_2 << " exp: " << exp(LL_average_2_2) << "\n";
							std::cerr << "LL_thisRead_cluster1" << ": " << LL_thisRead_cluster1 << " exp: " << exp(LL_thisRead_cluster1) << "\n";
							std::cerr << "LL_thisRead_cluster2" << ": " << LL_thisRead_cluster2 << " exp: " << exp(LL_thisRead_cluster2) << "\n";
							std::cerr << std::flush;	
						}						
					}
					
					assert((exp(LL_average_2_2) >= 0) && (exp(LL_average_2_2) <= 1));
					
					// if(! (abs(LL_average_2 - LL_average_2_2) < 1e-5))
					// {
						// std::cerr << "Cluster 1: " << LL_thisRead_cluster1 << " " << exp(LL_thisRead_cluster1) << "\n";
						// std::cerr << "Cluster 2: " << LL_thisRead_cluster2 << " " << exp(LL_thisRead_cluster2) << "\n";
						// std::cerr << LL_average_2 << "\n";
						// std::cerr << LL_average_2_2 << "\n" << std::flush;
					// }
					// assert(abs(LL_average_2 - LL_average_2_2) < 1e-5);

					int mismatches_cluster1 = mismatches_perCluster_perRead.at(clusterI1).at(readI);
					int mismatches_cluster2 = mismatches_perCluster_perRead.at(clusterI2).at(readI);

					mismatches_sum_averages += ((double)(mismatches_cluster1 + mismatches_cluster2) / 2.0);
					mismatches_sum_min += ((mismatches_cluster1 < mismatches_cluster2) ? mismatches_cluster1 : mismatches_cluster2);

					 pair_log_likelihood += LL_average_2_2;

					// if ( ((clusterI1 == 1160) && (clusterI2 == 1640)) || ((clusterI1 == 1127) && (clusterI2 == 1640)))
					// {
						// std::cout << "Cluster pair " << clusterI1 << " / " << clusterI2 << " read " << readID << ": " << LL_thisPositionSpecifier_cluster1 << " and " << LL_thisPositionSpecifier_cluster2 << ": " << LL_average_2 << "\n" << std::flush;
					// }
				}

				/*
				assert(! combineReadAndBaseLikelihoods);
				for(std::map<std::string, std::set<unsigned int> >::const_iterator readIDit = readID_noAlignment_2_readIndex.begin(); readIDit != readID_noAlignment_2_readIndex.end(); readIDit++)
				{
					const std::string readID = readIDit->first;
					unsigned int idx_for_readID = readID_noAlignment_2_idx.at(readIDit->first);
					double LL_thisRead_cluster1 = likelihoods_perCluster_perReads_allAlignments.at(clusterI1).at(idx_for_readID);
					double LL_thisRead_cluster2 = likelihoods_perCluster_perReads_allAlignments.at(clusterI2).at(idx_for_readID);
					double LL_average_2_2 = logAvg(LL_thisRead_cluster1, LL_thisRead_cluster2);

					pair_log_likelihood += LL_average_2_2;
				}
				
				*/
				

				LLs_completeReads_perThread.push_back(pair_log_likelihood);

				Mismatches_avg_perThread.push_back(mismatches_sum_averages);
				Mismatches_min_perThread.push_back(mismatches_sum_min);

				LLs_clusterIs_perThread.push_back(make_pair(clusterI1, clusterI2));

				// LLs_completeReads_indices_perThread.push_back(LLs_completeReads.size() - 1);

				double pair_log_likelihood_fromBases = 0;
				if(combineReadAndBaseLikelihoods)
				{
					assert(likelihoods_perCluster_perObservedBase.at(clusterI1).size() == likelihoods_perCluster_perObservedBase.at(clusterI2).size());
					for(unsigned int baseI = 0; baseI < likelihoods_perCluster_perObservedBase.at(clusterI1).size(); baseI++)
					{
						double LL_thisBase_cluster1 = likelihoods_perCluster_perObservedBase.at(clusterI1).at(baseI);
						double LL_thisBase_cluster2 = likelihoods_perCluster_perObservedBase.at(clusterI2).at(baseI);

						double LL_thisBase_avg = Utilities::logAvg(LL_thisBase_cluster1, LL_thisBase_cluster2);
						pair_log_likelihood_fromBases += LL_thisBase_avg;
					}
				}
				LLs_observedBases_perThread.push_back(pair_log_likelihood_fromBases);
			}
			
			#pragma omp critical
			{
				LLs_completeReads.insert(LLs_completeReads.end(), LLs_completeReads_perThread.begin(), LLs_completeReads_perThread.end());

				LLs_observedBases.insert(LLs_observedBases.end(), LLs_observedBases_perThread.begin(), LLs_observedBases_perThread.end());

				Mismatches_avg.insert(Mismatches_avg.end(), Mismatches_avg_perThread.begin(), Mismatches_avg_perThread.end());
				Mismatches_min.insert(Mismatches_min.end(), Mismatches_min_perThread.begin(), Mismatches_min_perThread.end());

				LLs_clusterIs.insert(LLs_clusterIs.end(), LLs_clusterIs_perThread.begin(), LLs_clusterIs_perThread.end());
				// LLs_completeReads_indices.insert(LLs_completeReads_indices.end(), LLs_completeReads_indices_perThread.begin(), LLs_completeReads_indices_perThread.end());
			}
		}

		std::vector<size_t> LLs_completeReads_indices;		
		for(size_t i = 0; i < LLs_completeReads.size(); i++)
		{
			LLs_completeReads_indices.push_back(i);
		}
		
		assert(LLs_observedBases.size() == LLs_completeReads.size());
		assert(Mismatches_avg.size() == LLs_completeReads.size());
		assert(LLs_completeReads.size() == LLs_completeReads.size());
			
		bool verbose = false;

		if(combineReadAndBaseLikelihoods)
		{
			std::cout << "\n\n" << Utilities::timestamp() << "\tCombine likelihoods...\n\n" << std::flush;
		
			std::vector<double> LLs_combined;

			convertLogVectorToP(LLs_completeReads);
			convertLogVectorToP(LLs_observedBases);

			normalizeVector(LLs_completeReads);
			normalizeVector(LLs_observedBases);

			for(unsigned int i = 0; i < LLs_completeReads.size(); i++)
			{
				assert(LLs_completeReads.at(i) >= 0);
				assert(LLs_completeReads.at(i) <= 1);
				assert(LLs_observedBases.at(i) >= 0);
				assert(LLs_observedBases.at(i) <= 1);

				double l_combined = 0.5 * LLs_completeReads.at(i) + 0.5 * LLs_observedBases.at(i);

				if(verbose)
					std::cout << "Alternative " << i << ": " <<  LLs_completeReads.at(i) << " x " << LLs_observedBases.at(i) << "\n" << std::flush;


				assert(l_combined >= 0);
				assert(l_combined <= 1);

				LLs_combined.push_back(log(l_combined));
			}

			assert(LLs_combined.size() == LLs_completeReads.size());
			LLs_completeReads = LLs_combined;
		}

		std::cout << "\n\n" << Utilities::timestamp() << "Sorting...\n\n" << std::flush;

		assert(LLs_completeReads.size() == Mismatches_avg.size());
		assert(LLs_completeReads.size() == LLs_completeReads_indices.size());

		std::sort(LLs_completeReads_indices.begin(), LLs_completeReads_indices.end(), [&](unsigned int a, unsigned int b) {
			if(!(a < LLs_completeReads.size()))
			{
				std::cerr << "a" << ": " << a << "\n";
				std::cerr << "LLs_completeReads.size()" << ": " << LLs_completeReads.size() << "\n" << std::flush;
			}
			assert(a < LLs_completeReads.size());
			assert(b < LLs_completeReads.size());
			assert(a < Mismatches_avg.size());
			assert(b < Mismatches_avg.size());
			
			
			if(LLs_completeReads.at(a) == LLs_completeReads.at(b))
			{
				return (Mismatches_avg.at(b) < Mismatches_avg.at(a));
			}
			else
			{
				return (LLs_completeReads.at(a) < LLs_completeReads.at(b));
			}
		});

		std::reverse(LLs_completeReads_indices.begin(), LLs_completeReads_indices.end());

		std::cout << "\n\n" << Utilities::timestamp() << "Sorting done!...\n\n" << std::flush;


		std::pair<double, unsigned int> maxPairI = Utilities::findVectorMax(LLs_completeReads);
		std::cout << Utilities::timestamp() << "Done. (One) maximum pair is " << maxPairI.second << " with LL = " << maxPairI.first << "\n" << std::flush;

		std::vector<double> LLs_normalized;
		double LL_max = maxPairI.first;
		double P_sum = 0;
		for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
		{
			double LL = LLs_completeReads.at(cI);
			double P = exp(LL - LL_max);
			P_sum += P;
		}
		if(P_sum > 0)
		{
			for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
			{
				double LL = LLs_completeReads.at(cI);
				double P = exp(LL - LL_max);
				double P_normalized = P / P_sum;
				if(!((P_normalized >= 0) && (P_normalized <= 1)))
				{
					std::cerr << "P_normalized: " << P_normalized << "\n";
					std::cerr << "P: " << P << "\n";
					std::cerr << "LL: " << LL << "\n";
					std::cerr << "LL_max: " << LL_max << "\n";
					
					std::cerr << std::flush;
				}
				assert(P_normalized >= 0);
				assert(P_normalized <= 1);
				LLs_normalized.push_back(P_normalized);
			}
		}
		else
		{
			for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
			{
				LLs_normalized.push_back(1.0/(double)LLs_clusterIs.size());
			}
		
		}

		std::ofstream allPairsStream;
		allPairsStream.open(outputFN_allPairs.c_str());
		assert(allPairsStream.is_open());
		allPairsStream << "ClusterID" << "\t" << "P" << "\t" << "LL" << "\t" << "Mismatches_avg" << "\n";

		std::vector<std::string> LLs_identifiers;
		std::map<int, double> clusterI_overAllPairs;
		for(unsigned int cII = 0; cII < LLs_completeReads_indices.size(); cII++)
		{
			unsigned int cI = LLs_completeReads_indices.at(cII);
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(cI);

			std::vector<std::string> cluster1_members(HLAtype_clusters.at(clusters.first).begin(), HLAtype_clusters.at(clusters.first).end());
			std::vector<std::string> cluster2_members(HLAtype_clusters.at(clusters.second).begin(), HLAtype_clusters.at(clusters.second).end());

			std::string id = Utilities::join(cluster1_members, ";") + "/" + Utilities::join(cluster2_members, ";");

			LLs_identifiers.push_back(id);

			allPairsStream << id << "\t" << LLs_normalized.at(cI) << "\t" << LLs_completeReads.at(cI) << "\t" << Mismatches_avg.at(cI) << "\n";

			if(clusterI_overAllPairs.count(clusters.first) == 0)
			{
				clusterI_overAllPairs[clusters.first] = 0;
			}
			clusterI_overAllPairs[clusters.first] += LLs_normalized.at(cI);

			if(clusters.second != clusters.first)
			{
				if(clusterI_overAllPairs.count(clusters.second) == 0)
				{
					clusterI_overAllPairs[clusters.second] = 0;
				}
				clusterI_overAllPairs[clusters.second] += LLs_normalized.at(cI);
			}
		}

		allPairsStream.close();

		std::pair<double, int> bestGuess_firstAllele = Utilities::findIntMapMax(clusterI_overAllPairs);
		std::string bestGuess_firstAllele_ID = Utilities::join(std::vector<std::string>(HLAtype_clusters.at(bestGuess_firstAllele.second).begin(), HLAtype_clusters.at(bestGuess_firstAllele.second).end()), ";");
		std::string bestGuess_firstAllele_oneType = *(HLAtype_clusters.at(bestGuess_firstAllele.second).begin());


		assert(bestGuess_firstAllele.first >= 0);
		assert(bestGuess_firstAllele.second >= 0);

		std::map<int, double> bestGuess_secondAllele_alternatives;
		std::map<int, double> bestGuess_secondAllele_alternatives_mismatches;

		for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
		{
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(cI);

			if((int)clusters.first == bestGuess_firstAllele.second)
			{
				assert(bestGuess_secondAllele_alternatives.count(clusters.second) == 0);
				bestGuess_secondAllele_alternatives[clusters.second] = LLs_normalized.at(cI);
				bestGuess_secondAllele_alternatives_mismatches[clusters.second] = Mismatches_min.at(cI);
			}
			else
			{
				if((int)clusters.second == bestGuess_firstAllele.second)
				{
					assert(bestGuess_secondAllele_alternatives.count(clusters.first) == 0);
					bestGuess_secondAllele_alternatives[clusters.first] = LLs_normalized.at(cI);
					bestGuess_secondAllele_alternatives_mismatches[clusters.first] = Mismatches_min.at(cI);
				}
			}
		}

		std::pair<double, int> oneBestGuess_secondAllele = Utilities::findIntMapMax(bestGuess_secondAllele_alternatives);
		assert(oneBestGuess_secondAllele.first >= 0);
		assert(oneBestGuess_secondAllele.first <= 1);

		std::map<int, double> mismatches_allBestGuessPairs;
		for(std::map<int, double>::iterator secondAlleleIt = bestGuess_secondAllele_alternatives.begin(); secondAlleleIt != bestGuess_secondAllele_alternatives.end(); secondAlleleIt++)
		{
			int cluster = secondAlleleIt->first;
			double LL_normalized = secondAlleleIt->second;
			if(LL_normalized == oneBestGuess_secondAllele.first)
			{
				mismatches_allBestGuessPairs[cluster] =  -1 * bestGuess_secondAllele_alternatives_mismatches.at(cluster);
			}
		}

		std::pair<double, int> bestGuess_secondAllele = Utilities::findIntMapMax(mismatches_allBestGuessPairs);
		std::string bestGuess_secondAllele_oneType = *(HLAtype_clusters.at(bestGuess_secondAllele.second).begin());
		std::string bestGuess_secondAllele_ID = Utilities::join(std::vector<std::string>(HLAtype_clusters.at(bestGuess_secondAllele.second).begin(), HLAtype_clusters.at(bestGuess_secondAllele.second).end()), ";");

		forReturn_starting_haplotype_1_vec.push_back(bestGuess_firstAllele_ID);
		forReturn_starting_haplotype_2_vec.push_back(bestGuess_secondAllele_ID);

		double locus_coverage = (double)HLATypeInference_thisLocus_bases_used / (double)thisLocus_totalColumns;;
		std::cout << "Locus " << locus << " " << HLATypeInference_thisLocus_bases_used << " bases across " << thisLocus_totalColumns << " columns utilized." << "\n";
		std::cout << "\tCoverage " << locus_coverage << "\n";
		
		std::vector<double> positionalCoverages;
		for(unsigned int pI = 0; pI < combined_exon_sequences_graphLevels.size(); pI++)
		{
			int graphLevel = combined_exon_sequences_graphLevels.at(pI);
			int individualExon = graphLevel_2_exonPosition_individualExon.at(graphLevel);
			int exonPosition = graphLevel_2_exonPosition_individualExonPosition.at(graphLevel);
			int coverage = pileUpPerPosition[individualExon][exonPosition].size();
			positionalCoverages.push_back(coverage);
		}
		assert(positionalCoverages.size() > 0);
		std::sort(positionalCoverages.begin(), positionalCoverages.end(), std::less<int>());
		if(positionalCoverages.size() > 1)
		{
			assert(positionalCoverages.at(0) <= positionalCoverages.at(1));
		}

		size_t allColumns_totalAlleles = 0;
		size_t allColumns_incompatibleAlleles = 0;
		std::vector<int> perColumn_totalAlleles;
		std::vector<int> perColumn_incomptibeAlleles;
		double average_perColumn_error_rate;
		double locus_minimumColumnP = -1;
		int n_columns_unaccounted_alleles = 0;

		double proportionkMersCovered_A1;
		double proportionkMersCovered_A2;
		// compute kMer coverages for both called alleles and statistics of deviant alleles
		{
			int clusterID_allele1 = HLAtype_2_clusterID.at(bestGuess_firstAllele_oneType);
			int clusterID_allele2 = HLAtype_2_clusterID.at(bestGuess_secondAllele_oneType);
			assert(clusterID_allele1 == bestGuess_firstAllele.second);
			assert(clusterID_allele2 == bestGuess_secondAllele.second);

			std::string combinedExonSequence_allele1 = combined_exon_sequences.at(bestGuess_firstAllele_oneType);
			std::string combinedExonSequence_allele2 = combined_exon_sequences.at(bestGuess_secondAllele_oneType);

			assert(combinedExonSequence_allele1.length() == combinedExonSequence_allele2.length());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combinedExonSequence_allele1.size());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combinedExonSequence_allele2.size());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combined_exon_sequences_graphLevels_individualExonPosition.size());

			std::vector<std::string> combinedExonSequence_allele1_byExon;
			std::vector<std::string> combinedExonSequence_allele2_byExon;

			for(unsigned int pI = 0; pI < combinedExonSequence_allele1.size(); pI++)
			{
				int thisPosition_individualExon = combined_exon_sequences_graphLevels_individualExon.at(pI);
				int thisPosition_individualExonPosition = combined_exon_sequences_graphLevels_individualExonPosition.at(pI);

				bool needNewExon = ((pI == 0) || (combined_exon_sequences_graphLevels_individualExon.at(pI) != combined_exon_sequences_graphLevels_individualExon.at(pI-1)));
				if(needNewExon)
				{
					combinedExonSequence_allele1_byExon.push_back("");
					combinedExonSequence_allele2_byExon.push_back("");
				}

				std::string underlyingAllele1 = combinedExonSequence_allele1.substr(pI, 1);
				std::string underlyingAllele2 = combinedExonSequence_allele2.substr(pI, 1);

				combinedExonSequence_allele1_byExon.back().push_back(underlyingAllele1.at(0));
				combinedExonSequence_allele2_byExon.back().push_back(underlyingAllele2.at(0));

				int totalAlleles = 0;
				int incompatibeAlleles = 0;
				for(auto piledUpAllele : pileUpPerPosition.at(thisPosition_individualExon).at(thisPosition_individualExonPosition))
				{
					totalAlleles++;
					if((piledUpAllele.genotype != underlyingAllele1) && (piledUpAllele.genotype != underlyingAllele2))
					{
						incompatibeAlleles++;
					}
				}

				allColumns_totalAlleles += totalAlleles;
				allColumns_incompatibleAlleles += incompatibeAlleles;

				perColumn_totalAlleles.push_back(totalAlleles);
				perColumn_incomptibeAlleles.push_back(incompatibeAlleles);

				if(perPosition_allele_counts_postFiltering.count(pI))
				{
					int totalCoverage_position = 0;
					for(auto alleleCount : perPosition_allele_counts_postFiltering.at(pI))
					{
						totalCoverage_position += alleleCount.second;
					}
					if(totalCoverage_position >= threshold_reportColumn_forPresenceOfUnaccountedAlleles_minCoverage)
					{
						for(auto alleleCount : perPosition_allele_counts_postFiltering.at(pI))
						{
							if(alleleCount.first == underlyingAllele1)
								continue;
							if(alleleCount.first == underlyingAllele2)
								continue;

							double aF = (double)alleleCount.second / (double)totalCoverage_position;
							if(aF >= threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction)
							{
								n_columns_unaccounted_alleles++;
							}
						}
					}
				}
			}

			auto calculcatekMerPresence = [&](std::vector<std::string> exons) -> double {
				int kMers_total = 0; 
				int kMers_present = 0;
				for(auto exonSeq : exons)
				{
					std::string exonSeqNoGaps = Utilities::removeGaps(exonSeq);
					std::vector<std::string> kMers = partitionStringIntokMers(exonSeqNoGaps, k_for_kMer_index);
					for(auto kMer : kMers)
					{
						kMers_total++;
						if(kMer.find("*") == std::string::npos)
						{
							std::string kMerKey = kMer_canonical_representation(kMer);						
							if(kMer_counts.count(kMerKey) && (kMer_counts.at(kMerKey) > 0))
							{
								kMers_present++;
							}
						}
					}
				}
				if(kMers_total == 0)
				{
					return -1;
					std::cerr << "kMers_total: " << kMers_total << "\n" << std::flush;
					for(unsigned int exonI = 0; exonI < exons.size(); exonI++)
					{
						std::cerr << "\texon " << exonI << ": " << exons.at(exonI) << "\n" << std::flush;
					}
					
				} 
				assert(kMers_total > 0);
				return (double) kMers_present / (double) kMers_total;
			};

			proportionkMersCovered_A1 = calculcatekMerPresence(combinedExonSequence_allele1_byExon);
			proportionkMersCovered_A2 = calculcatekMerPresence(combinedExonSequence_allele2_byExon);

						
			std::ofstream columnErrorRateStream;
			columnErrorRateStream.open(outputFN_columnError.c_str());
			assert(columnErrorRateStream.is_open());
			columnErrorRateStream << Utilities::join({"Column", "Coverage", "ExpectedIncompatible", "ObservedIncompatible", "p"}, "\t") << "\n";

			if(allColumns_totalAlleles > 0)
			{
				average_perColumn_error_rate = (double)allColumns_incompatibleAlleles / (double)allColumns_totalAlleles;
			}
			else
			{
				average_perColumn_error_rate = 0;
			}

			for(unsigned int columnI = 0; columnI < combinedExonSequence_allele1.size(); columnI++)
			{
				std::vector<std::string> outputFields;
				outputFields.push_back(Utilities::ItoStr(columnI));

				int columnCoverage = perColumn_totalAlleles.at(columnI);
				outputFields.push_back(Utilities::ItoStr(columnCoverage));

				double expectedError = average_perColumn_error_rate * columnCoverage;
				outputFields.push_back(Utilities::DtoStr(expectedError));

				int observedError = perColumn_incomptibeAlleles.at(columnI);
				outputFields.push_back(Utilities::ItoStr(observedError));

				double p = 1;
				if(observedError > expectedError)
				{
					std::vector<double> observed;
						observed.push_back(columnCoverage - observedError);
						observed.push_back(observedError);

					std::vector<double> expected;
						expected.push_back(columnCoverage - expectedError);
						expected.push_back(expectedError);

					p = simpleChiSq(observed, expected);

				}

				outputFields.push_back(Utilities::DtoStr(p));
				columnErrorRateStream << Utilities::join(outputFields, "\t") << "\n";
				
				if((locus_minimumColumnP < 0) || (p < locus_minimumColumnP))
				{
					locus_minimumColumnP = p;
				}
			}
		}



		int index_for_decile = (int)((double)positionalCoverages.size() / 10.0);
		double firstDecileCoverage = positionalCoverages.at(index_for_decile);
		double minimumCoverage = positionalCoverages.at(0);
		bestGuess_outputStream << locus << "\t" << 1 << "\t" << bestGuess_firstAllele_ID << "\t" << bestGuess_firstAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A1 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\n";
		bestGuess_outputStream << locus << "\t" << 2 << "\t" << bestGuess_secondAllele_ID << "\t" << oneBestGuess_secondAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A2 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\n" << std::flush;

		if(can_translateToG_locus(locus))
		{
			bool a1_perfectly;
			bool a2_perfectly;
			std::string bestGuess_firstAllele_ID_G = translate_allele_list_to_G_allele({HLAtype_clusters.at(bestGuess_firstAllele.second).begin(), HLAtype_clusters.at(bestGuess_firstAllele.second).end()}, a1_perfectly);
			std::string bestGuess_secondAllele_ID_G = translate_allele_list_to_G_allele({HLAtype_clusters.at(bestGuess_secondAllele.second).begin(), HLAtype_clusters.at(bestGuess_secondAllele.second).end()}, a2_perfectly);

			bestGuess_G_outputStream << locus << "\t" << 1 << "\t" << bestGuess_firstAllele_ID_G << "\t" << bestGuess_firstAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A1 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles <<  "\t" << a1_perfectly << "\n";
			bestGuess_G_outputStream << locus << "\t" << 2 << "\t" << bestGuess_secondAllele_ID_G << "\t" << oneBestGuess_secondAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A2 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\t" << a2_perfectly << "\n" << std::flush;
		}

		unsigned int maxPairPrint = (LLs_completeReads_indices.size() > 10) ? 10 : LLs_completeReads_indices.size();
		for(unsigned int LLi = 0; LLi < maxPairPrint; LLi++)
		{
			unsigned int pairIndex = LLs_completeReads_indices.at(LLi);
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(pairIndex);
			// std::cout << "#" << (LLi+1) << ": " << clusters.first << " / " << clusters.second << ": " << LLs.at(pairIndex) << " absolute and " << exp(LLs.at(pairIndex) - maxPairI.first) << " relative." << "\n" << std::flush;

			std::vector<std::string> cluster1_members(HLAtype_clusters.at(clusters.first).begin(), HLAtype_clusters.at(clusters.first).end());
			std::vector<std::string> cluster2_members(HLAtype_clusters.at(clusters.second).begin(), HLAtype_clusters.at(clusters.second).end());

			// std::cout << "\tcluster " << clusters.first << ": " << Utilities::join(cluster1_members, ", ") << "\n";
			// std::cout << "\tcluster " << clusters.second << ": " << Utilities::join(cluster2_members, ", ") << "\n" << std::flush;

			// std::cout << "\tMismatches " << Mismatches_avg.at(pairIndex) << " avg / " << Mismatches_min.at(pairIndex) << " min\n" << std::flush;

			// bool equal = (LLs.at(pairIndex) == LLs.at(LLs_indices.at(0)));
			// std::cout << "\tEqual: " << equal << "\n" << std::flush;


		}

		assert(LLs_completeReads.at(LLs_completeReads_indices.at(0)) == maxPairI.first);
		if(LLs_completeReads_indices.size() > 1)
		{
			assert(LLs_completeReads.at(LLs_completeReads_indices.at(0)) >= LLs_completeReads.at(LLs_completeReads_indices.at(1)));
		}
		

	}

	bestGuess_outputStream.close();

	forReturn_starting_haplotype_1 = Utilities::join(forReturn_starting_haplotype_1_vec, ",");
	forReturn_starting_haplotype_2 = Utilities::join(forReturn_starting_haplotype_2_vec, ",");

	std::string outputFN_parameters = outputDirectory + "/R1_parameters.txt";
	std::ofstream outputFN_parameters_outputStream;
	outputFN_parameters_outputStream.open(outputFN_parameters.c_str());
	assert(outputFN_parameters_outputStream.is_open());  
	outputFN_parameters_outputStream << "Loci" << " = " << forReturn_lociString << "\n";
	outputFN_parameters_outputStream << "restrictToFullHaplotypes" << " = " << restrictToFullHaplotypes << "\n";
	outputFN_parameters_outputStream << "veryConservativeReadLikelihoods" << " = " << veryConservativeReadLikelihoods << "\n";
	outputFN_parameters_outputStream.close();
	
	double coverage = (double)HLATypeInference_totalBases_used / (double) HLAtypeInference_totalColumns;
	std::cout << "Sample " << sampleName << " " << HLATypeInference_totalBases_used << " bases utilized (total) across" << HLAtypeInference_totalColumns <<  "columns\n" << std::flush;
	std::cout << "\tCoverage " << coverage << "\n" << std::flush;
}



void fillAS()
{
	if(codon2AS.size() == 0)
	{
		codon2AS["TAG"] = "End";
		codon2AS["CTT"] = "Leu";
		codon2AS["GCC"] = "Ala";
		codon2AS["GGA"] = "Gly";
		codon2AS["GTC"] = "Val";
		codon2AS["TGC"] = "Cys";
		codon2AS["AGT"] = "Ser";
		codon2AS["TGT"] = "Cys";
		codon2AS["TGA"] = "End";
		codon2AS["TCA"] = "Ser";
		codon2AS["CGA"] = "Arg";
		codon2AS["ATT"] = "Ile";
		codon2AS["TAT"] = "Tyr";
		codon2AS["ATC"] = "Ile";
		codon2AS["AAC"] = "Asn";
		codon2AS["AGC"] = "Ser";
		codon2AS["TAC"] = "Tyr";
		codon2AS["AAT"] = "Asn";
		codon2AS["TCG"] = "Ser";
		codon2AS["ACT"] = "Thr";
		codon2AS["ACA"] = "Thr";
		codon2AS["CAA"] = "Gln";
		codon2AS["GAC"] = "Asp";
		codon2AS["CCG"] = "Pro";
		codon2AS["CTG"] = "Leu";
		codon2AS["GGT"] = "Gly";
		codon2AS["GCA"] = "Ala";
		codon2AS["AAG"] = "Lys";
		codon2AS["GTG"] = "Val";
		codon2AS["TCC"] = "Ser";
		codon2AS["TTT"] = "Phe";
		codon2AS["CAC"] = "His";
		codon2AS["GTT"] = "Val";
		codon2AS["AGG"] = "Arg";
		codon2AS["CGT"] = "Arg";
		codon2AS["CAT"] = "His";
		codon2AS["CGG"] = "Arg";
		codon2AS["AGA"] = "Arg";
		codon2AS["ATA"] = "Ile";
		codon2AS["CCC"] = "Pro";
		codon2AS["GGG"] = "Gly";
		codon2AS["ACC"] = "Thr";
		codon2AS["TTA"] = "Leu";
		codon2AS["GAG"] = "Glu";
		codon2AS["CCA"] = "Pro";
		codon2AS["CTA"] = "Leu";
		codon2AS["GAT"] = "Asp";
		codon2AS["TCT"] = "Ser";
		codon2AS["TGG"] = "Trp";
		codon2AS["TTC"] = "Phe";
		codon2AS["CTC"] = "Leu";
		codon2AS["CGC"] = "Arg";
		codon2AS["TTG"] = "Leu";
		codon2AS["GCG"] = "Ala";
		codon2AS["TAA"] = "End";
		codon2AS["GGC"] = "Gly";
		codon2AS["CAG"] = "Gln";
		codon2AS["GCT"] = "Ala";
		codon2AS["GAA"] = "Glu";
		codon2AS["CCT"] = "Pro";
		codon2AS["ACG"] = "Thr";
		codon2AS["AAA"] = "Lys";
		codon2AS["ATG"] = "Met";
		codon2AS["GTA"] = "Val";
		vector<string> ASLeucodons;
		ASLeucodons.push_back("CTG");
		ASLeucodons.push_back("CTA");
		ASLeucodons.push_back("CTT");
		ASLeucodons.push_back("CTC");
		ASLeucodons.push_back("TTG");
		ASLeucodons.push_back("TTA");
		AS2codon["Leu"] = ASLeucodons;
		vector<string> ASLyscodons;
		ASLyscodons.push_back("AAG");
		ASLyscodons.push_back("AAA");
		AS2codon["Lys"] = ASLyscodons;
		vector<string> ASGlucodons;
		ASGlucodons.push_back("GAG");
		ASGlucodons.push_back("GAA");
		AS2codon["Glu"] = ASGlucodons;
		vector<string> ASGlncodons;
		ASGlncodons.push_back("CAG");
		ASGlncodons.push_back("CAA");
		AS2codon["Gln"] = ASGlncodons;
		vector<string> ASEndcodons;
		ASEndcodons.push_back("TGA");
		ASEndcodons.push_back("TAG");
		ASEndcodons.push_back("TAA");
		AS2codon["End"] = ASEndcodons;
		vector<string> ASHiscodons;
		ASHiscodons.push_back("CAT");
		ASHiscodons.push_back("CAC");
		AS2codon["His"] = ASHiscodons;
		vector<string> ASAsncodons;
		ASAsncodons.push_back("AAT");
		ASAsncodons.push_back("AAC");
		AS2codon["Asn"] = ASAsncodons;
		vector<string> ASSercodons;
		ASSercodons.push_back("AGT");
		ASSercodons.push_back("AGC");
		ASSercodons.push_back("TCG");
		ASSercodons.push_back("TCA");
		ASSercodons.push_back("TCT");
		ASSercodons.push_back("TCC");
		AS2codon["Ser"] = ASSercodons;
		vector<string> ASTrpcodons;
		ASTrpcodons.push_back("TGG");
		AS2codon["Trp"] = ASTrpcodons;
		vector<string> ASAlacodons;
		ASAlacodons.push_back("GCG");
		ASAlacodons.push_back("GCA");
		ASAlacodons.push_back("GCT");
		ASAlacodons.push_back("GCC");
		AS2codon["Ala"] = ASAlacodons;
		vector<string> ASArgcodons;
		ASArgcodons.push_back("AGG");
		ASArgcodons.push_back("AGA");
		ASArgcodons.push_back("CGG");
		ASArgcodons.push_back("CGA");
		ASArgcodons.push_back("CGT");
		ASArgcodons.push_back("CGC");
		AS2codon["Arg"] = ASArgcodons;
		vector<string> ASCyscodons;
		ASCyscodons.push_back("TGT");
		ASCyscodons.push_back("TGC");
		AS2codon["Cys"] = ASCyscodons;
		vector<string> ASGlycodons;
		ASGlycodons.push_back("GGG");
		ASGlycodons.push_back("GGA");
		ASGlycodons.push_back("GGT");
		ASGlycodons.push_back("GGC");
		AS2codon["Gly"] = ASGlycodons;
		vector<string> ASAspcodons;
		ASAspcodons.push_back("GAT");
		ASAspcodons.push_back("GAC");
		AS2codon["Asp"] = ASAspcodons;
		vector<string> ASPhecodons;
		ASPhecodons.push_back("TTT");
		ASPhecodons.push_back("TTC");
		AS2codon["Phe"] = ASPhecodons;
		vector<string> ASTyrcodons;
		ASTyrcodons.push_back("TAT");
		ASTyrcodons.push_back("TAC");
		AS2codon["Tyr"] = ASTyrcodons;
		vector<string> ASMetcodons;
		ASMetcodons.push_back("ATG");
		AS2codon["Met"] = ASMetcodons;
		vector<string> ASValcodons;
		ASValcodons.push_back("GTG");
		ASValcodons.push_back("GTA");
		ASValcodons.push_back("GTT");
		ASValcodons.push_back("GTC");
		AS2codon["Val"] = ASValcodons;
		vector<string> ASThrcodons;
		ASThrcodons.push_back("ACG");
		ASThrcodons.push_back("ACA");
		ASThrcodons.push_back("ACT");
		ASThrcodons.push_back("ACC");
		AS2codon["Thr"] = ASThrcodons;
		vector<string> ASProcodons;
		ASProcodons.push_back("CCG");
		ASProcodons.push_back("CCA");
		ASProcodons.push_back("CCT");
		ASProcodons.push_back("CCC");
		AS2codon["Pro"] = ASProcodons;
		vector<string> ASIlecodons;
		ASIlecodons.push_back("ATA");
		ASIlecodons.push_back("ATT");
		ASIlecodons.push_back("ATC");
		AS2codon["Ile"] = ASIlecodons;
	}
}

void _read_HLA_alleles_for_haplotypeInference(std::string graphDir, std::string locus, std::vector<int>& combined_sequences_graphLevels, std::vector<std::string>& combined_sequences_graphLevels_individualType, std::vector<int>& combined_sequences_graphLevels_individualTypeNumber, std::vector<int>& combined_sequences_graphLevels_individualPosition, std::vector<std::string>& combined_sequences_locusIDs, std::map<std::string, std::string>& combined_sequences, const std::map<std::string, unsigned int>& graphLocus_2_levels,  const std::vector<std::string>& padding_sequences_to_alleles, bool exonsWithStars)
{
	combined_sequences_graphLevels.clear();
	combined_sequences_graphLevels_individualType.clear();
	combined_sequences_graphLevels_individualTypeNumber.clear();
	combined_sequences_graphLevels_individualPosition.clear();
	combined_sequences_locusIDs.clear();
	combined_sequences.clear();

	// find intron and exon files belonging to locus
	std::vector<std::string> files_in_order;
	std::vector<std::string> files_in_order_type;
	std::vector<int> files_in_order_number;

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			std::string file_locus = split_by_underscore.at(0);



			if(file_locus == locus)
			{
				if((split_by_underscore.at(2) == "intron") || (split_by_underscore.at(2) == "exon"))
				{
					files_in_order.push_back(graphDir + "/" + line);
					files_in_order_type.push_back(split_by_underscore.at(2));
					files_in_order_number.push_back(Utilities::StrtoI(split_by_underscore.at(3)));
				}
				else
				{
					std::vector<std::string> third_split_by_colon = Utilities::split(split_by_underscore.at(2), ".");

					files_in_order.push_back(graphDir + "/" + line);
					files_in_order_type.push_back(third_split_by_colon.at(0));
					files_in_order_number.push_back(1);
				}
			}
		}
	}
	assert(files_in_order.size() > 0);

	std::cout << "read_HLA_alleles_for_haplotypeInference(..): Files read in for locus " << locus << "\n";
	std::cout << Utilities::join(files_in_order, ",   ") << "\n\n" << std::flush;

	if(exonsWithStars || 1)
	{
		// std::string takeTypesFromFile;
		// for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
		// {
			// if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
			// {
				// continue;
			// }
			// if(files_in_order_type.at(fileI) == "exon")
			// {
				// takeTypesFromFile = files_in_order.at(fileI);
				// break;
			// }
		// }
		// if(takeTypesFromFile.length() == 0)
		// {
			// for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
			// {
				// if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
				// {
					// continue;
				// }
				// if(files_in_order_type.at(fileI) == "intron")
				// {
					// takeTypesFromFile = files_in_order.at(fileI);
				// }
			// }
		// }

		// assert(takeTypesFromFile.length());

		// std::cout << "exonsWithStars is on, take types from file " << takeTypesFromFile << "\n" << std::flush;

		
		for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
		{
			if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
			{
				continue;
			}
				
			std::ifstream fileInputStream;
			fileInputStream.open(files_in_order.at(fileI).c_str());
			assert(fileInputStream.is_open());
			std::vector<std::string> file_lines;
			int n_line = 0;
			while(fileInputStream.good())
			{
				std::string line;
				std::getline(fileInputStream, line);
				Utilities::eraseNL(line);
				if(line.length())
				{
					if(n_line > 0)
					{
						std::vector<std::string> line_fields = Utilities::split(line, " ");
						std::string HLA_type = line_fields.at(0);
						combined_sequences[HLA_type] = "";
					}
				}
				n_line++;
			}
			fileInputStream.close();
		}

		std::cout << "exonsWithStars, total sequences " << combined_sequences.size() << "\n" << std::flush;
	}

	for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
	{
		std::cout << Utilities::timestamp() << "\tLocus" << locus << ", file " << fileI << "\n" << std::flush;
		std::string file = files_in_order.at(fileI);

		std::string type = files_in_order_type.at(fileI);
		int typeNumber = files_in_order_number.at(fileI);

		if(! Utilities::fileReadable(file))
		{
			std::cerr << "HLAHaplotypeInference(..): Locus " << locus << ", fileI " << fileI << ": Can't read file " << file << "\n";
		}
		assert(Utilities::fileReadable(file));

		std::ifstream fileInputStream;
		fileInputStream.open(file.c_str());
		assert(fileInputStream.is_open());
		std::vector<std::string> file_lines;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			if(line.length())
			{
				file_lines.push_back(line);
			}
		}
		fileInputStream.close();

		std::string firstLine = file_lines.at(0);
		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
		assert(firstLine_fields.at(0) == "IndividualID");

		std::vector<std::string> exon_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
		std::string first_graph_locusID = exon_level_names.front();
		std::string last_graph_locusID = exon_level_names.back();

		assert(graphLocus_2_levels.count(first_graph_locusID));
		assert(graphLocus_2_levels.count(last_graph_locusID));


		unsigned int first_graph_level = graphLocus_2_levels.at(first_graph_locusID);
		unsigned int last_graph_level = graphLocus_2_levels.at(last_graph_locusID);

		std::cout << Utilities::timestamp() << "\tLocus" << locus << ", fileI " << fileI << ": from " << first_graph_locusID << " (" << first_graph_level << ") to " << last_graph_locusID << " (" << last_graph_level << ").\n" << std::flush;

		assert(last_graph_level > first_graph_level);
		unsigned int expected_allele_length = last_graph_level - first_graph_level + 1;
		if(!(exon_level_names.size() == expected_allele_length))
		{
			std::cerr << "For locus " << locus << " fileI " << fileI << " (" << file << "), we have a problem with expected graph length.\n";
			std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
			std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
			std::cerr << std::flush;
		}
		assert(exon_level_names.size() == expected_allele_length);

		combined_sequences_locusIDs.insert(combined_sequences_locusIDs.end(), exon_level_names.begin(), exon_level_names.end());
		for(unsigned int lI = 0; lI < expected_allele_length; lI++)
		{
			unsigned int graphLevel = first_graph_level + lI;
			assert(graphLocus_2_levels.at(exon_level_names.at(lI)) == graphLevel);
			combined_sequences_graphLevels.push_back(graphLevel);
			combined_sequences_graphLevels_individualType.push_back(type);
			combined_sequences_graphLevels_individualTypeNumber.push_back(typeNumber);
			combined_sequences_graphLevels_individualPosition.push_back(lI);
		}


		if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
		{
			assert(file_lines.size() > 1);

			for(unsigned int sI = 0; sI < padding_sequences_to_alleles.size(); sI++)
			{
				int getLine = sI % (file_lines.size()-1) + 1;
				assert(getLine >= 1);
				assert(getLine < (int)file_lines.size());
				assert(file_lines.at(getLine).length());

				std::vector<std::string> line_fields = Utilities::split(file_lines.at(getLine), " ");
				assert(line_fields.size() == firstLine_fields.size());
				std::string HLA_type = padding_sequences_to_alleles.at(sI);
				std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

				std::string HLA_type_sequence = Utilities::join(line_alleles, "");

				if(exonsWithStars)
				{
					assert(combined_sequences.count(HLA_type));
				}
				// assert(combined_sequences.count(HLA_type));
				combined_sequences[HLA_type] += HLA_type_sequence;

				// std::cout << HLA_type << " " << files_in_order_type.at(fileI) << " add " << HLA_type_sequence.length() << "\n" << std::flush;

			}
		}
		else
		{
			int L = -1;
			std::set<std::string> saw_alleles;
			for(unsigned int lI = 1; lI < file_lines.size(); lI++)
			{
				assert(file_lines.at(lI).length());
				if(file_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					std::string HLA_type_sequence = Utilities::join(line_alleles, "");

						// assert(combined_sequences.count(HLA_type));
					if(exonsWithStars)
					{
						if(!combined_sequences.count(HLA_type))
						{
							std::cerr << "Don't have sequence for type " << HLA_type << "\n" << std::flush;
						}
						assert(combined_sequences.count(HLA_type));
					}
					combined_sequences[HLA_type] += HLA_type_sequence;
					saw_alleles.insert(HLA_type);

					if(L == -1)
					{
						L = HLA_type_sequence.length();
					}
					else
					{
						if(L != (int)HLA_type_sequence.length())
						{
							std::cerr << "Length error" << "\n";
							std::cerr << "\t" << "L" << ": " << L << "\n";
							std::cerr << "\t" << "HLA_type_sequence.length()" << ": " << HLA_type_sequence.length() << "\n";
							std::cerr << "\t" << "file" << ": " << file << "\n";
							std::cerr << "\t" << "HLA_type_sequence" << ": " << HLA_type_sequence << "\n";
							std::cerr << "\t" << "lI" << ": " << lI << "\n";
							std::cerr << std::flush;
						}
						assert(L == (int)HLA_type_sequence.length());
					}
				}
			}

			// if(exonsWithStars)
			// {
				// if(files_in_order_type.at(fileI) == "intron")
				// {
					if(L != -1)
					{
						for(std::map<std::string, std::string>::iterator typeIt = combined_sequences.begin(); typeIt != combined_sequences.end(); typeIt++)
						{
							std::string allele = typeIt->first;
							if(saw_alleles.count(allele) == 0)
							{
								std::string forAppend;
								forAppend.resize(L, '*');
								if(exonsWithStars)
								{
									assert(combined_sequences.count(allele));
								}
								combined_sequences[allele] += forAppend;
							}
						}
					}
				// }
			// }
		}
	}
	
	if(locus == "A")
	{
		if(!combined_sequences.count("A*33:24"))
		{
			std::cerr << "Type missing\nHave:\n";
			for(std::map<std::string, std::string>::iterator typeIt = combined_sequences.begin(); typeIt != combined_sequences.end(); typeIt++)
			{
				std::cerr << " - " << typeIt->first << "\n";
			}
			std::cerr << std::flush;
		}
		assert(combined_sequences.count("A*33:24"));
		
	}
}


void read_HLA_alleles_for_haplotypeInference(std::string graphDir, std::string locus, std::vector<int>& combined_sequences_graphLevels, std::vector<std::string>& combined_sequences_graphLevels_individualType, std::vector<int>& combined_sequences_graphLevels_individualTypeNumber, std::vector<int>& combined_sequences_graphLevels_individualPosition, std::vector<std::string>& combined_sequences_locusIDs, std::map<std::string, std::string>& combined_sequences, const std::map<std::string, unsigned int>& graphLocus_2_levels, bool removeIncompleteSequences, bool setIntronsToStars, bool setPaddingToStars)
{
	// removeIncompleteSequences: remove all sequences that contain stars in exons or introns, or all sequences that contain stars in exons (if setPaddingToStars)
	// setIntronsToStars: set all intronic sequences to stars
	// set all padding sequences to stars

	combined_sequences_graphLevels.clear();
	combined_sequences_graphLevels_individualType.clear();
	combined_sequences_graphLevels_individualTypeNumber.clear();
	combined_sequences_graphLevels_individualPosition.clear();
	combined_sequences_locusIDs.clear();
	combined_sequences.clear();

	// find intron and exon files belonging to locus
	std::vector<std::string> files_in_order;
	std::vector<std::string> files_in_order_type;
	std::vector<int> files_in_order_number;

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			std::string file_locus = split_by_underscore.at(0);



			if(file_locus == locus)
			{
				if((split_by_underscore.at(2) == "intron") || (split_by_underscore.at(2) == "exon"))
				{
					files_in_order.push_back(graphDir + "/" + line);
					files_in_order_type.push_back(split_by_underscore.at(2));
					files_in_order_number.push_back(Utilities::StrtoI(split_by_underscore.at(3)));
				}
				else
				{
					std::vector<std::string> third_split_by_colon = Utilities::split(split_by_underscore.at(2), ".");

					files_in_order.push_back(graphDir + "/" + line);
					files_in_order_type.push_back(third_split_by_colon.at(0));
					files_in_order_number.push_back(1);
				}
			}
		}
	}
	assert(files_in_order.size() > 0);

	std::cout << "read_HLA_alleles_for_haplotypeInference(..): Files read in for locus " << locus << "\n";
	std::cout << Utilities::join(files_in_order, ",   ") << "\n\n" << std::flush;

	std::set<std::string> sequenceIDs;
	for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
	{
		if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
		{
			continue;
		}

		std::ifstream fileInputStream;
		fileInputStream.open(files_in_order.at(fileI).c_str());
		assert(fileInputStream.is_open());
		std::vector<std::string> file_lines;
		int n_line = 0;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			if(line.length())
			{
				if(n_line > 0)
				{
					std::vector<std::string> line_fields = Utilities::split(line, " ");
					std::string HLA_type = line_fields.at(0);
					combined_sequences[HLA_type] = "";
					sequenceIDs.insert(HLA_type);
				}
			}
			n_line++;
		}
		fileInputStream.close();
	}

	std::set<std::string> incompleteSequences;

	for(unsigned int fileI = 0; fileI < files_in_order.size(); fileI++)
	{
		std::cout << Utilities::timestamp() << "\tLocus" << locus << ", file " << fileI << "\n" << std::flush;
		std::string file = files_in_order.at(fileI);

		std::string type = files_in_order_type.at(fileI);
		int typeNumber = files_in_order_number.at(fileI);

		if(! Utilities::fileReadable(file))
		{
			std::cerr << "HLAHaplotypeInference(..): Locus " << locus << ", fileI " << fileI << ": Can't read file " << file << "\n";
		}
		assert(Utilities::fileReadable(file));

		std::ifstream fileInputStream;
		fileInputStream.open(file.c_str());
		assert(fileInputStream.is_open());
		std::vector<std::string> file_lines;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			if(line.length())
			{
				file_lines.push_back(line);
			}
		}
		fileInputStream.close();

		std::string firstLine = file_lines.at(0);
		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
		assert(firstLine_fields.at(0) == "IndividualID");

		std::vector<std::string> exon_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
		std::string first_graph_locusID = exon_level_names.front();
		std::string last_graph_locusID = exon_level_names.back();

		if(!graphLocus_2_levels.count(first_graph_locusID))
		{
			std::cerr << "Don't have " << first_graph_locusID << " in set of " << graphLocus_2_levels.size() << " elements.\n" << std::flush;
		}
		assert(graphLocus_2_levels.count(first_graph_locusID));
		assert(graphLocus_2_levels.count(last_graph_locusID));

		unsigned int first_graph_level = graphLocus_2_levels.at(first_graph_locusID);
		unsigned int last_graph_level = graphLocus_2_levels.at(last_graph_locusID);

		std::cout << Utilities::timestamp() << "\tLocus" << locus << ", fileI " << fileI << ": from " << first_graph_locusID << " (" << first_graph_level << ") to " << last_graph_locusID << " (" << last_graph_level << ").\n" << std::flush;

		assert(last_graph_level > first_graph_level);
		unsigned int expected_allele_length = last_graph_level - first_graph_level + 1;
		if(!(exon_level_names.size() == expected_allele_length))
		{
			std::cerr << "For locus " << locus << " fileI " << fileI << " (" << file << "), we have a problem with expected graph length.\n";
			std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
			std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
			std::cerr << std::flush;
		}
		assert(exon_level_names.size() == expected_allele_length);

		combined_sequences_locusIDs.insert(combined_sequences_locusIDs.end(), exon_level_names.begin(), exon_level_names.end());
		for(unsigned int lI = 0; lI < expected_allele_length; lI++)
		{
			unsigned int graphLevel = first_graph_level + lI;
			assert(graphLocus_2_levels.at(exon_level_names.at(lI)) == graphLevel);
			combined_sequences_graphLevels.push_back(graphLevel);
			combined_sequences_graphLevels_individualType.push_back(type);
			combined_sequences_graphLevels_individualTypeNumber.push_back(typeNumber);
			combined_sequences_graphLevels_individualPosition.push_back(lI);
		}


		if((files_in_order_type.at(fileI) == "paddingLeft") || (files_in_order_type.at(fileI) == "paddingRight"))
		{
			assert(file_lines.size() > 1);
			std::vector<std::vector<std::string>> paddingSequences;
			for(unsigned int lI = 1; lI < file_lines.size(); lI++)
			{
				if(file_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					if(paddingSequences.size() > 0)
					{
						assert(line_alleles.size() == paddingSequences.at(0).size());
					}

					paddingSequences.push_back(line_alleles);
				}
			}

			if(paddingSequences.size() > 0)
			{
				int sI = 0;
				for(std::set<std::string>::iterator sequenceIDit = sequenceIDs.begin(); sequenceIDit != sequenceIDs.end(); sequenceIDit++)
				{
					std::string sequenceID = *sequenceIDit;
					sI++;
					int paddingSequenceI = sI % paddingSequences.size();
					assert(paddingSequenceI >= 0);
					assert(paddingSequenceI < (int)paddingSequences.size());
					std::vector<std::string> paddingSequence = paddingSequences.at(paddingSequenceI);

					if(setPaddingToStars)
					{
						unsigned int paddingLength = paddingSequence.size();
						paddingSequence.clear();
						paddingSequence.resize(paddingLength, "*");
						if(paddingLength > 0)
						{
							assert(paddingSequence.front() == "*");
							assert(paddingSequence.back() == "*");
						}  
					}

					combined_sequences.at(sequenceID) += Utilities::join(paddingSequence, "");
					
					if(sequenceID == "A*01:11N")
					{
						// std::cerr << "Length for " << sequenceID << ": " << combined_sequences.at(sequenceID).length() << "\n" << std::flush;
					}							
				}
			}
		}
		else
		{

			assert(((files_in_order_type.at(fileI) == "exon") || (files_in_order_type.at(fileI) == "intron")));

			std::map<std::string, std::string> local_alleles;
			for(unsigned int lI = 1; lI < file_lines.size(); lI++)
			{
				assert(file_lines.at(lI).length());
				if(file_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					std::string HLA_type_sequence = Utilities::join(line_alleles, "");

					if(local_alleles.size() > 0)
					{
						std::string arbitraryAllele = local_alleles.begin()->first;
						assert(local_alleles.at(arbitraryAllele).length() == HLA_type_sequence.length());
					}

					local_alleles[HLA_type] = HLA_type_sequence;
				}
			}

			if(local_alleles.size())
			{
				std::string arbitraryAllele = local_alleles.begin()->first;
				unsigned int alleleLength = local_alleles.at(arbitraryAllele).length();

				for(std::set<std::string>::iterator sequenceIDit = sequenceIDs.begin(); sequenceIDit != sequenceIDs.end(); sequenceIDit++)
				{
					std::string sequenceID = *sequenceIDit;

					std::string localAlleleSequence;

					if(((files_in_order_type.at(fileI) == "intron") && setIntronsToStars) || (local_alleles.count(sequenceID) == 0))
					{
						localAlleleSequence.resize(alleleLength, '*');
					}
					else
					{
						localAlleleSequence = local_alleles.at(sequenceID);
					}
					assert(localAlleleSequence.length() == alleleLength);

					bool sequence_is_incomplete = false;
					if(setIntronsToStars)
					{
						if(files_in_order_type.at(fileI) == "exon")
						{
							sequence_is_incomplete = (localAlleleSequence.find('*') != std::string::npos);
						}
					}
					else
					{
						sequence_is_incomplete = (localAlleleSequence.find('*') != std::string::npos);
					}

					if(sequence_is_incomplete)
					{
						incompleteSequences.insert(sequenceID);
					}
				
					combined_sequences.at(sequenceID) += localAlleleSequence;
					
					if(sequenceID == "A*01:11N")
					{
						// std::cerr << "Length for " << sequenceID << ": " << combined_sequences.at(sequenceID).length() << "\n" << std::flush;
					}					
				}
			}
		}
	}

	if(removeIncompleteSequences)
	{
		for(std::set<std::string>::iterator incompleteSequenceIt = incompleteSequences.begin(); incompleteSequenceIt != incompleteSequences.end(); incompleteSequenceIt++)
		{
			std::string sequenceID = *incompleteSequenceIt;
			combined_sequences.erase(sequenceID);
		}
	}

	for(std::map<std::string, std::string>::iterator sequenceIt = combined_sequences.begin(); sequenceIt != combined_sequences.end(); sequenceIt++)
	{
		std::string sequenceID = sequenceIt->first;
		if(sequenceIt->second.length() != combined_sequences_locusIDs.size())
		{
			std::cerr << "sequenceID" << ": " << sequenceID << "\n";
			std::cerr << "sequenceIt->second.length()" << ": " << sequenceIt->second.length() << "\n";
			std::cerr << "combined_sequences_locusIDs.size()" << ": " << combined_sequences_locusIDs.size() << "\n";
			// std::cerr << "combined_sequences->second" << ": " << sequenceIt->second << "\n";
			std::cerr << std::flush;
		}
		assert(sequenceIt->second.length() == combined_sequences_locusIDs.size());
		assert(sequenceIt->second.length() == combined_sequences_graphLevels.size());
	}

	if(! removeIncompleteSequences)
	{
		if(locus == "A")
		{
			if(!combined_sequences.count("A*33:24"))
			{
				std::cerr << "Type missing\nHave:\n";
				for(std::map<std::string, std::string>::iterator typeIt = combined_sequences.begin(); typeIt != combined_sequences.end(); typeIt++)
				{
					std::cerr << " - " << typeIt->first << "\n";
				}
				std::cerr << std::flush;
			}
			assert(combined_sequences.count("A*33:24"));

		}
	}
}

void read_HLA_alleles_for_6_8_digits(std::string graphDir, std::string locus, const std::vector<int>& combined_sequences_graphLevels, const std::map<std::string, unsigned int>& graphLocus_2_levels, std::map<std::string, std::vector<std::string>>& combined_sequences)
{
	combined_sequences.clear();


	std::vector<int> local_combined_sequences_graphLevels;
	std::vector<std::string> local_combined_sequences_graphLevels_individualType;
	std::vector<int> local_combined_sequences_graphLevels_individualTypeNumber;
	std::vector<int> local_combined_sequences_graphLevels_individualPosition;

	std::vector<std::string> local_combined_sequences_locusIDs;
	std::map<std::string, std::string> combined_incomplete_sequences_paddingSetToStars;
	std::vector<std::string> starting_haplotypes_combined_vec;

	// std::map<std::string, unsigned int> local_graphLocus_2_levels = graphLocus_2_levels;

	std::map<std::string, std::string> local_combined_sequences;

	read_HLA_alleles_for_haplotypeInference(
		graphDir,
		locus,
		local_combined_sequences_graphLevels,
		local_combined_sequences_graphLevels_individualType,
		local_combined_sequences_graphLevels_individualTypeNumber,
		local_combined_sequences_graphLevels_individualPosition,
		local_combined_sequences_locusIDs,
		local_combined_sequences,
		graphLocus_2_levels,
		false,
		false,
		true
	);

	assert(combined_sequences_graphLevels.size() == local_combined_sequences_graphLevels.size());
	// assert(graphLocus_2_levels.size() == local_graphLocus_2_levels.size());
	// for(std::map<std::string, unsigned int>::const_iterator graphLocusIt = graphLocus_2_levels.begin(); graphLocusIt != graphLocus_2_levels.end(); graphLocusIt++)
	// {
		// assert(local_graphLocus_2_levels.count(graphLocusIt->first));
		// assert(local_graphLocus_2_levels.at(graphLocusIt->first) == graphLocusIt->second);
	// }

	for(unsigned int i = 0; i < combined_sequences_graphLevels.size(); i++)
	{
		assert(combined_sequences_graphLevels.at(i) == local_combined_sequences_graphLevels.at(i));
	}

	std::vector<std::string> splitTestResult = Utilities::split("abc$", "");
	assert(splitTestResult.at(0) == "a");
	assert(splitTestResult.at(1) == "b");
	assert(splitTestResult.at(2) == "c");
	assert(splitTestResult.at(3) == "$");

	for(std::map<std::string, std::string>::iterator sequenceIt = local_combined_sequences.begin(); sequenceIt != local_combined_sequences.end(); sequenceIt++)
	{
		std::vector<std::string> sequence_parts = Utilities::split(sequenceIt->second, "");
		assert(sequence_parts.size() == sequenceIt->second.length());
		assert(sequence_parts.size() == combined_sequences_graphLevels.size());
		combined_sequences[sequenceIt->first] = sequence_parts;
	}
}

int compute_Hamming_distance(const std::string& completeTypeSequence, const std::string& partialTypeSequence, bool ignoreStars)
{
	int forReturn = 0;
	assert(completeTypeSequence.length() == partialTypeSequence.length());
	for(unsigned int i = 0; i < completeTypeSequence.length(); i++)
	{
		const char& c1 = completeTypeSequence.at(i);
		const char& c2 = partialTypeSequence.at(i);
		if(ignoreStars)
		{
			if((c1 == '*') || (c2 == '*'))
			{
				continue;
			}
		}

		if(c1 != c2)
		{
			forReturn++;
		}
	}
	assert(forReturn >= 0);
	assert(forReturn <= (int)completeTypeSequence.length());
	return forReturn;
}


void compute_weird_Edit_distance(const std::vector<std::string>& S1, const std::vector<std::string>& S2, int& comparisons, int& differences)
{
	comparisons = 0;
	differences = 0;
	assert(S1.size() == S2.size());
	for(unsigned int i = 0; i < S1.size(); i++)
	{
		const std::string& s1 = S1.at(i);
		const std::string& s2 = S2.at(i);

		int smaller_L = (s1.length() < s2.length()) ? s1.length() : s2.length();
		int greater_L = (s1.length() > s2.length()) ? s1.length() : s2.length();

		for(int cI = 0; cI < smaller_L; cI++)
		{
			char c1 = s1.at(cI);
			char c2 = s2.at(cI);

			if((c1 == '*') || (c2 == '*'))
			{
				continue;
			}

			comparisons++;

			if(c1 != c2)
			{
				differences++;
			}
		}

		if(smaller_L != greater_L)
		{
			int l_diff = greater_L - smaller_L;
			assert(l_diff > 0);
			comparisons += l_diff;
			differences += l_diff;
		}
	}
}



void oneReadAlignment_2_exonPositions_paired(seedAndExtend_return_local& alignment, oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, seedAndExtend_return_local& paired_alignment, oneRead& paired_read, int read_1_or_2, const std::vector<int>& combined_exon_sequences_graphLevels, const std::map<int, unsigned int>& graphLevel_2_exonPosition)
{
	int alignment_firstLevel = alignment.alignment_firstLevel();
	int alignment_lastLevel = alignment.alignment_lastLevel();

	double thisRead_fractionOK = alignmentFractionOK(alignment);
	double pairedRead_fractionOK = alignmentFractionOK(paired_alignment);

	double thisRead_WeightedCharactersOK = alignmentWeightedOKFraction(read, alignment);
	double pairedRead_WeightedCharactersOK = alignmentWeightedOKFraction(paired_read, paired_alignment);

	std::pair<seedAndExtend_return_local, seedAndExtend_return_local> alignedReadPair = make_pair(alignment, paired_alignment);
	double pairs_strands_OK = alignedReadPair_strandsValid(alignedReadPair);
	double pairs_strands_distance = alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair);

	// std::cout << "This alignment " << alignment_firstLevel << " - " << alignment_lastLevel << "\n";
	// std::cout << "\tvs combined exon " << combined_exon_sequences_graphLevels.front() << " - " << combined_exon_sequences_graphLevels.back() << "\n\n" << std::flush;

	if( ((alignment_firstLevel >= combined_exon_sequences_graphLevels.front()) && (alignment_firstLevel <= combined_exon_sequences_graphLevels.back())) ||
		((alignment_lastLevel >= combined_exon_sequences_graphLevels.front()) && (alignment_lastLevel <= combined_exon_sequences_graphLevels.back())) )
	{
		std::vector<oneExonPosition> readAlignment_exonPositions;

		int indexIntoOriginalReadData = -1;
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
			int graphLevel = alignment.graph_aligned_levels.at(cI);
			unsigned char alignmentQualityCharacter_thisPosition = alignment.mapQ_genomic_perPosition.at(cI);
			double alignmentQuality_thisPosition = Utilities::PhredToPCorrect(alignmentQualityCharacter_thisPosition);

			if(graphLevel == -1)
			{
				// insertion relative to the graph - we need to extend last character

				assert(graphCharacter == "_");
				assert(sequenceCharacter != "_");

				indexIntoOriginalReadData++;
				int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
				if(alignment.reverse)
				{
					indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
				}
				assert(indexIntoOriginalReadData_correctlyAligned >= 0);
				assert(indexIntoOriginalReadData_correctlyAligned <(int) read.sequence.length());

				std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
				if(alignment.reverse)
				{
					underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
				}
				assert(underlyingReadCharacter == sequenceCharacter);
				char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

				if(readAlignment_exonPositions.size() > 0)
				{
					readAlignment_exonPositions.back().genotype.append(sequenceCharacter);
					readAlignment_exonPositions.back().alignment_edgelabels.append(graphCharacter);
					readAlignment_exonPositions.back().qualities.push_back(qualityCharacter);
					if(!(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length()))
					{
						assert(readAlignment_exonPositions.back().genotype.length() == (readAlignment_exonPositions.back().qualities.length()+1));
						assert(readAlignment_exonPositions.back().genotype.at(0) == '_');
						readAlignment_exonPositions.back().genotype = readAlignment_exonPositions.back().genotype.substr(1);
						readAlignment_exonPositions.back().alignment_edgelabels = readAlignment_exonPositions.back().alignment_edgelabels.substr(1);
						assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());

						// std::cerr << "readAlignment_exonPositions.back().genotype.length()" << ": " << readAlignment_exonPositions.back().genotype.length() << "\n";
						// std::cerr << "readAlignment_exonPositions.back().qualities.length()" << ": " << readAlignment_exonPositions.back().qualities.length()<< "\n";

						// std::cerr << std::flush;
					}
					assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());
					assert(readAlignment_exonPositions.back().alignment_edgelabels.length() == readAlignment_exonPositions.back().genotype.length());
				}
			}
			else
			{
				if(sequenceCharacter != "_")
				{
					indexIntoOriginalReadData++;
					int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
					if(alignment.reverse)
					{
						indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
					}
					assert(indexIntoOriginalReadData_correctlyAligned >= 0);
					assert(indexIntoOriginalReadData_correctlyAligned < (int)read.sequence.length());

					std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
					if(alignment.reverse)
					{
						underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
					}
					assert(underlyingReadCharacter == sequenceCharacter);

					if(graphCharacter == "_")
					{

						assert(graphLevel != -1);

						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;

						readAlignment_exonPositions.push_back(thisPosition);

					}
					else
					{
						// two well-defined characters
						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}

				}
				else
				{
					assert(sequenceCharacter == "_");
					if(graphCharacter == "_")
					{
						assert(graphLevel != -1);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}
					else
					{
						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;


						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}
				}
			}
		}

		int alongReadMode = 0;
		int lastPositionInExon = -1;
		for(unsigned int posInAlignment = 0; posInAlignment < readAlignment_exonPositions.size(); posInAlignment++)
		{
			oneExonPosition& thisPosition = readAlignment_exonPositions.at(posInAlignment);
			assert(thisPosition.graphLevel != -1);

			// std::cout << alongReadMode << " " << thisPosition.graphLevel << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n" << std::flush;

			if(graphLevel_2_exonPosition.count(thisPosition.graphLevel))
			{
				// if(!((alongReadMode == 0) || (alongReadMode == 1)))
				// {
					// std::cerr << "alongReadMode" << ": " << alongReadMode << "\n";
					// std::cerr << "alignment.sequence_aligned.substr(cI, 1)" << ": " << alignment.sequence_aligned << "\n";
					// std::cerr << "alignment.graph_aligned.substr(cI, 1)" << ": " << alignment.graph_aligned << "\n";
					// std::cerr << "alignment.graph_aligned_levels)" << ": " << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), ", ") << "\n";
					// std::cerr << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
					// std::cerr << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.front()" << ": " << combined_exon_sequences_graphLevels.front() << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.back()" << ": " << combined_exon_sequences_graphLevels.back() << "\n";
					// for(unsigned int i = 0; i < alignment.graph_aligned_levels.size(); i++)
					// {
						// std::cout << "\t" << i << " " << alignment.graph_aligned_levels.at(i) << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n";
					// }
					// std::cerr << std::flush;
				// }
				if(alongReadMode == 2)
				{
					lastPositionInExon = -1;
				}
				// assert((alongReadMode == 0) || (alongReadMode == 1));


				thisPosition.positionInExon = graphLevel_2_exonPosition.at(thisPosition.graphLevel);
				assert((lastPositionInExon == -1) || ((int)thisPosition.positionInExon == ((int)lastPositionInExon + 1)));
				lastPositionInExon = thisPosition.positionInExon;
				alongReadMode = 1;

				ret_exonPositions.push_back(thisPosition);
			}
			else
			{
				if(alongReadMode == 1)
				{
					alongReadMode = 2;
				}
			}
		}

		// std::cout << "\n";
	}
}


void oneReadAlignment_2_exonPositions_unpaired(seedAndExtend_return_local& alignment, oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, const std::vector<int>& combined_exon_sequences_graphLevels, const std::map<int, unsigned int>& graphLevel_2_exonPosition)
{
	int alignment_firstLevel = alignment.alignment_firstLevel();
	int alignment_lastLevel = alignment.alignment_lastLevel();

	double thisRead_fractionOK = alignmentFractionOK(alignment);
	double pairedRead_fractionOK = -1;

	double thisRead_WeightedCharactersOK = alignmentWeightedOKFraction(read, alignment);
	double pairedRead_WeightedCharactersOK = -1;

	double pairs_strands_OK = 1;
	double pairs_strands_distance = -1;

	// std::cout << "This alignment " << alignment_firstLevel << " - " << alignment_lastLevel << "\n";
	// std::cout << "\tvs combined exon " << combined_exon_sequences_graphLevels.front() << " - " << combined_exon_sequences_graphLevels.back() << "\n\n" << std::flush;

	if( ((alignment_firstLevel >= combined_exon_sequences_graphLevels.front()) && (alignment_firstLevel <= combined_exon_sequences_graphLevels.back())) ||
		((alignment_lastLevel >= combined_exon_sequences_graphLevels.front()) && (alignment_lastLevel <= combined_exon_sequences_graphLevels.back())) )
	{
		std::vector<oneExonPosition> readAlignment_exonPositions;

		int indexIntoOriginalReadData = -1;
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
			int graphLevel = alignment.graph_aligned_levels.at(cI);
			unsigned char alignmentQualityCharacter_thisPosition = alignment.mapQ_genomic_perPosition.at(cI);
			double alignmentQuality_thisPosition = Utilities::PhredToPCorrect(alignmentQualityCharacter_thisPosition);

			if(graphLevel == -1)
			{
				// insertion relative to the graph - we need to extend last character

				assert(graphCharacter == "_");
				assert(sequenceCharacter != "_");

				indexIntoOriginalReadData++;
				int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
				if(alignment.reverse)
				{
					indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
				}
				assert(indexIntoOriginalReadData_correctlyAligned >= 0);
				assert(indexIntoOriginalReadData_correctlyAligned <(int) read.sequence.length());

				std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
				if(alignment.reverse)
				{
					underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
				}
				assert(underlyingReadCharacter == sequenceCharacter);
				char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

				if(readAlignment_exonPositions.size() > 0)
				{
					readAlignment_exonPositions.back().genotype.append(sequenceCharacter);
					readAlignment_exonPositions.back().alignment_edgelabels.append(graphCharacter);
					readAlignment_exonPositions.back().qualities.push_back(qualityCharacter);
					if(!(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length()))
					{
						assert(readAlignment_exonPositions.back().genotype.length() == (readAlignment_exonPositions.back().qualities.length()+1));
						assert(readAlignment_exonPositions.back().genotype.at(0) == '_');
						readAlignment_exonPositions.back().genotype = readAlignment_exonPositions.back().genotype.substr(1);
						readAlignment_exonPositions.back().alignment_edgelabels = readAlignment_exonPositions.back().alignment_edgelabels.substr(1);
						assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());

						// std::cerr << "readAlignment_exonPositions.back().genotype.length()" << ": " << readAlignment_exonPositions.back().genotype.length() << "\n";
						// std::cerr << "readAlignment_exonPositions.back().qualities.length()" << ": " << readAlignment_exonPositions.back().qualities.length()<< "\n";

						// std::cerr << std::flush;
					}
					assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());
					assert(readAlignment_exonPositions.back().alignment_edgelabels.length() == readAlignment_exonPositions.back().genotype.length());
				}
			}
			else
			{
				if(sequenceCharacter != "_")
				{
					indexIntoOriginalReadData++;
					int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
					if(alignment.reverse)
					{
						indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
					}
					assert(indexIntoOriginalReadData_correctlyAligned >= 0);
					assert(indexIntoOriginalReadData_correctlyAligned < (int)read.sequence.length());

					std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
					if(alignment.reverse)
					{
						underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
					}
					assert(underlyingReadCharacter == sequenceCharacter);

					if(graphCharacter == "_")
					{

						assert(graphLevel != -1);

						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						readAlignment_exonPositions.push_back(thisPosition);

					}
					else
					{
						// two well-defined characters
						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}

				}
				else
				{
					assert(sequenceCharacter == "_");
					if(graphCharacter == "_")
					{
						assert(graphLevel != -1);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}
					else
					{
						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ_genomic;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;


						thisPosition.read1_ID = read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}
				}
			}
		}

		int alongReadMode = 0;
		int lastPositionInExon = -1;
		for(unsigned int posInAlignment = 0; posInAlignment < readAlignment_exonPositions.size(); posInAlignment++)
		{
			oneExonPosition& thisPosition = readAlignment_exonPositions.at(posInAlignment);
			assert(thisPosition.graphLevel != -1);

			// std::cout << alongReadMode << " " << thisPosition.graphLevel << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n" << std::flush;

			if(graphLevel_2_exonPosition.count(thisPosition.graphLevel))
			{
				// if(!((alongReadMode == 0) || (alongReadMode == 1)))
				// {
					// std::cerr << "alongReadMode" << ": " << alongReadMode << "\n";
					// std::cerr << "alignment.sequence_aligned.substr(cI, 1)" << ": " << alignment.sequence_aligned << "\n";
					// std::cerr << "alignment.graph_aligned.substr(cI, 1)" << ": " << alignment.graph_aligned << "\n";
					// std::cerr << "alignment.graph_aligned_levels)" << ": " << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), ", ") << "\n";
					// std::cerr << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
					// std::cerr << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.front()" << ": " << combined_exon_sequences_graphLevels.front() << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.back()" << ": " << combined_exon_sequences_graphLevels.back() << "\n";
					// for(unsigned int i = 0; i < alignment.graph_aligned_levels.size(); i++)
					// {
						// std::cout << "\t" << i << " " << alignment.graph_aligned_levels.at(i) << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n";
					// }
					// std::cerr << std::flush;
				// }
				if(alongReadMode == 2)
				{
					lastPositionInExon = -1;
				}
				// assert((alongReadMode == 0) || (alongReadMode == 1));


				thisPosition.positionInExon = graphLevel_2_exonPosition.at(thisPosition.graphLevel);
				assert((lastPositionInExon == -1) || ((int)thisPosition.positionInExon == ((int)lastPositionInExon + 1)));
				lastPositionInExon = thisPosition.positionInExon;
				alongReadMode = 1;

				ret_exonPositions.push_back(thisPosition);
			}
			else
			{
				if(alongReadMode == 1)
				{
					alongReadMode = 2;
				}
			}
		}

		// std::cout << "\n";
	}
}


std::vector<oneExonPosition> removeDoublePositionsFromRead(const std::vector<oneExonPosition>& positions)
{
	std::map<int, std::vector<oneExonPosition>> positions_per_graphLevel;
	std::vector<oneExonPosition> forReturn;
	
	for(unsigned int i = 0; i < positions.size(); i++)
	{
		const oneExonPosition& p = positions.at(i);
		int graphLevel = p.graphLevel;
		positions_per_graphLevel[graphLevel].push_back(p);
	}

	auto getWorstQuality = [](const std::string& qualities) -> unsigned char {
		assert(qualities.size() > 0);
		unsigned char worstQ;
		for(unsigned int i = 0; i < qualities.size(); i++)
		{
			if((i == 0) || (qualities.at(i) < worstQ))
			{
				worstQ = qualities.at(i);
			}
		}
		return worstQ;
	};
	
	auto getBestExonPosition = [&](const std::vector<oneExonPosition>& alternatives) -> oneExonPosition {
		assert(alternatives.size() > 0);
		
		unsigned int bestI;
		unsigned char bestI_quality;
		
		for(unsigned int i = 0; i < alternatives.size(); i++)
		{
			const oneExonPosition& thisAlternative = alternatives.at(i);
			if(!(((thisAlternative.genotype == "_") || (thisAlternative.qualities.size() > 0))))
			{
				std::cerr << "thisAlternative.qualities.size() == 0" << "\n";
				std::cerr << "thisAlternative.graphLevel: " << thisAlternative.graphLevel << "\n";
				std::cerr << "thisAlternative.genotype: " << thisAlternative.genotype << "\n";
				std::cerr << std::flush;
				
			}
			assert((thisAlternative.genotype == "_") || (thisAlternative.qualities.size() > 0));
			unsigned char Q = (thisAlternative.genotype == "_") ? 0 : getWorstQuality(thisAlternative.qualities);
			if((i == 0) || (Q > bestI_quality))
			{
				bestI = i;
				bestI_quality = Q;
			}
		}
		
		return alternatives.at(bestI);
	};
	
	for(std::map<int, std::vector<oneExonPosition> >::iterator graphLevelIt = positions_per_graphLevel.begin(); graphLevelIt != positions_per_graphLevel.end(); graphLevelIt++)
	{
		oneExonPosition bestPositionAtLevel = getBestExonPosition(graphLevelIt->second);
		forReturn.push_back(bestPositionAtLevel);
	}
	
	// std::cout << "removeDoublePositionsFromRead(..): Removed " << (positions.size() - forReturn.size()) << " duplicated paired-end exon positions (from " << positions.size() << " to " << forReturn.size() << ").\n" << std::flush;
	
	return forReturn;
}


double simpleChiSq(std::vector<double> observed, std::vector<double> expected)
{
        double statistic = 0;
        assert(observed.size() == expected.size());

        double observed_sum = 0;
        for(unsigned int i = 0; i < observed.size(); i++)
        {
                if(observed.at(i) < 0)
                {
                        std::cerr << "Error: the " << i << "-th element of observed is < 0\n" << std::flush;
                }
                assert(observed.at(i) >= 0);
                observed_sum += observed.at(i);
        }
        assert(observed_sum > 0);

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                // observed.at(i) = observed.at(i)/observed_sum;
        }

        double expected_sum = 0;
        for(unsigned int i = 0; i < expected.size(); i++)
        {
                assert(expected.at(i) > 0);
                expected_sum += expected.at(i);
        }
        assert(expected_sum > 0);

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                // expected.at(i) = expected.at(i)/expected_sum;
        }

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                assert(expected.at(i) >= 0);
                double thisSummand = pow((observed.at(i) - expected.at(i)), 2)/expected.at(i);
                if(!(thisSummand >= 0))
                {
                        std::cerr << "thisSummand" << ": " << thisSummand << "\n";
                        std::cerr << "observed.at(i)" << ": " << observed.at(i) << "\n";
                        std::cerr << "expected.at(i)" << ": " << expected.at(i) << "\n";
                        std::cerr << "pow((observed.at(i) - expected.at(i)), 2)" << ": " << pow((observed.at(i) - expected.at(i)), 2) << "\n";
                        std::cerr << std::flush;
                }
                assert(thisSummand >= 0);
                statistic += thisSummand;
        }

        assert(statistic >= 0);
        // std::cout << "ChiSq statistic: " << statistic << "\n";

        int df = observed.size() - 1;
        assert(df > 0);


        boost::math::chi_squared chiSqdist(1);
        double pValue = 1 -  boost::math::cdf(chiSqdist, statistic);
        assert((pValue >= 0) && (pValue <= 1));
        return pValue;
}

bool can_translateToG_locus(std::string locus)
{
	read_G_alleles();
	assert(alleles_to_G.size() > 0);
	assert(locus.find("*") == std::string::npos);
	return (G_loci.count(locus));
}

std::string translate_allele_list_to_G_allele(const std::vector<std::string>& alleles, bool& ret_perfectly)
{
	read_G_alleles();
	assert(alleles_to_G.size() > 0);

	std::map<std::string, int> g_groups;
	for(auto a : alleles)
	{
		std::vector<std::string> locus_and_allele = Utilities::split(a, "*");
		if(locus_and_allele.size() != 2)
		{
			std::string oneAllele = alleles_to_G.begin()->first;
			throw std::runtime_error("Weird allele: "+a+"; have only things like "+oneAllele);
		}
		std::string locus_without_star = locus_and_allele.at(0);
		assert(can_translateToG_locus(locus_without_star));
		if(! alleles_to_G.count(a))
		{
			std::string oneAllele = alleles_to_G.begin()->first;			
			std::cerr << "Warning: Can't G-translate: " << a << "; have only things like " << oneAllele;	
			continue;
		}
		assert(alleles_to_G.count(a));
		std::string g_group = alleles_to_G.at(a);
		if(g_groups.count(g_group) == 0)
		{
			g_groups[g_group] = 0;
		}
		g_groups.at(g_group)++;
	}

	if(g_groups.size() == 0)
	{
		std::string forReturn = Utilities::join(alleles, ";");
		ret_perfectly = false;
		return forReturn;
	}
	assert(g_groups.size());
	
	if(g_groups.size() == 1)
	{
		ret_perfectly = true;
		return g_groups.begin()->first;
	}
	else
	{
		ret_perfectly = false;

		std::vector<std::string> possible_groups_sorted = Utilities::get_map_keys_sorted_by_value(g_groups);
		if(possible_groups_sorted.size() > 1)
		{
			assert(g_groups.at(possible_groups_sorted.at(0)) >= g_groups.at(possible_groups_sorted.at(1)));
		}

		return possible_groups_sorted.at(0);
	}
}

void read_G_alleles()
{
	if(alleles_to_G.size() == 0)
	{
		std::vector<std::string> _test_split = Utilities::split("A*;01:01:22;", ";");
		assert(_test_split.size() == 3);

		std::string g_filename = "hla_nom_g.txt";
		std::ifstream gStream;
		gStream.open(g_filename.c_str());
		if(! gStream.is_open())
		{
			throw std::runtime_error("Can't open file " + g_filename + " - are you executing me from the right directory?");
		}

		std::string line;
		while(gStream.good())
		{
			std::getline(gStream, line);
			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			if(line.substr(0, 1) == "#")
				continue;

			std::vector<std::string> components = Utilities::split(line, ";");
			std::string locus_with_star = components.front();
			assert(locus_with_star.back() == '*');

			std::string locus_without_star = locus_with_star.substr(0, locus_with_star.size() - 1);
			G_loci.insert(locus_without_star);

			std::string g_code;
			if(components.back() != "")
			{
				g_code = components.back();
			}
			else
			{
				assert(components.size() == 3);
				g_code = components.at(1);
			}

			g_code = locus_with_star + g_code;
			std::vector<std::string> differentAlleles = Utilities::split(components.at(1), "/");

			for(auto a : differentAlleles)
			{
				std::string allele_with_locus = locus_with_star + a;
				alleles_to_G[allele_with_locus] = g_code;				
			}
		}
	}
}
