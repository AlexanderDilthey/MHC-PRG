/*
 * HLAtypes.h
 *
 *  Created on: 20.06.2014
 *      Author: AlexanderDilthey
 */

#ifndef HLATYPES_H_
#define HLATYPES_H_

#include <string>
#include <vector>
#include <map>

void HLATypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName, bool restrictToFullHaplotypes, std::string& forReturn_lociString, std::string& forReturn_starting_haplotype_1, std::string& forReturn_starting_haplotype_2);
void HLAHaplotypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName, std::string loci_str, std::string starting_haplotype_1, std::string starting_haplotype_2);

void simulateHLAreads(std::string graphDir, int nIndividuals, bool exon23, bool perturbHaplotypes, bool readError, std::string outputDirectory, std::string qualityMatrixFile, int readLength, double insertSize_mean, double insertSize_sd, double haploidCoverage);
void simulateHLAreads_perturbHaplotype(std::vector<std::string>& haplotype);

void read_HLA_alleles_for_haplotypeInference(std::string graphDir, std::string locus, std::vector<int>& combined_sequences_graphLevels, std::vector<std::string>& combined_sequences_graphLevels_individualType, std::vector<int>& combined_sequences_graphLevels_individualTypeNumber, std::vector<int>& combined_sequences_graphLevels_individualPosition, std::vector<std::string>& combined_sequences_locusIDs, std::map<std::string, std::string>& combined_sequences, const std::map<std::string, unsigned int>& graphLocus_2_levels,  const std::vector<std::string>& padding_sequences_to_alleles, bool exonsWithStars);
void read_HLA_alleles_for_6_8_digits(std::string graphDir, std::string locus, const std::vector<int>& combined_sequences_graphLevels, const std::map<std::string, unsigned int>& graphLocus_2_levels, std::map<std::string, std::vector<std::string>>& combined_sequences);

int compute_Hamming_distance(const std::string& completeTypeSequence, const std::string& partialTypeSequence, bool ignoreStars);
void compute_weird_Edit_distance(const std::vector<std::string>& S1, const std::vector<std::string>& S2, int& comparisons, int& differences);

#endif /* HLATYPES_H_ */
