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

void HLATypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName);
void HLAHaplotypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName, std::string starting_haplotype_1, std::string starting_haplotype_2);

void simulateHLAreads(std::string graphDir, int nIndividuals, bool perturbHaplotypes, bool readError, std::string outputDirectory, std::string qualityMatrixFile, int readLength, double insertSize_mean, double insertSize_sd, double haploidCoverage);
void simulateHLAreads_perturbHaplotype(std::vector<std::string>& haplotype);

#endif /* HLATYPES_H_ */
