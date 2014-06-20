/*
 * HLAtypes.h
 *
 *  Created on: 20.06.2014
 *      Author: AlexanderDilthey
 */

#ifndef HLATYPES_H_
#define HLATYPES_H_

#include <string>

void HLATypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName);
void HLAHaplotypeInference(std::string alignedReads_file, std::string graphDir, std::string sampleName, std::string starting_haplotype_1, std::string starting_haplotype_2);


#endif /* HLATYPES_H_ */
