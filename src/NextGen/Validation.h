/*
 * Validation.h
 *
 *  Created on: 13.06.2013
 *      Author: AlexanderDilthey
 */

#ifndef VALIDATION_H_
#define VALIDATION_H_

#include <string>
#include <vector>
#include <set>
#include <map>
#include "../Graph/Graph.h"


using namespace std;

typedef std::vector< std::vector<std::string> > diploidGenomeString;

Graph* genomeString2Graph(diploidGenomeString gS, bool verbose = false);
diploidGenomeString VCF2GenomeString(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::string referenceGenomePath, std::vector<std::vector<int> >& ret_graph_referencePositions, bool ignoreVCF = false);

Graph* VCF2Graph(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::string referenceGenomePath);
void storeGenomeStringInFile(diploidGenomeString& gS, std::string filename);
diploidGenomeString readGenomeStringFromFile(std::string filename, bool ignorePuffer = false);
diploidGenomeString readGenomeStringFromChromotypesFile(std::string filename, int kMer_size, std::string graphDir, int chromotypes_startCoordinate, std::vector<int>& amendedChromotypes_genomicGraphLoci);

void validateGenomeString(diploidGenomeString& gS, std::string deBruijnGraph_fileName, int kMer_size);
void testValidation(std::string viterbi_diploid_gS,std::string deBruijnGraph, int k, std::string VCFfile, std::string referenceGenomeFile);
diploidGenomeString compressGenomeString(diploidGenomeString gS);

void validateChromotypesVsVCF(std::string chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width);
void validateAmendedChromotypesVsVCF(std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width);
void validateAllChromotypesVsVCF(std::string chromotypes_file, std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width, std::string outputDirectory, std::string graphDir);
void alignContigsToAllChromotypes(std::string chromotypes_file, std::string amended_chromotypes_file, int chromotypes_startCoordinate, int chromotypes_stopCoordinate, std::string VCFfile, int VCF_minRange, int VCF_maxRange, std::string referenceGenome, std::string deBruijnGraph, int kMer_size, int cortex_height, int cortex_width, std::string outputDir_contigs, std::string contigsFile_Fasta, std::string graphDir);

void vennDiagrams(std::vector<std::string> setNames, std::vector<std::set<std::string>*> kMers, std::vector<std::set<std::string>*> kMers_present, std::vector<std::map<std::string, double>* > kMer_optimalities, std::string outputFile);

void _addPufferChromotypes(diploidGenomeString& gS);

std::vector<int> getGenomicGraphLoci(std::string graphDir, int chromotypes_startCoordinate);
std::vector<std::string> readGraphLoci(std::string graphDir);
std::vector<int> graphLoci_2_PGFpositions(std::vector<std::string> graphLoci);
#endif /* VALIDATION_H_ */
