/*
 * GenotypePanel.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef GENOTYPEPANEL_H_
#define GENOTYPEPANEL_H_
#include <map>
#include <string>
#include <vector>
#include "../LocusCodeAllocation.h"
#include "HaplotypePanel.h"

extern string haplotype_number_separator;

class GenotypePanel {
protected:
	map<string, int> individualID_cache;
	bool preserveNAs;
	
public:
	GenotypePanel();
	GenotypePanel(LocusCodeAllocation lCODE);

	void readFromFile(string filename, string positions, bool doPreserveNAs);
	void writeToFile(string filename);
	basic_string<unsigned char> getIndividualGenotype(vector<string>& loci, string IndividualID);
	basic_string<unsigned char> getIndividualGenotypeHideHLA(vector<string>& loci, string IndividualID);
	
	HaplotypePanel sampleRandomHaplotypes(int samples);
	vector<string> getLoci();

	vector<string> IndividualIDs;
	map<string, basic_string<unsigned char> > GenotypesByLoci;
	map<string, int> LocusPositions;
	LocusCodeAllocation CODE;

};

#endif /* GENOTYPEPANEL_H_ */
