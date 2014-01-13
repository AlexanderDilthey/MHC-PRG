/*
 * HaplotypePanel.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef HAPLOTYPEPANEL_H_
#define HAPLOTYPEPANEL_H_
#include <map>
#include <string>
#include <vector>
#include "../LocusCodeAllocation.h"

extern string haplotype_field_separator;

using namespace std;

class HaplotypePanel {
protected:
	map<string, int> haplotypeID_cache;

public:
	HaplotypePanel();
	HaplotypePanel(LocusCodeAllocation lCODE);
	HaplotypePanel(LocusCodeAllocation lCODE, vector<string> lHaplotypeIDs, map<string, basic_string<unsigned char> > lHaplotypesByLoci, map<string, int> lLocusPositions);

	void readFromFile(string filename, string positions);
	void writeToFile(string filename);

	vector<string> getLoci();

	vector<string> HaplotypeIDs;
	map<string, basic_string<unsigned char> > HaplotypesByLoci;
	map<string, int> LocusPositions;
	map<string, string> LocusStrands;
	LocusCodeAllocation CODE;

	basic_string<unsigned char> getIndividualHaplotype(vector<string>& loci, string haplotypeID, bool pseudoTwoHaplotypes);

};

#endif /* HAPLOTYPEPANEL_H_ */
