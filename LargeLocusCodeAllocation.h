/*
 * LocusCodeAllocation.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef LARGELOCUSCODEALLOCATION_H_
#define LARGELOCUSCODEALLOCATION_H_

#include <string>
#include <map>
#include <vector>
#include <set>

using namespace std;

struct LargeMutateAwayTo {
	int targetAllele;
	double p;
};

extern string separatorForSerialization;

class LargeLocusCodeAllocation
{
protected:
	map<string, map<string, int> > coded_values;
	map<string, map<int, string> > coded_values_rev;
	map<string, set<int> > restrictedHLAcache;
	
public:
	LargeLocusCodeAllocation();
	virtual ~LargeLocusCodeAllocation();
	int doCode(string locus, string value);
	int invert(string locus, int emission);

	string deCode(string locus, int value);
	vector<int> getAlleles(string locus);

	// vector<LargeMutateAwayTo> mutateAway(string locus, int baseAllele, int mode, bool restrict2DHLAsingletons);

	bool knowCode(string locus, int code);
	bool knowAllele(string locus, string allele);
	vector<string> getLoci();
	
	set<int> getAlleles_2DrestrictOnSingletons (string locus);
	vector<int> getAlleles4DOnly(string locus);
	string alleleString2D(string allele);
	bool allele4D(string locus, int allele);
	bool locusIsHLA(string locus);

	vector<string> serializeIntoVector();
	void readFromVector(vector<string> lines);
};

#endif /* LARGELOCUSCODEALLOCATION_H_ */
