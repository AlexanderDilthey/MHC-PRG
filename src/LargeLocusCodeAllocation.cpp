/*
 * LargeLocusCodeAllocation.cpp
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#include "LargeLocusCodeAllocation.h"
#include "MHC-PRG.h"
#include "Utilities.h"
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <set>
#include <algorithm>
#include <iostream>

string separatorForSerialization = "|||";

LargeLocusCodeAllocation::LargeLocusCodeAllocation() {

}

LargeLocusCodeAllocation::~LargeLocusCodeAllocation() {
}

std::string LargeLocusCodeAllocation::deCode(std::string locus, int value)
{
	assert(coded_values_rev.count(locus) > 0);
	if(coded_values_rev[locus].count(value) > 0)
	{
		return coded_values_rev[locus][value];
	}
	else
	{
		cerr << "Non-assigned allele for "+locus+ " value: "+Utilities::ItoStr((unsigned int)value);
		assert(coded_values_rev[locus].count(value) > 0);
		return 0;
	}
}

int LargeLocusCodeAllocation::invert(string locus, int emission)
{
	string nucleotide = deCode(locus, emission);
	string invertedNucleotide;
	if(nucleotide == "A")
	{
		invertedNucleotide = "T";
	}
	else if(nucleotide == "C")
	{
		invertedNucleotide = "G";
	}
	else if(nucleotide == "G")
	{
		invertedNucleotide = "C";
	}
	else if(nucleotide == "T")
	{
		invertedNucleotide = "A";
	}
	else
	{
		errEx("Cannot invert this nucleotide: "+nucleotide);
		assert(1==0);
	}
	return doCode(locus, invertedNucleotide);
}
int LargeLocusCodeAllocation::doCode(std::string locus, std::string value)
{
	if((value == "?") || (value == "????"))
	{
		return '0';
	}

	if(coded_values[locus].count(value) > 0)
	{
		return coded_values[locus][value];
	}
	else
	{
		int start = '0';
		int elements_in = coded_values[locus].size();
		int new_index = (int) start+elements_in+1;

		coded_values[locus][value] = new_index;
		coded_values_rev[locus][new_index] = value;

		return new_index;
	}
}

vector<string> LargeLocusCodeAllocation::getLoci()
{
	vector<string> loci;
	for(map<string, map<string, int> >::iterator StringCodes = coded_values.begin(); StringCodes != coded_values.end(); StringCodes++)
	{
		string locus = StringCodes->first;
		loci.push_back(locus);
	}
	return loci;	
}

bool LargeLocusCodeAllocation::knowCode(string locus, int code)
{
	return (coded_values_rev[locus].count(code) > 0);
}

bool LargeLocusCodeAllocation::knowAllele(string locus, string allele)
{
	return (coded_values[locus].count(allele) > 0);
}


vector<int> LargeLocusCodeAllocation::getAlleles(string locus)
{
	vector<int> alleles;
	for(map<int, string>::iterator CodeAllele = coded_values_rev[locus].begin(); CodeAllele != coded_values_rev[locus].end(); CodeAllele++)
	{
		int allele = CodeAllele->first;
		assert(allele != '0');
		alleles.push_back(allele);
	}
	return alleles;
}

/*
vector<LargeMutateAwayTo> LargeLocusCodeAllocation::mutateAway(string locus, int baseAllele, int mode, bool restrict2DHLAsingletons)
{
//	*
//	 * The current behaviour of this function is best described as follows:
//	 * 	- assume that we are using it on a non-HLA locus. Then, the probability to observe the same allele as the baseAllele
//	 * 	  is given by the mode setting, and if we observe missing data, it is split evenly among all possible (usually
//	 * 	  SNP) alleles.
//	 * 	- now, think about HLA alleles. In general, the function does the same as for SNPs, except for one special
//	 * 	  situation: mode == 3 (i.e. we are doing HLA type inference) && restrict2DHLAsingletons == true
//	 * 	  (i.e. we only want, wherever possible, 4-digit alleles in the inference results). In this case, the missing
//	 * 	  data probabilities are manipulated to exclude all unnecessary 2-digit alleles, and the probabilities are split
//	 * 	  evenly for missing data (remember that loci that inference should be made for are marked as missing data).
//	 *

	vector<int> locus_alleles = getAlleles(locus);
	int num_of_alleles = locus_alleles.size();
	bool locus_is_HLA = locusIsHLA(locus);
	
	assert(baseAllele != '0');

	double emission_same_allele = 0.0;
	// todo check later if required


	if(mode == 0)
	{
		emission_same_allele = 1 - CONFIG.mutation_build_graph;
		assert(restrict2DHLAsingletons == false);
	}
	else if(mode == 1)
	{
		emission_same_allele = 1 - CONFIG.mutation_sample_graph;
		assert(restrict2DHLAsingletons == false);
	}
	else if(mode == 2)
	{
		emission_same_allele = 1 - CONFIG.mutation_inference;

	}
	else
	{
		assert(string("Unknown mode for mutateAway") == string(""));
	}
	
	if(locus_is_HLA)
	{
		double mutation = 1 - emission_same_allele;
		mutation = mutation / CONFIG.HLA_mutation_scaling;
		emission_same_allele = 1 - mutation;
	}


	vector<LargeMutateAwayTo> forReturn;

	if(! restrict2DHLAsingletons)
	{
		LargeMutateAwayTo missingData;
		missingData.targetAllele = '0';
		missingData.p = (double)1/(double)num_of_alleles;
		forReturn.push_back(missingData);
	}
	else
	{
		if(locus_is_HLA)
		{
			set<int> locus_alleles_2Drestricted = getAlleles_2DrestrictOnSingletons(locus);
			int num_of_alleles_2Drestricted = locus_alleles_2Drestricted.size();
			
			if(locus_alleles_2Drestricted.count(baseAllele) > 0)
			{
				LargeMutateAwayTo missingData;
				missingData.targetAllele = '0';
				missingData.p = (double)1/(double)num_of_alleles_2Drestricted;
				forReturn.push_back(missingData);
			}
			else
			{
				LargeMutateAwayTo missingData;
				missingData.targetAllele = '0';
				missingData.p = 0;
				forReturn.push_back(missingData);
			}
		}
		else
		{
			LargeMutateAwayTo missingData;
			missingData.targetAllele = '0';
			missingData.p = (double)1/(double)num_of_alleles;
			forReturn.push_back(missingData);
		}

	}

	for(unsigned int i = 0; i < locus_alleles.size(); i++)
	{
		int c = locus_alleles.at(i);
		LargeMutateAwayTo mut;
		mut.targetAllele = c;
		if(c == baseAllele)
		{
			if(num_of_alleles == 1)
			{
				mut.p = 1;
			}
			else
			{
				mut.p = emission_same_allele;
			}
		}
		else
		{
			mut.p = (1 - emission_same_allele) / ((double)num_of_alleles-1);
		}
		forReturn.push_back(mut);

	}

	return forReturn;
}

*/

bool LargeLocusCodeAllocation::locusIsHLA(string locus)
{
	if(locus.length() < 3)
		return false;

	return ((locus.substr(0, 3) == "HLA") || (locus.substr(0, 3) == "KIR"));
}

bool LargeLocusCodeAllocation::allele4D(string locus, int allele)
{
	assert(locusIsHLA(locus) == true);
	string actualType = deCode(locus, allele);
	assert(actualType.length() == 4);
	return (actualType.substr(2,2) != "00");
}

string LargeLocusCodeAllocation::alleleString2D(string allele)
{
	assert(allele.length() == 4);
	return allele.substr(0,2)+"00";
}

vector<int> LargeLocusCodeAllocation::getAlleles4DOnly(string locus)
{
	assert(locusIsHLA(locus) == true);
	vector<int> all_alleles = getAlleles(locus);
	vector<int> alleles_4D;
	for(unsigned int i = 0; i < all_alleles.size(); i++)
	{
		if(allele4D(locus, all_alleles.at(i)))
		{
			alleles_4D.push_back(all_alleles.at(i));
		}
	}
	return alleles_4D;
}

set<int> LargeLocusCodeAllocation::getAlleles_2DrestrictOnSingletons (string locus)
{
	if(restrictedHLAcache.count(locus) > 0)
	{
		return restrictedHLAcache[locus];
	}
	else
	{
		assert(locusIsHLA(locus) == true);
	
		set<string> alleles_2D_implied_by_4D;
		vector<int> alleles_4D = getAlleles4DOnly(locus);
		for(unsigned int i = 0; i < alleles_4D.size(); i++)
		{
			string allele_4D_string = deCode(locus, alleles_4D.at(i));
			string allele_2D_string = alleleString2D(allele_4D_string);
			alleles_2D_implied_by_4D.insert(allele_2D_string);
		}
	
		set<int> forReturn;
		vector<int> alleles_all = getAlleles(locus);
		for(unsigned int i = 0; i < alleles_all.size(); i++)
		{
			int allele = alleles_all.at(i);
			if(allele4D(locus, allele))
			{
				forReturn.insert(allele);
			}
			else
			{
				string alleleString = deCode(locus, allele);
				assert(alleleString.length() == 4);
				if(alleles_2D_implied_by_4D.count(alleleString) == 0)
				{
					forReturn.insert(allele);
				}
			}
		}
		
		restrictedHLAcache[locus] = forReturn;
		
		return forReturn;
	}
}

vector<string> LargeLocusCodeAllocation::serializeIntoVector()
{
	vector<string> forReturn;
	for(map<string, map<string, int> >::iterator locusIt = coded_values.begin(); locusIt != coded_values.end(); locusIt++)
	{
		string locus = locusIt->first;
		for(map<string, int>::iterator codeIt = coded_values[locus].begin(); codeIt != coded_values[locus].end(); codeIt++)
		{
			string originalAllele = codeIt->first;
			int codedAllele = codeIt->second;

			vector<string> lineVector;
			lineVector.push_back(locus);
			lineVector.push_back(originalAllele);
			lineVector.push_back(Utilities::ItoStr((unsigned int)codedAllele));

			forReturn.push_back(boost::algorithm::join(lineVector, separatorForSerialization));
		}
	}

	return forReturn;
}

void LargeLocusCodeAllocation::readFromVector(vector<string> lines)
{
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		string line = lines.at(i);
		vector<string> lineFields;
		boost::iter_split(lineFields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));

		if(lineFields.size() != 3)
		{
			errEx("Cannot read CODE from line, expect 3 fields, got" + Utilities::ItoStr(lineFields.size())+ "! Line: "+line);
		}
		string locus = lineFields.at(0);
		string originalAllele = lineFields.at(1);
		string integerString = lineFields.at(2);
		int codedChar = Utilities::StrtoI(integerString);
		if(codedChar < 0)
		{
			errEx("Weird codedChar value: cannot convert back! "+integerString);
		}
		int codedAllele = (int) codedChar;

		coded_values[locus][originalAllele] = codedAllele;
		coded_values_rev[locus][codedAllele] = originalAllele;

//		if(locus == "L106379")
//		{
//			cout << locus << " -" << codedAllele << "-==" << integerString << " stands for " << originalAllele << "\n";
//		}
	}
}

