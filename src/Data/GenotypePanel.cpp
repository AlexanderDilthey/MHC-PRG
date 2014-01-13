/*
 * GenotypePanel.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "GenotypePanel.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "../MHC-PRG.h"
#include "HaplotypePanel.h"
#include "../Utilities.h"

string haplotype_number_separator = "-x-";

GenotypePanel::GenotypePanel() {
	LocusCodeAllocation lCODE;
	CODE = lCODE;
	preserveNAs = false;
}

GenotypePanel::GenotypePanel(LocusCodeAllocation lCODE)
{
	CODE = lCODE;
}

basic_string<unsigned char> GenotypePanel::getIndividualGenotype(vector<string>& loci, string IndividualID)
{
	basic_string<unsigned char> forReturn;
	
	if(individualID_cache.count(IndividualID) == 0)
	{
		int found_id = -1;
		for(unsigned int i = 0; i < IndividualIDs.size(); i++)
		{
			if(IndividualIDs.at(i) == IndividualID)
			{
				found_id = i;
			}
		}
		assert(found_id != -1);
		individualID_cache[IndividualID] = found_id;
	}
	
	assert(IndividualIDs.at(individualID_cache[IndividualID]) == IndividualID);
	
	for(unsigned int l = 0; l < loci.size(); l++)
	{
		string locus = loci.at(l);
		if(GenotypesByLoci.count(locus) == 0)
		{
			cout << "Problem: requested locus " << locus << ", but not present in genotype panel -- add as missing data!\n";
		}
		assert(GenotypesByLoci.count(locus) > 0);
		basic_string<unsigned char> toAppend = GenotypesByLoci[locus].substr(individualID_cache[IndividualID]*2, 2);
		forReturn.append(toAppend);
	}
	
	return forReturn;
}

basic_string<unsigned char> GenotypePanel::getIndividualGenotypeHideHLA(vector<string>& loci, string IndividualID)
{
	basic_string<unsigned char> forReturn = getIndividualGenotype(loci, IndividualID);
	for(unsigned int l = 0; l < loci.size(); l++)
	{
		string locus = loci.at(l);
		if(CODE.locusIsHLA(locus))
		{
			forReturn.at(2*l) = '0';
			forReturn.at(2*l+1) = '0';
		}
	}
	return forReturn;
}

void GenotypePanel::readFromFile(string filename, string positions, bool doPreserveNAs)
{
	preserveNAs = doPreserveNAs;
	
	ifstream positionsStr;
	positionsStr.open (positions.c_str(), ios::in);
	if(positionsStr.is_open())
	{
		string line;
		while(positionsStr.good())
		{
			getline (positionsStr, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			vector<string> fields = Utilities::split(line, ' ');
			if(fields.size() != 2)
				errEx("Strange format in positions file. Expect two fields: field name and position (integer), separated by a whitespace. Problem occured in line: "+line);

			string locusID = fields.at(0);
			string posStr = fields.at(1);
			int pos = Utilities::StrtoI(posStr);

			LocusPositions[locusID] = pos;
		}
		positionsStr.close();
	}
	else
	{
		errEx("Cannot open positions file: "+positions);
	}

	ifstream genotypes;
	genotypes.open (filename.c_str(), ios::in);

	int line_counter = -1;
	if (genotypes.is_open())
	{
		string line;
		if(! genotypes.good())
		{
			errEx("Read error: first line!");
		}
		getline (genotypes, line);
		Utilities::eraseNL(line);

		vector<string> header_fields = Utilities::split(line, ' ');

		if(header_fields.at(0) != "IndividualID")
			errEx("File error: no IndividualID!");
		if(header_fields.at(1) == "Chromosome")
			errEx("File error: there should not be a chromosome field!");

		vector<string> locus_IDs = header_fields;
		locus_IDs.erase(locus_IDs.begin(), locus_IDs.begin()+1);

		while(genotypes.good())
		{
			line_counter++;

			getline (genotypes, line);
			Utilities::eraseNL(line);
			vector<string> line_fields = Utilities::split(line, ' ');

			if(line_fields.size() == 0)
				continue;

			if(line_fields.size() != header_fields.size())
				errEx("Wrong field size: " + Utilities::ItoStr(line_fields.size()) + " vs expected " + Utilities::ItoStr(header_fields.size()) + " in line "+Utilities::ItoStr(line_counter));

			string ID = line_fields.at(0);
			IndividualIDs.push_back(ID);

			for(unsigned int i = 1; i < line_fields.size(); i++)
			{
				string locusID = locus_IDs.at(i-1);
				string genotype = line_fields.at(i);
				vector<string> genotypes = Utilities::split(genotype, '/');
				if(genotypes.size() != 2)
				{
					errEx("Wrong field length - genotype specified as "+genotype);
				}
				unsigned char thisG1 = CODE.doCode(locusID, genotypes.at(0));
				unsigned char thisG2 = CODE.doCode(locusID, genotypes.at(1));
				if((thisG1 < thisG2) || (preserveNAs && (ID.substr(0, 2) == "NA")))
				{
					GenotypesByLoci[locusID].push_back(thisG1);
					GenotypesByLoci[locusID].push_back(thisG2);
				}
				else
				{
					GenotypesByLoci[locusID].push_back(thisG2);
					GenotypesByLoci[locusID].push_back(thisG1);					
				}
			}
		}
		genotypes.close();
	}
	else
	{
		errEx("Could not open file: "+filename);
	}

}

void GenotypePanel::writeToFile(string filename)
{
	ofstream genotypes;
	genotypes.open (filename.c_str(), ios::out | ios::trunc);
	string sep = " ";
	if (genotypes.is_open())
	{
		vector<string> fieldNames = getLoci();

		genotypes << "IndividualID" << sep << Utilities::join(fieldNames, " ") << "\n";

		int genotypes_counter = -1;
		for(vector<string>::iterator IndividualID = IndividualIDs.begin(); IndividualID != IndividualIDs.end(); IndividualID++)
		{
			genotypes_counter++;

			vector<string> fields_to_print;
			fields_to_print.push_back(*IndividualID);
			for(vector<string>::iterator fieldName = fieldNames.begin(); fieldName != fieldNames.end(); fieldName++)
			{
				unsigned char value1 = GenotypesByLoci[*fieldName][2*genotypes_counter];
				unsigned char value2 = GenotypesByLoci[*fieldName][2*genotypes_counter+1];

				string decoded1 = CODE.deCode(*fieldName, value1);
				string decoded2 = CODE.deCode(*fieldName, value2);

				fields_to_print.push_back(decoded1+"/"+decoded2);

			}

			genotypes << Utilities::join(fields_to_print, " ") << "\n";
		}

		genotypes.close();
	}
	else
	{
		errEx("Could not open file: "+filename);
	}
}

vector<string> GenotypePanel::getLoci()
{
	vector<string> fieldNames;
	for(map<string, basic_string<unsigned char> >::iterator iter = GenotypesByLoci.begin(); iter != GenotypesByLoci.end(); ++iter)
	{
		fieldNames.push_back(iter->first);
	}
	return fieldNames;
}

HaplotypePanel GenotypePanel::sampleRandomHaplotypes(int samples)
{
	vector<string> haplotypeIDs;
	map<string, basic_string<unsigned char> > HaplotypesByLoci;

	//TODO activate at later point
	// srand ( time(NULL) );

	vector<string> loci_in_order = getLoci();

	int individualID_counter = -1;
	for(vector<string>::iterator IndividualID = IndividualIDs.begin(); IndividualID != IndividualIDs.end(); IndividualID++)
	{
		individualID_counter++;

		for(int sampleN = 1; sampleN <= samples; sampleN++)
		{
			string ID1 = *IndividualID+haplotype_number_separator+Utilities::ItoStr(sampleN)+haplotype_field_separator+"1";
			string ID2 = *IndividualID+haplotype_number_separator+Utilities::ItoStr(sampleN)+haplotype_field_separator+"2";

			haplotypeIDs.push_back(ID1);
			haplotypeIDs.push_back(ID2);

			for(vector<string>::iterator locusID = loci_in_order.begin(); locusID != loci_in_order.end(); locusID++)
			{
				unsigned char v1 = GenotypesByLoci[*locusID][2*individualID_counter];
				unsigned char v2 = GenotypesByLoci[*locusID][2*individualID_counter+1];

				if(v1 == v2)
				{
					HaplotypesByLoci[*locusID].push_back(v1);
					HaplotypesByLoci[*locusID].push_back(v2);
				}
				else
				{
					if(preserveNAs && (IndividualID->substr(0, 2) == "NA"))
					{
						HaplotypesByLoci[*locusID].push_back(v1);
						HaplotypesByLoci[*locusID].push_back(v2);						
					}
					else
					{
						float r = (float)rand()/(float)RAND_MAX;
						if(r <= 0.5)
						{
							HaplotypesByLoci[*locusID].push_back(v1);
							HaplotypesByLoci[*locusID].push_back(v2);
						}
						else
						{
							HaplotypesByLoci[*locusID].push_back(v2);
							HaplotypesByLoci[*locusID].push_back(v1);
						}						
					}

				}
			}
		}
	}

	HaplotypePanel hp(CODE, haplotypeIDs, HaplotypesByLoci, LocusPositions);
	return hp;
}

