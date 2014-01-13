/*
 * Utilities.cpp
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#include "Utilities.h"

#include <iostream>
#include <sstream>
#include "MHC-PRG.h"
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <exception>
#include <stdexcept>

// todo at some point activate...
// unsigned int globalRandRSeed = time(NULL);
unsigned int globalRandRSeed = 0;

Utilities::Utilities() {
}

Utilities::~Utilities() {
}

std::pair<double, int> Utilities::findIntMapMax(std::map<int, double>& m)
{
	assert(m.size() > 0);
	double max;
	unsigned iMax;
	for(std::map<int, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		if((mIt == m.begin()) || (mIt->second > max))
		{
			max = mIt->second;
			iMax = mIt->first;
		}
	}

	return std::pair<double, unsigned int>(max, iMax);
}


std::pair<double, int> Utilities::findIntMapMaxP_nonCritical(std::map<int, double>& m, unsigned int* thisSeed)
{
	assert(m.size() > 0);
	double max;
	std::vector<unsigned int> iMax;
	for(std::map<int, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		if((mIt == m.begin()) || (mIt->second > max))
		{
			max = mIt->second;
			iMax.clear();
			iMax.push_back(mIt->first);
		}
	}

	if(iMax.size() > 1)
	{
		int selectedI = randomNumber_nonCritical(iMax.size() - 1, thisSeed);
		assert((selectedI >= 0) && (selectedI < iMax.size()));
		return std::pair<double, unsigned int>(max, iMax.at(selectedI));
	}
	else
	{
		return std::pair<double, unsigned int>(max, iMax.at(0));
	}
}

bool Utilities::extractBit(unsigned int number, unsigned int bit)
{
	unsigned int bitmask = pow((unsigned int)2, bit);
	bool bitValue =  ((number & bitmask) == bitmask);
	return bitValue;
}

std::pair<double, unsigned int> Utilities::findVectorMax(std::vector<double>& v)
{
	assert(v.size() > 0);
	double max;
	unsigned iMax;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMax = i;
			max = v.at(i);
		}
	}

	return std::pair<double, unsigned int>(max, iMax);
}


std::pair<double, unsigned int> Utilities::findVectorMaxP(std::vector<double>& v)
{
	assert(v.size() > 0);
	double max;
	// unsigned iMax;
	std::vector<unsigned int> iMaxs;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMaxs.clear();
			iMaxs.push_back(i);
			max = v.at(i);
		}
	}

	unsigned int selectedI;
	if(iMaxs.size() == 1)
	{
		selectedI = iMaxs.at(0);
	}
	else
	{
		selectedI = iMaxs.at(randomNumber(iMaxs.size() - 1));
	}

	return std::pair<double, unsigned int>(max, selectedI);
}

std::pair<double, unsigned int> Utilities::findVectorMaxP_nonCritical(std::vector<double>& v, unsigned int* thisSeed)
{
	assert(v.size() > 0);
	double max;
	// unsigned iMax;
	std::vector<unsigned int> iMaxs;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMaxs.clear();
			iMaxs.push_back(i);
			max = v.at(i);
		}
	}

	unsigned int selectedI;
	if(iMaxs.size() == 1)
	{
		selectedI = iMaxs.at(0);
	}
	else
	{
		selectedI = iMaxs.at(randomNumber_nonCritical(iMaxs.size() - 1, thisSeed));
	}

	return std::pair<double, unsigned int>(max, selectedI);
}

vector<string> Utilities::ItoStr(vector<int> i)
{
	vector<string> forReturn;
	for(int I = 0; I < (int)i.size(); I++)
	{
		forReturn.push_back(ItoStr(i.at(I)));
	}
	return forReturn;
}

double Utilities::randomDouble()
{
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	return f;
}


int Utilities::chooseFromVector(vector<double>& v)
{
	//TODO activate at later point
	// srand ( time(NULL) );
	
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	double resolution = 1 / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	
	double sum = 0.0;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		sum += v.at(i);
	}
	assert(sum != 0);
	
	int problematic_items = 0;
	double problematic_proportion = 0.0;
	double running_sum = 0.0;
	int selected_i = -1;
	
	assert(v.size() > 0);
	for(unsigned int i = 0; i < v.size(); i++)
	{
		double proportion = v.at(i)/sum;
		running_sum += proportion;
		if(proportion < resolution)
		{
			problematic_items++;
			problematic_proportion += proportion;
		}
		
		if(f <= running_sum)
		{
			selected_i = i;
			break;
		}
		if(i == (v.size() - 1))
		{
			selected_i = i;
			cout << "WARNING!!\n";
			cout << "selected_i: " << selected_i << "\n";
			cout << "f: " << f << "\n";
			cout << "running_sum: " << running_sum << "\n";
			cout << "resolution: " << resolution << "\n";
			cout << "problematic_items: " << problematic_items << "\n";
			cout << "sum: " << sum << "\n";

			double deb_running_sum = 0.00;
			for(unsigned int k = 0; k < v.size(); k++)
			{
				double deb_proportion = v.at(k)/sum;
				cout << "Vector v element " << k << "\n";
				cout << "\tv.at(k): " << v.at(k) << "\n";
				cout << "\tdeb_proportion: "<< deb_proportion << "\n";
				cout << "\tdeb_running_sum before: "<< deb_running_sum << "\n";

				deb_running_sum += deb_proportion;

				cout << "\tdeb_running_sum after: "<< deb_running_sum << "\n";


			}

		}
	}
	
	assert(selected_i != -1);
	assert(problematic_proportion < 0.25);
	
	return selected_i;
}

vector<string> Utilities::split(const string &s, char delim, vector<string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

vector<string> Utilities::split(string input, string delimiter)
{
	vector<string> output;
	if(input.length() == 0)
	{
		return output;
	}

	if(input.find(delimiter) == string::npos)
	{
		output.push_back(input);
	}
	else
	{
		int s = 0;
		int p = input.find(delimiter);

		do {
			output.push_back(input.substr(s, p - s));
			s = p + delimiter.size();
			p = input.find(delimiter, s);
		} while (p != (int)string::npos);
		output.push_back(input.substr(s));
	}

	return output;
}


vector<string> Utilities::split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}

string Utilities::ItoStr(int i)
{
	std::stringstream sstm;
	sstm << i;
	return sstm.str();
}


string Utilities::PtoStr(void* p)
{
	std::stringstream sstm;
	sstm << p;
	return sstm.str();
}


string Utilities::DtoStr(double d)
{
	std::stringstream sstm;
	sstm << d;
	return sstm.str();
}

string Utilities::join(vector<string> parts, string delim)
{
	if(parts.size() == 0)
		return "";

	string ret = parts.at(0);

	for(unsigned int i = 1; i < parts.size(); i++)
	{
		ret.append(delim);
		ret.append(parts.at(i));
	}

	return ret;
}

void Utilities::eraseNL(string& s)
{
	if (!s.empty() && s[s.length()-1] == '\r') {
	    s.erase(s.length()-1);
	}
	if (!s.empty() && s[s.length()-1] == '\n') {
	    s.erase(s.length()-1);
	}
}

long long Utilities::StrtoLongLong(string s)
{
	  stringstream ss(s);
	  long long i;
	  ss >> i;
	  return i;
}

int Utilities::StrtoI(string s)
{
	  stringstream ss(s);
	  int i;
	  ss >> i;
	  return i;
}

double Utilities::StrtoD(string s)
{
	  stringstream ss(s);
	  double d;
	  ss >> d;
	  return d;
}

bool Utilities::StrtoB(string s)
{
	  stringstream ss(s);
	  bool b;
	  ss >> b;
	  return b;
}



std::string Utilities::timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    char* timeCString = asctime( localtime(&ltime) );
    std::string forReturn(timeCString);
    return " [ "+ forReturn.substr(0, forReturn.length() - 1)+" ] ";
}

map<string, string> Utilities::readFASTA(std::string file, bool fullIdentifier)
{
	map<string, string> forReturn;

	ifstream FASTAstream;

	FASTAstream.open(file.c_str());
	if(! FASTAstream.is_open())
	{
		throw std::runtime_error("readFASTA(): Cannot open file "+file);
	}
	
	while(FASTAstream.good())
	{
		std::string line;
		size_t lineCounter = 0;

		std::string currentSequenceIdentifier;
		while(FASTAstream.good())
		{
			std::getline(FASTAstream, line);
			Utilities::eraseNL(line);

			lineCounter++;

			if(line.substr(0, 1) == ">")
			{
				std::string ident = line.substr(1);
				if(! fullIdentifier)
				{
					for(int i = 0; i < ident.size(); i++)
					{
						if(ident.at(i) == ' ')
						{
							ident = ident.substr(0, i);
							break;
						}
					}
				}
				currentSequenceIdentifier = ident;
				assert(forReturn.count(ident) == 0);
			}
			else
			{
				forReturn[currentSequenceIdentifier] += line;
			}
		}
	}

	return forReturn;
}

std::string Utilities::repeatString(std::string s, int repeatNumber)
{
	std::string forReturn;
	assert(repeatNumber >= 0);
	for(int i = 1; i <= repeatNumber; i++)
	{
		forReturn.append(s);
	}
	return forReturn;
}


std::pair<std::string,std::vector<int>> Utilities::modifySequence(std::string sequence, std::vector<int> positionOrigin, double mutationFrequence, double insertionFrequence, int insertionMaxLength, double deletionFrequence, int deletionMaxLength)
{
	assert(sequence.size() == positionOrigin.size());

	assert((mutationFrequence >= 0) && (mutationFrequence <= 1));
	assert((insertionFrequence >= 0) && (insertionFrequence <= 1));
	assert((deletionFrequence >= 0) && (deletionFrequence <= 1));
	assert(insertionMaxLength > 1);
	assert(deletionMaxLength > 1);

	std::pair<std::string,std::vector<int>> forReturn;
	for(int i = 0; i < (int)sequence.length(); i++)
	{
		if(randomDouble() <= mutationFrequence)
		{
			char newNuc = randomNucleotide();
			forReturn.first.push_back(newNuc);
			if(newNuc == sequence.at(i))
			{
				forReturn.second.push_back(positionOrigin.at(i));
			}
			else
			{
				forReturn.second.push_back(-1);
			}
		}
		else if(randomDouble() <= deletionFrequence)
		{
			int deletionLength = randomNumber(deletionMaxLength - 1);
			for(int dI = 0; dI < deletionLength; dI++)
			{
				forReturn.first.push_back('_');
				forReturn.second.push_back(-1);
			}
			i += deletionLength;
		}
		else if(randomDouble() <= insertionFrequence)
		{
			int insertionLength = randomNumber(insertionMaxLength - 1);
			forReturn.first.append(generateRandomSequence(insertionLength));
			for(int dI = 0; dI < insertionLength; dI++)
			{
				forReturn.second.push_back(-1);
			}
		}
		else
		{
			forReturn.first.push_back(sequence.at(i));
			forReturn.second.push_back(positionOrigin.at(i));

		}
	}

	assert(forReturn.first.size() == forReturn.second.size());

	return forReturn;
}

std::string Utilities::generateRandomSequenceWithGaps(int length, double gapFrequency)
{
	std::string forReturn;
	forReturn.resize(length);
	for(int i = 0; i < length; i++)
	{
		if(randomDouble() <= gapFrequency)
		{
			forReturn.at(i) = '_';
		}
		else
		{
			forReturn.at(i) = randomNucleotide();
		}
	}
	return forReturn;
}

std::string Utilities::generateRandomSequence(int length)
{
	std::string forReturn;
	forReturn.resize(length);
	for(int i = 0; i < length; i++)
	{
		forReturn.at(i) = randomNucleotide();
	}
	return forReturn;
}

int Utilities::randomNumber(int max)
{
	int n;
	#pragma omp critical
	{
		// n = rand() % (max + 1);
		n = rand_r(&globalRandRSeed) % (max + 1);
	}
	assert((n >= 0) && (n <= max));
	return n;
}

int Utilities::randomNumber_nonCritical(int max, unsigned int* thisSeed)
{
	int n = rand_r(thisSeed) % (max + 1);
	assert((n >= 0) && (n <= max));
	return n;
}

char Utilities::randomNucleotide()
{
	char nucleotides[4] = {'A', 'C', 'G', 'T'};
	int n = rand() % 4;
	assert((n >= 0) && (n <= 3));
	return nucleotides[n];
}


void Utilities::writeStatus(std::string statusFile, int status)
{
	std::ofstream statusStream;
	statusStream.open(statusFile.c_str(), std::ios::out);
	if(! statusStream.is_open())
	{
		throw std::runtime_error("Cannot open file for writing: "+statusFile);
	}
	statusStream << status << "\n";
	statusStream.close();
}


int Utilities::readStatus(std::string statusFile)
{
	std::ifstream statusStream;
	statusStream.open(statusFile.c_str(), std::ios::in);
	if(! statusStream.is_open())
	{
		return 0;
	}
	assert(statusStream.good());
	std::string line;
	getline (statusStream, line);
	Utilities::eraseNL(line);
	if(line.length() == 0)
	{
		return 0;
	}
	int status = Utilities::StrtoI(line);
	assert(status >= 0);

	return status;
}






