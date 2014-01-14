/*
 * Utilities.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>

using namespace std;

class Utilities {
public:
	Utilities();
	virtual ~Utilities();

	static map<string, string> readFASTA(std::string file, bool fullIdentifier = false);

	static vector<string> split(const string &s, char delim, vector<string> &elems);
	static vector<string> split(const string &s, char delim);
	static string ItoStr(int i);
	static string DtoStr(double d);
	static string PtoStr(void* p);

	static int StrtoI(string s);

	static vector<string> ItoStr(vector<int> i);
	static double StrtoD(string s);
	static bool StrtoB(string s);
	static long long StrtoLongLong(string s);

	static string join(vector<string> parts, string delim);
	static void eraseNL(string& s);
	static int chooseFromVector(vector<double>& v);

	static vector<string> split(string input, string delimiter);


	static std::string timestamp();

	static std::pair<double, unsigned int> findVectorMax(std::vector<double>& v);
	static std::pair<double, int> findIntMapMax(std::map<int, double>& m);
	static std::pair<double, int> findIntMapMaxP_nonCritical(std::map<int, double>& m, unsigned int* thisSeed);
	static std::pair<double, unsigned int> findVectorMaxP(std::vector<double>& v);
	static std::pair<double, unsigned int> findVectorMaxP_nonCritical(std::vector<double>& v, unsigned int* thisSeed);

	static std::string generateRandomSequence(int length);
	static std::string generateRandomSequenceWithGaps(int length, double gapFrequency = 0.1);
	static std::pair<std::string,std::vector<int>> modifySequence(std::string sequence, std::vector<int> positionOrigin, double mutationFrequence = 0.1, double insertionFrequence = 0.05, int insertionMaxLength = 4, double deletionFrequence = 0.1, int deletionMaxLength = 5);

	static char randomNucleotide();
	static int randomNumber(int max);
	static int randomNumber_nonCritical(int max, unsigned int* thisSeed);

	static double randomDouble();

	static std::string repeatString(std::string s, int repeatNumber);

	static bool extractBit(unsigned int number, unsigned int bit);

	static int readStatus(std::string statusFile);
	static void writeStatus(std::string statusFile, int status);


};

extern unsigned int globalRandRSeed;
#endif /* UTILITIES_H_ */

