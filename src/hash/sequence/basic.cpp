#include "basic.h"
#include <iostream>
#include <string>
#include <exception>
#include <stdexcept>
#include <assert.h>

char reverse_char_nucleotide(char c)
{
    switch (c)
    {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'N':
			return 'N';			
		case 'a':
			return 't';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		case 't':
			return 'a';
		case 'n':
			return 'n';			
		default:
			throw std::runtime_error("Nucleotide not existing!");
    }
}

std::string seq_reverse_complement(std::string sequence)
{
	int length = sequence.size();
	std::string forReturn;
	forReturn.resize(length);
    for(int k=0; k < length; k++)
    {
        forReturn[k] = reverse_char_nucleotide(sequence.at(length-k-1));
    }
    return forReturn;
}

Nucleotide char_to_binary_nucleotide(char c)
{
    switch (c)
    {
    case 'A':
        return Adenine;
    case 'C':
        return Cytosine;
    case 'G':
        return Guanine;
    case 'T':
        return Thymine;
    case 'a':
        return Adenine;
    case 'c':
        return Cytosine;
    case 'g':
        return Guanine;
    case 't':
        return Thymine;
    default:
        return Undefined;
    }
}

Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
    switch (n)
    {
    case Adenine:
        return Thymine;
    case Cytosine:
        return Guanine;
    case Guanine:
        return Cytosine;
    case Thymine:
        return Adenine;
    default:
        throw std::runtime_error("Calling reverse_binary_nucleotide on non-existent nucleotide");
    }
}

char binary_nucleotide_to_char(Nucleotide n)
{
    switch (n)
    {
    case Adenine:
        return 'A';
    case Cytosine:
        return 'C';
    case Guanine:
        return 'G';
    case Thymine:
        return 'T';
    default:
    	throw std::runtime_error("Non existent binary nucleotide");
    }
}

int find_first_valid_kMer_position(std::string sequence, int k)
{
	return find_first_valid_kMer_position(sequence, k, 0);
}

int find_first_valid_kMer_position(std::string sequence, int k, int startAtPosition)
{
	//std::cout << "Enter find_first_valid_kMer_position: " << sequence.length() << " " << k << " " << startAtPosition << "\n" << std::flush;
	if(sequence.length() < (k+startAtPosition))
	{
		std::cout << "\treturn -1\n" << std::flush;
		return -1;
	}

	int currentStart = startAtPosition;
	do {
		std::string candidate = sequence.substr(currentStart, k);
		assert((int)candidate.length() == k);

		int lastInvalidPosition = -1;
		for(int i = 0; i < k; i++)
		{
			if(char_to_binary_nucleotide(candidate.at(i)) == Undefined)
			{
				lastInvalidPosition = i;
			}

			if(lastInvalidPosition == -1)
			{
				//std::cout << "\treturn " << currentStart << "\n" << std::flush;

				return currentStart;
			}
			else
			{
				currentStart += (lastInvalidPosition + 1);
			}
		}
	} while((currentStart + k) <= sequence.length());

	//std::cout << "\treturn -1\n" << std::flush;
	return -1;
}

std::vector<std::string> partitionStringIntokMers(std::string str, int k)
{
	std::vector<std::string> forReturn;
	if((int)str.length() >= k)
	{
		for(int i = 0; i <= (str.length() - k); i++)
		{
			std::string kMer = str.substr(i, k);
			assert((int)kMer.length() == k);
			forReturn.push_back(kMer);
		}
	}
	return forReturn;
}

