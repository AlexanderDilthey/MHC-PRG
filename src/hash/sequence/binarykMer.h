/*
 * binarykMer.h
 *
 *  Created on: 29.11.2012
 *      Author: AlexanderDilthey
 */

#ifndef BINARYKMER_H_
#define BINARYKMER_H_

#include <inttypes.h>
#include <array>
#include <vector>
#include <string>
#include <exception>
#include <stdexcept>
#include "basic.h"
#include <assert.h>

extern uint32_t hashlittle( const void *key, size_t length, uint32_t initval);

typedef uint64_t bitfield_of_64bits;

template<int bitFieldMultiplicity, int kMerSize>
class binarykMer
{
public:
	std::array<bitfield_of_64bits, bitFieldMultiplicity> kMerBinaryRepresentation;

	binarykMer()
	{
		if(bitFieldMultiplicity == 1)
		{
			assert((kMerSize >= 1) && (kMerSize <= 31));
		}
		else if(bitFieldMultiplicity == 2)
		{
			assert((kMerSize >= 32) && (kMerSize <= 63));
		}
		else
		{
			throw new std::runtime_error("Unsupported bitFieldMultiplicity (you just need to extend the check range to fix this!)");
		}

		for(int i = 0; i < bitFieldMultiplicity; i++)
		{
			kMerBinaryRepresentation.at(i) = 0;
		}
	}

	binarykMer(std::array<bitfield_of_64bits, bitFieldMultiplicity>& binaryTemplate)
	{
		for(int i = 0; i < bitFieldMultiplicity; i++)
		{
			kMerBinaryRepresentation.at(i) = binaryTemplate.at(i);
		}
	}

	binarykMer(const binarykMer<bitFieldMultiplicity, kMerSize>& that)
	{
		for(int i = 0; i < bitFieldMultiplicity; i++)
		{
			kMerBinaryRepresentation[i] = that.kMerBinaryRepresentation[i];
		}
	}

	friend void swap (binarykMer<bitFieldMultiplicity, kMerSize>& first, binarykMer<bitFieldMultiplicity, kMerSize>& second)
	{
	    using std::swap;
	    swap(first.kMerBinaryRepresentation, second.kMerBinaryRepresentation);
	}

	binarykMer<bitFieldMultiplicity, kMerSize>& operator=(binarykMer<bitFieldMultiplicity, kMerSize> other)
	{
	    swap(*this, other);
	    return *this;
	}


	void loadSeq(std::string seq);
	std::string getSeq();

	void bitOps_right_shift_one_base();
	void bitOps_left_shift_one_base();
	void bitOps_left_shift_one_base_and_insert_new_base_at_right_end(Nucleotide n);

	bitfield_of_64bits& at(int i);

	friend bool operator==(const binarykMer<bitFieldMultiplicity, kMerSize>& kMer1, const binarykMer<bitFieldMultiplicity, kMerSize>& kMer2)
	{
	    bool identical = true;
	    for (int i=0; i < bitFieldMultiplicity; i++)
	    {
	        if (kMer1.kMerBinaryRepresentation.at(i) != kMer2.kMerBinaryRepresentation[i])
	        {
	        	identical = false;
	            break;
	        }
	    }
	    return identical;
	}

	friend bool operator<=(const binarykMer<bitFieldMultiplicity, kMerSize>& left, const binarykMer<bitFieldMultiplicity, kMerSize>& right)
	{
	    bool left_is_less_than_right = false;

	    //need the following to work out which bits to ignore
	    int number_of_bitfields_fully_used = kMerSize/32;

	    int i;

	    //start at most significant end
	    // this would break if we had number_of_bitfields_fully_used==NUMBER_OF_BITFIELDS_IN_BINARY_KMER. But we can never have that as k is always odd.
	    for (i=bitFieldMultiplicity-number_of_bitfields_fully_used-1; i<bitFieldMultiplicity ; i++)
	    {
	        if (left.kMerBinaryRepresentation.at(i) < right.kMerBinaryRepresentation.at(i))
	        {
	            left_is_less_than_right=true;
	            break;
	        }
	        else if (left.kMerBinaryRepresentation.at(i) > right.kMerBinaryRepresentation.at(i))
	        {
	            left_is_less_than_right=false;
	            break;
	        }
	    }

	    return left_is_less_than_right;
	}

	friend bool operator!=(const binarykMer<bitFieldMultiplicity, kMerSize>& kMer1, const binarykMer<bitFieldMultiplicity, kMerSize>& kMer2)
	{
		return (! operator==(kMer1, kMer2));
	}

	Nucleotide binary_kmer_get_first_nucleotide()
	{
	    int number_of_bits_in_most_sig_bitfield = 2 * (kMerSize % 32);
	    return (Nucleotide)((kMerBinaryRepresentation.at(0) >> (number_of_bits_in_most_sig_bitfield - 2)) & 0x3);
	}

	Nucleotide binary_kmer_get_last_nucleotide()
	{
		return (Nucleotide)(kMerBinaryRepresentation.at(bitFieldMultiplicity-1) & 0x3);
	}


	size_t hash_value(int number_buckets)
	{
		size_t hashval = hashlittle(kMerBinaryRepresentation.data(), bitFieldMultiplicity * sizeof(bitfield_of_64bits), 10);
	    return (hashval & (number_buckets-1));
	}

	void reverse_complement()
	{
		std::array<bitfield_of_64bits, bitFieldMultiplicity> kMerCopy = kMerBinaryRepresentation;

	    // Complement
	    int i;
	    for(i = 0; i < bitFieldMultiplicity; i++)
	    {
	    	kMerCopy.at(i) = ~kMerBinaryRepresentation[i];
	    }

	    // Reverse
	    // Loop over full bitfields
	    int j, k;
	    for(i = bitFieldMultiplicity-1; i > 0; i--)
	    {
	        for(j = 0; j < 32; j++)
	        {
	            // Shift destination left
	            for(k = 0; k < bitFieldMultiplicity-1; k++)
	            {
	                kMerBinaryRepresentation.at(k) <<= 2;
	                kMerBinaryRepresentation.at(k) |= (kMerBinaryRepresentation.at(k+1) >> 62);
	            }

	            kMerBinaryRepresentation.at(bitFieldMultiplicity-1) <<= 2;

	            // Append new base
	            kMerBinaryRepresentation.at(bitFieldMultiplicity-1) |= kMerCopy.at(i) & 0x3;

	            // Shift source right
	            kMerCopy.at(i) >>= 2;
	        }
	    }

	    // Do remaining bases in last bitfield [0]
	    int top_bases = kMerSize % 32;

	    for(i = 0; i < top_bases; i++)
	    {
	        // Shift destination left
	        for(k = 0; k < bitFieldMultiplicity-1; k++)
	        {
	        	kMerBinaryRepresentation.at(k) <<= 2;
	        	kMerBinaryRepresentation.at(k) |= (kMerBinaryRepresentation.at(k+1) >> 62);
	        }

	        kMerBinaryRepresentation.at(bitFieldMultiplicity-1) <<= 2;

	        // Append new base
	        kMerBinaryRepresentation.at(bitFieldMultiplicity-1) |= kMerCopy.at(0) & 0x3;

	        // Shift source right
	        kMerCopy.at(0) >>= 2;
	    }

	    // Mask top word (we didn't zero at the beginning!)
	    short top_bits = 2 * top_bases; // bits in top word
	    kMerBinaryRepresentation.at(0) &= (~(uint64_t)0 >> (64 - top_bits));
	}

	binarykMer<bitFieldMultiplicity, kMerSize> element_get_key(bool& was_inverted)
	{
	    // Get first and last nucleotides
	    Nucleotide first = binary_kmer_get_first_nucleotide();
	    Nucleotide last  = binary_kmer_get_last_nucleotide();
	    Nucleotide rev_last = (Nucleotide)((int)~last & 0x3);

	    binarykMer<bitFieldMultiplicity, kMerSize> forReturn(*this);

	    if(first < rev_last)
	    {
	    	// nothing
	    	was_inverted = false;
	    }
	    else if(first > rev_last)
	    {
	        forReturn.reverse_complement();
	        was_inverted = true;
	    }
	    else
	    {
	        // Don't know which is going to be correct
	        // This will happen 1 in 4 times
	        forReturn.reverse_complement();

	        if(*this <= forReturn)
	        {
	        	forReturn =  binarykMer<bitFieldMultiplicity, kMerSize>(*this);
	        	was_inverted = false;
	        }
	        else
	        {
	        	was_inverted = true;
	        }
	    }

	    return forReturn;
	}

	binarykMer<bitFieldMultiplicity, kMerSize> element_get_key()
	{
		bool ignore;
	    return element_get_key(ignore);
	}
};




template<int bitFieldMultiplicity, int kMerSize>
void binarykMer<bitFieldMultiplicity, kMerSize>::loadSeq(std::string seq)
{
	if((int)seq.size() != kMerSize)
	{
		throw std::runtime_error("Length of seq != kMerSize");
	}

    for(int j=0; j < kMerSize; j++)
    {
    	bitOps_left_shift_one_base_and_insert_new_base_at_right_end(char_to_binary_nucleotide(seq[j]));
    }
}

template<int bitFieldMultiplicity, int kMerSize>
std::string binarykMer<bitFieldMultiplicity, kMerSize>::getSeq()
{
	binarykMer<bitFieldMultiplicity, kMerSize> localkMer (*this);

	std::string forReturn;
	forReturn.resize(kMerSize);

	int mask = 3; // 0000011 mask used to extract the two least significative bits

	for(int j=kMerSize-1; j>=0; j--)  //start from the back of the sequence
	{
		//get translation for the two least significant bits
		forReturn.at(j) =  binary_nucleotide_to_char((Nucleotide)(localkMer.at(bitFieldMultiplicity-1) & mask));
		localkMer.bitOps_right_shift_one_base(); //note this is a local copy internal to this function - not altering the original BinaryKmer
	}

	return forReturn;
}

template<int bitFieldMultiplicity, int kMerSize>
void binarykMer<bitFieldMultiplicity, kMerSize>::bitOps_right_shift_one_base()
{
	for(int i = bitFieldMultiplicity-1; i > 0; i--)
	{
		kMerBinaryRepresentation.at(i) >>= 2;
		kMerBinaryRepresentation.at(i) |= (kMerBinaryRepresentation.at(i-1) << 62); // & 0x3
	}
	kMerBinaryRepresentation.at(0) >>= 2;
}

template<int bitFieldMultiplicity, int kMerSize>
void binarykMer<bitFieldMultiplicity, kMerSize>::bitOps_left_shift_one_base()
{
    int top_word = bitFieldMultiplicity - (kMerSize+31)/32;

    for(int i = top_word; i < bitFieldMultiplicity-1; i++)
    {
    	kMerBinaryRepresentation.at(i) <<= 2;
        kMerBinaryRepresentation.at(i) |= (kMerBinaryRepresentation.at(i+1) >> 62); // & 0x3
    }

    kMerBinaryRepresentation.at(bitFieldMultiplicity-1) <<= 2;

    // Mask top word
    short top_bits = 2 * (kMerSize % 32); // bits in top word
    kMerBinaryRepresentation.at(top_word) &= (~ (bitfield_of_64bits)0 >> (64 - top_bits));
}

template<int bitFieldMultiplicity, int kMerSize>
void binarykMer<bitFieldMultiplicity, kMerSize>::bitOps_left_shift_one_base_and_insert_new_base_at_right_end(Nucleotide n)
{
	bitOps_left_shift_one_base();
	kMerBinaryRepresentation.at(bitFieldMultiplicity-1) |= n;
}

template<int bitFieldMultiplicity, int kMerSize>
bitfield_of_64bits& binarykMer<bitFieldMultiplicity, kMerSize>::at(int i)
{
	if(i < 0)
	{
		throw std::runtime_error("Indices < 0 not supported!");
	}
	if(i >= bitFieldMultiplicity)
	{
		throw std::runtime_error("Requested index >= bitFieldMultiplicity");
	}
	return kMerBinaryRepresentation.at(i);
}


/*

template<int bitFieldMultiplicity>
class binarykMerHandler {
protected:
	int kMerSize;
public:
	binarykMerHandler(int k);
	void binary_kmer_initialise_to_zero(typename binarykMer<bitFieldMultiplicity>::type& kMer);
	void binary_kmer_assignment_operator(typename binarykMer<bitFieldMultiplicity>::type& kMerTarget, typename binarykMer<bitFieldMultiplicity>::type& kMerSource);
	void binary_kmer_comparison_operator(typename binarykMer<bitFieldMultiplicity>::type& kMer1, typename binarykMer<bitFieldMultiplicity>::type& kMer2);
	void seq_to_binary_kmer(std::string seq, typename binarykMer<bitFieldMultiplicity>::type& forReturn);
	std::string binary_kmer_to_seq(typename binarykMer<bitFieldMultiplicity>::type& kMer);
	void binary_kmer_right_shift_one_base(typename binarykMer<bitFieldMultiplicity>::type& kMer);
	void binary_kmer_left_shift_one_base(typename binarykMer<bitFieldMultiplicity>::type& kMer);
	void binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(typename binarykMer<bitFieldMultiplicity>::type& kMer, Nucleotide n);

};

template<int bitFieldMultiplicity>
binarykMerHandler<bitFieldMultiplicity>::binarykMerHandler(int k)
{
	if(bitFieldMultiplicity == 1)
	{
		if(! ((k >= 1) && (k <= 31)))
		{
			throw std::runtime_error("Wrong kMer size for bitFieldMultiplicity == 1");
		}
	}
	else if(bitFieldMultiplicity == 1)
	{
		if(! ((k >= 32) && (k <= 63)))
		{
			throw std::runtime_error("Wrong kMer size for bitFieldMultiplicity == 2");
		}
	}
	else
	{
		throw std::runtime_error("Cannot deal with bitFieldMultiplicity");
	}

	kMerSize = k;
}


template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_initialise_to_zero(typename binarykMer<bitFieldMultiplicity>::type& kMer)
{
    for (int i=0; i< bitFieldMultiplicity; i++)
    {
    	kMer[i] = (bitfield_of_64bits) 0;
    }
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_assignment_operator(typename binarykMer<bitFieldMultiplicity>::type& kMerTarget, typename binarykMer<bitFieldMultiplicity>::type& kMerSource)
{
    for (int i=0; i < bitFieldMultiplicity; i++)
    {
        kMerTarget[i]=kMerSource[i];
    }
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_comparison_operator(typename binarykMer<bitFieldMultiplicity>::type& kMer1, typename binarykMer<bitFieldMultiplicity>::type& kMer2)
{
    bool identical = true;
    for (int i=0; i < bitFieldMultiplicity; i++)
    {
        if (kMer1[i]!=kMer2[i])
        {
        	identical = false;
            break;
        }
    }
    return identical;
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::seq_to_binary_kmer(std::string seq, typename binarykMer<bitFieldMultiplicity>::type& forReturn)
{
	if((int)seq.size() != kMerSize)
	{
		throw std::runtime_error("Length of seq != kMerSize");
	}

    binary_kmer_initialise_to_zero(forReturn);

    for(int j=0; j < kMerSize; j++)
    {
        binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(forReturn, char_to_binary_nucleotide(seq[j]));
    }
}

template<int bitFieldMultiplicity>
std::string binarykMerHandler<bitFieldMultiplicity>::binary_kmer_to_seq(typename binarykMer<bitFieldMultiplicity>::type& kMer)
{
	typename binarykMer<bitFieldMultiplicity>::type local_bkmer;
	binary_kmer_assignment_operator(local_bkmer, kMer);

	std::string forReturn;
	forReturn.resize(kMerSize);

	int mask = 3; // 0000011 mask used to extract the two least significative bits

	for(int j=kMerSize-1; j>=0; j--)  //start from the back of the sequence
	{
		//get translation for the two least significant bits
		forReturn.at(j) =  binary_nucleotide_to_char((Nucleotide)(local_bkmer[bitFieldMultiplicity-1] & mask));
		binary_kmer_right_shift_one_base(local_bkmer); //note this is a local copy internal to this function - not altering the original BinaryKmer
	}

	return forReturn;
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_right_shift_one_base(typename binarykMer<bitFieldMultiplicity>::type& kMer)
{
	for(int i = bitFieldMultiplicity-1; i > 0; i--)
	{
		kMer[i] >>= 2;
		kMer[i] |= (kMer[i-1] << 62); // & 0x3
	}
	kMer[0] >>= 2;
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_left_shift_one_base(typename binarykMer<bitFieldMultiplicity>::type& kMer)
{
    int top_word = bitFieldMultiplicity - (kMerSize+31)/32;

    for(int i = top_word; i < bitFieldMultiplicity-1; i++)
    {
        kMer[i] <<= 2;
        kMer[i] |= (kMer[i+1] >> 62); // & 0x3
    }

    kMer[bitFieldMultiplicity-1] <<= 2;

    // Mask top word
    short top_bits = 2 * (kMerSize % 32); // bits in top word
    kMer[top_word] &= (~ (bitfield_of_64bits)0 >> (64 - top_bits));
}

template<int bitFieldMultiplicity>
void binarykMerHandler<bitFieldMultiplicity>::binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(typename binarykMer<bitFieldMultiplicity>::type& kMer, Nucleotide n)
{
    binary_kmer_left_shift_one_base(kMer);
    kMer[bitFieldMultiplicity-1] |= n;
}

*/

#endif /* BINARYKMER_H_ */
