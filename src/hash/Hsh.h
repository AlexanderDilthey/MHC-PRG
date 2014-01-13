/*
 * Hsh.h
 *
 *  Created on: 29.11.2012
 *      Author: AlexanderDilthey
 *
 *  Interface and some code adapted from / inspired by Zamin Iqbal / cortex_var.
 *
 */

#ifndef HSH_H_
#define HSH_H_

#include <cstdlib>
#include <string>
#include <exception>
#include <stdexcept>
#include <vector>
#include <assert.h>
#include <iostream>

#include "sequence/binarykMer.h"

uint32_t hashlittle( const void *key, size_t length, uint32_t initval);

struct flags_struct {
	bool assigned;
	bool marked;
};

template<int m, int k, class ElementType >
class hashTableElement
{
public:
	binarykMer<m, k> kMer;
	flags_struct flags;
	ElementType element;

	hashTableElement()
	{
		flags.assigned = false;
		flags.marked = false;
	}
};

template<int m, int k, class ElementType>
class Hsh {
	int kmer_size;
	size_t number_buckets;
	int bucket_size;
	std::vector<short> next_element; //keeps index of the next free element in bucket
	std::vector<long long> collisions;
	size_t unique_kmers;
	int max_rehash_tries;

	std::vector< hashTableElement<m, k, ElementType> > table;

public:
	Hsh(int log2_buckets, long long onebucket_size);

	ElementType* hash_table_find_or_insert(binarykMer<m, k>& key, bool& found);
	ElementType* hash_table_insert(binarykMer<m, k>& key);
	ElementType* hash_table_find(binarykMer<m, k>& key);
	hashTableElement<m, k, ElementType >* hash_table_find_fullElement(binarykMer<m, k>& key);

	bool hash_table_find_in_bucket(binarykMer<m, k>& key, size_t& current_pos, bool& overflow, int rehash);

	size_t getTableSize()
	{
		return table.size();
	}

	hashTableElement<m, k, ElementType>* getTableElement(size_t i)
	{
		return &(table.at(i));
	}

	hashTableElement<m, k, ElementType>* quicklyGetRandomElement()
	{
		unsigned int maxI = table.size();
		unsigned int randomIStart = rand() % maxI;
		assert(randomIStart >= 0);
		assert(randomIStart < maxI);
		bool abort = false;
		unsigned int currentI = randomIStart;
		while(! abort)
		{
			if(table.at(currentI).flags.assigned == true)
			{
				return &(table.at(currentI));
			}
			else
			{
				currentI++;
				if(currentI >= maxI)
				{
					currentI = 0;
				}
				if(currentI == randomIStart)
				{
					abort = true;
				}
			}
		}

		return 0;
	}
};

template<int m, int k, class ElementType>
Hsh<m, k, ElementType>::Hsh(int log2_buckets, long long onebucket_size)
{
	kmer_size = k;
	number_buckets = (size_t) 1 << log2_buckets;
	bucket_size = onebucket_size;
	max_rehash_tries = 100;

	collisions.resize(max_rehash_tries);
	unique_kmers = 0;

	table.resize(number_buckets * bucket_size);
	next_element.resize(number_buckets, 0);
}

template<int m, int k, class ElementType>
hashTableElement<m, k, ElementType >* Hsh<m, k, ElementType>::hash_table_find_fullElement(binarykMer<m, k>& key)
{
    int rehash = 0;
    bool overflow;
    size_t current_pos;
    bool found;

    hashTableElement<m, k, ElementType >* forReturn = 0;

	do
	{
		found = hash_table_find_in_bucket(key, current_pos, overflow, rehash);
		if (found) //then we know overflow is false - this is checked in find_in_bucket
		{
			forReturn =  &(table.at(current_pos));
		}
		else if (overflow)
		{
			rehash++;
			if (rehash > max_rehash_tries)
			{
				//fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
				throw std::runtime_error("Dear user - you have not allocated enough memory to contain your sequence data. Either allocate more memory (have you done your calculations right? have you allowed for sequencing errors?), or threshold more harshly on quality score, and try again. Aborting mission.\n");
			}
		}
	}
	while(overflow);

    return forReturn;
}

template<int m, int k, class ElementType>
ElementType* Hsh<m, k, ElementType>::hash_table_find(binarykMer<m, k>& key)
{
    int rehash = 0;
    bool overflow;
    size_t current_pos;
    bool found;

    ElementType* forReturn = 0;

	do
	{
		found = hash_table_find_in_bucket(key, current_pos, overflow, rehash);
		if (found) //then we know overflow is false - this is checked in find_in_bucket
		{
			forReturn =  &(table.at(current_pos).element);
		}
		else if (overflow)
		{
			rehash++;
			if (rehash > max_rehash_tries)
			{
				//fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
				throw std::runtime_error("Dear user - you have not allocated enough memory to contain your sequence data. Either allocate more memory (have you done your calculations right? have you allowed for sequencing errors?), or threshold more harshly on quality score, and try again. Aborting mission.\n");
			}
		}
	}
	while(overflow);

    return forReturn;
}

template<int m, int k, class ElementType>
ElementType* Hsh<m, k, ElementType>::hash_table_find_or_insert(binarykMer<m, k>& key, bool& found)
{
    int rehash = 0;
    bool overflow;

    size_t current_pos;

    do
    {

        found = hash_table_find_in_bucket(key, current_pos, overflow, rehash);


        if (! found)
        {
            if (! overflow) //it is definitely nowhere in the hashtable, so free to insert
            {
                //sanity check
            	assert(table.at(current_pos).flags.assigned == false);

            	table.at(current_pos).kMer = key;
            	table.at(current_pos).flags.assigned = true;

                unique_kmers++;

            }
            else
            {
                //overflow -> rehashing

                rehash++;
                if (rehash > max_rehash_tries)
                {
                    throw std::runtime_error("Dear user - you have not allocated enough memory to contain your sequence data. Either allocate more memory (have you done your calculations right? have you allowed for sequencing errors?), or threshold more harshly on quality score, and try again. Aborting mission.\n");
                }
            }
        }
        else //it is found
        {

        }
    }
    while (overflow);


    collisions[rehash]++;

    ElementType* forReturn = &(table.at(current_pos).element);

    // std::cout << "\t\t\t" << "found key " << key.getSeq() <<  " at address " << forReturn << "\n" << std::flush;

    return forReturn;
}

template<int m, int k, class ElementType>
ElementType* Hsh<m, k, ElementType>::hash_table_insert(binarykMer<m, k>& key)
{
	ElementType* forReturn;
	int rehash = 0;
	bool inserted = false;

    do
    {
        //add the rehash to the final bitfield in the BinaryKmer

    	binarykMer<m, k> kMer_with_rehash = key;
    	kMer_with_rehash.at(m - 1) =	kMer_with_rehash.at(m - 1) + (bitfield_of_64bits) rehash;

        size_t hashval = kMer_with_rehash.hash_value(number_buckets);

        if (next_element.at(hashval) < bucket_size)
        {
            //can insert element
            size_t  current_pos   = (size_t) hashval * bucket_size + (size_t) next_element.at(hashval);   //position in hash table
            assert(table.at(current_pos).flags.assigned == false);
            table.at(current_pos).kMer = key;
            table.at(current_pos).flags.assigned = true;
            unique_kmers++;
            next_element.at(hashval)++;
            forReturn = &(table.at(current_pos).element);
            inserted=true;
        }
        else
        {
            //rehash
            rehash++;
            if (rehash > max_rehash_tries)
            {
                //fprintf(stderr,"too much rehashing!! Reserve more memory.  Rehash=%d\n", rehash);
                throw std::runtime_error("Dear user - you have not allocated enough memory to contain your sequence data. Either allocate more memory (have you done your calculations right? have you allowed for sequencing errors?), or threshold more harshly on quality score, and try again. Aborting mission.\n");
            }
        }

    }
    while (! inserted);

	return forReturn;
}

// Lookup for key in bucket defined by the hash value.
// If key is in bucket, returns true and the position of the key/element in current_pos.
// If key is not in bucket, and bucket is not full, returns the next available position in current_pos (and overflow is returned as false)
// If key is not in bucket, and bucket is full, returns overflow=true
template<int m, int k, class ElementType>
bool Hsh<m, k, ElementType>::hash_table_find_in_bucket(binarykMer<m, k>& key, size_t& current_pos, bool& overflow, int rehash)
{

	binarykMer<m, k> kMer_with_rehash = key;
	kMer_with_rehash.at(m - 1) = (kMer_with_rehash.at(m - 1) + (bitfield_of_64bits) rehash);

    size_t hashval = kMer_with_rehash.hash_value(number_buckets);

	bool found = false;
	overflow = false;

	size_t i = 0; //position in bucket
	size_t bucket_start_pos = (size_t) hashval * bucket_size;

    while(
    		(i<bucket_size) &&   // still within the bucket
            (table.at(bucket_start_pos+i).flags.assigned == true)  && // not yet reached an empty space
            (!found)
         )
    {
		if (key == table.at(bucket_start_pos+i).kMer)
		{
			found = true;
		}
		else
		{
			i++;
		}
    }

    current_pos = bucket_start_pos + i;
    if (i == bucket_size)
    {
        overflow = true;
    }

    assert(! (found && overflow));

    // I do not understand this - if found, we cannot be overflown and vice versa. But if not found,
    // we can either be overflown or not. Have therefore substituted with above.
    // assert( (! found) || (! overflow));

    return found;
}




#endif /* HSH_H_ */
