/*
 * DeBruijnGraph.h
 *
 *  Created on: 04.12.2012
 *      Author: AlexanderDilthey
 *
 *  Interface and some code adapted from / inspired by Zamin Iqbal / cortex_var.
 *
 */

#ifndef DEBRUIJNGRAPH_H_
#define DEBRUIJNGRAPH_H_

#include <cstdio>
#include <cstring>
#include <string>
#include <array>
#include <assert.h>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <set>
#include <omp.h>
#include <map>
#include <cstdlib>
#include <queue>

#include "../Hsh.h"
#include "../../test.h"

#include "DeBruijnElement.h"

extern const int BINVERSION;
extern const int MAX_LEN_SAMPLE_NAME;

template<int m, int k>
class superNodeSpecifier
{
public:
	binarykMer<m, k> begin_kMer_key;					// key of first kMer
	binarykMer<m, k> end_kMer_key;						// key of last kMer
	bool begin_reverse;								// Is the first kMer (begin_kMer) of supernode inverted relative to the key of the first kMer?  Not valid for 1-kmer supernodes!!!
	bool end_reverse;								// Is the last kMer (end_kMer) of supernode inverted relative to the key of the last kMer?  Not valid for 1-kmer supernodes!!!
	signed char begin_direction_intoSupernode;		// if we are at begin_kmer, walking direction to enter supernode. Relative to key of begin_kMer. Not valid for 1-kmer supernodes!!!
	signed char end_direction_intoSupernode;		// if we are at end_kMer, what direction to walk to enter supernode? Relative to key of end_kMer. Not valid for 1-kmer supernodes!!!
	std::string sequence;							// sequence of the supernode - might be removed later
	int length;										// length of the supernode in kMers
	bool circular;
	superNodeSpecifier(){length = 0; begin_reverse = false; end_reverse = false; begin_direction_intoSupernode = 0; end_direction_intoSupernode = 0; circular = false;}

	std::string getSequence()
	{
		return sequence;
	}
};

typedef enum {
  EValid                        = 0,
  ECannotReadMagicNumber        = 1,
  ECanReadMagicNumberButIsWrong = 2,
  ECannotReadBinversion         = 3,
  EValidButOldBinVersion        = 4,
  EInvalidBinversion            = 5,
  ECannotReadKmer               = 6,
  EWrongKmer                    = 7,
  ECannotReadNumBitfields       = 8,
  EWrongNumberBitfields         = 9,
  ECannotReadNumColours         = 10,
  EBadColours                   = 11,
  EFailedToReadReadLensAndCovgs = 12,
  ECannotReadEndOfHeaderMagicNumber = 13,
  EFailedToReadSampleIds        =14,
  EFailedToReadSampleIdsSeemsTooLong = 15,
  EFailedToReadSeqErrRates      =16,
  EFailedToReadErrorCleaningInfo=17,
  ECanReadEndOfHeaderMagicNumberButIsWrong = 18,
  EGarbage = 19,
  EBinaryHasTooManyColoursGivenFirstColour=20, //ie starting at colour 10, loading 100 colours but 110>NUMBER_OF_COLOURS-1
} BinaryHeaderErrorCode;

class ErrorCleaning
{
public:
  bool tip_clipping;
  bool remv_low_cov_sups;
  bool remv_low_cov_nodes;
  int len_name_of_graph_against_which_was_cleaned;

  int remv_low_cov_sups_thresh;
  int remv_low_cov_nodes_thresh;

  bool cleaned_against_another_graph;
  std::string name_of_graph_against_which_was_cleaned;

  ErrorCleaning()
  {
	  tip_clipping = false;
	  remv_low_cov_sups = false;
	  remv_low_cov_nodes = false;
	  len_name_of_graph_against_which_was_cleaned = -1;
	  remv_low_cov_sups_thresh = -1;
	  remv_low_cov_nodes_thresh = -1;
	  cleaned_against_another_graph = false;
  }
};


template <int colours>
class GraphInfo
{
public:
  std::array<int, colours> sample_id_lens;
  std::array<std::string, colours> sample_id;
  std::array<uint64_t, colours> total_sequence;
  std::array<uint32_t, colours> mean_read_length;
  std::array<long double, colours> seq_err;
  std::array<ErrorCleaning, colours> cleaning;

  GraphInfo()
  {
	  sample_id_lens.fill(0);
	  total_sequence.fill(0);
	  mean_read_length.fill(0);
	  seq_err.fill(0);
  }
};

template <int colours>
class BinaryHeaderInfo
{
public:
  int version;
  int kmer_size;
  int number_of_bitfields;
  int number_of_colours;
  GraphInfo<colours> ginfo;

  BinaryHeaderInfo()
  {
	  version = 0;
	  kmer_size = 0;
	  number_of_bitfields = 0;
	  number_of_colours = 0;
  }
};

class supernodeNeighborDistance
{
public:
	std::map<int, int> fromLastkMer;
	std::map<int, int> fromFirstkMer;
};

template<int m, int k, int colours>
class DeBruijnGraph {
protected:
	Hsh<m, k, DeBruijnElement<colours> >* hash;
	std::vector< superNodeSpecifier<m, k> > graph_supernode_representation;

	std::map<unsigned int, std::set<int> > graph_supernode_neighbors;

	// the following structure will hold supernode neighbors and the associated distances
	// first unsigned int: 		supernode ID
	// second supernodeNeighborDistance: contains distance information to the neighbors, stratified by whether we
	//	start from the last or first kMer of the supernode. The first integer of the maps is really signed and negative
	//		if we enter the neighboring supernode via the last kMer.

	std::map<unsigned int, supernodeNeighborDistance> graph_supernode_neighborsWithDistance;

	int supernode_neighbor_distance;

	bool read_next_error_cleaning_object(FILE* fp, ErrorCleaning& cl);
	bool db_node_read_multicolour_binary(FILE* fp, binarykMer<m, k>& kMer, DeBruijnElement<colours>& node, int num_colours_in_binary, int binversion_in_binheader);
	bool check_binary_compatibility(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into);
	bool get_binversion6_extra_data(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into);
	bool get_read_lengths_and_total_seqs_from_header(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into);

	unsigned int threads;

	bool paranoid;

public:
	DeBruijnGraph(int log2_rows, int bucketsPerRow)
	{
		hash = new Hsh<m, k, DeBruijnElement<colours> >(log2_rows, bucketsPerRow);
		supernode_neighbor_distance = -1;
		threads = 1;
		paranoid = true;
	}

	Hsh<m, k, DeBruijnElement<colours> >* getHash()
	{
		return hash;
	}

	// return the supernode neighbors for a given supernode
	// the results set has signed integers!
	std::set<int> supernodeGetNeighbors(unsigned int supernodeIndex)
	{
		assert(supernodeIndex >= 1);
		return graph_supernode_neighbors[supernodeIndex-1];
	}

	~DeBruijnGraph()
	{
		delete(hash);
	}

	unsigned int getSupernodeTableSize()
	{
		return graph_supernode_representation.size();
	}

	superNodeSpecifier<m, k> getSupernode(unsigned int i)
	{
		assert(i >= 1);
		assert((i - 1) < graph_supernode_representation.size());
		return graph_supernode_representation.at(i-1);
	}

	long long loadMultiColourBinary(std::string filename)
	{
	    //printf("Load this binary - %s\n", filename);
	    FILE* fp_bin = fopen(filename.c_str(), "r");
	    long long  seq_length = 0;

	    if (fp_bin == NULL)
	    {
	        throw std::runtime_error("load_multicolour_binary_from_filename_into_graph cannot open file "+filename);
	    }

	    int count=0;

	    BinaryHeaderErrorCode ecode = EValid;
	    BinaryHeaderInfo<colours> binfo;
	    if(! check_binary_compatibility(fp_bin, binfo, ecode, 0))
	    {
	    	throw std::runtime_error("Cannot load binary "+filename+": binary problem detected!");
	    }

	    int num_cols_in_loaded_binary = binfo.number_of_colours;

	    //always reads the multicol binary into successive colours starting from 0 - assumes the hash table is empty prior to this

	    binarykMer<m, k> kMer;
	    DeBruijnElement<colours> node;

	    while (db_node_read_multicolour_binary(fp_bin, kMer, node, num_cols_in_loaded_binary, binfo.version))
	    {
	        count++;

			
	        DeBruijnElement<colours>* element_insert_position = hash->hash_table_insert(kMer);
	        *element_insert_position = node;
	        seq_length += k;
			node.nullify();
	    }

	    fclose(fp_bin);
	    return seq_length;
	}

	void checkGraphIntegrity()
	{
		long long global_kmers_present = 0;
		long long global_forward_edges_present = 0;
		long long global_backward_edges_present = 0;
		long long global_nonAssigned_cells = 0;
		
		omp_set_num_threads(threads);		

		long long max_i = hash->getTableSize() - 1;
		long long chunk_size = max_i / threads;
	
		#pragma omp parallel
		{
			assert(omp_get_num_threads() == (int)threads);
			unsigned int thisThread = omp_get_thread_num();
			long long firstPair = thisThread * chunk_size;
			long long lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (threads-1)) && (lastPair < max_i))
			{
				lastPair = max_i;
			}		
			
			long long local_kmers_present = 0;
			long long local_forward_edges_present = 0;
			long long local_backward_edges_present = 0;		
			long long local_nonAssigned_cells = 0;
			
			for(long long i = firstPair; i <= lastPair; i++)
			{
				hashTableElement<m, k, DeBruijnElement<colours>>* thisEntry = hash->getTableElement(i);
				if(thisEntry->flags.assigned == true)
				{
					local_kmers_present++;
					std::set<Nucleotide> forward_extensions = thisEntry->element.getEdges(false);
					local_forward_edges_present += forward_extensions.size();
					
					std::set<Nucleotide> backward_extensions = thisEntry->element.getEdges(true);
					local_backward_edges_present += backward_extensions.size();

					for(std::set<Nucleotide>::iterator nIt = forward_extensions.begin(); nIt != forward_extensions.end(); nIt++)
					{
						binarykMer<m, k> kMerCopy = thisEntry->kMer;
						kMerCopy.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
						binarykMer<m, k> kMerKey = kMerCopy.element_get_key();
						hashTableElement<m, k, DeBruijnElement<colours>>* nextFullElement = hash->hash_table_find_fullElement(kMerKey);
						
						if(nextFullElement == 0)
						{
							std::cout << "Graph integrity error - kMer as specified by edge not found!\n";
							std::cout << "Original kMer key: " << thisEntry->kMer.getSeq() << "\n";
							std::cout << "Reverse complement: no" << "\n";
							std::cout << "Edge: " << binary_nucleotide_to_char(*nIt) << "\n";
							std::cout << "Resulting shifted kMer: " << kMerCopy.getSeq() << "\n";
							std::cout << "Key of shifted kMer:: " << kMerKey.getSeq() << "\n" << std::flush;
						}
						
						assert(nextFullElement != 0);	
					}
					
					for(std::set<Nucleotide>::iterator nIt = backward_extensions.begin(); nIt != backward_extensions.end(); nIt++)
					{
						binarykMer<m, k> kMerCopy = thisEntry->kMer;
						kMerCopy.reverse_complement();
						kMerCopy.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
						binarykMer<m, k> kMerKey = kMerCopy.element_get_key();
						hashTableElement<m, k, DeBruijnElement<colours>>* nextFullElement = hash->hash_table_find_fullElement(kMerKey);
						
						if(nextFullElement == 0)
						{
							std::cout << "Graph integrity error - kMer as specified by edge not found!\n";
							std::cout << "Original kMer key: " << thisEntry->kMer.getSeq() << "\n";
							std::cout << "Reverse complement: yes" << "\n";
							std::cout << "Edge: " << binary_nucleotide_to_char(*nIt) << "\n";
							std::cout << "Resulting shifted kMer: " << kMerCopy.getSeq() << "\n";
							std::cout << "Key of shifted kMer:: " << kMerKey.getSeq() << "\n" << std::flush;
						}
						
						assert(nextFullElement != 0);	
					}		
				}
				else
				{
					local_nonAssigned_cells++;
				}
			}
			
			#pragma omp critical
			{
				global_kmers_present += local_kmers_present;
				global_forward_edges_present += local_forward_edges_present;
				global_backward_edges_present += local_backward_edges_present;			
				global_nonAssigned_cells += local_nonAssigned_cells;
			}
		}
		
		std::cout << "Graph integrity checked and OK. We have " << global_kmers_present << " Mers in total, with " << global_forward_edges_present << " forward edges and " << global_backward_edges_present << " backward edges\n" << std::flush;		
	}
	
	long long totalCoverage(int colour = 0)
	{
		long long tC = 0;

		omp_set_num_threads(threads);

		long long max_i = hash->getTableSize() - 1;
		long long chunk_size = max_i / threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == (int)threads);
			long long thisThread = omp_get_thread_num();
			long long firstPair = thisThread * chunk_size;
			long long lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (threads-1)) && (lastPair < max_i))
			{
				lastPair = max_i;
			}

			long long tC_thisThread = 0;

			std::vector< superNodeSpecifier<m, k> > local_found_supernodes;
			std::set< binarykMer<m, k>* > local_supernodeStartsorEnds;
			std::set< binarykMer<m, k>* > local_supernodeAlready;

			for(long long i = firstPair; i <= lastPair; i++)
			{
				hashTableElement<m, k, DeBruijnElement<colours>>* thisEntry = hash->getTableElement(i);
				tC_thisThread += thisEntry->element.coverage.at(colour);
			}

			#pragma omp critical
			{
				tC += tC_thisThread;
			}
		}

		return tC;
	}


	uint32_t kMer_getCoverage(std::string kMerSeq, int colour = 0)
	{
		assert(kMerSeq.length() == k);
		binarykMer<m, k> key;
		key.loadSeq(kMerSeq);
		binarykMer<m, k> queryKey = key.element_get_key();
		DeBruijnElement<colours>* element = hash->hash_table_find(queryKey);
		if(element == 0)
		{
			return 0;
		}
		else
		{
			return element->coverage.at(colour);
		}
	}

	bool kMerinGraph(std::string kMerSeq)
	{
		assert(kMerSeq.length() == k);
		binarykMer<m, k> key;
		key.loadSeq(kMerSeq);   
		return kMerinGraph(key);
	}

	bool kMerinGraph(binarykMer<m, k>& key)
	{
		bool egal;
		binarykMer<m, k> queryKey = key.element_get_key(egal);
		//std::cout << key.getSeq() << " " << queryKey.getSeq() << "\n" << std::flush;
		return (hash->hash_table_find(queryKey) != 0);
	}

	std::vector<binarykMer<m, k>> walkOneStep(binarykMer<m, k>& key, bool reverse)
	{
		return walkOneStep(0, key, reverse);
	}


	std::vector<binarykMer<m, k>> walkOneStep(int colour, binarykMer<m, k>& kMer, bool reverse)
	{
		assert(colour >= 0);
		assert(colour < colours);

		bool key_was_inverted;
		binarykMer<m, k> key = kMer.element_get_key(key_was_inverted);
		DeBruijnElement<colours>* thisElement = hash->hash_table_find(key);

		std::vector<binarykMer<m, k>> forReturn;

		if(reverse == false)
		{
			if(key_was_inverted == false)
			{
				std::set<Nucleotide> extensions = thisElement->getEdges(colour, false);
				for(std::set<Nucleotide>::iterator nIt = extensions.begin(); nIt != extensions.end(); nIt++)
				{
					binarykMer<m, k> newkMer(key);
					newkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
					forReturn.push_back(newkMer);
				}
			}
			else if(key_was_inverted == true)
			{
				std::set<Nucleotide> extensions = thisElement->getEdges(colour, true);
				for(std::set<Nucleotide>::iterator nIt = extensions.begin(); nIt != extensions.end(); nIt++)
				{
					binarykMer<m, k> newkMer(kMer); // kMer is the inverted key
					newkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
					forReturn.push_back(newkMer);
				}
			}
		}
		else if(reverse == true)
		{
			if(key_was_inverted == false)
			{
				key.reverse_complement();
				std::set<Nucleotide> extensions = thisElement->getEdges(colour, true);
				for(std::set<Nucleotide>::iterator nIt = extensions.begin(); nIt != extensions.end(); nIt++)
				{
					binarykMer<m, k> newkMer(key);
					newkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
					newkMer.reverse_complement();
					forReturn.push_back(newkMer);
				}
			}
			else if(key_was_inverted == true)
			{
				std::set<Nucleotide> extensions = thisElement->getEdges(colour, false);
				for(std::set<Nucleotide>::iterator nIt = extensions.begin(); nIt != extensions.end(); nIt++)
				{
					binarykMer<m, k> newkMer(key);
					newkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
					newkMer.reverse_complement();
					forReturn.push_back(newkMer);
				}
			}
		}

		assert(thisElement != 0);

		return forReturn;
	}

	void loadReadChunk(std::vector<readPair>& reads)
	{
		for(unsigned int i = 0; i < reads.size(); i++)
		{
			loadSequence(reads.at(i).r1);
			loadSequence(reads.at(i).r2);
		}
	}

	void loadSequence(std::string sequence)
	{
		_loadSequence(sequence);
		_loadSequence(seq_reverse_complement(sequence));
	}



	void _loadSequence(std::string sequence)
	{
		assert(sequence.length() >= k);

		int currentFirstPos = find_first_valid_kMer_position(sequence, k);
		while(currentFirstPos != -1)
		{
			// std::cout << "currentFirstPos: " << currentFirstPos << "\n" << std::flush;

			std::string firstkMerSeq = sequence.substr(currentFirstPos, k);
			assert(firstkMerSeq.length() == k);

			binarykMer<m, k> runningkMer;
			runningkMer.loadSeq(firstkMerSeq);
			bool key_is_inverted;
			binarykMer<m, k> runningkMerKey = runningkMer.element_get_key(key_is_inverted);
			bool ignore;
			DeBruijnElement<colours>* runningElement = hash->hash_table_find_or_insert(runningkMerKey, ignore);
			// std::cout << "Inserted kMer " << runningkMer.getSeq() << "\n" << std::flush;

			//std::cout << "Address runningElement: " << runningElement << "\n" << std::flush;

			// currentFirstPos + k is the position of the character coming after the current kMer
			while((currentFirstPos + k) < (int)sequence.length())
			{
				char nextCharacter = sequence.at(currentFirstPos + k);
				Nucleotide nextCharacterAsNucleotide = char_to_binary_nucleotide(nextCharacter);

				if(nextCharacterAsNucleotide == Undefined)
				{
					currentFirstPos = find_first_valid_kMer_position(sequence, k, currentFirstPos + k + 1);
					break;
				}
				else
				{

					// std::cout << "Address runningElement: " << runningElement << "\n" << std::flush;

					// std::cout << "Setting edge for character " << nextCharacter << " for kMer " << runningkMer.getSeq() << " with key " << runningkMerKey.getSeq() << " which is inverted? " << key_is_inverted << " runningElement address: " << (void*) &runningElement << "\n" << std::flush;

					//std::cout << "\ncurrent kMer: " << runningkMer.getSeq() << "\n";
					//std::cout << "\t current associated element address: " << runningElement << "\n";

					runningElement->setEdge(nextCharacterAsNucleotide, key_is_inverted);
					runningkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(nextCharacterAsNucleotide);
					runningkMerKey = runningkMer.element_get_key(key_is_inverted);
					runningElement = hash->hash_table_find_or_insert(runningkMerKey, ignore);

					// std::cout << "Inserted kMer " << runningkMer.getSeq() << "\n" << std::flush;

					//std::cout << "\nnext kMer: " << runningkMer.getSeq() << "\n";
					//std::cout << "\tnext kMer key: " << runningkMerKey.getSeq() << "\n";

					//std::cout << "\tnext kMer associated element address: " << runningElement << "\n" << std::flush;

					currentFirstPos++;
				}
			}
			if((currentFirstPos + k) == (int)sequence.length())
			{
				currentFirstPos = -1;
			}
		}
	}

	bool checkSequenceForPresence(std::string sequence)
	{
		if(sequence.length() < k)
		{
			std::cout << "sequence.length() < k" << "\n" << std::flush;
			return false;
		}

		int currentFirstPos = find_first_valid_kMer_position(sequence, k);
		if(currentFirstPos != 0)
		{
			std::cout << "currentFirstPos != 0" << "\n" << std::flush;
			return false;
		}

		std::string firstkMerSeq = sequence.substr(currentFirstPos, k);
		assert(firstkMerSeq.length() == k);

		binarykMer<m, k> runningkMer;
		runningkMer.loadSeq(firstkMerSeq);
		bool key_is_inverted;
		binarykMer<m, k> runningkMerKey = runningkMer.element_get_key(key_is_inverted);
		DeBruijnElement<colours>* runningElement = hash->hash_table_find(runningkMerKey);
		if(runningElement == 0)
		{
			std::cout << "I) runningElement == 0" << "\n" << std::flush;
			return false;
		}

		// currentFirstPos + k is the position of the character coming after the current kMer
		while((currentFirstPos + k) < (int)sequence.length())
		{
			char nextCharacter = sequence.at(currentFirstPos + k);
			Nucleotide nextCharacterAsNucleotide = char_to_binary_nucleotide(nextCharacter);
			if(nextCharacterAsNucleotide == Undefined)
			{
				//std::cout << "nextCharacterAsNucleotide == Undefined" << "\n" << std::flush;
				return false;
			}
			else
			{
				bool foundEdge = runningElement->checkEdge(nextCharacterAsNucleotide, key_is_inverted);

				// std::cout << "Checking edge for character " << nextCharacter << " for kMer " << runningkMer.getSeq() << " which is inverted? " << key_is_inverted << " runningElement address: " << runningElement << "\n" << std::flush;

				if(! foundEdge)
				{
					// std::cout << "! foundEdge" << "\n" << std::flush;
					return false;
				}

				runningkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(nextCharacterAsNucleotide);
				runningkMerKey = runningkMer.element_get_key(key_is_inverted);
				runningElement = hash->hash_table_find(runningkMerKey);

				if(runningElement == 0)
				{
					// std::cout << "II) runningElement == 0" << "\n" << std::flush;
					return false;
				}

				currentFirstPos++;
			}
		}

		return true;
	}

	// this function determines whether a kMer is the beginning or the end of a supernode
	// cannot detect entirely circular stuff.
	// we apply the following criteria:
	// - if a kMer has 0 outgoing edges in either direction, it has to be beginning or end
	//		- if there is one direction with exactly one outgoing edge, and there is only one reverse edge
	//        from the targeted edge, we can extend the supernode that way
	// - if a kMer has more than one outgoing edge in either direction, it is beginning or end
	//		- if there is one direction with exactly one outgoing edge, and there is only one reverse edge
	//        from the targeted edge, we can extend the supernode that way
	// - if a kMer has exactly one edge in both directions, we need to find out whether the targeted edges are
	//   ends of supernodes by their own criteria (that is, if there is more than one reverse edge from the targeted node,
	//   they have to be beginnings/ends of supernodes by themselves)
	//
	// We return true/false indicating whether we are at the end of a supernode, and
	// and integer that specifies possible move directions. The integer is *always relative to the kMer key*!
	//   0: can move neither or either direction
	//   -1: can move along the reverse complement edges
	//  1: can move along the non-complemented edges

	bool kMerCircularConnection(binarykMer<m, k>& kMerKey)
	{
		bool forReturn = false;

		DeBruijnElement<colours>* thisElement = hash->hash_table_find(kMerKey);
		assert(thisElement != 0);
		int num_outgoing_edges = thisElement->numEdges(false);
		int num_incoming_edges = thisElement->numEdges(true);

		if(num_incoming_edges == 1)
		{
			std::vector<binarykMer<m, k>> previousMer = walkOneStep(kMerKey, true);
			assert(previousMer.size() == 1);

			if(previousMer.at(0).element_get_key() == kMerKey)
			{
				forReturn = true;
			}
		}

		if(num_outgoing_edges == 1)
		{
			std::vector<binarykMer<m, k>> nextkMer = walkOneStep(kMerKey, false);
			assert(nextkMer.size() == 1);

			if(nextkMer.at(0).element_get_key() == kMerKey)
			{
				forReturn = true;
			}
		}

		return forReturn;

	}

	bool kMerIsBeginningOrEndOfSupernode(binarykMer<m, k>& kMer, int& extensionDirection)
	{
		bool key_was_inverted;
		binarykMer<m, k> key = kMer.element_get_key(key_was_inverted);
		DeBruijnElement<colours>* thisElement = hash->hash_table_find(key);
		assert(thisElement != 0);

		int num_outgoing_edges = thisElement->numEdges(false);
		int num_incoming_edges = thisElement->numEdges(true);


		bool canExtendForward = true;
		bool canExtendBackward = true;

		bool forReturn = false;

		if(num_incoming_edges == 0)
		{
			canExtendBackward = false;
			forReturn = true;
		}
		else if(num_incoming_edges > 1)
		{
			canExtendBackward = false;
			forReturn = true;
		}
		else if(num_incoming_edges == 1)
		{
			std::vector<binarykMer<m, k>> previousMer = walkOneStep(key, true);
			assert(previousMer.size() == 1);

			if(previousMer.at(0).element_get_key() == key)
			{
				forReturn = true;
				canExtendBackward = false;
				canExtendForward = false;
			}
			else
			{
				bool previous_key_was_inverted;
				binarykMer<m, k> previousKey = previousMer.at(0).element_get_key(previous_key_was_inverted);
				DeBruijnElement<colours>* previousElement = hash->hash_table_find(previousKey);
				assert(previousElement != 0);

				int num_previous_outgoing_edges = previousElement->numEdges(previous_key_was_inverted);
				assert(num_previous_outgoing_edges >= 1);
				if((num_previous_outgoing_edges > 1) || (kMerCircularConnection(previousKey)))
				{
					forReturn = true;
					canExtendBackward = false;
				}
			}
		}

		if(num_outgoing_edges == 0)
		{
			canExtendForward = false;
			forReturn = true;
		}
		else if(num_outgoing_edges > 1)
		{
			canExtendForward = false;
			forReturn = true;
		}
		else if (num_outgoing_edges == 1)
		{
			std::vector<binarykMer<m, k>> nextkMer = walkOneStep(key, false);
			assert(nextkMer.size() == 1);

			if(nextkMer.at(0).element_get_key() == key)
			{
				forReturn = true;
				canExtendForward = false;
				canExtendBackward = false;
			}
			else
			{
				bool next_key_was_inverted;
				binarykMer<m, k> nextKey = nextkMer.at(0).element_get_key(next_key_was_inverted);
				DeBruijnElement<colours>* nextElement = hash->hash_table_find(nextKey);
				assert(nextElement != 0);

				int num_next_incoming_edges = nextElement->numEdges(! next_key_was_inverted);
				assert(num_next_incoming_edges >= 1);
				if((num_next_incoming_edges > 1) || (kMerCircularConnection(nextKey)))
				{
					forReturn = true;
					canExtendForward = false;
				}
			}
		}

		if(forReturn == false)
		{
			assert(canExtendForward && canExtendBackward);
			extensionDirection = 0;
		}
		else
		{
			assert(! (canExtendForward && canExtendBackward));


			if(canExtendForward)
			{
				extensionDirection = 1;
			}
			if(canExtendBackward)
			{
				extensionDirection = -1;
			}
			if((! canExtendForward) && (! canExtendBackward))
			{
				extensionDirection = 0;
			}
		}

		return forReturn;

	}

	bool testkMerWalkValidity(std::vector<binarykMer<m, k>>& kMers)
	{
		std::string prefix = "";
		for(unsigned int i = 0; i < kMers.size(); i++)
		{
			std::string kMerSequence = kMers.at(i).getSeq();
			if(prefix != "")
			{
				if(kMerSequence.substr(0, k - 1) != prefix)
				{
					return false;
				}
			}
			prefix = kMerSequence.substr(1);
			assert(prefix.size() == (k - 1));

			bool inverted_key;
			binarykMer<m, k> kMerKey = kMers.at(i).element_get_key(inverted_key);

			DeBruijnElement<colours>* kMerElement = hash->hash_table_find(kMerKey);
			if(i != (kMers.size() - 1))
			{
				std::string nextMerSequence = kMers.at(i+1).getSeq();
				char edgeCharacter = nextMerSequence.at(k - 1);
				if( ! kMerElement->checkEdge(char_to_binary_nucleotide(edgeCharacter), inverted_key))
				{
					return false;
				}
			}

		}
		return true;
	}

	// generates a random walk through the graph. maximumLength specifies the number of steps.
	std::vector< std::vector<binarykMer<m, k> >> DFlimitedSearch(binarykMer<m, k> startkMer, int depth)
	{
		bool startInverted;
		binarykMer<m, k> currentKey = startkMer.element_get_key(startInverted);

		bool currentlyInverted = startInverted;

		class extensionTip
		{
		public:
			std::vector<binarykMer<m, k>> path;
			binarykMer<m, k> tip;
		};

		std::queue<extensionTip> search;
		std::vector< std::vector<binarykMer<m, k>>> results;

		extensionTip start;
		start.tip = startkMer;
		start.path.reserve(depth);
		search.push(start);

		while(search.size() > 0)
		{
			extensionTip extension = search.front();
			search.pop();

			if((int)extension.path.size() >= depth)
			{
				results.push_back(extension.path);
			}
			else
			{
				bool kMerInverted = false;
				binarykMer<m, k> kMerKey = extension.tip.element_get_key(kMerInverted);

				hashTableElement<m, k, DeBruijnElement<colours>>* currentFullElement = hash->hash_table_find_fullElement(kMerKey);
				assert(currentFullElement != 0);

				std::set<Nucleotide> extensions = currentFullElement->element.getEdges(kMerInverted);
				for(std::set<Nucleotide>::iterator nIt = extensions.begin(); nIt != extensions.end(); nIt++)
				{
					extensionTip newExtension = extension;
					newExtension.tip.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(*nIt);
					newExtension.path.push_back(newExtension.tip);
					search.push(newExtension);
				}

				if((extensions.size() == 0) && (extension.path.size() > 0))
				{
					results.push_back(extension.path);
				}
			}
		}

		return results;
	}

	// generates a random walk through the graph. maximumLength specifies the number of steps.
	std::string randomWalk(int maximumLength, bool startOnlyAtSupernode = false)
	{
		binarykMer<m, k> startkMer = quickRandomkMer(startOnlyAtSupernode);
		std::string forReturn = startkMer.getSeq();
		unsigned int startInverted = rand() % 2;
		bool currentlyInverted = false;
		assert((startInverted == 0) || (startInverted == 1));
		if(startOnlyAtSupernode)
		{
			unsigned int supernodeID = supernodeIndexForkMer(startkMer);
			superNodeSpecifier<m, k> supernode = getSupernode(supernodeID);
			assert(supernode.begin_kMer_key == startkMer);
			if(supernode.length == 1)
			{
				// startInverted = startInverted;
			}
			else
			{
				if(supernode.begin_direction_intoSupernode == 1)
				{
					startInverted = false;
				}
				else
				{
					assert(supernode.begin_direction_intoSupernode == -1);
					startInverted = true;
				}
			}
		}

		if(startInverted)
		{
			startkMer.reverse_complement();
			forReturn = seq_reverse_complement(forReturn);
			currentlyInverted = true;
		}

		bool _inv;
		binarykMer<m, k> currentKey = startkMer.element_get_key(_inv);
		assert(_inv == currentlyInverted);

		bool madeExtension;
		do {
			hashTableElement<m, k, DeBruijnElement<colours>>* currentFullElement = hash->hash_table_find_fullElement(currentKey);
			assert(currentFullElement != 0);

			std::set<Nucleotide> extensions = currentFullElement->element.getEdges(currentlyInverted);
			madeExtension = false;
			if(extensions.size() > 0)
			{
				madeExtension = true;

				std::vector<Nucleotide> extensions_vec(extensions.begin(), extensions.end());
				assert(extensions_vec.size() == extensions.size());

				unsigned int selectExtensionIndex = rand() % extensions_vec.size();
				assert(selectExtensionIndex >= 0);
				assert(selectExtensionIndex < extensions_vec.size());

				Nucleotide n = extensions_vec.at(selectExtensionIndex);

				forReturn.push_back(binary_nucleotide_to_char(n));

				binarykMer<m, k> newkMer(currentKey);
				if(currentlyInverted)
				{
					newkMer.reverse_complement();
				}
				newkMer.bitOps_left_shift_one_base_and_insert_new_base_at_right_end(n);

				currentKey = newkMer.element_get_key(currentlyInverted);
			}
		} while(madeExtension && ((forReturn.length() - k + 1) < maximumLength));


		return forReturn;
	}

	// quick random KMer selection - note - this might not truly be random!
	binarykMer<m, k> quickRandomkMer(bool supernodeStartOnly)
	{
		hashTableElement<m, k, DeBruijnElement<colours>>* randomElement = hash->quicklyGetRandomElement();
		assert(randomElement != 0);
		if(supernodeStartOnly)
		{
			unsigned int supernodeID = supernodeIndexForkMer(randomElement->kMer);
			return getSupernode(supernodeID).begin_kMer_key;
		}
		return randomElement->kMer;
	}



};


template<int m, int k, int colours>
bool DeBruijnGraph<m, k, colours>::read_next_error_cleaning_object(FILE* fp, ErrorCleaning& cl)
{
    int read;
    bool no_problem=true;

    signed char tip_clipping;
    signed char remv_low_cov_sups;
    signed char remv_low_cov_nodes;
    signed char cleaned_against_another_graph;


    read = fread(&tip_clipping, sizeof(signed char), 1, fp);

    if (read==0)
    {
        no_problem=false;
    }
    if (no_problem==true)
    {
        read = fread(&remv_low_cov_sups, sizeof(signed char), 1, fp);
        if (read==0)
        {
            no_problem=false;
        }
    }

    if (no_problem==true)
    {
        read = fread(&remv_low_cov_nodes, sizeof(signed char), 1, fp);
        if (read==0)
        {
            no_problem=false;
        }
    }
    if (no_problem==true)
    {
        read = fread(&cleaned_against_another_graph, sizeof(signed char), 1, fp);
        if (read==0)
        {
            no_problem=false;
        }
    }

    if(no_problem)
    {
    	cl.tip_clipping = (bool)tip_clipping;
    	cl.remv_low_cov_nodes = (bool) remv_low_cov_nodes;
    	cl.remv_low_cov_sups = (bool) remv_low_cov_sups;
    	cl.cleaned_against_another_graph = (bool) cleaned_against_another_graph;
    }

    if (no_problem==true)
    {
        read = fread(&(cl.remv_low_cov_sups_thresh),sizeof(int),1,fp);
        if (read==0)
        {
            no_problem=false;
        }
    }

    if (no_problem==true)
    {
        read = fread(&(cl.remv_low_cov_nodes_thresh),sizeof(int),1,fp);
        if (read==0)
        {
            no_problem=false;
        }
    }

    if (no_problem==true)
    {
        read = fread(&(cl.len_name_of_graph_against_which_was_cleaned),sizeof(int),1,fp);
        if (read==0)
        {
            no_problem=false;
        }
    }

    if(cl.len_name_of_graph_against_which_was_cleaned == -1)
    {
		cl.name_of_graph_against_which_was_cleaned = "";
    }
    else
    {

		if (no_problem==true)
		{
			char tmp_name[cl.len_name_of_graph_against_which_was_cleaned];

			read = fread(tmp_name, sizeof(char), cl.len_name_of_graph_against_which_was_cleaned, fp);

			if (read==0)
			{
				no_problem=false;
				cl.name_of_graph_against_which_was_cleaned = std::string(tmp_name);

			}
		}
    }


    return no_problem;

}

template<int m, int k, int colours>
bool DeBruijnGraph<m, k, colours>::get_read_lengths_and_total_seqs_from_header(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into)
{
    int read;
    int i;
    bool no_problem=true;

    for (i=first_colour_loading_into; ( i < (first_colour_loading_into+binfo.number_of_colours)) && (no_problem == true); i++)
    {
        int mean_read_len=0;
        read = fread(&mean_read_len, sizeof(int),1,fp);
        if (read==0)
        {
            no_problem=false;
            ecode= EFailedToReadReadLensAndCovgs;
        }
        else
        {
            binfo.ginfo.mean_read_length.at(i) += mean_read_len;
        }
    }
    if (no_problem==true)
    {
        for (i=first_colour_loading_into; (i< (first_colour_loading_into + binfo.number_of_colours) ) && (no_problem==true); i++)
        {
            long long tot=0;
            read = fread(&tot, sizeof(long long), 1, fp);
            if (read==0)
            {
                no_problem=false;
                ecode= EFailedToReadReadLensAndCovgs;
            }
            else
            {
                binfo.ginfo.total_sequence.at(i) +=tot;
            }
        }
    }

    return no_problem;
}

template<int m, int k, int colours>
bool DeBruijnGraph<m, k, colours>::get_binversion6_extra_data(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into)
{
    int read;
    int i;
    bool no_problem = true;

    //first get sample information
    for (i=first_colour_loading_into; (i<(first_colour_loading_into+binfo.number_of_colours)) && (no_problem==true); i++)
    {
        //get the length of the sample id name
        read = fread(&(binfo.ginfo.sample_id_lens.at(i)), sizeof(int), 1, fp);

        if (read==0)
        {
            no_problem = false;
            ecode = EFailedToReadSampleIds;
        }
        else
        {
            //now get the actual sample id for colour i
            if (binfo.ginfo.sample_id_lens.at(i) < MAX_LEN_SAMPLE_NAME)
            {
                char tmp_name[binfo.ginfo.sample_id_lens.at(i)+1];

                memset(tmp_name, 0, binfo.ginfo.sample_id_lens.at(i)+1);

                read = fread(tmp_name, sizeof(char), binfo.ginfo.sample_id_lens.at(i), fp);

                if (read==0)
                {
                    no_problem = false;
                    ecode = EFailedToReadSampleIds;
                }
                else
                {
                	binfo.ginfo.sample_id.at(i) = std::string(tmp_name);
                }

            }
            else
            {
                no_problem=false;
                ecode = EFailedToReadSampleIdsSeemsTooLong;

            }
        }
    }

    //now get the sequencing error rate
    for (i=first_colour_loading_into; (i<(first_colour_loading_into+binfo.number_of_colours)) && (no_problem==true); i++)
    {
    	long double seqTemp;

    	read = fread(&seqTemp, sizeof(long double), 1, fp);

        if (read == 0)
        {
            no_problem=false;
            throw std::runtime_error("i is %d and num colours is %d, but read is zero  - problem reading binary header. Contact Alex -- and yes, there is no substitution!\n");
            ecode = EFailedToReadSeqErrRates;
        }

        binfo.ginfo.seq_err.at(i) = seqTemp;

    }

    //now get the error cleaning information for each colour
    for (i=first_colour_loading_into; (i< (first_colour_loading_into+binfo.number_of_colours) ) && (no_problem==true); i++)
    {
        no_problem = read_next_error_cleaning_object(fp, binfo.ginfo.cleaning.at(i));
        if (no_problem==false)
        {
           ecode = EFailedToReadErrorCleaningInfo;
        }

    }

    return no_problem;
}


//return true if signature is readable, checks binversion, number of bitfields, magic number.
//does not check kmer is compatible with number of bitfields, leaves that to caller.
template<int m, int k, int colours>
bool DeBruijnGraph<m, k, colours>::check_binary_compatibility(FILE* fp, BinaryHeaderInfo<colours>& binfo, BinaryHeaderErrorCode& ecode, int first_colour_loading_into)
{
    int read;
    char magic_number[6];

    ecode = EValid;
    read = fread(magic_number, sizeof(char), 6, fp);

    if (read>0)
    {
        if (       magic_number[0]=='C' &&
                   magic_number[1]=='O' &&
                   magic_number[2]=='R' &&
                   magic_number[3]=='T' &&
                   magic_number[4]=='E' &&
                   magic_number[5]=='X' )
        {


            read = fread(&(binfo.version),sizeof(int),1,fp);
            if (read>0)
            {
                //can read version
                if ((binfo.version >= 4) && (binfo.version <= BINVERSION) )
                {
                    //version is good
                    read = fread(&(binfo.kmer_size), sizeof(int), 1, fp);
                    if (read>0)
                    {
                    	if(binfo.kmer_size == k)
                    	{
							//can read bitfields
							read = fread(&(binfo.number_of_bitfields), sizeof(int), 1, fp);

							if (binfo.number_of_bitfields == m)
							{
								//bitfields are good


								read = fread(&(binfo.number_of_colours), sizeof(int), 1, fp);

								if (read > 0)
								{
									//can read colours

									if (binfo.number_of_colours <= colours)
								{
										//colours are good

										//ok, the basic information looks OK
										// get extra information. In all cases, return false and an error code if anything looks bad.
										bool no_problem=true;


										if ((first_colour_loading_into + binfo.number_of_colours) > (colours - 1))
										{
											ecode = EBinaryHasTooManyColoursGivenFirstColour;
											no_problem=false;
										}

										if (( binfo.version == 5 ) || ( binfo.version == 4 ))//legacy
										{

											no_problem = get_read_lengths_and_total_seqs_from_header(fp, binfo, ecode, first_colour_loading_into);
										}
										else if (binfo.version==6)
										{

											no_problem = get_read_lengths_and_total_seqs_from_header(fp, binfo, ecode, first_colour_loading_into);

											if (no_problem==true)
											{
												//get the extra stuff for binary version 6

												no_problem = get_binversion6_extra_data(fp, binfo, ecode, first_colour_loading_into);

											}
										}
										else
										{
											ecode = EInvalidBinversion;
											no_problem = false;
										}

										if (no_problem == true)
										{

											//only thing remaining to check is the end of header magic number
											int read;
											char magic_number[6];
											magic_number[0]='\0';
											magic_number[1]='\0';
											magic_number[2]='\0';
											magic_number[3]='\0';
											magic_number[4]='\0';
											magic_number[5]='\0';
											read = fread(magic_number, sizeof(char), 6, fp);
											if (read==0)
											{
												ecode = ECannotReadEndOfHeaderMagicNumber;
												no_problem = false;
											}
											else if (
												magic_number[0]=='C' &&
												magic_number[1]=='O' &&
												magic_number[2]=='R' &&
												magic_number[3]=='T' &&
												magic_number[4]=='E' &&
												magic_number[5]=='X' )
											{
												//all good.
											}
											else
											{
												no_problem = false;
												ecode = ECanReadEndOfHeaderMagicNumberButIsWrong;
											}
										}

										return no_problem;
									}
									else
									{
										//colours bad
										ecode = EBadColours;
										return false;
									}
								}
								else
								{
									//cant read colours
									ecode = ECannotReadNumColours;
									return false;
								}
							}//bitfields are good
							else
							{
								//bitfields are bad
								ecode =  EWrongNumberBitfields;
								return false;
							}
                    	}
                    	else
                    	{
                    		ecode = EWrongKmer;
                    		return false;
                    	}
                    }//can read bitfields
                    else
                    {
                        //cannot read bitfields
                        ecode = ECannotReadNumBitfields;
                        return false;
                    }
                }//version is good
                else
                {
                    //version is bad
                    ecode = EInvalidBinversion;
                    return false;
                }
            }//can read version
            else
            {
                //cannot read version
                ecode =  ECannotReadBinversion;
                return false;
            }
        }//magic number good
        else
        {
            ecode=  ECanReadMagicNumberButIsWrong;
            return false;
        }
    }
    else
    {
        ecode = ECannotReadMagicNumber;
        return false;
    }
}

template<int m, int k, int colours>
bool DeBruijnGraph<m, k, colours>::db_node_read_multicolour_binary(FILE* fp, binarykMer<m, k>& kMer, DeBruijnElement<colours>& node, int num_colours_in_binary, int binversion_in_binheader)
{
	assert ((num_colours_in_binary <= colours) || (num_colours_in_binary > 0));

	if (binversion_in_binheader == 4)//legacy
	{
		int covg_reading_from_binary[num_colours_in_binary];
		Edges    individual_edges_reading_from_binary[num_colours_in_binary];
		int read;

		read = fread(kMer.kMerBinaryRepresentation.data(), sizeof(bitfield_of_64bits)*m, 1, fp);

		if (read>0)
		{

			read = fread(covg_reading_from_binary, sizeof(int), num_colours_in_binary, fp);
			if (read==0)
			{
				throw std::runtime_error("error with input file - failed to read covg in db_node_read_multicolour_binary");
			}

			read = fread(individual_edges_reading_from_binary, sizeof(Edges), num_colours_in_binary, fp);
			if (read==0)
			{
				throw std::runtime_error("error with input file - failed to read Edges in db_node_read_multicolour_binary\n");
			}
		}
		else
		{
			return false;
		}

		int i;
		for (i=0; i< num_colours_in_binary; i++)
		{
			node.coverage.at(i) = covg_reading_from_binary[i];
			node.individual_edges.at(i) = individual_edges_reading_from_binary[i];
		}
	}
	else//higher than 4
	{
		Covg covg_reading_from_binary[num_colours_in_binary];
		Edges individual_edges_reading_from_binary[num_colours_in_binary];

		int read;

		read = fread(kMer.kMerBinaryRepresentation.data(), sizeof(bitfield_of_64bits)*m, 1, fp);

		if (read > 0)
		{

			read = fread(covg_reading_from_binary, sizeof(uint32_t), num_colours_in_binary, fp);
			if (read==0)
			{
				throw std::runtime_error("Failed to read covg in db_node_read_multicolour_binary\n");
			}

			read = fread(individual_edges_reading_from_binary, sizeof(Edges), num_colours_in_binary, fp);
			if (read==0)
			{
				throw std::runtime_error("Failed to read Edges in db_node_read_multicolour_binary\n");
			}
		}
		else
		{
			return false;
		}

		int i;
		for (i=0; i < num_colours_in_binary; i++)
		{
			node.coverage.at(i)         = covg_reading_from_binary[i];
			node.individual_edges.at(i) = individual_edges_reading_from_binary[i];
		}
	}

	return true;
}


#endif /* DEBRUIJNGRAPH_H_ */
