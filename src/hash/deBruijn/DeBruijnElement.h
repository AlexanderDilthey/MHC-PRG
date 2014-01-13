/*
 * DeBruijnElement.h
 *
 *  Created on: 29.11.2012
 *      Author: AlexanderDilthey
 */

#ifndef DEBRUIJNELEMENT_H_
#define DEBRUIJNELEMENT_H_

#include <inttypes.h>
#include <array>
#include <set>
#include <assert.h>
#include "../sequence/basic.h"

typedef char Edges;

typedef uint32_t Covg;

typedef enum
{
	unassigned   = 0,
	none         = 1,
	visited      = 2,
	pruned       = 3,
	exists_in_reference = 4,
	visited_and_exists_in_reference = 5,
	to_be_dumped = 6, //to be dumped as binary
	read_start_forward = 7,//used when removing duplicate reads
	read_start_reverse = 8,//used when removing duplicate reads
	read_start_forward_and_reverse = 9,//used when removing duplicate reads
	ignore_this_node = 10,
	in_desired_genotype = 11,
	special_visited = 12,
	special_none = 13,
	special_pruned = 14
} NodeStatus;


template<int NUMBER_OF_COLOURS>
class DeBruijnElement {
public:
	  std::array<Covg, NUMBER_OF_COLOURS> coverage;
	  std::array<Edges, NUMBER_OF_COLOURS> individual_edges;


	  DeBruijnElement();

	  void nullify()
	  {
		  coverage.fill(0);
		  individual_edges.fill(0);
	  }


	  void setEdge(int colour, Nucleotide n, bool reverse)
	  {
		  assert(colour >= 0);
		  assert(colour < NUMBER_OF_COLOURS);

		  char edge = 1 << n;
		  if (reverse == true)
		  {
			  edge <<= 4; //move to next nibble
		  }

		  individual_edges.at(colour) |= edge;
	  }

	  void setEdge(Nucleotide n, bool reverse)
	  {
		  setEdge(0, n, reverse);
	  }

	  int numEdges(int colour, bool reverse)
	  {
		  return checkEdge(colour, Adenine, reverse)+checkEdge(colour, Cytosine, reverse)+checkEdge(colour, Guanine, reverse)+checkEdge(colour, Thymine, reverse);
	  }

	  int numEdges(bool reverse)
	  {
		  return numEdges(0, reverse);
	  }

	  std::set<Nucleotide> getEdges(bool reverse)
	  {
		  return getEdges(0, reverse);
	  }

	  std::set<Nucleotide> getEdges(int colour, bool reverse)
	  {
		  Edges eInfo = individual_edges.at(colour);
		  if(reverse == true)
		  {
			  eInfo >>= 4;
		  }

		  std::set<Nucleotide> forReturn;
		  char mask = 1;
		  for(int i = 0; i <= 3; i++)
		  {
			  if((eInfo & mask) != 0)
			  {
				  forReturn.insert((Nucleotide)i);
			  }
			  mask = mask << 1;
		  }

		  return forReturn;
	  }

	  bool checkEdge(int colour, Nucleotide n, bool reverse)
	  {
		  assert(colour >= 0);
		  assert(colour < NUMBER_OF_COLOURS);

		  char edge = 1 << n;
		  if (reverse == true)
		  {
			  edge <<= 4; //move to next nibble
		  }

		  return ((individual_edges.at(colour) & edge) != 0);
	  }



	  bool checkEdge(Nucleotide n, bool reverse)
	  {
		  return checkEdge(0, n, reverse);
	  }

	  void removeEdge(int colour, Nucleotide n, bool reverse)
	  {
		  assert(colour >= 0);
		  assert(colour < NUMBER_OF_COLOURS);

		  char edge = 1 << n;

		  if (reverse == true)
		  {
			  edge <<= 4;
		  }

		  // toggle 1->0 0->1
		  // xor with all 1's, ie 00010000 -> 11101111
		  edge ^= (unsigned char) 0xFF;

		  // reset one edge
		  individual_edges.at(colour) &= edge;
	  }

	  void removeEdge(Nucleotide n, bool reverse)
	  {
		  removeEdge(0, n, reverse);
	  }
};

template<int NUMBER_OF_COLOURS>
DeBruijnElement<NUMBER_OF_COLOURS>::DeBruijnElement()
{
	for(int i = 0; i < NUMBER_OF_COLOURS; i++)
	{
		coverage[i] = 0;
		individual_edges[i] = 0;
	}
}

#endif /* DEBRUIJNELEMENT_H_ */
