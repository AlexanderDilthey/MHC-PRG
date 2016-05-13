/*
 * basic.h
 *
 *  Created on: 29.11.2012
 *      Author: AlexanderDilthey
 *  Interface and some code adapted from / inspired by Zamin Iqbal / cortex_var.*
 */

#ifndef BASIC_H_
#define BASIC_H_

#include <string>
#include <vector>

typedef enum
{
  forward = 0,
  reverse = 1
} Orientation;


typedef enum
{
    Adenine   = 0,
    Cytosine  = 1,
    Guanine   = 2,
    Thymine   = 3,
    Undefined = 4,
} Nucleotide;

char reverse_char_nucleotide(char c);
std::string seq_reverse_complement(std::string sequence);

Nucleotide char_to_binary_nucleotide(char c);
Nucleotide reverse_binary_nucleotide(Nucleotide n);
char binary_nucleotide_to_char(Nucleotide n);

int find_first_valid_kMer_position(std::string sequence, int k);
int find_first_valid_kMer_position(std::string sequence, int k, int startAtPosition);

std::vector<std::string> partitionStringIntokMers(std::string str, int k);

std::string kMer_canonical_representation(const std::string& kMer);
std::string kMer_canonical_representation(const std::string& kMer, bool& was_inverted);


#endif /* BASIC_H_ */
