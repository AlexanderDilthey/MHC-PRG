/*
 * NextGen.cpp
 *
 *  Created on: 5 Jan 2012
 *      Author: dilthey
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include "assert.h"
#include "NextGen.h"
#include <stdlib.h>
#include "../Data/HaplotypePanel.h"
#include "../Data/GenotypePanel.h"
#include "../Graph/Graph.h"
#include "../Graph/Node.h"
#include "../Graph/HMM.h"
#include "../MHC-PRG.h"
#include "../Utilities.h"
#include "../Graph/LargeGraph.h"
#include "../Graph/MultiGraph.h"
#include "../Graph/AlphaHMM.h"
#include <iomanip>

vector<kMerInfo> forwardScanRec(Node* currentNode, int depth, int realdepth, int limit, int firstEdgeGap);
map<string, string> readHLAalleleAlignments(string file, vector<int>& fourNumbers);
string alleleAndLocusIdentifier(string locus, string allele);
int pointerStrintLength = -1;

Graph* variationGraph(string input_panel, string positions_file, bool wantPGFprotection)
{
	cout << "Read haplotypes panel...\n" << flush;

	HaplotypePanel hp;
	hp.readFromFile(input_panel, positions_file);

	cout << "\tdone...\n" << flush;

	Graph* g = new Graph();
	g->CODE = hp.CODE;

	Node* n0 = new Node();
	n0->level = 0;
	n0->terminal = false;
	g->registerNode(n0, 0);

	vector<int> realHaplotypeIndices;
	map<int, set<int> > snpHaplotypeIndices;

	map<Node*, set<int> > NodeToHaplotype;

	int last_real_haplotype_index = -1;
	int pgf_haplotype_id = -1;
	for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
	{
		string haplotypeID = hp.HaplotypeIDs.at(haplotypeI);
		if(haplotypeID.substr(0,4) == "SNPs")
		{
			assert(last_real_haplotype_index != -1);
			snpHaplotypeIndices[last_real_haplotype_index].insert(haplotypeI);
		}
		else
		{
			realHaplotypeIndices.push_back(haplotypeI);
			NodeToHaplotype[n0].insert(haplotypeI);
			last_real_haplotype_index = haplotypeI;
		}

		if((haplotypeID == "pgf_____1") || ((haplotypeID.length() > 10) && (haplotypeID.substr(0, 10) == "pgfallele_")))
		{
			pgf_haplotype_id = haplotypeI;
			cout << "PGF haplotype index: " << pgf_haplotype_id << "\n";
		}
	}

	if(wantPGFprotection)
	{
		if(pgf_haplotype_id == -1)
		{
			cout << "Haplotype IDs:\n";
			for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
			{
				string haplotypeID = hp.HaplotypeIDs.at(haplotypeI);
				cout << haplotypeID << "\n";
			}
		}
		assert(pgf_haplotype_id != -1);
	}

	vector<string> loci = hp.getLoci();
	positionsSorter sortClass;
	sortClass.p = &(hp.LocusPositions);
	sort(loci.begin(), loci.end(), sortClass);

	string lastLocus = loci.at(loci.size()-1);
	int lastLocusPosition = hp.LocusPositions[lastLocus];

	string pufferLocusID = "END_PUFFER";
	int pufferLocusLocation = lastLocusPosition+1;
	loci.push_back(pufferLocusID);
	hp.LocusPositions[pufferLocusID] = pufferLocusLocation;
	hp.LocusStrands[pufferLocusID] = "+";
	unsigned char puffer_symbol = hp.CODE.doCode(pufferLocusID, "*");
	g->CODE = hp.CODE;
	for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
	{
		hp.HaplotypesByLoci[pufferLocusID].push_back(puffer_symbol);
	}

	for(int level = 0; level < (int)loci.size(); level++)
	{

		if(level > 5)
		{
			//exit(0);
		}
		
		if((level % 100000) == 0)
		{
		}

		if((level % 1000) == 0)
		{
			cout << "\rLevel " << level << flush;
		}

		string locusID = loci.at(level);
		int star_symbol = g->CODE.doCode(locusID, "*");

		// cout << "\n";
		// cout << "\n\nLevel " << level << " locusID " << locusID << flush;

		map<Edge*, set<int> > looseEdges;
		map<Edge*, set<Edge*> > snpEdges;
		set<Edge*> snpEdges_to_catch;

		set<int> seen_haplotypes;

		bool level_protected_PGF = false;
		Edge* pgfEdge = 0;

		
		for(map< Node*, set<int> >::iterator nIt = NodeToHaplotype.begin(); nIt != NodeToHaplotype.end(); nIt++)
		{
			// cout << "\tNode " << nIt->first << "\n";

			map<unsigned char, set<int> > emission2Haplotypes;
			map<unsigned char, set<unsigned char> > emissions_possible_SNPs;

			for(set<int>::iterator hIt = nIt->second.begin(); hIt != nIt->second.end(); hIt++)
			{
				int haplotype_index = *hIt;

				assert(hp.HaplotypesByLoci.count(locusID) > 0);

				unsigned char haplo_symbol = hp.HaplotypesByLoci[locusID].at(haplotype_index);
				emission2Haplotypes[haplo_symbol].insert(haplotype_index);

				for(set<int>::iterator snpHaplotypesIndicesIt = snpHaplotypeIndices[haplotype_index].begin(); snpHaplotypesIndicesIt != snpHaplotypeIndices[haplotype_index].end(); snpHaplotypesIndicesIt++)
				{
					unsigned char haplo_symbol_possibleSNP = hp.HaplotypesByLoci[locusID].at(*snpHaplotypesIndicesIt);

					if(haplo_symbol_possibleSNP != star_symbol)
					{
						emissions_possible_SNPs[haplo_symbol].insert(haplo_symbol_possibleSNP);
						assert(haplo_symbol != star_symbol);
					}
				}

				seen_haplotypes.insert(haplotype_index);

				// cout << "\t\t attached haplotype  " << hp.HaplotypeIDs.at(haplotype_index) << "char: " << haplo_symbol << " real symbol: " << g->CODE.deCode(locusID, haplo_symbol) << "\n";
			}

			/*
			if(emission2Haplotypes.count(star_symbol) > 0)
			{
				if(emission2Haplotypes.size() == 2)
				{
					set<int> haplosWithStar = emission2Haplotypes[star_symbol];
					emission2Haplotypes.erase(star_symbol);
					unsigned char remainingNonStar = (emission2Haplotypes.begin()->first);
					emission2Haplotypes[remainingNonStar].insert(haplosWithStar.begin(), haplosWithStar.end());
				}
			}
			*/

			for(map<unsigned char, set<int> >::iterator charIt = emission2Haplotypes.begin(); charIt != emission2Haplotypes.end(); charIt++)
			{
				assert(charIt->second.size() > 0);
				unsigned char symbol = charIt->first;

				//cout << "level " << level << " node " << nIt->first << " char " << symbol << " value " << g->CODE.deCode(locusID, symbol) << "\n";
				Edge* newE = new Edge();
						
				newE->count = 1;
				newE->emission = symbol;
				newE->locus_id = locusID;
				if((locusID.substr(0,2) == "rs") || (locusID.substr(0,3) == "HLA"))
				{
					newE->label = locusID+"*"+hp.CODE.deCode(locusID, symbol);
				}

				newE->pgf_protect = false;
			
				if((pgf_haplotype_id != -1) && (charIt->second.count(pgf_haplotype_id) > 0))
				{
					newE->pgf_protect = true;
					level_protected_PGF = true;

					/*
					if(level > 0)
					{
						bool found_preceding_PGF = false;
						Node* fromNode = nIt->first;
						for(set<Edge*>::iterator InEdgeIt = fromNode->Incoming_Edges.begin(); InEdgeIt != fromNode->Incoming_Edges.end(); InEdgeIt++)
						{
							Edge *preEdge = *InEdgeIt;
							if((bool)preEdge->pgf_protect)
							{
								found_preceding_PGF = true;
							}
						}
						assert(found_preceding_PGF);
					}
					*/

					assert(pgfEdge == 0);
					pgfEdge = newE;
				}

				g->registerEdge(newE);
				newE->From = nIt->first;
				nIt->first->Outgoing_Edges.insert(newE);

				looseEdges[newE] = charIt->second;

				if(emissions_possible_SNPs.count(symbol) > 0)
				{
					assert(symbol != star_symbol);
					for(set<unsigned char>::iterator possibleSNPit = emissions_possible_SNPs[symbol].begin(); possibleSNPit != emissions_possible_SNPs[symbol].end(); possibleSNPit++ )
					{
						unsigned char alternativeEdgeSymbol = *possibleSNPit;
						if(alternativeEdgeSymbol != symbol)
						{
							Edge* newE_SNP = new Edge();
										
							newE_SNP->count = 1;
							newE_SNP->emission = alternativeEdgeSymbol;
							newE_SNP->locus_id = locusID;
							newE_SNP->pgf_protect = false;
							
							if((locusID.substr(0,2) == "rs") || (locusID.substr(0,3) == "HLA"))
							{
								newE_SNP->label = locusID+"*"+hp.CODE.deCode(locusID, alternativeEdgeSymbol);
							}
							g->registerEdge(newE_SNP);
							newE_SNP->From = nIt->first;
							nIt->first->Outgoing_Edges.insert(newE_SNP);

							snpEdges[newE].insert(newE_SNP);

							snpEdges_to_catch.insert(newE_SNP);
						}
					}
				}
			}
			
			
		}

		// PGF protection
		assert((! wantPGFprotection) || level_protected_PGF);

		vector< set<int> > groupingHaplotypes;
		vector< set<Edge*> > groupingEdges;
		
		for(map<Edge*, set<int> >::iterator looseEdgeIt = looseEdges.begin(); looseEdgeIt != looseEdges.end(); looseEdgeIt++)
		{
			assert(looseEdgeIt->second.size() > 0);
			groupingHaplotypes.push_back(looseEdgeIt->second);
			set<Edge*> edgeSet;
			edgeSet.insert(looseEdgeIt->first);
			groupingEdges.push_back(edgeSet);
		}

		vector< vector<int> > pairsToCheck;
		for(int i = 0; i < (int)groupingHaplotypes.size(); i++)
		{
			for(int j = i+1; j < (int)groupingHaplotypes.size(); j++)
			{
				vector<int> p;
				p.push_back(i);
				p.push_back(j);
				pairsToCheck.push_back(p);
			}
		}
		set<int> groupsDeleted;

		int want_suffix_length = 20;

		while(pairsToCheck.size() != 0)
		{
			vector<int> thisPair = pairsToCheck.at(0);
			assert(thisPair.size() == 2);

			pairsToCheck.erase(pairsToCheck.begin());

			int g1 = thisPair.at(0);
			int g2 = thisPair.at(1);
			if((groupsDeleted.count(g1) > 0) || (groupsDeleted.count(g2) > 0))
			{
				continue;
			}
			else
			{

				bool joinPair = true;

				int local_need_suffix_length = want_suffix_length;
				int got_suffix_length = 0;

				while((joinPair == true) && (want_suffix_length != got_suffix_length))
				{
					//cout << "want_suffix_length: " << want_suffix_length << " local_need_suffix_length: " << local_need_suffix_length << " got_suffix_length: " << got_suffix_length << "\n";

					map< basic_string<unsigned char>, int> suffixes_g1;
					map< basic_string<unsigned char>, int> suffixes_g2;
					set< basic_string<unsigned char > > suffixes;

					if(local_need_suffix_length > (want_suffix_length*10))
					{
						joinPair = false;
						continue;
					}

					if((level + local_need_suffix_length) > (loci.size() - 1))
					{
						joinPair = false;
					}
					else
					{
						vector<string> takeSuffixLoci;

						for(int lI = 1; lI <= local_need_suffix_length; lI++)
						{
							int suffixL = level + lI;
							takeSuffixLoci.push_back(loci.at(suffixL));
						}

						vector<unsigned char> takeSuffixLociStars;
						for(vector<string>::iterator locusIt = takeSuffixLoci.begin(); locusIt != takeSuffixLoci.end(); locusIt++)
						{
							string locus = *locusIt;
							unsigned char starSymbol = g->CODE.doCode(locus, "*");
							takeSuffixLociStars.push_back(starSymbol);
						}

						for(set<int>::iterator hIt = groupingHaplotypes.at(g1).begin(); hIt != groupingHaplotypes.at(g1).end(); hIt++)
						{
							// normal haplotype

							basic_string<unsigned char> suffix = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*hIt), false);
							assert((int)suffix.length() == local_need_suffix_length);
							if(suffixes_g1.count(suffix) == 0)
							{
								suffixes_g1[suffix] = 0;
							}
							suffixes_g1[suffix]++;
							suffixes.insert(suffix);

							// SNPs

							for(set<int>::iterator snpHaplotypeIt = snpHaplotypeIndices[*hIt].begin(); snpHaplotypeIt != snpHaplotypeIndices[*hIt].end(); snpHaplotypeIt++)
							{
								basic_string<unsigned char> suffix_2 = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*snpHaplotypeIt), false);
								assert((int)suffix_2.length() == local_need_suffix_length);
								for(int pos = 0; pos < (int)suffix_2.length(); pos++)
								{
									if(suffix_2.at(pos) == takeSuffixLociStars.at(pos))
									{
										suffix_2.at(pos) = suffix.at(pos);
									}
								}
								if(suffixes_g1.count(suffix_2) == 0)
								{
									suffixes_g1[suffix_2] = 0;
								}
								suffixes_g1[suffix_2]++;
								suffixes.insert(suffix_2);
							}

						}

						for(set<int>::iterator hIt = groupingHaplotypes.at(g2).begin(); hIt != groupingHaplotypes.at(g2).end(); hIt++)
						{
							// normal haplotypes

							basic_string<unsigned char> suffix = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*hIt), false);
							assert((int)suffix.length() == local_need_suffix_length);

							if(suffixes_g2.count(suffix) == 0)
							{
								suffixes_g2[suffix] = 0;
							}
							suffixes_g2[suffix]++;
							suffixes.insert(suffix);

							// SNPs

							for(set<int>::iterator snpHaplotypeIt = snpHaplotypeIndices[*hIt].begin(); snpHaplotypeIt != snpHaplotypeIndices[*hIt].end(); snpHaplotypeIt++)
							{
								basic_string<unsigned char> suffix_2 = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*snpHaplotypeIt), false);
								assert((int)suffix_2.length() == local_need_suffix_length);

								for(int pos = 0; pos < (int)suffix_2.length(); pos++)
								{
									if(suffix_2.at(pos) == takeSuffixLociStars.at(pos))
									{
										suffix_2.at(pos) = suffix.at(pos);
									}
								}

								if(suffixes_g2.count(suffix_2) == 0)
								{
									suffixes_g2[suffix_2] = 0;
								}
								suffixes_g2[suffix_2]++;
								suffixes.insert(suffix_2);
							}
						}

						got_suffix_length =  local_need_suffix_length;


						bool one_node_all_stars = false;
						for(set< basic_string<unsigned char> >::iterator suffixIt = suffixes.begin(); suffixIt != suffixes.end(); suffixIt++)
						{
							assert((int)suffixIt->size() == local_need_suffix_length);
							bool all_stars = true;
							for(int i = 0; i < local_need_suffix_length; i++)
							{
								unsigned char suffixC = suffixIt->at(i);
								string suffixCharLocus = takeSuffixLoci.at(i);
								string decodedSuffix = g->CODE.deCode(suffixCharLocus, suffixC);
								all_stars = (all_stars && (decodedSuffix == "*"));
							}

							if(all_stars)
							{
								if(((suffixes_g1.count(*suffixIt) > 0) && (suffixes_g1.size() == 1)) || ((suffixes_g2.count(*suffixIt) > 0) && (suffixes_g2.size() == 1)))
								{
									one_node_all_stars = true;
								}
							}
						}


						for(set< basic_string<unsigned char> >::iterator suffixIt = suffixes.begin(); suffixIt != suffixes.end(); suffixIt++)
						{
							int gapCounter = 0;
							assert((int)suffixIt->size() == local_need_suffix_length);
							bool all_stars = true;
							for(int i = 0; i < local_need_suffix_length; i++)
							{
								unsigned char suffixC = suffixIt->at(i);
								string suffixCharLocus = takeSuffixLoci.at(i);
								string decodedSuffix = g->CODE.deCode(suffixCharLocus, suffixC);
								all_stars = (all_stars && (decodedSuffix == "*"));

								if(decodedSuffix == "_")
								{
									gapCounter++;
									if(i == 0)
									{
										joinPair = false;
									}
								}
							}

							if((local_need_suffix_length - gapCounter) < got_suffix_length)
							{
								got_suffix_length = (local_need_suffix_length - gapCounter);
							}

							if(! one_node_all_stars)
							{
								if(! all_stars)
								{
									if(!((suffixes_g1.count(*suffixIt) > 0) && (suffixes_g2.count(*suffixIt) > 0)))
									{
										joinPair = false;
									}
								}
							}
						}
						assert(got_suffix_length <= local_need_suffix_length);
						local_need_suffix_length = want_suffix_length + (local_need_suffix_length - got_suffix_length);
					}
				}

				if(joinPair)
				{
					set<int> newGroupHaplos;
					set<Edge*> newGroupEdges;

					newGroupHaplos.insert(groupingHaplotypes.at(g1).begin(), groupingHaplotypes.at(g1).end());
					newGroupHaplos.insert(groupingHaplotypes.at(g2).begin(), groupingHaplotypes.at(g2).end());
					groupingHaplotypes.at(g1).clear();
					groupingHaplotypes.at(g2).clear();

					newGroupEdges.insert(groupingEdges.at(g1).begin(), groupingEdges.at(g1).end());
					newGroupEdges.insert(groupingEdges.at(g2).begin(), groupingEdges.at(g2).end());
					groupingEdges.at(g1).clear();
					groupingEdges.at(g2).clear();

					groupingHaplotypes.push_back(newGroupHaplos);
					groupingEdges.push_back(newGroupEdges);
					assert(groupingHaplotypes.size() == groupingEdges.size());

					groupsDeleted.insert(g1);
					groupsDeleted.insert(g2);

					int index_of_new_element = groupingHaplotypes.size()-1;
					assert(groupingHaplotypes.at(index_of_new_element).size() > 0);
					for(int i = 0; i < index_of_new_element; i++)
					{
						if(groupsDeleted.count(i) > 0)
						{
							continue;
						}
						else
						{
							vector<int> newP;
							newP.push_back(i);
							newP.push_back(index_of_new_element);
							pairsToCheck.push_back(newP);
						}
					}
				}
			}
		}

		for(unsigned int haplotypeII = 0; haplotypeII < realHaplotypeIndices.size(); haplotypeII++)
		{
			int haplotypeI = realHaplotypeIndices.at(haplotypeII);
			if(!(seen_haplotypes.count(haplotypeI) > 0))
			{
				cout << "level: " << level << ", haplotypeI: " << haplotypeI << ", seen_haplotypes.size(): " << seen_haplotypes.size() << ", hp.HaplotypeIDs.size(): " << hp.HaplotypeIDs.size() << "\n";
				cout << flush;
			}
			assert(seen_haplotypes.count(haplotypeI) > 0);
		}

		NodeToHaplotype.clear();
		
		int addedNodeCounter = 0;
		for(int i = 0; i < (int)groupingEdges.size(); i++)
		{
			if(groupsDeleted.count(i) > 0)
			{
				continue;
			}
			else
			{
				addedNodeCounter++;

				Node* newN = new Node();
				if(level == (loci.size()-1))
				{
					newN->terminal = true;
				}
				else
				{
					newN->terminal = false;
				}

				newN->level = level + 1;
				g->registerNode(newN, newN->level);

				for(set<Edge*>::iterator eIt = groupingEdges.at(i).begin(); eIt != groupingEdges.at(i).end(); eIt++)
				{
					newN->Incoming_Edges.insert(*eIt);
					(*eIt)->To = newN;

					if(snpEdges.count(*eIt) > 0)
					{
						for(set<Edge*>::iterator snpEdgeIt = snpEdges[*eIt].begin(); snpEdgeIt != snpEdges[*eIt].end(); snpEdgeIt++)
						{
							Edge* snpEdge = *snpEdgeIt;
							newN->Incoming_Edges.insert(snpEdge);
							snpEdge->To = newN;
							snpEdges_to_catch.erase(snpEdge);
						}
					}
				}


				assert(groupingHaplotypes.size() > 0);
				NodeToHaplotype[newN] = groupingHaplotypes.at(i);
				// cout << "New node memory address: " << newN << "\n";
				for(set<int>::iterator hIndex = groupingHaplotypes.at(i).begin(); hIndex != groupingHaplotypes.at(i).end(); hIndex++)
				{
					// cout << "\t associated haplotype index : " << *hIndex << "\n";
				}
				// cout << flush;
			}
		}

		assert(snpEdges_to_catch.size() == 0);

	}

	g->removeStarPaths();
	
	if(wantPGFprotection)
	{
		for(unsigned int level = 0; level < (g->NodesPerLevel.size()-1); level++)
		{
			bool foundPGF = false;
			for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(level).begin(); nodeIt != g->NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;

				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					foundPGF = (foundPGF || (bool)e->pgf_protect);
					if((level > 0) && (level < (g->NodesPerLevel.size()-2)))
					{
						if((bool)e->pgf_protect)
						{
							bool found_preceding_PGF = false;

							for(set<Edge*>::iterator beforeIt = e->From->Incoming_Edges.begin(); beforeIt != e->From->Incoming_Edges.end(); beforeIt++)
							{
								if((bool)(*beforeIt)->pgf_protect)
								{
									found_preceding_PGF = true;
									break;
								}
							}
							assert(found_preceding_PGF);

							bool found_following_PGF = false;						
							for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
							{
								if((bool)(*afterIt)->pgf_protect)
								{
									found_following_PGF = true;
									break;
								}
							}

							if(! found_following_PGF)
							{
								cout << "Level " << level << " of " << (g->NodesPerLevel.size()-2) << ", edge " << e << " pgf protect: " << (bool)e->pgf_protect << " starting at level " << e->From->level << "\n";
								for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
								{
									cout << "\t following " << (*afterIt) << ": " << (bool)(*afterIt)->pgf_protect << "\n" << flush;
									assert( e->To->level == (level + 1));
								}
								
								cout << "We now look at all edges from the this level " << (level) << ":\n";
								for(set<Node*>::iterator nodeIt2 = g->NodesPerLevel.at(level).begin(); nodeIt2 != g->NodesPerLevel.at(level).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}
								
								cout << "We now look at all edges from the next level " << (level+1) << ":\n";
								for(set<Node*>::iterator nodeIt2 = g->NodesPerLevel.at(level+1).begin(); nodeIt2 != g->NodesPerLevel.at(level+1).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}
							}
							
							assert(found_following_PGF);
						}
					}
				}
			}
			assert(foundPGF);
		}
	}

	for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(g->NodesPerLevel.size()-1).begin(); nodeIt != g->NodesPerLevel.at(g->NodesPerLevel.size()-1).end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			string emission = g->CODE.deCode(edge->locus_id, edge->emission);
			assert(emission == "*");
		}
	}
	
	for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(g->NodesPerLevel.size()-1).begin(); nodeIt != g->NodesPerLevel.at(g->NodesPerLevel.size()-1).end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			string emission = g->CODE.deCode(edge->locus_id, edge->emission);
			assert(emission != "_");
		}
	}
	
	cout << "\n\nBuilding nucleotide graph done!\n\n" << flush;



	return g;

}

Graph* HLAalignmentGraph(string hla_allele_dir, string locus)
{
	string expected_alignment_file = hla_allele_dir+"/"+locus+" short alleles.txt";
	vector<int> fourNumbers;
	map<string, string> alleles_alignment = readHLAalleleAlignments(expected_alignment_file, fourNumbers);

	cout << "Build HLA allele graph...\n";
	assert(fourNumbers.size()==2);

	Graph* g = new Graph();
	Node* n0 = new Node();
	n0->level = 0;
	n0->terminal = false;
	g->registerNode(n0, 0);

	for(map<string, string>::iterator allelesIt = alleles_alignment.begin(); allelesIt != alleles_alignment.end(); allelesIt++)
	{
		string allele = allelesIt->first;
		string alleleSequence = allelesIt->second;

		Node* toAppendTo = n0;
		for(int i = 0; i < (int)alleleSequence.length(); i++)
		{
			string locusID;
			if(i == 0)
			{
				locusID = "chr6:"+Utilities::ItoStr(fourNumbers.at(0));
			}
			else if(i == (alleleSequence.length()-1))
			{
				locusID = "chr6:"+Utilities::ItoStr(fourNumbers.at(1));
			}
			else
			{
				locusID = "A"+Utilities::ItoStr(i);
			}

			string c = alleleSequence.substr(i, 1);

			Edge* newE = new Edge();
			newE->count = 1;
			newE->emission = g->CODE.doCode(locusID, c);
			newE->locus_id = locusID;
			g->registerEdge(newE);
			if(i == 0)
			{
				newE->label = "HLA"+allele;
			}

			toAppendTo->Outgoing_Edges.insert(newE);
			newE->From = toAppendTo;

			Node* newN = new Node();
			if(i == (int)(alleleSequence.npos - 1))
			{
				newN->terminal = true;
			}
			else
			{
				newN->terminal = false;
			}
			newN->level = i+1;
			g->registerNode(newN, i+1);

			newN->Incoming_Edges.insert(newE);
			newE->To = newN;

			toAppendTo = newN;

		}
	}

	cout << "\t\t .. done!\n" << flush;

	return g;
}

MultiGraph* simplifyAccordingToCoverage(MultiGraph* mG, map<string, long long> estimatedEmissions, set<int> protectLevels)
{
	set<string> utilizekMers;
	set<string> ignorekMers;
	kMerPositionInfo kMersInGraphInfo;
	set<string> kMersInGraph;
	kMerUniquenessInfo uniqueMers;

	kMersInGraphInfo = mG->getkMerPositions();
	kMersInGraph = kMersInGraphInfo.kMers;
	uniqueMers = mG->kMerUniqueness();

	// find out which kMers we want to keep
	for(set<string>::iterator kMerIt = kMersInGraph.begin(); kMerIt != kMersInGraph.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		bool usekMer = true;

		// non-unique within graph
		if(uniqueMers.levelUniqueKMers.count(kMer) == 0)
		{
			usekMer = false;
		}

		assert(kMersInGraphInfo.kMerEdgeNum.count(kMer) > 0);
		if(kMersInGraphInfo.kMerEdgeNum[kMer] != 1)
		{
			usekMer = false;
		}

		// non-unique within genome
		assert(estimatedEmissions.count(kMer) > 0);
		if(estimatedEmissions[kMer] == -1)
		{
			usekMer = false;
		}

		if(usekMer)
		{
			utilizekMers.insert(kMer);
		}
		else
		{
			ignorekMers.insert(kMer);
		}
	}

	cout << "Graph simplification: out of " << kMersInGraph.size() << " kMers, we are now using " << utilizekMers.size() << "\n";
	cout << "\t\t Number of protected levels: " << protectLevels.size() << "\n";

	set<Edge*> deleteEdge;
	map<Edge*, set<int> > deleteEdgeCausalLevels;
	set<Node*> deleteNode;

	vector<int> remainingEdgesPerLevel;
	vector<vector<Edge*> > edgesFromLevel;

	int levels = mG->NodesPerLevel.size();
	remainingEdgesPerLevel.resize(levels - 1);
	edgesFromLevel.resize(levels - 1);

	for(int level = 0; level < levels - 1; level++)
	{
		set<Node*> nodesAtLevel = mG->NodesPerLevel.at(level);
		set<Edge*> edgesFromLevel_set;

		for(set<Node*>::iterator nIt = nodesAtLevel.begin(); nIt != nodesAtLevel.end(); nIt++)
		{
			Node* n = *nIt;
			for(set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				edgesFromLevel_set.insert(*eIt);
			}
		}

		vector<Edge*> edgesFromLevel_thisLevel(edgesFromLevel_set.begin(), edgesFromLevel_set.end());

		edgesFromLevel.at(level) = edgesFromLevel_thisLevel;
		remainingEdgesPerLevel.at(level) = edgesFromLevel_thisLevel.size();
	}

	for(int level = 0; level < levels - 1; level++)
	{
		if(protectLevels.count(level) > 0)
		{
			continue;
		}

		int edgeIMax = edgesFromLevel.at(level).size();

		if(edgeIMax < 10)
		{
			continue;
		}

		for(int edgeI = 0; edgeI < edgeIMax; edgeI++)
		{
			Edge* e = edgesFromLevel.at(level).at(edgeI);

			if(deleteEdge.count(e) > 0)
			{
				continue;
			}

			map<int, int> edgeEmissions = e->multiEmission;

			int kMer_present = 0;
			int kMer_total = 0;

			int kMer_present_all = 0;
			int kMer_total_all = 0;

			for(map<int, int>::iterator emissionIt = edgeEmissions.begin(); emissionIt != edgeEmissions.end(); emissionIt++)
			{
				string kMer = mG->CODE.deCode(e->locus_id, emissionIt->first);
				if((kMer == "*") || (kMer == "_"))
				{
					continue;
				}

				assert(estimatedEmissions.count(kMer) > 0);
				if((estimatedEmissions[kMer] > 0) || (estimatedEmissions[kMer] == -1))
				{
					kMer_present_all++;
				}
				kMer_total_all++;

				if(utilizekMers.count(kMer) > 0)
				{
					assert(estimatedEmissions[kMer] != -1);
					if(estimatedEmissions[kMer] > 0)
					{
						kMer_present++;
					}
					kMer_total++;
				}
			}

			double fraction_missing = 0;
			if(kMer_total > 0)
			{
				fraction_missing = 1 - (double)kMer_present/(double)kMer_total;
			}

			double fraction_missing_all = 0;
			if(kMer_total_all > 0)
			{
				fraction_missing_all = 1 - (double)kMer_present_all/(double)kMer_total_all;
			}

			if(edgeIMax > 66050)
			{
				//cout << "level " << level << " edge " << edgeI << " misssing: " << fraction_missing << " missing_all: " << fraction_missing_all << " total kMers " << kMer_total << " total kMers all: " << kMer_total_all << "\n";
			}
			if(( (fraction_missing_all >= 0.5) || ((kMer_total > 5) && (fraction_missing >= 0.5)) || ((kMer_total_all > 5) && (fraction_missing_all >= 0.5))) && (protectLevels.count(level) == 0))
			{
				// cout << "Level " << level << " edge " << edgeI << "/" << edgeIMax << " missing: " << fraction_missing << " ( ie " << kMer_present << " out of " <<  kMer_total << " kmers are present)\n";

				if(! (bool)e->pgf_protect)
				{
					assert(deleteEdge.count(e) == 0);
					deleteEdge.insert(e);
					deleteEdgeCausalLevels[e].insert(level);

					remainingEdgesPerLevel.at(level)--;

					// cout << "\t remaining after: " << remainingEdgesPerLevel.at(level) << "\n";

					if(remainingEdgesPerLevel.at(level) == 0)
					{
						assert(1 == 0); // this should not happen because of PGF protection

						set<int> nextRoundProtectLevel(protectLevels.begin(), protectLevels.end());
						for(int levelEdgeI = 0; levelEdgeI < (int)edgesFromLevel.at(level).size(); levelEdgeI++)
						{
							Edge* deletedEdge = edgesFromLevel.at(level).at(levelEdgeI);
							assert(deleteEdge.count(deletedEdge) > 0);
							assert(deleteEdgeCausalLevels.count(deletedEdge) > 0);
							set<int> causallyRelatedLevels = deleteEdgeCausalLevels[deletedEdge];
							nextRoundProtectLevel.insert(causallyRelatedLevels.begin(), causallyRelatedLevels.end());
						}

						// cout << "\t level " << level << " 0 edges, go recursive! [1]\n";

						return simplifyAccordingToCoverage(mG, estimatedEmissions, nextRoundProtectLevel);
					}

					set<Node*> nodeDeletionCheck;
					nodeDeletionCheck.insert(e->From);
					nodeDeletionCheck.insert(e->To);


					while(nodeDeletionCheck.size() > 0)
					{
						Node* thisNode = *(nodeDeletionCheck.begin());

						bool allIncomingDeleted = (thisNode->Incoming_Edges.size() > 0);
						for(set<Edge*>::iterator eIt = thisNode->Incoming_Edges.begin(); eIt != thisNode->Incoming_Edges.end(); eIt++)
						{
							Edge* checkEdge = *eIt;
							allIncomingDeleted = (allIncomingDeleted && (deleteEdge.count(checkEdge) > 0));
						}

						bool allOutgoingDeleted = (thisNode->Outgoing_Edges.size() > 0);
						for(set<Edge*>::iterator eIt = thisNode->Outgoing_Edges.begin(); eIt != thisNode->Outgoing_Edges.end(); eIt++)
						{
							Edge* checkEdge = *eIt;
							allOutgoingDeleted = (allOutgoingDeleted && (deleteEdge.count(checkEdge) > 0));
						}

						// cout << "Check a node at level " << thisNode->level << "\n";
						// cout << "\t allIncomingDeleted: " << allIncomingDeleted << " // allOutgoingDeleted: " << allOutgoingDeleted << "\n";
						if(allIncomingDeleted || allOutgoingDeleted)
						{

							deleteNode.insert(thisNode);

							set<Edge*> fromNodeDeleteEdges;

							set<int> deletionCausality;
							for(set<Edge*>::iterator eIt = thisNode->Incoming_Edges.begin(); eIt != thisNode->Incoming_Edges.end(); eIt++)
							{
								fromNodeDeleteEdges.insert(*eIt);
								if(deleteEdgeCausalLevels.count(*eIt) > 0)
								{
									deletionCausality.insert(deleteEdgeCausalLevels[*eIt].begin(), deleteEdgeCausalLevels[*eIt].end());
								}
							}

							for(set<Edge*>::iterator eIt = thisNode->Outgoing_Edges.begin(); eIt != thisNode->Outgoing_Edges.end(); eIt++)
							{
								fromNodeDeleteEdges.insert(*eIt);
								if(deleteEdgeCausalLevels.count(*eIt) > 0)
								{
									deletionCausality.insert(deleteEdgeCausalLevels[*eIt].begin(), deleteEdgeCausalLevels[*eIt].end());
								}
							}

							for(set<Edge*>::iterator eIt = fromNodeDeleteEdges.begin(); eIt != fromNodeDeleteEdges.end(); eIt++)
							{
								if(deleteEdge.count(*eIt) == 0)
								{
									// cout << "delete edge at level " << (*eIt)->From->level << "\n";

									deleteEdge.insert(*eIt);
									deleteEdgeCausalLevels[*eIt].insert(deletionCausality.begin(), deletionCausality.end());
									int edgeLevel = (*eIt)->From->level;
									remainingEdgesPerLevel.at(edgeLevel)--;

									if(remainingEdgesPerLevel.at(level) == 0)
									{
										set<int> nextRoundProtectLevel(protectLevels.begin(), protectLevels.end());
										for(int levelEdgeI = 0; levelEdgeI < (int)edgesFromLevel.at(level).size(); levelEdgeI++)
										{
											Edge* deletedEdge = edgesFromLevel.at(level).at(levelEdgeI);
											assert(deleteEdge.count(deletedEdge) > 0);
											assert(deleteEdgeCausalLevels.count(deletedEdge) > 0);
											set<int> causallyRelatedLevels = deleteEdgeCausalLevels[deletedEdge];
											nextRoundProtectLevel.insert(causallyRelatedLevels.begin(), causallyRelatedLevels.end());
										}

										return simplifyAccordingToCoverage(mG, estimatedEmissions, nextRoundProtectLevel);
									}

									if(deleteNode.count((*eIt)->To) == 0)
									{
										nodeDeletionCheck.insert((*eIt)->To);
									}

									if(deleteNode.count((*eIt)->From) == 0)
									{
										nodeDeletionCheck.insert((*eIt)->From);
									}
								}
							}
						}

						nodeDeletionCheck.erase(thisNode);
					}
				}

			}
		}
	}

	cout << "Will kick out " << deleteEdge.size() << " of " << mG->Edges.size() << " edges \n";

	set<int> zeroEdgeLevels;
	for(int level = 0; level < levels - 1; level++)
	{
		assert(edgesFromLevel.at(level).size() > 0);
		bool allDeleted = true;
		for(int edgeI = 0; edgeI < (int)edgesFromLevel.at(level).size(); edgeI++)
		{
			Edge* e = edgesFromLevel.at(level).at(edgeI);
			allDeleted = (allDeleted && (deleteEdge.count(e) > 0));
		}
		if(allDeleted)
		{
			zeroEdgeLevels.insert(level);
		}
	}

	if(protectLevels.size() > 0)
	{
		assert(zeroEdgeLevels.size() == 0);
	}

	if(zeroEdgeLevels.size() > 0)
	{
		return simplifyAccordingToCoverage(mG, estimatedEmissions, zeroEdgeLevels);
	}

	MultiGraph* mGsimple = new MultiGraph();
	mGsimple->CODE = mG->CODE;
	mGsimple->kMerSize = mG->kMerSize;
	mGsimple->underlyingLargeGraph = mG->underlyingLargeGraph;

	Node* multiN0 = new Node();
	multiN0->level = 0;
	multiN0->terminal = false;
	mGsimple->registerNode(multiN0, 0);

	map<Node*, Node*> oldNodeToNewNode;
	oldNodeToNewNode[*(mG->NodesPerLevel.at(0).begin())] = multiN0;
	assert(deleteNode.count(*(mG->NodesPerLevel.at(0).begin())) == 0);

	int max_remaining_edges = 0;
	int max_remaining_edges_before = 0;

	for(int level = 0; level < levels - 1; level++)
	{
		set<Node*> nodesAtLevel = mG->NodesPerLevel.at(level);

		if((int)edgesFromLevel.at(level).size() > max_remaining_edges_before)
		{
			max_remaining_edges_before = edgesFromLevel.at(level).size();
		}

		int level_remaining_edges = 0;

		for(set<Node*>::iterator nIt = nodesAtLevel.begin(); nIt != nodesAtLevel.end(); nIt++)
		{
			Node* n = *nIt;

			int remainingEdges = 0;

			for(set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				if(deleteEdge.count(e) == 0)
				{
					remainingEdges++;
				}
			}

			level_remaining_edges += remainingEdges;

			if(remainingEdges > 0)
			{
				for(set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
				{
					Edge* oldE = *eIt;

					if(deleteEdge.count(oldE) == 0)
					{
						assert(oldNodeToNewNode.count(oldE->From) > 0);

						Edge* newE = new Edge();
						newE->From = oldNodeToNewNode[oldE->From];
						oldNodeToNewNode[oldE->From]->Outgoing_Edges.insert(newE);
						newE->locus_id = oldE->locus_id;
						newE->label = oldE->label;
						newE->multiEmission = oldE->multiEmission;
						newE->count = oldE->count;
						newE->multiEmission_full = oldE->multiEmission_full;
						newE->_largeGraphEdge = oldE->_largeGraphEdge;

						if(oldNodeToNewNode.count(oldE->To) == 0)
						{
							assert(deleteNode.count(oldE->To) == 0);

							Node* newToNode = new Node();
							newToNode->terminal = oldE->To->terminal;
							newToNode->level =  oldE->To->level;
							mGsimple->registerNode(newToNode, newToNode->level);
							oldNodeToNewNode[oldE->To] = newToNode;
						}

						newE->To = oldNodeToNewNode[oldE->To];
						oldNodeToNewNode[oldE->To]->Incoming_Edges.insert(newE);
						mGsimple->registerEdge(newE);
					}
				}
			}
			else
			{
				assert(deleteNode.count(n) > 0);
			}


			if(level_remaining_edges > max_remaining_edges)
			{
				max_remaining_edges = level_remaining_edges;
			}

		}
	}

	mGsimple->checkConsistency(false);

	cout << "\nSimplification done. Maximum edges per level is now " << max_remaining_edges << ", whereas it was " << max_remaining_edges_before << " before.\n\n" << flush;
	assert(max_remaining_edges < 5000);

	return mGsimple;
}


MultiGraph* multiBeautifyForAlpha2(LargeGraph* g, string kMerCountsGenomePath, bool quiet, bool pgf_protect)
{
	MultiGraph* mG = new MultiGraph();
	Node* multiN0 = new Node();
	multiN0->level = 0;
	multiN0->terminal = false;
	mG->registerNode(multiN0, 0);
	mG->underlyingLargeGraph = g;
	mG->kMerSize = g->kMerSize;

	// scan levels of conventional graph, determine levels
	if(! quiet)
		cout << "Determine graph segmentation...\n";

	assert(g->NodesPerLevel.at(0).size() == 1);

	vector< vector<int> > levelsInMultiGraph;
	int conventionalLevels = g->NodesPerLevel.size();

	vector<bool> levelCarryStartingGap;
	set<int> levelWantEndSegment;
	levelCarryStartingGap.push_back(false);

	for(int level = 1; level < conventionalLevels; level++)
	{
		bool hasStartingGap = false;

		for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(level).begin(); nodeIt != g->NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* n = *nodeIt;

			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				Edge* e = *edgeIt;
				string kMer = g->CODE.deCode(e->locus_id, e->largeEmission);
				if(kMer == "_")
				{
					bool thisEdgehasStartingGap = true;
					if((n->Incoming_Edges.size() == 1) && (n->Outgoing_Edges.size() == 1))
					{
						Edge* previousEdge = *(n->Incoming_Edges.begin());
						string previousKMer = g->CODE.deCode(previousEdge->locus_id, previousEdge->largeEmission);
						if(previousKMer == "_")
						{
							thisEdgehasStartingGap = false;
						}
					}

					hasStartingGap = (hasStartingGap || thisEdgehasStartingGap);
				}
			}
		}

		levelCarryStartingGap.push_back(hasStartingGap);
		if(hasStartingGap)
		{
			int breakLevel = level - g->kMerSize;
			if(breakLevel > 0)
			{
				levelWantEndSegment.insert(breakLevel);
			}
		}
	}

	for(int level = 0; level < (conventionalLevels-1); level++)
	{
		if((! quiet) && ((level % 100000) == 0))
			cout << "\r" << level << "/" << conventionalLevels;

		int startMultiSegment = level;
		int	stopMultiSegment = level+1;
		bool checkSegmentExtension = true;
		while(checkSegmentExtension)
		{
			bool cannotExtend = false;

			if(levelWantEndSegment.count(stopMultiSegment) > 0)
			{
				cannotExtend = true;
			}
			else
			{
				for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(stopMultiSegment).begin(); nodeIt != g->NodesPerLevel.at(stopMultiSegment).end(); nodeIt++)
				{
					Node* n = *nodeIt;
					if((n->Incoming_Edges.size() !=1) || (n->Outgoing_Edges.size() != 1))
					{
						cannotExtend = true;
					}

					for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
					{
						Edge* e = *edgeIt;
						string kMer = g->CODE.deCode(e->locus_id, e->largeEmission);
						if(kMer == "_")
						{

						}
					}

				}
			}
			if(cannotExtend)
			{
				checkSegmentExtension = false;
			}
			else
			{
				stopMultiSegment++;
			}
		}

		vector<int> segmentSpec;
		segmentSpec.push_back(startMultiSegment);
		segmentSpec.push_back(stopMultiSegment);
		levelsInMultiGraph.push_back(segmentSpec);

		level = stopMultiSegment - 1;
	}

	set<Node*> processedNodes;
	set<Edge*> processedEdges;

	map<Node*, Node*> oldNodeToNewNode;
	oldNodeToNewNode[*(g->NodesPerLevel.at(0).begin())] = multiN0;
	processedNodes.insert(*(g->NodesPerLevel.at(0).begin()));

	int multiGraphLevels = levelsInMultiGraph.size();
	for(int multiGraphLevel = 0; multiGraphLevel < multiGraphLevels; multiGraphLevel++)
	{
		if((! quiet) && ((multiGraphLevel % 10000) == 0))
			cout << "\r" << multiGraphLevel << "/" << multiGraphLevels << flush;

		assert(levelsInMultiGraph.at(multiGraphLevel).size()==2);

		int startLevelNormalGraph = levelsInMultiGraph.at(multiGraphLevel).at(0);
		int stopLevelNormalGraph = levelsInMultiGraph.at(multiGraphLevel).at(1);

		string newLocusID = "MultiLevel"+Utilities::ItoStr(multiGraphLevel);
		for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(startLevelNormalGraph).begin(); nodeIt != g->NodesPerLevel.at(startLevelNormalGraph).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			assert(oldNodeToNewNode.count(n) > 0);
			Node* multiGraphN = oldNodeToNewNode[n];

			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				Edge* e = *edgeIt;

				bool edgeAllPGF = (bool)e->pgf_protect;

				vector<Edge*> processEdges;
				processEdges.push_back(e);

				Edge* workingEdge = e;
				while((int)workingEdge->To->level != stopLevelNormalGraph)
				{
					Node* jumpOver = workingEdge->To;
					processedNodes.insert(jumpOver);
					assert(jumpOver->Outgoing_Edges.size()==1);
					Edge* nextEdge = *(jumpOver->Outgoing_Edges.begin());
					edgeAllPGF = (edgeAllPGF && (bool)nextEdge->pgf_protect);
					processEdges.push_back(nextEdge);
					workingEdge = nextEdge;
				}

				map<int, int> multiEmission;
				vector<int> emissionsInOrder;
				vector<string> labels;
				for(int eI = 0; eI < (int)processEdges.size(); eI++)
				{
					Edge* processEdge = processEdges.at(eI);
					processedEdges.insert(processEdge);
					int emissionCode = mG->CODE.doCode(newLocusID, g->CODE.deCode(processEdge->locus_id, processEdge->largeEmission));
					if(multiEmission.count(emissionCode) == 0)
						multiEmission[emissionCode] = 0;

					multiEmission[emissionCode]++;

					if(processEdge->label != "")
					{
						labels.push_back(processEdge->label);
					}

					emissionsInOrder.push_back(emissionCode);
				}

				Node* edgeFinalNode = workingEdge->To;

				Edge* multiEdge = new Edge();
				multiEdge->From = multiGraphN;
				multiGraphN->Outgoing_Edges.insert(multiEdge);
				multiEdge->locus_id = newLocusID;
				multiEdge->label = Utilities::join(labels, "||");
				multiEdge->multiEmission = multiEmission;
				multiEdge->count = e->count;
				multiEdge->multiEmission_full = emissionsInOrder;
				multiEdge->_largeGraphEdge = e;
				multiEdge->pgf_protect = (bool)edgeAllPGF;

				if(oldNodeToNewNode.count(edgeFinalNode) == 0)
				{
					Node* multiNode = new Node();
					multiNode->terminal = edgeFinalNode->terminal;
					multiNode->level = multiGraphLevel+1;
					mG->registerNode(multiNode, multiGraphLevel+1);
					oldNodeToNewNode[edgeFinalNode] = multiNode;
				}

				multiEdge->To = oldNodeToNewNode[edgeFinalNode];
				oldNodeToNewNode[edgeFinalNode]->Incoming_Edges.insert(multiEdge);
				mG->registerEdge(multiEdge);

				processedNodes.insert(edgeFinalNode);
			}
		}
	}

	if(!quiet)
		cout << "\ndone.\n\n";

	for(set<Edge*>::iterator edgeIt = g->Edges.begin(); edgeIt != g->Edges.end(); edgeIt++)
	{
		assert(processedEdges.count(*edgeIt) > 0);
	}
	for(set<Node*>::iterator nodeIt = g->Nodes.begin(); nodeIt != g->Nodes.end(); nodeIt++)
	{
		assert(processedNodes.count(*nodeIt) > 0);
	}

	mG->checkConsistency(false);

	if(pgf_protect)
	{
		for(unsigned int level = 0; level < (mG->NodesPerLevel.size()-1); level++)
		{
			bool foundPGF = false;
			for(set<Node*>::iterator nodeIt = mG->NodesPerLevel.at(level).begin(); nodeIt != mG->NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;

				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					foundPGF = (foundPGF || (bool)e->pgf_protect);
					if((level > 0) && (level < (mG->NodesPerLevel.size()-2)))
					{
						if((bool)e->pgf_protect)
						{
							bool found_preceding_PGF = false;
							bool found_following_PGF = false;

							for(set<Edge*>::iterator beforeIt = e->From->Incoming_Edges.begin(); beforeIt != e->From->Incoming_Edges.end(); beforeIt++)
							{
								if((bool)(*beforeIt)->pgf_protect)
								{
									found_preceding_PGF = true;
									break;
								}
							}
							assert(found_preceding_PGF);

							for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
							{
								if((bool)(*afterIt)->pgf_protect)
								{
									found_following_PGF = true;
									break;
								}
							}
							assert(found_following_PGF);
						}
					}
				}
			}
			assert(foundPGF);
		}
	}


	return mG;
}


MultiGraph* multiBeautifyForAlpha(LargeGraph* g, string kMerCountsGenomePath, bool quiet)
{
	MultiGraph* mG = new MultiGraph();
	Node* multiN0 = new Node();
	multiN0->level = 0;
	multiN0->terminal = false;
	mG->registerNode(multiN0, 0);
	mG->underlyingLargeGraph = g;
	mG->kMerSize = g->kMerSize;

	// scan levels of conventional graph, determine levels
	if(! quiet)
		cout << "Determine graph segmentation...\n";

	assert(g->NodesPerLevel.at(0).size() == 1);

	vector< vector<int> > firstLevelsInMultiGraph;
	int conventionalLevels = g->NodesPerLevel.size();

	for(int level = 0; level < (conventionalLevels-1); level++)
	{
		if((! quiet) && ((level % 100000) == 0))
			cout << "\r" << level << "/" << conventionalLevels;

		int startMultiSegment = level;
		int	stopMultiSegment = level+1;
		bool checkSegmentExtension = true;
		while(checkSegmentExtension)
		{
			bool cannotExtend = false;
			for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(stopMultiSegment).begin(); nodeIt != g->NodesPerLevel.at(stopMultiSegment).end(); nodeIt++)
			{
				Node* n = *nodeIt;
				if((n->Incoming_Edges.size() !=1) || (n->Outgoing_Edges.size() != 1))
				{
					cannotExtend = true;
				}
			}
			if(cannotExtend)
			{
				checkSegmentExtension = false;
			}
			else
			{
				stopMultiSegment++;
			}
		}

		vector<int> segmentSpec;
		segmentSpec.push_back(startMultiSegment);
		segmentSpec.push_back(stopMultiSegment);
		firstLevelsInMultiGraph.push_back(segmentSpec);

		level = stopMultiSegment - 1;
	}



	vector< vector<int> > levelsInMultiGraph;
	vector<int> minRealEdges;
	vector<int> nominalArcLengths;

	int firstMultiGraphLevels = firstLevelsInMultiGraph.size();
	for(int multiGraphLevel = 0; multiGraphLevel < firstMultiGraphLevels; multiGraphLevel++)
	{
		assert(firstLevelsInMultiGraph.at(multiGraphLevel).size()==2);

		int startLevelNormalGraph = firstLevelsInMultiGraph.at(multiGraphLevel).at(0);
		int stopLevelNormalGraph = firstLevelsInMultiGraph.at(multiGraphLevel).at(1);

		int nominalArcLength = stopLevelNormalGraph - startLevelNormalGraph + 1;
		nominalArcLengths.push_back(nominalArcLength);

		int minArcRealEdges = -1;

		for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(startLevelNormalGraph).begin(); nodeIt != g->NodesPerLevel.at(startLevelNormalGraph).end(); nodeIt++)
		{
			Node* n = *nodeIt;

			int arcRealEdges = 0;
			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				Edge* e = *edgeIt;

				vector<Edge*> processEdges;
				Edge* workingEdge = e;
				while((int)workingEdge->To->level != stopLevelNormalGraph)
				{
					Node* jumpOver = workingEdge->To;
					assert(jumpOver->Outgoing_Edges.size()==1);
					Edge* nextEdge = *(jumpOver->Outgoing_Edges.begin());
					processEdges.push_back(nextEdge);
					workingEdge = nextEdge;
				}

				for(int eI = 0; eI < (int)processEdges.size(); eI++)
				{
					Edge* processEdge = processEdges.at(eI);
					string kMer = g->CODE.deCode(processEdge->locus_id, processEdge->largeEmission);
					if((kMer != "_") && (kMer != "*"))
					{
						arcRealEdges++;
					}
				}
			}

			if((minArcRealEdges == -1) || (arcRealEdges < minArcRealEdges))
			{
				minArcRealEdges = arcRealEdges;
			}
		}

		minRealEdges.push_back(minArcRealEdges);
	}

	assert((int)minRealEdges.size() == firstMultiGraphLevels);
	assert((int)nominalArcLengths.size() == firstMultiGraphLevels);

	/*
	for(int multiGraphLevel = 1; multiGraphLevel < firstMultiGraphLevels; multiGraphLevel++)
	{
		int startLevelNormalGraph = (int)firstLevelsInMultiGraph.at(multiGraphLevel).at(0);
		int stopLevelNormalGraph = (int)firstLevelsInMultiGraph.at(multiGraphLevel).at(1);

		if((minRealEdges.at(multiGraphLevel) < g->kMerSize) && (nominalArcLengths.at(multiGraphLevel-1) > g->kMerSize))
		{

		}
	}
	*/

	set<Node*> processedNodes;
	set<Edge*> processedEdges;

	map<Node*, Node*> oldNodeToNewNode;
	oldNodeToNewNode[*(g->NodesPerLevel.at(0).begin())] = multiN0;
	processedNodes.insert(*(g->NodesPerLevel.at(0).begin()));

	int multiGraphLevels = levelsInMultiGraph.size();
	for(int multiGraphLevel = 0; multiGraphLevel < multiGraphLevels; multiGraphLevel++)
	{
		if((! quiet) && ((multiGraphLevel % 10000) == 0))
			cout << "\r" << multiGraphLevel << "/" << multiGraphLevels << flush;

		assert(levelsInMultiGraph.at(multiGraphLevel).size()==2);

		int startLevelNormalGraph = levelsInMultiGraph.at(multiGraphLevel).at(0);
		int stopLevelNormalGraph = levelsInMultiGraph.at(multiGraphLevel).at(1);

		string newLocusID = "MultiLevel"+Utilities::ItoStr(multiGraphLevel);
		for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(startLevelNormalGraph).begin(); nodeIt != g->NodesPerLevel.at(startLevelNormalGraph).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			assert(oldNodeToNewNode.count(n) > 0);
			Node* multiGraphN = oldNodeToNewNode[n];

			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				Edge* e = *edgeIt;

				vector<Edge*> processEdges;
				processEdges.push_back(e);

				Edge* workingEdge = e;
				while((int)workingEdge->To->level != stopLevelNormalGraph)
				{
					Node* jumpOver = workingEdge->To;
					processedNodes.insert(jumpOver);
					assert(jumpOver->Outgoing_Edges.size()==1);
					Edge* nextEdge = *(jumpOver->Outgoing_Edges.begin());
					processEdges.push_back(nextEdge);
					workingEdge = nextEdge;
				}

				map<int, int> multiEmission;
				vector<int> emissionsInOrder;
				vector<string> labels;
				for(int eI = 0; eI < (int)processEdges.size(); eI++)
				{
					Edge* processEdge = processEdges.at(eI);
					processedEdges.insert(processEdge);
					int emissionCode = mG->CODE.doCode(newLocusID, g->CODE.deCode(processEdge->locus_id, processEdge->largeEmission));
					if(multiEmission.count(emissionCode) == 0)
						multiEmission[emissionCode] = 0;

					multiEmission[emissionCode]++;

					if(processEdge->label != "")
					{
						labels.push_back(processEdge->label);
					}

					emissionsInOrder.push_back(emissionCode);
				}

				Node* edgeFinalNode = workingEdge->To;

				Edge* multiEdge = new Edge();
				multiEdge->From = multiGraphN;
				multiGraphN->Outgoing_Edges.insert(multiEdge);
				multiEdge->locus_id = newLocusID;
				multiEdge->label = Utilities::join(labels, "||");
				multiEdge->multiEmission = multiEmission;
				multiEdge->count = e->count;
				multiEdge->multiEmission_full = emissionsInOrder;
				multiEdge->_largeGraphEdge = e;

				if(oldNodeToNewNode.count(edgeFinalNode) == 0)
				{
					Node* multiNode = new Node();
					multiNode->terminal = edgeFinalNode->terminal;
					multiNode->level = multiGraphLevel+1;
					mG->registerNode(multiNode, multiGraphLevel+1);
					oldNodeToNewNode[edgeFinalNode] = multiNode;
				}

				multiEdge->To = oldNodeToNewNode[edgeFinalNode];
				oldNodeToNewNode[edgeFinalNode]->Incoming_Edges.insert(multiEdge);
				mG->registerEdge(multiEdge);

				processedNodes.insert(edgeFinalNode);
			}
		}
	}

	if(!quiet)
		cout << "\ndone.\n\n";

	for(set<Edge*>::iterator edgeIt = g->Edges.begin(); edgeIt != g->Edges.end(); edgeIt++)
	{
		assert(processedEdges.count(*edgeIt) > 0);
	}
	for(set<Node*>::iterator nodeIt = g->Nodes.begin(); nodeIt != g->Nodes.end(); nodeIt++)
	{
		assert(processedNodes.count(*nodeIt) > 0);
	}

	mG->checkConsistency(false);

	return mG;
}

MultiGraph* multiBeautify(LargeGraph* g, string kMerCountsGenomePath, bool quiet)
{
	MultiGraph* mG = new MultiGraph();
	Node* multiN0 = new Node();
	multiN0->level = 0;
	multiN0->terminal = false;
	mG->registerNode(multiN0, 0);

	mG->kMerSize = g->kMerSize;

	// read in kMer count in reference genome minus area covered by graph

	map<string, int> kMerCountsInGenomeExceptGraph;
	ifstream kMerCountsGenomeFile;
	kMerCountsGenomeFile.open (kMerCountsGenomePath.c_str(), ios::in);
	int lineCounter = 0;
	if(kMerCountsGenomeFile.is_open())
	{
		string line;
		while(kMerCountsGenomeFile.good())
		{
			lineCounter++;
			getline (kMerCountsGenomeFile, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			vector<string> fields = Utilities::split(line, ' ');
			if(fields.size() != 2)
				errEx("Strange format in kMerCountsInGenomeExceptGraph file. Expect two fields (kMer and count). Problem occured in line: "+Utilities::ItoStr(lineCounter));

			assert((int)fields.at(0).length() == g->kMerSize);

			kMerCountsInGenomeExceptGraph[fields.at(0)] = Utilities::StrtoI(fields.at(1));
		}
		kMerCountsGenomeFile.close();
	}
	else
	{
		errEx("Cannot open kMer counts file: "+kMerCountsGenomePath);
	}

	// scan levels of conventional graph, determine levels
	if(! quiet)
		cout << "Determine graph segmentation...\n";
	vector< vector<int> > levelsInMultiGraph;
	int conventionalLevels = g->NodesPerLevel.size();
	bool continueStretchNextLevel = false;
	for(int level = 0; level < conventionalLevels; level++)
	{
		if(! quiet)
			cout << "\r" << level << "/" << conventionalLevels;

		bool nodeCondition = true;
		set<Node*> seenTargetNodes;
		for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(level).begin(); nodeIt != g->NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			if(n->Outgoing_Edges.size() != 1)
			{
				nodeCondition = false;
			}
			else
			{
				for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					Node* targetNode = e->To;

					if(seenTargetNodes.count(targetNode))
					{
						nodeCondition = false;
					}
					else
					{
						seenTargetNodes.insert(targetNode);
					}
				}
			}
		}

		if(nodeCondition && continueStretchNextLevel)
		{
			levelsInMultiGraph.back().push_back(level);
		}
		else
		{
			vector<int> newLast;
			newLast.push_back(level);
			levelsInMultiGraph.push_back(newLast);
		}

		continueStretchNextLevel = nodeCondition;
	}

	if(! quiet)
		cout << "\nDone. Found " << levelsInMultiGraph.size() << " segments.\n\n";

	// find out in which segments which kMers occur
	if(! quiet)
		cout << "Examine kMer segment-uniqueness\n\n";

	int multiGraphLevels = levelsInMultiGraph.size();
	set<int> ConventionalLevelNonUniqueKMers;
	map<string, set<int> > kMerToSegment;
	map<string, set<int> > kMerToOriginalLevel;
	set<string> problematicKMers;
	for(int multiGraphLevel = 0; multiGraphLevel < multiGraphLevels; multiGraphLevel++)
	{
		int attachedNormalLevelsIndex = levelsInMultiGraph.at(multiGraphLevel).size();
		for(int levelI = 0; levelI < attachedNormalLevelsIndex; levelI++)
		{
			int oldGraphLevel = levelsInMultiGraph.at(multiGraphLevel).at(levelI);
			for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(oldGraphLevel).begin(); nodeIt != g->NodesPerLevel.at(oldGraphLevel).end(); nodeIt++)
			{
				Node* oldNode = *nodeIt;
				for(set<Edge*>::iterator edgeIt = oldNode->Outgoing_Edges.begin(); edgeIt != oldNode->Outgoing_Edges.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					string kMer = g->CODE.deCode(e->locus_id, e->largeEmission);
					kMerToSegment[kMer].insert(multiGraphLevel);
					kMerToOriginalLevel[kMer].insert(oldGraphLevel);
				}
			}
		}
	}



	for(map<string, set<int> >::iterator kMerIt = kMerToSegment.begin(); kMerIt != kMerToSegment.end(); kMerIt++)
	{
		string kMer = kMerIt->first;
		if(!(((kMer == "_")) || (kMer == "*") || (kMerCountsInGenomeExceptGraph.count(kMerIt->first) > 0)))
		{
			cerr << "Data for kMer " << kMerIt->first << " missing from kMer input file " << kMerCountsGenomePath << "\n\n";
		}
		assert(((kMer == "_")) || (kMer == "*") || (kMerCountsInGenomeExceptGraph.count(kMerIt->first) > 0));
		if(kMerIt->second.size() > 1)
		{
			assert(kMerToOriginalLevel[kMer].size() > 1);
		}

		if((kMer != "_") && ((kMerToOriginalLevel[kMer].size() > 1) || (kMerCountsInGenomeExceptGraph[kMerIt->first] > 0)))
		{
			problematicKMers.insert(kMer);
		}

		if(kMer == "*")
		{
			problematicKMers.insert(kMer);
		}
	}

	for(map<string, set<int> >::iterator kMerIt = kMerToOriginalLevel.begin(); kMerIt != kMerToOriginalLevel.end(); kMerIt++)
	{
		string kMer = kMerIt->first;
		if(problematicKMers.count(kMer) > 0)
		{
			for(set<int>::iterator levelIt = kMerIt->second.begin(); levelIt != kMerIt->second.end(); levelIt++)
			{
				int associatedLevel = *levelIt;
				ConventionalLevelNonUniqueKMers.insert(associatedLevel);
			}
		}
	}

	// create new multiGraph
	if(! quiet)
		cout << "Create multigraph...\n" << flush;

	set<Node*> processedNodes;
	set<Edge*> processedEdges;
	map<Node*, Node*> oldNodeToNewNode;
	assert(g->NodesPerLevel.at(0).size() == 1);
	oldNodeToNewNode[*(g->NodesPerLevel.at(0).begin())] = multiN0;
	for(int multiGraphLevel = 0; multiGraphLevel < multiGraphLevels; multiGraphLevel++)
	{
		if(! quiet)
			cout << "\r" << multiGraphLevel << "/" << multiGraphLevels << flush;

		int attachedNormalLevelsIndex = levelsInMultiGraph.at(multiGraphLevel).size();
		if(attachedNormalLevelsIndex == 1)
		{
			int level = levelsInMultiGraph.at(multiGraphLevel).at(0);

			for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(level).begin(); nodeIt != g->NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* oldNode = *nodeIt;
				assert(oldNodeToNewNode.count(oldNode) > 0);
				Node* fromNodeMultiGraph = oldNodeToNewNode[oldNode];

				for(set<Edge*>::iterator edgeIt = oldNode->Outgoing_Edges.begin(); edgeIt != oldNode->Outgoing_Edges.end(); edgeIt++)
				{
					Edge* oldEdge = *edgeIt;

					Edge* multiEdge = new Edge();
					multiEdge->From = fromNodeMultiGraph;
					fromNodeMultiGraph->Outgoing_Edges.insert(multiEdge);
					multiEdge->count = oldEdge->count;
					multiEdge->locus_id = oldEdge->locus_id;
					multiEdge->label = oldEdge->label;

					int emissionCode = mG->CODE.doCode(oldEdge->locus_id, g->CODE.deCode(oldEdge->locus_id, oldEdge->largeEmission));

					if(ConventionalLevelNonUniqueKMers.count(level) == 1)
					{
						multiEdge->multiEmission[emissionCode] = 1;
					}
					else
					{
						multiEdge->multiEmission_full.push_back(emissionCode);
					}

					Node* oldEdgeOldTargetNode = oldEdge->To;
					if(oldNodeToNewNode.count(oldEdgeOldTargetNode) == 0)
					{
						Node* multiNode = new Node();
						multiNode->terminal = oldEdge->To->terminal;
						multiNode->level = multiGraphLevel+1;
						mG->registerNode(multiNode, multiGraphLevel+1);
						oldNodeToNewNode[oldEdgeOldTargetNode] = multiNode;
					}

					multiEdge->To = oldNodeToNewNode[oldEdgeOldTargetNode];
					oldNodeToNewNode[oldEdgeOldTargetNode]->Incoming_Edges.insert(multiEdge);

					mG->registerEdge(multiEdge);

					processedEdges.insert(oldEdge);

				}

				processedNodes.insert(oldNode);
			}

		}
		else
		{
			int startLevel = levelsInMultiGraph.at(multiGraphLevel).at(0);
			for(set<Node*>::iterator nodeIt = g->NodesPerLevel.at(startLevel).begin(); nodeIt != g->NodesPerLevel.at(startLevel).end(); nodeIt++)
			{
				Node* oldNode = *nodeIt;

				assert(oldNodeToNewNode.count(oldNode) > 0);
				Node* fromNodeMultiGraph = oldNodeToNewNode[oldNode];

				assert(oldNode->Outgoing_Edges.size() == 1);
				Edge* firstEdge = *(oldNode->Outgoing_Edges.begin());

				Edge* multiEdge = new Edge();
				string multiEdgeLocusID = "MultiLevel"+Utilities::ItoStr(multiGraphLevel);
				multiEdge->From = fromNodeMultiGraph;
				fromNodeMultiGraph->Outgoing_Edges.insert(multiEdge);
				multiEdge->count = firstEdge->count;
				multiEdge->locus_id = multiEdgeLocusID;

				Node* workingNode = oldNode;
				string newLabel = "";
				int addedLabels = 0;

				for(int levelI = 0; levelI < attachedNormalLevelsIndex; levelI++)
				{
					int oldGraphLevel = levelsInMultiGraph.at(multiGraphLevel).at(levelI);
					if(levelI != 0)
					{
						if(!(oldGraphLevel == (levelsInMultiGraph.at(multiGraphLevel).at(levelI-1)+1)))
						{
							cout << "oldGraphLevel: " << oldGraphLevel << " multiGraphLevel: " << multiGraphLevel << " levelsInMultiGraph.at(multiGraphLevel).at(levelI-1): " << levelsInMultiGraph.at(multiGraphLevel).at(levelI-1) << "\n";
						}
						assert(oldGraphLevel == (levelsInMultiGraph.at(multiGraphLevel).at(levelI-1)+1));
					}
					assert((int)workingNode->level == oldGraphLevel);
					assert(workingNode->Outgoing_Edges.size() == 1);
					Edge* attachedEdge = *(workingNode->Outgoing_Edges.begin());

					int emissionCode = mG->CODE.doCode(multiEdgeLocusID, g->CODE.deCode(attachedEdge->locus_id, attachedEdge->largeEmission));
					newLabel = newLabel + attachedEdge->label;
					if(attachedEdge->label != "")
					{
						addedLabels++;
					}

					if(ConventionalLevelNonUniqueKMers.count(oldGraphLevel) == 1)
					{
						if(multiEdge->multiEmission.count(emissionCode) == 0)
							multiEdge->multiEmission[emissionCode] = 0;
						multiEdge->multiEmission[emissionCode]++;

						//multiEdge->multiEmission[emissionCode] = 1;
					}
					else
					{
						multiEdge->multiEmission_full.push_back(emissionCode);
					}

					//multiEdge->multiEmission_full.push_back(emissionCode);

					workingNode = attachedEdge->To;
					processedEdges.insert(attachedEdge);
					processedNodes.insert(workingNode);

				}

				if(addedLabels > 1)
				{
					cout << "addedLabels > 1: " << addedLabels << "\n";
					cout << newLabel << "\n";
				}

				assert(addedLabels <= 1);
				multiEdge->label = newLabel;

				Node* oldEdgeOldTargetNode = workingNode;
				if(oldNodeToNewNode.count(oldEdgeOldTargetNode) == 0)
				{
					Node* multiNode = new Node();
					multiNode->terminal = oldEdgeOldTargetNode->terminal;
					multiNode->level = multiGraphLevel+1;
					mG->registerNode(multiNode, multiGraphLevel+1);
					oldNodeToNewNode[oldEdgeOldTargetNode] = multiNode;
				}

				multiEdge->To = oldNodeToNewNode[oldEdgeOldTargetNode];
				oldNodeToNewNode[oldEdgeOldTargetNode]->Incoming_Edges.insert(multiEdge);

				mG->registerEdge(multiEdge);
				processedNodes.insert(oldNode);
			}
		}
	}

	if(!quiet)
		cout << "\ndone.\n\n";


	for(set<Edge*>::iterator edgeIt = g->Edges.begin(); edgeIt != g->Edges.end(); edgeIt++)
	{
		assert(processedEdges.count(*edgeIt) > 0);
	}
	for(set<Node*>::iterator nodeIt = g->Nodes.begin(); nodeIt != g->Nodes.end(); nodeIt++)
	{
		assert(processedNodes.count(*nodeIt) > 0);
	}

	mG->checkConsistency(false);

	return mG;
}

LargeGraph* kMerify(Graph* g, bool quiet, int kMerSize, bool wantPGFprotection)
{
	LargeGraph* kMerGraph = new LargeGraph();

	if(! quiet)
		cout << "Initialize...\n";

	Node* KMerN0 = new Node();
	KMerN0->level = 0;
	KMerN0->terminal = false;
	kMerGraph->registerNode(KMerN0, 0);

	Node* classicalN0 = *(g->NodesPerLevel.at(0).begin());

	for(set<Edge*>::iterator edgeIt = classicalN0->Outgoing_Edges.begin(); edgeIt != classicalN0->Outgoing_Edges.end(); edgeIt++)
	{
		//Edge* e = *edgeIt;
		//assert((2 == 0) || g->CODE.deCode(e->locus_id, e->emission) != "_");
	}

	int levels = g->NodesPerLevel.size();
	vector<string> loci = g->getAssignedLoci();

	map<Node*, vector< kMerAtNode > > originalGraphAssignedKMers; // TODO Koennte theoretisch geloescht werden fuer vorhergehende Levels - 1
	if(! quiet)
		cout << "\tdone.\n";

	std::ostringstream nullPointerString;
	nullPointerString << setw(15) << (void const *) kMerGraph;
	pointerStrintLength = nullPointerString.str().length();

	/*
	for(int level = 4949500; level <= (levels - 1 - kMerSize); level++)
	{
		cout << "Node scan level " << level << "\n" << flush;
		for(set<Node*>::iterator nIt = g->NodesPerLevel.at(level).begin(); nIt != g->NodesPerLevel.at(level).end(); nIt++)
		{
			Node* n = *nIt;
			vector<kMerInfo> startFromN0_GAP = forwardScan(n, 1, 0);
		}
		cout << "\tdone\n" << flush;
	}
	*/

	set<Node*> nodesLastLevel = g->NodesPerLevel[levels-1];
	string lastLocusID = loci.at(loci.size()-1);


	cout << "Last " << (levels-1) << " level emissions:\n";

	for(set<Node*>::iterator nodeIt = nodesLastLevel.begin(); nodeIt != nodesLastLevel.end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			assert(lastLocusID == edge->locus_id);
			string emission = g->CODE.deCode(lastLocusID, edge->emission);
			cout << emission << "\n";
		}
	}
	
	cout << "\n\n" << flush;
	
	for(set<Node*>::iterator nodeIt = nodesLastLevel.begin(); nodeIt != nodesLastLevel.end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			assert(lastLocusID == edge->locus_id);
			string emission = g->CODE.deCode(lastLocusID, edge->emission);
			assert(emission != "_");
		}
	}

	for(int level = 0; level <= (levels - 1 - kMerSize); level++)
	{

		int min_span_level = -1;
		int max_span_level = -1;

		if((! quiet) && ((level % 100) == 0))
		{
			cout << "Level " << level << "/" << (levels - 1 - kMerSize) << " of original graph\n";
			Node* firstN = *(g->NodesPerLevel.at(level).begin());
			if(firstN->Outgoing_Edges.size() > 0)
			{
				Edge* e = *(firstN->Outgoing_Edges.begin());
				cout << "\t\t" << e->locus_id << "\t\t" << e->label << "\n";
			}
		}
		string locusID = loci.at(level);

		map<Node*, map<string, vector<Edge*> > > edgeTargetCache;
		map<Edge*, kMerInfo> newEdgeNodeInfos;

		if(level == 0)
		{
			// Level 0: Forward-Scan nach kMers der Laenge x

			vector<kMerInfo> startFromN0 = forwardScan(classicalN0, kMerSize, -1);

			if(! quiet)
				cout << "\tFound " << startFromN0.size() << " initial kMers with no gap \n";

			for(vector<kMerInfo>::iterator kMerIt = startFromN0.begin(); kMerIt != startFromN0.end(); kMerIt++)
			{
				// Fuer alle gefundenen kMers:
				//	- erzeuge eine Edge
				//	- lege einen Eintrag an im Cache, der ueber die zu erzeugenden Nodes bestimmt
				//		dabei verwenden wir die Target-Node und die Memory-Adressen der letzten x-1 Edges als Index

				string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
				assert((int)kMer_string.size() == kMerSize);

				string newLocusID = "L"+Utilities::ItoStr(level);

				// cout << kMer_string << "\n" << flush;

				Edge* kMerEdge = new Edge();



				kMerEdge->From = KMerN0;
				KMerN0->Outgoing_Edges.insert(kMerEdge);
				kMerEdge->count = kMerIt->p;
				kMerEdge->pgf_protect = (bool)kMerIt->allPGF;

				assert(g->CODE.deCode(kMerIt->traverseEdges.at(0)->locus_id, kMerIt->traverseEdges.at(0)->emission) != "_");

				kMerEdge->largeEmission = kMerGraph->CODE.doCode(newLocusID, kMer_string);
				kMerIt->gapEdge = false;

				kMerEdge->locus_id = newLocusID;
				kMerGraph->registerEdge(kMerEdge);

				if(kMerIt->traverseEdges.at(0)->label != "")
				{
					kMerEdge->label = kMerIt->traverseEdges.at(0)->label;
				}

				for(int lI = 0; lI < kMerIt->traverseEdges.size(); lI++)
				{
					Edge* traversedEdge = kMerIt->traverseEdges.at(lI);
					if(g->CODE.deCode(traversedEdge->locus_id, traversedEdge->emission) != "_")
					{
						kMerEdge->levelsNucleotideGraph.push_back(traversedEdge->From->level);
					}
				}
				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);

				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;
			}

			// Level 0: Forward-Scan nach kMers der Laenge x-1 mit Gaps vorne

			vector<kMerInfo> startFromN0_GAP = forwardScan(classicalN0, kMerSize-1, 1);

			if(! quiet)
				cout << "\tFound " << startFromN0_GAP.size() << " initial kMers with gap \n";

			for(vector<kMerInfo>::iterator kMerIt = startFromN0_GAP.begin(); kMerIt != startFromN0_GAP.end(); kMerIt++)
			{
				// Fuer alle gefundenen kMers:
				//	- erzeuge eine Edge
				//	- lege einen Eintrag an im Cache, der ueber die zu erzeugenden Nodes bestimmt
				//		dabei verwenden wir die Target-Node und die Memory-Adressen der letzten x-1 Edges als Index

				string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
				assert((int)kMer_string.size() == kMerSize-1);

				string newLocusID = "L"+Utilities::ItoStr(level);

				// cout << kMer_string << "\n" << flush;

				Edge* kMerEdge = new Edge();

				kMerEdge->From = KMerN0;
				KMerN0->Outgoing_Edges.insert(kMerEdge);
				kMerEdge->count = kMerIt->p;
				kMerEdge->pgf_protect = (bool)kMerIt->allPGF;

				assert(g->CODE.deCode(kMerIt->traverseEdges.at(0)->locus_id, kMerIt->traverseEdges.at(0)->emission) == "_");

				kMerEdge->largeEmission = kMerGraph->CODE.doCode(newLocusID, "_");
				kMerIt->gapEdge = true;


				kMerEdge->locus_id = newLocusID;
				kMerGraph->registerEdge(kMerEdge);

				if(kMerIt->traverseEdges.at(0)->label != "")
				{
					kMerEdge->label = kMerIt->traverseEdges.at(0)->label;
				}

				for(int lI = 0; lI < kMerIt->traverseEdges.size(); lI++)
				{
					Edge* traversedEdge = kMerIt->traverseEdges.at(lI);
					if(g->CODE.deCode(traversedEdge->locus_id, traversedEdge->emission) != "_")
					{
						kMerEdge->levelsNucleotideGraph.push_back(traversedEdge->From->level);
					}
				}

				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);

				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;
			}

		}
		else
		{
			// alle Nodes dieses Levels
			vector<Node*> originalNodes = vector<Node*>(g->NodesPerLevel.at(level).begin(), g->NodesPerLevel.at(level).end());
			int node_c = 0;
			int node_c_total = originalNodes.size();

			for(vector<Node*>::iterator originalNodeIt = originalNodes.begin(); originalNodeIt != originalNodes.end(); originalNodeIt++)
			{
				node_c++;

				// potentielle Kanten (attachte KMers) durchgehen
				Node* originalNode = *originalNodeIt;
				assert(originalNode->Incoming_Edges.size() > 0);

				if(!(originalGraphAssignedKMers.count(originalNode) > 0))
				{
					cout << "\t Level " << level << ", search for originalNode " << originalNode << " information with result " << (originalGraphAssignedKMers.count(originalNode) > 0) << "\n";
					Edge* edgesToPrevious = *(originalNode->Incoming_Edges.begin());
					Node* previousNode = edgesToPrevious->From;

					vector<kMerInfo> startFromLastNode = forwardScan(previousNode, kMerSize, 0);

					cout << "\t Found " << startFromLastNode.size() << " paths which should traverse our dubious node, taking edge " << edgesToPrevious << "\n";

					for(vector<kMerInfo>::iterator kMerIt = startFromLastNode.begin(); kMerIt != startFromLastNode.end(); kMerIt++)
					{
						string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
						cout << "\t\t kMer_string " << kMer_string << "\n";
						cout << "\t\t traversed Edges" << kMerIt->traverseEdges_string << "\n";
					}

				}

				assert(originalGraphAssignedKMers.count(originalNode) > 0);

				int attachedKMer_c = 0;
				int attachedKMer_c_total = originalGraphAssignedKMers[originalNode].size();
				for(vector< kMerAtNode >::iterator kMerBasisEdgeIt = originalGraphAssignedKMers[originalNode].begin();  kMerBasisEdgeIt != originalGraphAssignedKMers[originalNode].end(); kMerBasisEdgeIt++)
				{
					attachedKMer_c++;
					if((! quiet) && ((level % 100) == 0))
					{
						cout << "\r                                                            ";
						cout << "\r\t Node " << node_c << "/" << node_c_total << ", attached kMer " << attachedKMer_c << "/" << attachedKMer_c_total << flush;
					}
					// fuer jede attachte Kante neue Kanten erstellen, die sich von der letzten Node
					// im Original-Graph ergeben

					kMerAtNode BasisForNewEdges = *kMerBasisEdgeIt;

					assert(BasisForNewEdges.traverseEdges.at(0)->From == originalNode);

					// Wenn das erste Symbol des attachten kMers kein Gap ist, fuegen wir frohlich neue Kanten hinzu
					if(g->CODE.deCode(BasisForNewEdges.traverseEdges.at(0)->locus_id, BasisForNewEdges.traverseEdges.at(0)->emission) != "_")
					{
						Node* lastNodeInOriginalGraph = BasisForNewEdges.traverseEdges.back()->To;

						// some debug information

						if(min_span_level == -1)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}
						if(max_span_level == -1)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level > max_span_level)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level < min_span_level)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}

						vector<kMerInfo> startFromLastNode = forwardScan(lastNodeInOriginalGraph, 1, 0);


						if(startFromLastNode.size() == 0)
						{
							if(! quiet)
								cout << " -- reached end of graph" << "\n";

							//assert((int)BasisForNewEdges.traverseEdges.back()->To->level == (int)levels);

							// Wir sind am Ende des normalen Graphen angelangt, haben aber noch zu wenige Positionen im kmer-Graph
							// Das kommt durch Gaps in den letzten X Symbolen
							// Wir fuegen Gap-KMers hinzu

							// Wenn das erste Symbol dieser potentiellen Kante ein Gap ist, muessen wir nur eine Gap-Edge ans Ende hinzufuegen
							// neues kMerInfo-Objekt, GAP!

							kMerInfo newNodeKMerInfo;
							newNodeKMerInfo.gapEdge = true;
							newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
							newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
							newNodeKMerInfo.p = 1;
							newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
							newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;
							newNodeKMerInfo.allPGF = BasisForNewEdges.allPGF;

							string newLocusID = "L"+Utilities::ItoStr(level);
							string kMer_string = "_";

							Edge* kMerEdge = new Edge();
							kMerEdge->From = BasisForNewEdges.lastNewNode;
							BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

							kMerEdge->count = 1;
							kMerEdge->largeEmission = kMerGraph->CODE.doCode(newLocusID, kMer_string);
							kMerEdge->locus_id = newLocusID;
							kMerEdge->pgf_protect = (bool)newNodeKMerInfo.allPGF;

							kMerGraph->registerEdge(kMerEdge);

							if(newNodeKMerInfo.traverseEdges.at(0)->label != "")
							{
								kMerEdge->label = newNodeKMerInfo.traverseEdges.at(0)->label;
							}

							//vector<Edge*> traverseEdges_m1 = newNodeKMerInfo.traverseEdges;
							//traverseEdges_m1.erase (traverseEdges_m1.begin(),traverseEdges_m1.begin()+1);

							assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
							string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

							assert(g->Nodes.count(newNodeKMerInfo.traverseEdges.back()->To) > 0);

							edgeTargetCache[newNodeKMerInfo.traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
							newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;
						}
						else
						{
							if((! quiet) && ((level % 100) == 0))
								cout << " -- " << startFromLastNode.size() << " expansions" << "\n";

							for(vector<kMerInfo>::iterator attachKMerIt = startFromLastNode.begin(); attachKMerIt != startFromLastNode.end(); attachKMerIt++)
							{

									assert(attachKMerIt->kMer_coded.size() == 1);
									assert(attachKMerIt->kMer_deCoded.size() == 1);

									// neues kMerInfo-Objekt
									kMerInfo newNodeKMerInfo;
									newNodeKMerInfo.gapEdge = false;
									newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
									newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
									newNodeKMerInfo.p = attachKMerIt->p;
									newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
									newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;
									assert(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1));


									newNodeKMerInfo.kMer_coded.push_back(attachKMerIt->kMer_coded.at(0));
									newNodeKMerInfo.kMer_deCoded.push_back(attachKMerIt->kMer_deCoded.at(0));
									newNodeKMerInfo.traverseEdges.insert(newNodeKMerInfo.traverseEdges.end(), attachKMerIt->traverseEdges.begin(), attachKMerIt->traverseEdges.end());
									newNodeKMerInfo.traverseEdges_string = newNodeKMerInfo.traverseEdges_string.append(attachKMerIt->traverseEdges_string);
									newNodeKMerInfo.allPGF = (BasisForNewEdges.allPGF && attachKMerIt->allPGF);

									string newLocusID = "L"+Utilities::ItoStr(level);
									string kMer_string = Utilities::join(newNodeKMerInfo.kMer_deCoded, "");
									assert((int)kMer_string.size() == kMerSize);

									Edge* kMerEdge = new Edge();
									kMerEdge->From = BasisForNewEdges.lastNewNode;
									BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

									kMerEdge->count = attachKMerIt->p;
									kMerEdge->largeEmission = kMerGraph->CODE.doCode(newLocusID, kMer_string);
									kMerEdge->locus_id = newLocusID;
									kMerEdge->pgf_protect = newNodeKMerInfo.allPGF;

									kMerGraph->registerEdge(kMerEdge);

									if(newNodeKMerInfo.traverseEdges.at(0)->label != "")
									{
										kMerEdge->label = newNodeKMerInfo.traverseEdges.at(0)->label;
									}


									for(int lI = 0; lI < attachKMerIt->traverseEdges.size(); lI++)
									{
										Edge* traversedEdge = attachKMerIt->traverseEdges.at(lI);
										if(g->CODE.deCode(traversedEdge->locus_id, traversedEdge->emission) != "_")
										{
											kMerEdge->levelsNucleotideGraph.push_back(traversedEdge->From->level);
										}
									}


									assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
									string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);


									edgeTargetCache[attachKMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
									newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;
							}
						}
					}
					else
					{

						Node* lastNodeInOriginalGraph = BasisForNewEdges.traverseEdges.back()->To;

						if(min_span_level == -1)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}
						if(max_span_level == -1)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level > max_span_level)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level < min_span_level)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}

						// Wenn das erste Symbol dieser potentiellen Kante ein Gap ist, muessen wir nur eine Gap-Edge ans Ende hinzufuegen
						// neues kMerInfo-Objekt, GAP!

						if((! quiet) && ((level % 100) == 0))
							cout << " -- gap length " << BasisForNewEdges.traverseEdges.size() << "\n";

						kMerInfo newNodeKMerInfo;
						newNodeKMerInfo.gapEdge = true;
						newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
						newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
						newNodeKMerInfo.p = 1;
						newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
						newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;

						string newLocusID = "L"+Utilities::ItoStr(level);
						string kMer_string = "_";

						Edge* kMerEdge = new Edge();
						kMerEdge->From = BasisForNewEdges.lastNewNode;
						BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

						kMerEdge->count = 1;
						kMerEdge->largeEmission = kMerGraph->CODE.doCode(newLocusID, kMer_string);
						kMerEdge->locus_id = newLocusID;
						kMerEdge->pgf_protect = BasisForNewEdges.allPGF;

						kMerGraph->registerEdge(kMerEdge);

						if(newNodeKMerInfo.traverseEdges.at(0)->label != "")
						{
							kMerEdge->label = newNodeKMerInfo.traverseEdges.at(0)->label;
						}

						//vector<Edge*> traverseEdges_m1(newNodeKMerInfo.traverseEdges.begin()+1, newNodeKMerInfo.traverseEdges.end());
						//traverseEdges_m1.erase (traverseEdges_m1.begin(),traverseEdges_m1.begin()+1);

						assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
						string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

						edgeTargetCache[newNodeKMerInfo.traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
						newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;
					}
				}
			}
		}

		originalGraphAssignedKMers.clear();

		int newNodeC = 0;
		int targetEdgeC = 0;

		for(map<Node*, map<string, vector<Edge*> > >::iterator targetNodeIt = edgeTargetCache.begin(); targetNodeIt != edgeTargetCache.end(); targetNodeIt++)
		{
			for( map< string, vector<Edge*> >::iterator targetStringIt = targetNodeIt->second.begin(); targetStringIt != targetNodeIt->second.end(); targetStringIt++)
			{
				// In der inneren Schleife haben wir alle Edges, die zur selben Node im original Graph fuehren und die
				// das auf den letzten x-1 Edges auf demselben Weg tun -- i.e. die in dieselbe Node fuehren sollen,
				// weil die darauffolgenden kMers identisch sind

				newNodeC++;

				string targetString = targetStringIt->first;
				vector<Edge*> attachedEdges = targetStringIt->second;
				assert(attachedEdges.size()>0);

				// Wir generieren eine Node und haengen sie an alle Kanten
				Node* kMerTargetNode = new Node();
				kMerTargetNode->level = level+1;
				kMerTargetNode->terminal = (level == (levels - 1 - kMerSize)) ? true : false; // TODO korrekt?
				kMerGraph->registerNode(kMerTargetNode, level+1);

				int num_of_PGFprotected = 0;
				for(unsigned int idx = 0; idx < attachedEdges.size(); idx++)
				{
					attachedEdges.at(idx)->To = kMerTargetNode;
					kMerTargetNode->Incoming_Edges.insert(attachedEdges.at(idx));
					targetEdgeC++;
					if(attachedEdges.at(idx)->pgf_protect)
					{
						num_of_PGFprotected++;
					}
				}
				assert((num_of_PGFprotected == 0) || (num_of_PGFprotected == 1));

				// Aus der ersten Edge, die wir haben, generieren wir die Info fuer zukuenftige kMer-Kanten,
				// die wir an die 2. Node des originals Graphs heften

				Edge* firstEdge = attachedEdges.at(0);
				assert(newEdgeNodeInfos.count(firstEdge) > 0);
				kMerInfo firstEdgeKMerInfo = newEdgeNodeInfos[firstEdge];
				assert((int)firstEdgeKMerInfo.kMer_deCoded.size() <= kMerSize);

				vector<unsigned char> km1Mer_coded = firstEdgeKMerInfo.kMer_coded;
				vector<string> km1Mer_deCoded  = firstEdgeKMerInfo.kMer_deCoded;

				//vector<Edge*> traverseEdges_m1 = firstEdgeKMerInfo.traverseEdges;
				vector<Edge*> traverseEdges_m1(firstEdgeKMerInfo.traverseEdges.begin()+1, firstEdgeKMerInfo.traverseEdges.end());
				string traverseEdges_m1_string = firstEdgeKMerInfo.traverseEdges_string.substr(pointerStrintLength);
				assert(traverseEdges_m1_string.length() == (traverseEdges_m1.size()*pointerStrintLength));

				if(firstEdgeKMerInfo.gapEdge == false)
				{
					km1Mer_coded.erase (km1Mer_coded.begin(),km1Mer_coded.begin()+1);
					km1Mer_deCoded.erase (km1Mer_deCoded.begin(),km1Mer_deCoded.begin()+1);
				}
				else
				{
					//assert(attachedEdges.size()==1);
				}

				//traverseEdges_m1.erase (traverseEdges_m1.begin(),traverseEdges_m1.begin()+1);

				kMerAtNode infoForNode2;
				infoForNode2.km1Mer_coded = km1Mer_coded;
				infoForNode2.km1Mer_deCoded = km1Mer_deCoded;
				infoForNode2.traverseEdges = traverseEdges_m1;
				infoForNode2.lastNewNode = kMerTargetNode;
				infoForNode2.traverseEdges_string = traverseEdges_m1_string;
				infoForNode2.allPGF = (num_of_PGFprotected > 0) ? true : false;
				
				Node* originalGraphNode2 = firstEdgeKMerInfo.traverseEdges.at(0)->To;
				assert(g->Nodes.count(originalGraphNode2) > 0);

				if((! quiet) && ((level % 100) == 0) && (level > 6000000))
				{
					//cout << "\tFound originalGraphNode2 " << originalGraphNode2 << " with " << attachedEdges.size() << " attached edges!\n";
				}

				originalGraphAssignedKMers[originalGraphNode2].push_back(infoForNode2);

				/*
				for(unsigned int idx = 0; idx < attachedEdges.size(); idx++)
				{
					Edge* e = attachedEdges.at(idx);
					Node* shouldbeNode2 = newEdgeNodeInfos[e].traverseEdges.at(0)->To;
					assert(shouldbeNode2 == originalGraphNode2);
				}
				*/

			}
		}

		/*
		if(level < (levels - 1 - kMerSize))
		{
			vector<Node*> originalNodesNextLevel = vector<Node*>(g->NodesPerLevel.at(level+1).begin(), g->NodesPerLevel.at(level+1).end());

			bool foundAll = true;
			for(vector<Node*>::iterator originalNodeIt = originalNodesNextLevel.begin(); originalNodeIt != originalNodesNextLevel.end(); originalNodeIt++)
			{
				Node* originalNode = *originalNodeIt;
				foundAll = foundAll && (originalGraphAssignedKMers.count(originalNode) > 0);
				if((! quiet) && ((level % 100) == 0) && (level > 6000000))
				{
					// cout << "\t Level " << level << ", search for node at next level " << originalNode << " with result " << (originalGraphAssignedKMers.count(originalNode) > 0) << "\n";
				}
			}
			assert(foundAll = true);
		}
		*/

		assert(min_span_level <= max_span_level);

		string newLocusID = "L"+Utilities::ItoStr(level);

		if((! quiet) && ((level % 100) == 0))
		{
			cout << Utilities::timestamp() << "\t\t added " << newNodeC << " nodes, " << targetEdgeC << " edges, and have " << kMerGraph->CODE.getAlleles(newLocusID).size() << " encoded kMers \n";
			cout << "\t\t\t target node span original graph: " << min_span_level << " - " << max_span_level << "\n";
		}
	}

	kMerGraph->kMerSize = kMerSize;

	for(set<Edge*>::iterator edgeIt = kMerGraph->Edges.begin(); edgeIt != kMerGraph->Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;
		string kMer = kMerGraph->CODE.deCode(e->locus_id, e->largeEmission);
		e->largeEmission2 = e->largeEmission;

		bool foundStar = false;
		if(kMer != "_")
		{
			assert((int)kMer.size() == kMerSize);
			for(int i = 0; i < kMerSize; i++)
			{
				if(kMer.at(i) == '*')
				{
					foundStar = true;
					break;
				}
			}
			if(foundStar)
			{
				e->largeEmission = kMerGraph->CODE.doCode(e->locus_id, "*");
			}
		}
	}

	kMerGraph->removeStarPaths();

	if(wantPGFprotection)
	{
		for(unsigned int level = 0; level < (kMerGraph->NodesPerLevel.size()-1); level++)
		{
			bool foundPGF = false;
			for(set<Node*>::iterator nodeIt = kMerGraph->NodesPerLevel.at(level).begin(); nodeIt != kMerGraph->NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;

				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					foundPGF = (foundPGF || (bool)e->pgf_protect);
					if((level > 0) && (level < (kMerGraph->NodesPerLevel.size()-2)))
					{
						if((bool)e->pgf_protect)
						{
							bool found_preceding_PGF = false;
							bool found_following_PGF = false;
							for(set<Edge*>::iterator beforeIt = e->From->Incoming_Edges.begin(); beforeIt != e->From->Incoming_Edges.end(); beforeIt++)
							{
								if((bool)(*beforeIt)->pgf_protect)
								{
									found_preceding_PGF = true;
									break;
								}
							}
							if(! found_preceding_PGF)
							{
								cout << "Level " << level << " of " << (kMerGraph->NodesPerLevel.size()-2) << ", edge " << e << " pgf protect: " << (bool)e->pgf_protect << " starting at level " << e->From->level << "\n";
								for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
								{
									cout << "\t following " << (*afterIt) << ": " << (bool)(*afterIt)->pgf_protect << "\n" << flush;
									assert( e->To->level == (level + 1));
								}
								
								cout << "We now look at all edges from the last level " << (level-1) << ":\n";
								for(set<Node*>::iterator nodeIt2 = kMerGraph->NodesPerLevel.at(level-1).begin(); nodeIt2 != kMerGraph->NodesPerLevel.at(level-1).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}		
								
								cout << "We now look at all edges from the this level " << (level) << ":\n";
								for(set<Node*>::iterator nodeIt2 = kMerGraph->NodesPerLevel.at(level).begin(); nodeIt2 != kMerGraph->NodesPerLevel.at(level).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}
								
								cout << "We now look at all edges from the next level " << (level+1) << ":\n";
								for(set<Node*>::iterator nodeIt2 = kMerGraph->NodesPerLevel.at(level+1).begin(); nodeIt2 != kMerGraph->NodesPerLevel.at(level+1).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}						
							}
							
							assert(found_preceding_PGF);

							for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
							{
								if((bool)(*afterIt)->pgf_protect)
								{
									found_following_PGF = true;
									break;
								}
							}
							assert(found_following_PGF);
						}
					}
				}
			}
			assert(foundPGF);
		}
	}

	kMerGraph->checkConsistency(false);
	return kMerGraph;
}

vector<kMerInfo> forwardScan(Node* start, int limit, int firstEdgeGap)
{
	if(pointerStrintLength == -1)
	{
		std::ostringstream nullPointerString;
		nullPointerString << setw(15) << (void const *) 0;
		pointerStrintLength = nullPointerString.str().length();
	}

	vector<kMerInfo> forReturn = forwardScanRec(start, 0, 0, limit, firstEdgeGap);
	for(vector<kMerInfo>::iterator kMerIt = forReturn.begin(); kMerIt != forReturn.end(); kMerIt++)
	{
		assert((int)(kMerIt->kMer_coded.size()) == limit);
		assert((int)(kMerIt->kMer_deCoded.size()) == limit);
		reverse(kMerIt->kMer_coded.begin(), kMerIt->kMer_coded.end());
		reverse(kMerIt->kMer_deCoded.begin(), kMerIt->kMer_deCoded.end());
		reverse(kMerIt->traverseEdges.begin(), kMerIt->traverseEdges.end());

		std::ostringstream o_traverseEdges_string;
		bool allPGF = true;
		for(int i = 0; i < (int)kMerIt->traverseEdges.size(); i++)
		{
			o_traverseEdges_string << setw(15) << (void const * ) kMerIt->traverseEdges.at(i);
			allPGF = (allPGF && (bool)kMerIt->traverseEdges.at(i)->pgf_protect);
		}
		kMerIt->allPGF = allPGF;

		string traverseEdges_string = o_traverseEdges_string.str();
		if(traverseEdges_string.length() != (kMerIt->traverseEdges.size()*pointerStrintLength))
		{
			cout << "One pointer as string: " << setw(15) << (void const * ) kMerIt->traverseEdges.at(0) << "\n";
			cout << "kMerIt->traverseEdges.size(): " << kMerIt->traverseEdges.size() << "\n";
			cout << "pointerStrintLength: " << pointerStrintLength << "\n";
			cout << "traverseEdges_string.length(): " << traverseEdges_string.length() << "\n";
			cout << "traverseEdges_string: " << traverseEdges_string << "\n";
		}
		assert(traverseEdges_string.length() == (kMerIt->traverseEdges.size()*pointerStrintLength));


		kMerIt->traverseEdges_string = traverseEdges_string;

	}
	return forReturn;
}

vector<kMerInfo> forwardScanRec(Node* currentNode, int depth, int realdepth, int limit, int firstEdgeGap)
{
	vector<kMerInfo> returnMers;

	if(limit == depth)
	{
		kMerInfo i;
		i.p = 1;
		i.gapEdge = false;
		returnMers.push_back(i);
		return returnMers;
	}
	else
	{
		for(set<Edge*>::iterator eIt = currentNode->Outgoing_Edges.begin(); eIt != currentNode->Outgoing_Edges.end(); eIt++)
		{
			Edge* e = *eIt;
			Node* targetNode = e->To;

			// TODO maybe we have such nodes, then we need to ignore them
			assert(currentNode->Sum_Outgoing() != 0);
			double edgeP = e->count/currentNode->Sum_Outgoing();
			assert(edgeP >= 0);
			assert(edgeP <= 1);

			vector<kMerInfo> localReturn;
			bool neverTrue = false;
			if(realdepth == 0)
			{
				if(firstEdgeGap == 1)
				{
					if(currentNode->g->CODE.deCode(e->locus_id, e->emission) != "_")
					{
						neverTrue = true;
						continue;
					}
				}
				else if (firstEdgeGap == -1)
				{
					if(currentNode->g->CODE.deCode(e->locus_id, e->emission) == "_")
					{
						neverTrue = true;
						continue;
					}
				}
			}
			assert(neverTrue == false);

			if(currentNode->g->CODE.deCode(e->locus_id, e->emission) != "_")
			{
				localReturn = forwardScanRec(targetNode, depth+1, realdepth + 1, limit, firstEdgeGap);
				for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
				{
					kMerIt->kMer_coded.push_back(e->emission);
					kMerIt->kMer_deCoded.push_back(currentNode->g->CODE.deCode(e->locus_id, e->emission));

					kMerIt->p = kMerIt->p*edgeP;
					kMerIt->traverseEdges.push_back(e);
				}
			}
			else
			{
				// this is a simple look-forward to check whether we can non-recursively extend
				if((targetNode->Outgoing_Edges.size()==1) && (currentNode->g->CODE.deCode((*targetNode->Outgoing_Edges.begin())->locus_id, (*targetNode->Outgoing_Edges.begin())->emission) == "_"))
				{
					vector<Edge*> traversedEdges;
					Node* currentTargetNode = targetNode;
					traversedEdges.push_back(e);
					while((currentTargetNode->Outgoing_Edges.size()==1) && (currentNode->g->CODE.deCode((*currentTargetNode->Outgoing_Edges.begin())->locus_id, (*currentTargetNode->Outgoing_Edges.begin())->emission) == "_"))
					{
						traversedEdges.push_back((*currentTargetNode->Outgoing_Edges.begin()));
						currentTargetNode = (*currentTargetNode->Outgoing_Edges.begin())->To;
					}
					reverse(traversedEdges.begin(),traversedEdges.end());

					localReturn = forwardScanRec(currentTargetNode, depth, realdepth + 1, limit, firstEdgeGap);
					for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
					{
						kMerIt->p = kMerIt->p*edgeP;
						kMerIt->traverseEdges.insert(kMerIt->traverseEdges.end(), traversedEdges.begin(), traversedEdges.end());
					}
				}
				else
				{
					localReturn = forwardScanRec(targetNode, depth, realdepth + 1, limit, firstEdgeGap);
					for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
					{
						kMerIt->p = kMerIt->p*edgeP;
						kMerIt->traverseEdges.push_back(e);
					}
				}
			}

			returnMers.insert(returnMers.end(), localReturn.begin(), localReturn.end());

		}

		return returnMers;
	}
}

string alleleAndLocusIdentifier(string locus, string allele)
{
	map<string, string> locusTranslation;
	locusTranslation["HLAA"] = "A";
	locusTranslation["HLAB"] = "B";
	locusTranslation["HLAC"] = "C";
	locusTranslation["HLADQA"] = "DQA1";
	locusTranslation["HLADQB"] = "DQB1";
	locusTranslation["HLADRB"] = "DRB1";
	if(!(locusTranslation.count(locus) > 0))
	{
		errEx("Nextgen.cpp, alleleAndLocusIdentifier: No locus translation table for locus "+locus);
	}
	assert(locusTranslation.count(locus) > 0);
	assert(allele.length()==4);

	string newAllele = allele.substr(0,2)+":"+allele.substr(2,2);
	return locusTranslation[locus]+"*"+newAllele;
}

map<string, string> readHLAalleleAlignments(string file, vector<int>& fourNumbers)
{
	map<string, string> alignments;

	ifstream alignmentsFile;
	alignmentsFile.open (file.c_str(), ios::in);
	int lineCounter = 0;
	if(alignmentsFile.is_open())
	{
		string line;
		while(alignmentsFile.good())
		{
			lineCounter++;
			getline (alignmentsFile, line);
			if(lineCounter > 4)
			{
				Utilities::eraseNL(line);

				if(line.length() == 0)
					continue;

				vector<string> fields = Utilities::split(line, ' ');
				if(fields.size() != 2)
					errEx("Strange format in alignments file. Expect two fields (from line 5 onwards): allele name and alignment string, separated by a space. Problem occured in line: "+line);

				alignments[fields.at(0)] = fields.at(1);
			}
			else
			{
				Utilities::eraseNL(line);
				if(line.length() == 0)
					continue;
				vector<string> fields = Utilities::split(line, ' ');
				if(fields.size() != 2)
					errEx("Strange format in alignments file. Expect two fields (from line 5 onwards): allele name and alignment string, separated by a space. Problem occured in line: "+line);
				fourNumbers.push_back(Utilities::StrtoI(fields.at(1)));
			}
		}
		alignmentsFile.close();
	}
	else
	{
		errEx("Cannot open HLA allele alignment file: "+file);
	}

	return alignments;
}


