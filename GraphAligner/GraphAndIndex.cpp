/*
 * GraphAndIndex.cpp
 *
 *  Created on: 18.06.2013
 *      Author: AlexanderDilthey
 */

#include <assert.h>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "GraphAndIndex.h"

#include "../NextGen/NextGen.h"
#include "../Utilities.h"
#include "../hash/sequence/basic.h"

#include <queue>

using namespace std;

GraphAndIndex::GraphAndIndex(Graph* graph, int k) {
	kMerSize = k;
	g = graph;
	Index();
}

std::vector<kMerChain> GraphAndIndex::findChains(std::string sequence)
{
	std::vector<kMerChain> forReturn = findChainsPlusStrand(sequence);
	return forReturn;

	// todo activate later!

	std::string sequence_reverse = seq_reverse_complement(sequence);
	std::vector<kMerChain> chains_reverse = findChainsPlusStrand(sequence_reverse);
	for(unsigned int cI = 0; cI < chains_reverse.size(); cI++)
	{
		kMerChain kC = chains_reverse.at(cI);
		kC.sequence_reverse = true;
		forReturn.push_back(kC);
	}
	return forReturn;
}

std::vector<kMerChain> GraphAndIndex::findChainsPlusStrand(std::string sequence)
{
	std::vector<kMerChain> runningChains;

	bool superquiet = false;

	if(! superquiet)
		std::cout  << Utilities::timestamp() << "GraphAndIndex::findChainsPlusStrand(..): Enter function.\n" << std::flush;

	std::vector<std::string> kMers = partitionStringIntokMers(sequence, kMerSize);

	bool verbose = false;

	if(! superquiet)
		std::cout  << Utilities::timestamp() << "GraphAndIndex::findChainsPlusStrand(..): Process " << kMers.size() << " kMers.\n" << std::flush;

	size_t positions_unknown_kMers = 0;

	for(unsigned int kI = 0; kI < kMers.size(); kI++)
	{
		if((! superquiet) && ((kI % 10000) == 0))
			std::cout << "\r" << "kMer " << kI << " / " << kMers.size() << std::flush;

		std::string& kMer = kMers.at(kI);
		std::vector<std::pair<int, int> > kMer_positions = queryIndex(kMer);

		if(kMer_positions.size() == 0)
		{
			positions_unknown_kMers++;
		}

		if(verbose)
			std::cout << "kMer " << kI << " [" << kMer << "], found positions: " << kMer_positions.size() << ", investigate " << runningChains.size() << " existing chains.\n";

		int pos_in_sequence_begin = kI;
		int pos_in_sequence_end = kI+kMerSize-1;

		std::vector<bool> kMer_positions_covered;
		kMer_positions_covered.resize(kMer_positions.size(), false);

		unsigned int alreadyExistingChains = runningChains.size();
		for(unsigned int chainI = 0; chainI < alreadyExistingChains; chainI++)
		{
			if(verbose)
				std::cout << "\tChain" << chainI << "\n";

			kMerChain& existingChain = runningChains.at(chainI);

			std::vector<int> compatiblePosExtensions;
			for(unsigned int posI = 0; posI < kMer_positions.size(); posI++)
			{
				std::pair<int, int>& thisPosition = kMer_positions.at(posI);

				if(verbose)
					std::cout << "\t\tPosition alternative " << posI << ": " << thisPosition.first << " - " << thisPosition.second <<  "\n";


				int diff_in_starting_positions_graph = thisPosition.first - existingChain.lastAttachedkMer_graph_level;
				int diff_in_starting_positions_seq = kI - existingChain.lastAttachedkMer_level;
				int diff_in_stopping_positions_graph = thisPosition.second - existingChain.graph_lastLevel;

				if(verbose)
					std::cout << "\t\t\tVs existing this chain" << chainI << ", diff_in_starting_positions_graph: " << diff_in_starting_positions_graph << ", diff_in_starting_positions_seq: " << diff_in_starting_positions_seq << ", diff_in_stopping_positions_graph: " << diff_in_stopping_positions_graph << "\n";

				//if((diff_in_starting_positions > 0) && (diff_in_starting_positions < startingPositionTolerance))
				if(
						((diff_in_starting_positions_graph == 1) && (diff_in_starting_positions_seq == 1) && (diff_in_stopping_positions_graph >= 1)) ||
						((diff_in_starting_positions_graph > 1) && (diff_in_starting_positions_seq == 1) && (diff_in_stopping_positions_graph == 1))
				)
				{
					compatiblePosExtensions.push_back(posI);

					if(verbose)
						std::cout << "\t\t\t\tcompatible!\n";
				}
			}

			if(verbose)
				std::cout << "\t\tFound positions for attaching to this chain: " << compatiblePosExtensions.size() << "\n";

			if(compatiblePosExtensions.size() == 1)
			{
				int compatiblePositionI = compatiblePosExtensions.at(0);
				std::pair<int, int>& thisPosition = kMer_positions.at(compatiblePositionI);

				kMerChain& existingChain = runningChains.at(chainI);

				existingChain.graph_lastLevel = thisPosition.second;
				existingChain.sequence_end = pos_in_sequence_end;
				existingChain.lastAttachedkMer_level = kI;
				existingChain.lastAttachedkMer_graph_level = thisPosition.first;

				kMer_positions_covered.at(compatiblePositionI) = true;
			}
			else if(compatiblePosExtensions.size() > 1)
			{
				kMerChain originalChain = existingChain;

				for(unsigned int compatiblePositionII = 0; compatiblePositionII < compatiblePosExtensions.size(); compatiblePositionII++)
				{
					int compatiblePositionI = compatiblePosExtensions.at(compatiblePositionII);
					std::pair<int, int>& thisPosition = kMer_positions.at(compatiblePositionI);

					kMerChain newChain = originalChain;
					newChain.graph_lastLevel = thisPosition.second;
					newChain.sequence_end = pos_in_sequence_end;
					newChain.lastAttachedkMer_level = kI;
					newChain.lastAttachedkMer_graph_level = thisPosition.first;

					if(compatiblePositionII == 0)
					{
						runningChains.at(chainI) = newChain;
					}
					else
					{
						runningChains.push_back(newChain);
					}

					kMer_positions_covered.at(compatiblePositionI) = true;
				}
			}
		}

		for(unsigned int posI = 0; posI < kMer_positions.size(); posI++)
		{
			if(! kMer_positions_covered.at(posI))
			{
				std::pair<int, int>& thisPosition = kMer_positions.at(posI);

				if(verbose)
					std::cout << "\tPosition alternative " << posI << " still uncovered, produce new chain!\n";

				kMerChain newChain;
				newChain.sequence_begin = pos_in_sequence_begin;
				newChain.sequence_end = pos_in_sequence_end;
				newChain.graph_firstLevel = thisPosition.first;
				newChain.graph_lastLevel = thisPosition.second;
				newChain.sequence_reverse = false;
				newChain.matchedkMers = 1;
				newChain.lastAttachedkMer_level = kI;
				newChain.lastAttachedkMer_graph_level = thisPosition.first;
				runningChains.push_back(newChain);
			}
		}
	}

	if(! superquiet)
		std::cout << "\n";

	for(unsigned int cI = 0; cI < runningChains.size(); cI++)
	{
		kMerChain& kC = runningChains.at(cI);
		assert(kC.sequence_begin >= 0);
		assert(kC.sequence_end < (int)sequence.length());
	}

	if(! superquiet)
		std::cout  << Utilities::timestamp() << "GraphAndIndex::findChainsPlusStrand(..): Done!\n" << std::flush;


//	std::cout << "findChainsPlusStrand(..): Found " << runningChains.size() << " chains for input sequence of " << sequence << " characters; " << positions_unknown_kMers << " positions with unknown kMer.\n";

	return runningChains;
}

void GraphAndIndex::Index()
{
	LargeGraph* kMerGraph = new LargeGraph();
	bool quiet = true;
	bool superQuiet = false;

	int nonQuiet_levelModulo = 1;

	std::ostringstream nullPointerString;
	nullPointerString << setw(15) << (void const *) kMerGraph;
	int pointerStrintLength = nullPointerString.str().length();

	int levels = g->NodesPerLevel.size();
	vector<string> loci = g->getAssignedLoci();

	set<Node*> nodesLastLevel = g->NodesPerLevel[levels-1];
	string lastLocusID = loci.at(loci.size()-1);
	Node* classicalN0 = *(g->NodesPerLevel.at(0).begin());

	std::map<std::string, std::map<std::string, std::pair<int, int> > > _kMerPositions;
	auto kMer_bookkeeping = [&](std::string kMer, int firstLevel, int lastLevel) -> void {
		assert((int)kMer.length() == this->kMerSize);
		assert(firstLevel >= 0);
		assert(lastLevel >= 0);
		assert(lastLevel > firstLevel);

		std::string secondKey = Utilities::ItoStr(firstLevel)+"--"+Utilities::ItoStr(lastLevel);
		std::pair<int, int> p(firstLevel, lastLevel);
		_kMerPositions[kMer][secondKey] = p;

		if(! quiet)
		{
			std::cout << "+ log kMer " << kMer << " from " << firstLevel << " to " << lastLevel << "\n" << std::flush;
		}
	};

	// make sure that last graph level has no gaps and that all edges are correctly labelled
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

	map<Node*, vector< kMerAtNode > > originalGraphAssignedKMers;

	std::set<Node*> generated_nodes;
	std::set<Edge*> generated_edges;

	for(int level = 0; level <= (levels - 1 - kMerSize); level++)
	{
		if(! superQuiet)
		{
			if((level % 10000) == 0)
			{
				std::cout << "\r" << level << " / " << (levels - 1 - kMerSize) << std::flush;
			}
		}

		int min_span_level = -1;
		int max_span_level = -1;

		if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
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
				assert((int)kMerIt->kMer_deCoded.size() == kMerSize);

				int firstLevel = kMerIt->traverseEdges.at(0)->From->level;
				int lastLevel = kMerIt->traverseEdges.at(kMerIt->traverseEdges.size()-1)->To->level;
				kMer_bookkeeping(kMer_string, firstLevel, lastLevel);

				Edge* kMerEdge = new Edge();
				assert(g->CODE.deCode(kMerIt->traverseEdges.at(0)->locus_id, kMerIt->traverseEdges.at(0)->emission) != "_");

				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);
				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;

				generated_edges.insert(kMerEdge);
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
				assert(kMerIt->kMer_deCoded.size() == (kMerSize -1));
				assert((int)kMer_string.size() == kMerSize-1);

				string newLocusID = "L"+Utilities::ItoStr(level);

				// cout << kMer_string << "\n" << flush;

				Edge* kMerEdge = new Edge();
				kMerIt->gapEdge = true;

				assert(g->CODE.deCode(kMerIt->traverseEdges.at(0)->locus_id, kMerIt->traverseEdges.at(0)->emission) == "_");

				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);

				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;

				generated_edges.insert(kMerEdge);
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

				// only if next assertion fails...
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
					if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
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

							string kMer_string = "_";

							Edge* kMerEdge = new Edge();

							BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

							assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
							string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

							assert(g->Nodes.count(newNodeKMerInfo.traverseEdges.back()->To) > 0);

							edgeTargetCache[newNodeKMerInfo.traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
							newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;

							generated_edges.insert(kMerEdge);

						}
						else
						{
							if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
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

									if(!(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1)))
									{
										std::cerr << "!(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1))\n";
										std::cerr << "newNodeKMerInfo.kMer_deCoded.size(): " << newNodeKMerInfo.kMer_deCoded.size() << "\n";
										std::cerr << "(kMerSize-1): " << (kMerSize-1) << "\n" << std::flush;
										std::cerr << "newNodeKMerInfo.kMer_deCoded: " << Utilities::join(newNodeKMerInfo.kMer_deCoded, " ") << "\n" << std::flush;

									}

									assert(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1));


									newNodeKMerInfo.kMer_coded.push_back(attachKMerIt->kMer_coded.at(0));
									newNodeKMerInfo.kMer_deCoded.push_back(attachKMerIt->kMer_deCoded.at(0));
									newNodeKMerInfo.traverseEdges.insert(newNodeKMerInfo.traverseEdges.end(), attachKMerIt->traverseEdges.begin(), attachKMerIt->traverseEdges.end());
									newNodeKMerInfo.traverseEdges_string = newNodeKMerInfo.traverseEdges_string.append(attachKMerIt->traverseEdges_string);
									newNodeKMerInfo.allPGF = false; //(BasisForNewEdges.allPGF && attachKMerIt->allPGF);

									string kMer_string = Utilities::join(newNodeKMerInfo.kMer_deCoded, "");
									assert((int)kMer_string.size() == kMerSize);

									int firstLevel = newNodeKMerInfo.traverseEdges.at(0)->From->level;
									int lastLevel = newNodeKMerInfo.traverseEdges.at(newNodeKMerInfo.traverseEdges.size()-1)->To->level;
									kMer_bookkeeping(kMer_string, firstLevel, lastLevel);

									Edge* kMerEdge = new Edge();
									BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

									for(int lI = 0; lI < (int)attachKMerIt->traverseEdges.size(); lI++)
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

									generated_edges.insert(kMerEdge);
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

						if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
							cout << " -- gap length " << BasisForNewEdges.traverseEdges.size() << "\n";

						kMerInfo newNodeKMerInfo;
						newNodeKMerInfo.gapEdge = true;
						newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
						newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
						newNodeKMerInfo.p = 1;
						newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
						newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;

						string kMer_string = "_";

						Edge* kMerEdge = new Edge();
						BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

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
				generated_nodes.insert(kMerTargetNode);
				kMerTargetNode->level = level+1;
				kMerTargetNode->terminal = (level == (levels - 1 - kMerSize)) ? true : false; // TODO korrekt?

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


				kMerAtNode infoForNode2;
				infoForNode2.km1Mer_coded = km1Mer_coded;
				infoForNode2.km1Mer_deCoded = km1Mer_deCoded;
				infoForNode2.traverseEdges = traverseEdges_m1;
				infoForNode2.lastNewNode = kMerTargetNode;
				infoForNode2.traverseEdges_string = traverseEdges_m1_string;

				Node* originalGraphNode2 = firstEdgeKMerInfo.traverseEdges.at(0)->To;
				assert(g->Nodes.count(originalGraphNode2) > 0);

				if((! quiet) && ((level % nonQuiet_levelModulo) == 0) && (level > 6000000))
				{
					//cout << "\tFound originalGraphNode2 " << originalGraphNode2 << " with " << attachedEdges.size() << " attached edges!\n";
				}

				if(! quiet)
				{
					std::cout << "Attach to node " << originalGraphNode2 << " at level " << originalGraphNode2->level << " new kMerAtNode object!\n";
					std::cout << "\tkm1Mer_deCoded: " << Utilities::join(infoForNode2.km1Mer_deCoded, " ") << "\n";
					std::cout << "\ttraverseEdges.size(): " << infoForNode2.traverseEdges.size() << "\n";
					std::cout << "\tinfoForNode2.traverseEdges.back()->To->level: " << infoForNode2.traverseEdges.back()->To->level << "\n" << std::flush;
					std::cout << "\tfirstEdgeKMerInfo.gapEdge: " << firstEdgeKMerInfo.gapEdge << "\n" << std::flush;

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

		assert(min_span_level <= max_span_level);

		string newLocusID = "L"+Utilities::ItoStr(level);
		if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
		{
			cout << "\t\t added " << newNodeC << " nodes, " << targetEdgeC << " edges, and have " << kMerGraph->CODE.getAlleles(newLocusID).size() << " encoded kMers \n";
			cout << "\t\t\t target node span original graph: " << min_span_level << " - " << max_span_level << "\n";
		}
	}

	if(! superQuiet)
	{
		std::cout << "\n" << std::flush;
	}


	for(std::set<Node*>::iterator nIt = generated_nodes.begin(); nIt != generated_nodes.end(); nIt++)
	{
		delete(*nIt);
	}

	for(std::map<std::string, std::map<std::string, std::pair<int, int> > >::iterator kMerIt = _kMerPositions.begin(); kMerIt != _kMerPositions.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		std::map<std::string, std::pair<int, int> >& thisMerPositions = _kMerPositions.at(kMer);
		for(std::map<std::string, std::pair<int, int> >::iterator kMerPosIt = thisMerPositions.begin(); kMerPosIt != thisMerPositions.end(); kMerPosIt++)
		{
			kMerPositions[kMer].push_back(kMerPosIt->second);
		}
	}
}

void GraphAndIndex::printIndex()
{
	for(std::map<std::string, std::vector< std::pair<int, int> > >::iterator kMerPosIt = kMerPositions.begin(); kMerPosIt != kMerPositions.end(); kMerPosIt++)
	{
		std::string kMer = kMerPosIt->first;

		std::cout << "kMer " << kMer << "\n";

		std::vector< std::pair<int, int> >& kMerPos = kMerPositions.at(kMer);
		for(unsigned int i = 0; i < kMerPos.size(); i++)
		{
			std::cout << " - " << kMerPos.at(i).first << " " << kMerPos.at(i).second << "\n";
		}
	}

	std::cout << std::flush;
}

std::vector< std::pair<int, int> > GraphAndIndex::queryIndex(std::string kMer)
{
	if(kMerPositions.count(kMer) == 0)
	{
		return std::vector< std::pair<int, int> >();
	}
	else
	{
		return kMerPositions.at(kMer);
	}
}

void _printChains(std::vector<kMerChain>& chains)
{
	std::cout << "Chains: " << chains.size() << "\n\n";
	for(unsigned int cI = 0; cI < chains.size(); cI++)
	{
		std::cout << "\t" << "chain #" << cI << "\n\n";
		std::cout << "\t\t" << "graph_firstLevel: " << chains.at(cI).graph_firstLevel << "\n";
		std::cout << "\t\t" << "graph_lastLevel: " << chains.at(cI).graph_lastLevel << "\n";
		std::cout << "\t\t" << "sequence_begin: " << chains.at(cI).sequence_begin << "\n";
		std::cout << "\t\t" << "sequence_end: " << chains.at(cI).sequence_end << "\n\n";
		std::cout << "\t\t" << "sequence_reverse: " << chains.at(cI).sequence_reverse << "\n\n";
		std::cout << "\t\t" << "matched kMers: " << chains.at(cI).matchedkMers << "\n\n";

	}
	std::cout << std::flush;
}
