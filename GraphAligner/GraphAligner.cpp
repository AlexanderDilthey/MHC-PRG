/*
 * GraphAligner.cpp
 *
 *  Created on: 18.06.2013
 *      Author: AlexanderDilthey
 */

#include "GraphAligner.h"

#include "../NextGen/Validation.h"
#include "GraphAndIndex.h"

#include <assert.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <limits>
#include <list>
#include <omp.h>

#include "VirtualNW.h"

#include "../Utilities.h"
#include "../Graph/Edge.h"

#include "AlignerTests.h"

int threads = 40;

GraphAligner::GraphAligner(Graph* graph, int k) : g(graph), kMerSize(k), gI(g, k) {
	S_match = 2;
	S_mismatch = -5;
	S_gap = -2;
	S_graphGap = 0;

	S_openGap = -4;
	S_extendGap = -2;

	S_longGap_open = 0;
	S_longGap_perGap = -2;

	// enumerate states at levels
	unsigned int levels = g->NodesPerLevel.size();
	nodesPerLevel_ordered.resize(levels);
	nodesPerLevel_ordered_rev.resize(levels);
	for(unsigned int levelI = 0; levelI < levels; levelI++)
	{
		nodesPerLevel_ordered.at(levelI) = std::vector<Node*>(g->NodesPerLevel.at(levelI).begin(), g->NodesPerLevel.at(levelI).end());
		for(unsigned int nodeI = 0; nodeI < nodesPerLevel_ordered.at(levelI).size(); nodeI++)
		{
			nodesPerLevel_ordered_rev.at(levelI)[nodesPerLevel_ordered.at(levelI).at(nodeI)] = nodeI;
		}
	}
}



std::vector<NWPath*> GraphAligner::findVirtualNWExactSeeds(std::string& sequence, int sequenceStart, int sequenceStop, int graph_firstLevel, int graph_lastLevel, int minimumInexactMatchLength)
{
	assert(graph_firstLevel < graph_lastLevel);
	assert(sequenceStart <= sequenceStop);
	assert(sequenceStart >= 0);
	assert(sequenceStop < (int)sequence.length());

	if(minimumInexactMatchLength < 0)
	{
		minimumInexactMatchLength = kMerSize;
	}

	class candidatesPaths {
	public:
		Node* currentTip;
		int currentSequencePosition;
		std::vector<Edge*> edges;
		std::string aligned_sequence;
		std::vector<std::vector<int>> edge_coordinates;
	};

	std::vector<candidatesPaths> candidates;
	std::vector<candidatesPaths> matches;

	std::set<Node*> startingNodes = g->NodesPerLevel.at(graph_firstLevel);
	for(std::set<Node*>::iterator startNodeIt = startingNodes.begin(); startNodeIt != startingNodes.end(); startNodeIt++)
	{
		candidatesPaths p;
		p.currentTip = *startNodeIt;
		p.currentSequencePosition = sequenceStart;
		std::vector<int> coordinates;
		coordinates.push_back(graph_firstLevel);
		coordinates.push_back(sequenceStart);
		coordinates.push_back(nodesPerLevel_ordered_rev.at(graph_firstLevel).at(*startNodeIt));
		p.edge_coordinates.push_back(coordinates);
		candidates.push_back(p);
	}

	assert(candidates.size() > 0);

	std::vector<candidatesPaths> fallbackMatches;
	while(candidates.size() > 0)
	{
		// std::cout << candidates.size() << "\n" << std::flush;

		candidatesPaths c = candidates.at(candidates.size() - 1);
		unsigned int s_before = candidates.size();
		candidates.erase(candidates.end() - 1);
		assert((int)candidates.size() == (s_before - 1));
		bool couldExtend = false;

		if(!(((int)c.currentTip->level >= graph_firstLevel) && ((int)c.currentTip->level < graph_lastLevel)))
		{
			continue;
		}

		assert(c.currentSequencePosition < (int)sequence.length());
		std::string sequenceCharacter = sequence.substr(c.currentSequencePosition, 1);

		std::set<Edge*> edges = c.currentTip->Outgoing_Edges;
		std::set<Edge*> compatibleEdges;

 		for(std::set<Edge*>::iterator eIt = edges.begin(); eIt != edges.end(); eIt++)
		{
			Edge* e = *eIt;
			std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
			assert(edgeEmission.length() == 1);

			if((edgeEmission == sequenceCharacter) || (edgeEmission == "_"))
			{
				compatibleEdges.insert(e);
			}
		}

		for(std::set<Edge*>::iterator eIt = compatibleEdges.begin(); eIt != compatibleEdges.end(); eIt++)
		{
			Edge* compatibleEdge = *eIt;

			int nextLevel = c.currentTip->level + 1;
			std::string edgeEmission = g->CODE.deCode(compatibleEdge->locus_id, compatibleEdge->emission);
			int nextSequencePosition = c.currentSequencePosition + ((edgeEmission == "_") ? 0 : 1);

			std::vector<int> move_coordinates;
			move_coordinates.push_back(nextLevel);
			move_coordinates.push_back(nextSequencePosition);
			move_coordinates.push_back(nodesPerLevel_ordered_rev.at(nextLevel).at(compatibleEdge->To));

			if((nextSequencePosition == (int)(sequenceStop+1)) && (nextLevel <= graph_lastLevel))
			{
				candidatesPaths resultC = c;
				resultC.currentTip = compatibleEdge->To;
				resultC.edges.push_back(compatibleEdge);
				resultC.aligned_sequence += edgeEmission;
				resultC.edge_coordinates.push_back(move_coordinates);

				matches.push_back(resultC);
				couldExtend = true;
			}
			else
			{
				if(nextLevel < graph_lastLevel)
				{
					if(nextSequencePosition < (int)sequence.length())
					{
						candidatesPaths newC = c;
						newC.currentSequencePosition = nextSequencePosition;
						newC.currentTip = compatibleEdge->To;
						newC.edges.push_back(compatibleEdge);
						newC.aligned_sequence += edgeEmission;
						newC.edge_coordinates.push_back(move_coordinates);

						candidates.push_back(newC);
						couldExtend = true;
					}
				}
				else
				{
					//
				}
			}
		}

		if(! couldExtend)
		{
			if(((c.currentSequencePosition - sequenceStart) >= kMerSize) && (fallbackMatches.size() == 0))
			{
				fallbackMatches.push_back(c);
			}
		}
	}

	if(matches.size() == 0)
	{
		if(fallbackMatches.size() == 0)
		{
			std::cerr << "No fallback matches!?!\n";
			std::cerr << "\t" << "sequenceStart: " << sequenceStart << "\n";
			std::cerr << "\t" << "sequenceStop: " << sequenceStop << "\n";
			std::cerr << "\t" << "graph_firstLevel: " << graph_firstLevel << "\n";
			std::cerr << "\t" << "graph_lastLevel: " << graph_lastLevel << "\n";
			std::cerr << "\t" << "minimumInexactMatchLength: " << minimumInexactMatchLength << "\n";
			std::cerr << "\t" << "sequence: " << sequence << "\n";
			std::cerr << std::flush;
		}
		assert(fallbackMatches.size() > 0);
		matches.push_back(fallbackMatches.at(0));
	}

	assert(matches.size() >= 1);

	std::set<std::string> capturedNodes;
	std::vector<NWPath*> forReturn;

	for(unsigned int cI = 0; cI < matches.size(); cI++)
	{
		std::string nodesString = Utilities::PtoStr(matches.at(cI).edges.at(0)->From)+"__"+Utilities::PtoStr(matches.at(cI).edges.at(matches.at(cI).edges.size() - 1)->To);
		if(capturedNodes.count(nodesString) == 0)
		{
			NWPath* path = new NWPath();
			for(unsigned int nodeI = 0; nodeI <= (matches.at(cI).edge_coordinates.size() - 2); nodeI++)
			{
				std::vector<int> c1 = matches.at(cI).edge_coordinates.at(nodeI);
				std::vector<int> c2 = matches.at(cI).edge_coordinates.at(nodeI+1);

				int firstOrLast = 0;
				if(nodeI == 0)
				{
					firstOrLast = -1;
				}
				if(nodeI == (matches.at(cI).edge_coordinates.size() - 2))
				{
					firstOrLast = 1;
				}

				Edge* graphEdge = matches.at(cI).edges.at(nodeI);
				assert(g->Edges.count(graphEdge) > 0);
				path->createAndAddEdge(c1.at(0), c1.at(1), c1.at(2), c2.at(0), c2.at(1), c2.at(2), graphEdge, firstOrLast, firstOrLast);
			}
			forReturn.push_back(path);

			capturedNodes.insert(nodesString);
		}
	}

	return forReturn;
}



std::vector<std::pair<std::vector<Edge*>, std::string>> GraphAligner::findEdgePathForSequence(std::string sequence, int graph_firstLevel, int graph_lastLevel)
{
	// This function might not perform well for very long chains with very complicated graphs - essentially, in the current implementation,
	// I search for edge paths first and decide then which ones are redundant (same start- and end-node). This could be improved
	// by some dynamic programming approach.

	assert(graph_firstLevel < graph_lastLevel);

	class candidatesPaths {
	public:
		Node* currentTip;
		int currentSequencePosition;
		std::vector<Edge*> edges;
		std::string aligned_sequence;
	};

	std::vector<candidatesPaths> candidates;
	std::vector<candidatesPaths> matches;

	std::set<Node*> startingNodes = g->NodesPerLevel.at(graph_firstLevel);
	for(std::set<Node*>::iterator startNodeIt = startingNodes.begin(); startNodeIt != startingNodes.end(); startNodeIt++)
	{
		candidatesPaths p;
		p.currentTip = *startNodeIt;
		p.currentSequencePosition = 0;
		candidates.push_back(p);
	}

	while(candidates.size() > 0)
	{
		// std::cout << candidates.size() << "\n" << std::flush;

		candidatesPaths c = candidates.at(candidates.size() - 1);
		unsigned int s_before = candidates.size();
		candidates.erase(candidates.end() - 1);
		assert((int)candidates.size() == (s_before - 1));

		assert(((int)c.currentTip->level >= graph_firstLevel) && ((int)c.currentTip->level < graph_lastLevel));
		assert(c.currentSequencePosition < (int)sequence.length());
		std::string sequenceCharacter = sequence.substr(c.currentSequencePosition, 1);

		std::set<Edge*> edges = c.currentTip->Outgoing_Edges;
		std::set<Edge*> compatibleEdges;

 		for(std::set<Edge*>::iterator eIt = edges.begin(); eIt != edges.end(); eIt++)
		{
			Edge* e = *eIt;
			std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
			assert(edgeEmission.length() == 1);

			if((edgeEmission == sequenceCharacter) || (edgeEmission == "_"))
			{
				compatibleEdges.insert(e);
			}
		}

		for(std::set<Edge*>::iterator eIt = compatibleEdges.begin(); eIt != compatibleEdges.end(); eIt++)
		{
			Edge* compatibleEdge = *eIt;

			int nextLevel = c.currentTip->level + 1;
			std::string edgeEmission = g->CODE.deCode(compatibleEdge->locus_id, compatibleEdge->emission);
			int nextSequencePosition = c.currentSequencePosition + ((edgeEmission == "_") ? 0 : 1);

			if((nextSequencePosition == (int)sequence.length()) && (nextLevel <= graph_lastLevel))
			{
				candidatesPaths resultC = c;
				resultC.currentTip = compatibleEdge->To;
				resultC.edges.push_back(compatibleEdge);
				resultC.aligned_sequence += edgeEmission;

				matches.push_back(resultC);
			}
			else
			{
				if(nextLevel <= graph_lastLevel)
				{
					if(nextSequencePosition < (int)sequence.length())
					{
						candidatesPaths newC = c;
						newC.currentSequencePosition = nextSequencePosition;
						newC.currentTip = compatibleEdge->To;
						newC.edges.push_back(compatibleEdge);
						newC.aligned_sequence += edgeEmission;

						candidates.push_back(newC);
					}
				}
				else
				{
					//
				}
			}
		}
	}

	assert(matches.size() >= 1);

	std::set<std::string> capturedNodes;
	std::vector<std::pair<std::vector<Edge*>, std::string>> forReturn;
	for(unsigned int cI = 0; cI < matches.size(); cI++)
	{
		std::string nodesString = Utilities::PtoStr(matches.at(cI).edges.at(0)->From)+"__"+Utilities::PtoStr(matches.at(cI).edges.at(matches.at(cI).edges.size() - 1)->To);
		if(capturedNodes.count(nodesString) == 0)
		{
			std::pair<std::vector<Edge*>, std::string> returnElement;
			returnElement.first = matches.at(cI).edges;
			returnElement.second = matches.at(cI).aligned_sequence;
			forReturn.push_back(returnElement);

			capturedNodes.insert(nodesString);
		}
	}

	return forReturn;
}


int GraphAligner::countMatchesInSequence(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence)
{
	assert(reconstructed_graph.length() == reconstructed_sequence.length());
	assert(reconstructed_graph_levels.size() == reconstructed_graph.length());

	int matches = 0;

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::string G = reconstructed_graph.substr(i, 1);

		if(S != "_")
		{
			if(S == G)
			{
				matches++;
			}
		}
	}

	return matches;
}

std::vector<std::pair<int, Edge*> > GraphAligner::_graph_get_previous_z_values_and_edges(int x, int z)
{
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x > 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> incomingEdges = thisZ->Incoming_Edges;
	set<Node*> nodesPreviousLevel;
	assert(incomingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* fromNode = e->From;
		int z_for_fromNode = nodesPerLevel_ordered_rev.at(x-1).at(fromNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_fromNode, e));
	}

	return forReturn;
}

std::vector<std::pair<int, Edge*> > GraphAligner::_graph_get_next_z_values_and_edges(int x, int z)
{
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x >= 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> outgoingEdges = thisZ->Outgoing_Edges;
	set<Node*> nodesNextLevel;
	assert(outgoingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = outgoingEdges.begin(); eIt != outgoingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* toNode = e->To;
		int z_for_toNode = nodesPerLevel_ordered_rev.at(x+1).at(toNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_toNode, e));
	}

	return forReturn;
}

std::vector<int> GraphAligner::_graph_get_previous_z_values(int x, int z)
{
	std::vector<std::pair<int, Edge*> > previousZs = _graph_get_previous_z_values_and_edges(x, z);

	std::set<int> forReturn;
	for(unsigned int i = 0; i < previousZs.size(); i++)
	{
		int previousZ = previousZs.at(i).first;
		forReturn.insert(previousZ);
	}

	return std::vector<int>(forReturn.begin(), forReturn.end());
}

bool alignedReadPair_strandsValid(std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& p)
{
	if ((p.first.alignment_firstLevel() != -1) && (p.second.alignment_firstLevel() != -1) && (p.first.reverse != p.second.reverse))
	{
		if(! p.first.reverse)
		{
			if(p.first.alignment_firstLevel() < p.second.alignment_firstLevel())
			{
				return true;
			}
			else
			{
				return false;
			}	
		}
		else
		{
			if(p.first.alignment_lastLevel() > p.second.alignment_lastLevel())
			{
				return true;
			}
			else
			{
				return false;
			}			
		}
	}
	else
	{
		return false;
	}
}

int alignedReadPair_pairsDistanceInGraphLevels(std::pair<seedAndExtend_return_local, seedAndExtend_return_local>& p)
{
	// std::cerr << "alignedReadPair_pairsDistanceInGraphLevels(..)\n";
	// std::cerr << "\t" << "first reverse: " << p.first.reverse << "\n";
	// std::cerr << "\t" << "first coordinates: " << p.first.alignment_firstLevel() << " - " << p.first.alignment_lastLevel() << "\n";
	// std::cerr << "\t" << "second reverse: " << p.second.reverse << "\n";
	// std::cerr << "\t" << "second coordinates: " << p.second.alignment_firstLevel() << " - " << p.second.alignment_lastLevel() << "\n";
	
	if(p.first.alignment_firstLevel() < p.second.alignment_firstLevel())
	{
		int D = (p.second.alignment_firstLevel() - p.first.alignment_lastLevel());
		// std::cerr << "\t" << "first in front -- distance " << D << "\n";

		// if(alignedReadPair_strandsValid(p) && (! p.first.reverse))
		// {
			// std::cerr << "\t" << "OK" << "\n";
		// }
		// else
		// {
			// std::cerr << "\t" << "WARNING!" << "\n";
		// }
		// std::cerr << std::flush;
		return D;

	}
	else
	{
		assert(p.first.alignment_firstLevel() >= p.second.alignment_firstLevel());
		int D = (p.first.alignment_firstLevel() - p.second.alignment_lastLevel());
		// std::cerr << "\t" << "second in front -- distance " << D << "\n";

		// if(alignedRead Pair_strandsValid(p) && (p.first.reverse))
		// {
			// std::cerr << "\t" << "OK" << "\n";
		// }
		// else
		// {
			// std::cerr << "\t" << "WARNING!" << "\n";
		// }
		// std::cerr << std::flush;
		return D;
	}
}

