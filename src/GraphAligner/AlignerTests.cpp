/*
 * AlignerTests.cpp
 *
 *  Created on: 27.06.2013
 *      Author: AlexanderDilthey
 */

#include "AlignerTests.h"

#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <exception>
#include <stdexcept>

#include "../Graph/Graph.h"
#include "../Utilities.h"
#include "../NextGen/Validation.h"

#include "GraphAligner.h"
#include "GraphAlignerAffine.h"
#include "GraphAlignernonAffine.h"
#include "GraphAlignerendsFree.h"

int lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret = 0;
int lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret = 0;
int lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret = 0;
int lastCall_sampleString_forSeedAndExtend_trunk_end_sequence_ret = 0;

double minusInfinity = -1 * numeric_limits<double>::max();


void sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(Graph* g, int kMerSize, int S_match, int S_mismatch, int S_gapOpen, int S_gapExtend, std::string& sequenceLabels_ret, std::string& underlyingEdgeLabels_ret, std::vector<int>& underlyingEdges_levels_ret, int& string_begin_ret, int& string_end_ret, int minTrunkLength = 10, int maxTrunkLength = 15, int minTrunks = 3, int maxTrunks = 5, int minSideTrunkExtensions = 1, int maxSideTrunkExtensions = 3, int minIntermediatePenalty = 3, int maxIntermediatePenalty = 7, int betweenTrunkGapMin = 40, int betweenTrunkGapMax = 80, bool returnWholeGenomeString = true)
{
	bool verbose = false;
	if(verbose)
		std::cout << "\n-------------\nEnter sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(..).\n" << std::flush;

	int levels = g->NodesPerLevel.size();
	assert(levels > 150);
	assert(minTrunkLength >= kMerSize);
	assert(maxTrunkLength >= minTrunkLength);

	int trunkExtensions_size = kMerSize - 1;
	assert((trunkExtensions_size * S_match) > maxIntermediatePenalty);

	int trunks = minTrunks + Utilities::randomNumber(maxTrunks - minTrunks);
	int estimatedMaxLength = trunks * (maxTrunkLength + maxSideTrunkExtensions * (trunkExtensions_size + maxIntermediatePenalty) * 1.5);
	std::vector<int> bigGapLengths;
	for(int trunkI = 1; trunkI < trunks; trunkI++)
	{
		bigGapLengths.push_back(betweenTrunkGapMin + Utilities::randomNumber(betweenTrunkGapMax - betweenTrunkGapMin));
		estimatedMaxLength += bigGapLengths.at(trunkI - 1);
	}

	bool error = false;
	auto sampleCharactersFromNode = [&](Node* n, int characters, std::vector<Edge*>& traversedEdges, std::string& edgeLabels) {
		traversedEdges.clear();
		edgeLabels.clear();

		assert(n != 0);
		Node* currentNode = n;
		int have_characters = 0;
		while(((have_characters != characters)) && (! error))
		{
			assert(currentNode != 0);

			if(!( currentNode->level <= (levels - 2)))
			{
				error = true;
				break;
			}

			std::set<Edge*> availableEdges = currentNode->Outgoing_Edges;
			assert(availableEdges.size() > 0);
			std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

			int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
			Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);

			assert(selectedEdge != 0);
			traversedEdges.push_back(selectedEdge);
			string emission = g->CODE.deCode(selectedEdge->locus_id, selectedEdge->emission);
			assert(emission.length() == 1);

			edgeLabels.append(emission);
			currentNode = selectedEdge->To;

			if(emission != "_")
			{
				have_characters++;
			}
		}
	};


	auto selectEdge = [&](Node* currentNode, int direction) -> Edge* {
		assert(currentNode != 0);

		assert((direction == -1) || (direction == 1));
		std::set<Edge*> availableEdges;
		if(direction == 1)
		{
			availableEdges = currentNode->Outgoing_Edges;
		}
		else
		{
			availableEdges = currentNode->Incoming_Edges;
		}

		std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

		if(availableEdges.size() == 0)
		{
			return 0;
		}

		int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
		Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);
		return selectedEdge;
	};

	std::vector<Edge*> total_traversedEdges;

	// std::string total_edgeLabels;
	std::string total_sequenceLabels;

	if(verbose)
		std::cout << "Sample genome with " << trunks << " trunks.\n" << std::flush;

	if((estimatedMaxLength + 20) >= levels)
	{
		error = true;
	}
	else
	{

		assert((levels - 20 - estimatedMaxLength) > 0);
		int firstCharacterStart = 10 + Utilities::randomNumber(levels - 20 - estimatedMaxLength);

		std::set<Node*> nodesFirstLevel = g->NodesPerLevel.at(firstCharacterStart);
		std::vector<Node*> nodesFirstLevel_vec(nodesFirstLevel.begin(), nodesFirstLevel.end());
		Node* currentNode = nodesFirstLevel_vec.at(Utilities::randomNumber(nodesFirstLevel_vec.size() - 1));

		for(int trunkI = 0; trunkI < trunks; trunkI++)
		{
			int want_extensions =  minSideTrunkExtensions + Utilities::randomNumber(maxSideTrunkExtensions - minSideTrunkExtensions);
			int want_trunk_characters = minTrunkLength + Utilities::randomNumber(maxTrunkLength - minTrunkLength);
			if(verbose)
			{
				std::cout << "Sample trunk of length " << want_trunk_characters << ", with " << want_extensions << " either side. Begin level in graph of first extension: " << currentNode->level << "\n" << std::flush;
			}

			for(unsigned int aroundExtensionI = 0; aroundExtensionI <= 1; aroundExtensionI++)
			{
				for(int eI = 0; eI < want_extensions; eI++)
				{
					if(verbose)
						std::cout << "\t extension" << eI << " on side " << aroundExtensionI << ". Begin in graph: " << currentNode->level  << " / " << levels << "\n" << std::flush;

					// constant bit, before trunk

					if((!error) && (aroundExtensionI == 0))
					{
						std::vector<Edge*> traversed_edges;
						std::string traversedEdgeLabels;

						sampleCharactersFromNode(currentNode, trunkExtensions_size, traversed_edges, traversedEdgeLabels);
						currentNode = traversed_edges.back()->To;
						// total_edgeLabels.append(traversedEdgeLabels);
						total_traversedEdges.insert(total_traversedEdges.end(), traversed_edges.begin(), traversed_edges.end());
						total_sequenceLabels.append(traversedEdgeLabels);
					}

					if(error)
						break;

					// variable bit
					int proposed_extension_penalty = minIntermediatePenalty + Utilities::randomNumber(maxIntermediatePenalty - minIntermediatePenalty);

					bool extension_is_gap = false;
					bool gapExtension_gapsGraph = false;
					if(proposed_extension_penalty >= -1 * (S_gapOpen + S_gapExtend) )
					{
						if(Utilities::randomNumber(1) == 1)
						{
							extension_is_gap = true;
							if(Utilities::randomNumber(1) == 1)
							{
								gapExtension_gapsGraph = true;
							}
						}
					}

					bool in_affine_gap = false;
					int accumulated_penalty = 0;
					do {

						if(gapExtension_gapsGraph)
						{
							assert(extension_is_gap);
							in_affine_gap = false;
							int number_of_gaps = (proposed_extension_penalty - (-1) * S_gapOpen) / ( (-1) * S_gapExtend );
							assert(number_of_gaps > 0);

							for(int i = 0; i < number_of_gaps; i++)
							{
								total_sequenceLabels.push_back(Utilities::randomNucleotide());
								total_traversedEdges.push_back(0);
							}

							int realized_penalty = S_gapOpen + S_gapExtend * number_of_gaps;
							accumulated_penalty += ( (-1) * realized_penalty);
							break;
						}
						else
						{
							Edge* nextEdge = selectEdge(currentNode, 1);
							if(nextEdge == 0)
							{
								error = true;
								break;
							}
							else
							{
								string edgeEmission = g->CODE.deCode(nextEdge->locus_id, nextEdge->emission);
								assert(edgeEmission.length() == 1);

								if(extension_is_gap)
								{
									if(edgeEmission == "_")
									{
										in_affine_gap = false;
										total_sequenceLabels.push_back('_');
										total_traversedEdges.push_back(nextEdge);
										currentNode = nextEdge->To;
									}
									else
									{
										int thisStep_penalty;

										if(in_affine_gap)
										{
											thisStep_penalty = S_gapExtend;
										}
										else
										{
											thisStep_penalty = S_gapOpen + S_gapExtend;
											in_affine_gap = true;
										}

										int penalty_with_thisStep_added = accumulated_penalty + (-1)*thisStep_penalty;

										if(penalty_with_thisStep_added <= proposed_extension_penalty)
										{

											total_sequenceLabels.push_back('_');
											total_traversedEdges.push_back(nextEdge);
											accumulated_penalty = penalty_with_thisStep_added;
											currentNode = nextEdge->To;
										}
										else
										{
											break;
										}
									}
								}
								else
								{
									int thisStep_penalty = S_mismatch;
									int penalty_with_thisStep_added = accumulated_penalty + (-1)*thisStep_penalty;
									if(penalty_with_thisStep_added <= proposed_extension_penalty)
									{
										total_sequenceLabels.push_back(Utilities::randomNucleotide());
										total_traversedEdges.push_back(nextEdge);
										accumulated_penalty= penalty_with_thisStep_added;
										currentNode = nextEdge->To;
									}
									else
									{
										break;
									}
								}
							}
						}
					} while((! error) && (accumulated_penalty < proposed_extension_penalty));

					// constant bit, after trunk
					if((! error) && (aroundExtensionI == 1))
					{
						std::vector<Edge*> traversed_edges;
						std::string traversedEdgeLabels;

						sampleCharactersFromNode(currentNode, trunkExtensions_size, traversed_edges, traversedEdgeLabels);
						currentNode = traversed_edges.back()->To;
						// total_edgeLabels.append(traversedEdgeLabels);
						total_traversedEdges.insert(total_traversedEdges.end(), traversed_edges.begin(), traversed_edges.end());
						total_sequenceLabels.append(traversedEdgeLabels);
					}

				}

				if((! error) && (aroundExtensionI == 0))
				{

					if(verbose)
						std::cout << "\t trunk trunk. Begin in graph: " << currentNode->level << "\n" << std::flush;

					// sample trunk
					std::vector<Edge*> traversed_edges;
					std::string traversedEdgeLabels;

					sampleCharactersFromNode(currentNode, want_trunk_characters, traversed_edges, traversedEdgeLabels);
					currentNode = traversed_edges.back()->To;
					// total_edgeLabels.append(traversedEdgeLabels);
					total_traversedEdges.insert(total_traversedEdges.end(), traversed_edges.begin(), traversed_edges.end());
					total_sequenceLabels.append(traversedEdgeLabels);

				}
			}

			if((!error) && (trunkI != (trunks - 1)))
			{
				// make nice gap!
				int wantGapLength = (int)bigGapLengths.at(trunkI);

				if(verbose)
					std::cout << "Introduce gap of " << wantGapLength << " characters.\n" << std::flush;

				std::vector<Edge*> traversed_edges;
				std::string traversedEdgeLabels;

				sampleCharactersFromNode(currentNode, wantGapLength, traversed_edges, traversedEdgeLabels);

				currentNode = traversed_edges.back()->To;
				// total_edgeLabels.append(traversedEdgeLabels);
				total_traversedEdges.insert(total_traversedEdges.end(), traversed_edges.begin(), traversed_edges.end());

				std::string sequenceLabels_gaps;
				sequenceLabels_gaps.resize(traversedEdgeLabels.size(), '_');
				total_sequenceLabels.append(sequenceLabels_gaps);
			}
		}
	}

	if(!error)
	{
		assert(total_sequenceLabels.size() == total_traversedEdges.size());
		assert(total_sequenceLabels.size() > 0);
		int spannedLevels_begin = total_traversedEdges.at(0)->From->level;
		int spannedLevels_end = total_traversedEdges.back()->To->level;
		int spannedLevels = spannedLevels_end - spannedLevels_begin;
		assert(spannedLevels > 0);

		std::string total_realUnderlyingEdgeLabels;
		std::vector<int> total_underlyingGraphLevels;
		for(unsigned int eI = 0; eI < total_traversedEdges.size(); eI++)
		{
			Edge* e = total_traversedEdges.at(eI);
			if(e == 0)
			{
				total_realUnderlyingEdgeLabels.append("_");
				total_underlyingGraphLevels.push_back(-1);
			}
			else
			{
				string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
				assert(edgeEmission.length() == 1);
				total_realUnderlyingEdgeLabels.append(edgeEmission);

				int level = e->From->level;
				if(eI > 0)
				{
					assert(total_underlyingGraphLevels.size() > 0);
					assert((total_underlyingGraphLevels.at(total_underlyingGraphLevels.size() - 1) == -1) || (level == (total_underlyingGraphLevels.at(total_underlyingGraphLevels.size() - 1) + 1)));
				}

				total_underlyingGraphLevels.push_back(level);
			}
		}

		assert(total_underlyingGraphLevels.size() == total_sequenceLabels.length());
		assert(total_underlyingGraphLevels.size() == total_realUnderlyingEdgeLabels.length());

		for(int staticCapI = 1; staticCapI <= trunkExtensions_size; staticCapI++)
		{
			std::string char_l_from_emission = total_sequenceLabels.substr(staticCapI - 1, 1);
			std::string char_l_from_graph = total_realUnderlyingEdgeLabels.substr(staticCapI - 1, 1);
			assert((char_l_from_graph == "_") || (char_l_from_emission == char_l_from_graph));

			std::string char_r_from_emission = total_sequenceLabels.substr(total_sequenceLabels.size() - staticCapI, 1);
			std::string char_r_from_graph = total_realUnderlyingEdgeLabels.substr(total_realUnderlyingEdgeLabels.size() - staticCapI, 1);
			assert((char_r_from_graph == "_") || (char_r_from_emission == char_r_from_graph));
		}

		if(returnWholeGenomeString)
		{
			Node* targetNode_left = total_traversedEdges.at(0)->From;

			std::string leftMissing_emission;
			std::vector<int> leftMissing_underlyingGraphLevels;
			std::string leftMissing_underlyingEdgeLabels;

			Node* targetNode_right = total_traversedEdges.back()->To;

			std::string rightMissing_emission;
			std::vector<int> rightMissing_underlyingGraphLevels;
			std::string rightMissing_underlyingEdgeLabels;


			class WGfindExtensionPointer {
			public:
				double Score;
				std::string emission;
				std::string edgeLabels;
				std::vector<int> underlyingGraphLevels;
				WGfindExtensionPointer(){
					Score = 0;
				}
			};

			std::map<Node*, WGfindExtensionPointer> pathsToNodes_left;
			std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_endInAffineGap;

			std::map<Node*, WGfindExtensionPointer> pathsToNodes_right;
			std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_endInAffineGap;

			std::set<Node*> nodesFirstLevel = g->NodesPerLevel.at(0);
			for(std::set<Node*>::iterator nodeIt = nodesFirstLevel.begin(); nodeIt != nodesFirstLevel.end(); nodeIt++)
			{
				Node* n = *nodeIt;
				pathsToNodes_left[n] = WGfindExtensionPointer();
				pathsToNodes_left_endInAffineGap[n] = WGfindExtensionPointer();
				pathsToNodes_left_endInAffineGap[n].Score = minusInfinity;
			}

			std::set<Node*> nodesLastLevel = g->NodesPerLevel.at(levels - 1);
			for(std::set<Node*>::iterator nodeIt = nodesLastLevel.begin(); nodeIt != nodesLastLevel.end(); nodeIt++)
			{
				Node* n = *nodeIt;
				pathsToNodes_right[n] = WGfindExtensionPointer();
				pathsToNodes_right_endInAffineGap[n] = WGfindExtensionPointer();
				pathsToNodes_right_endInAffineGap[n].Score = minusInfinity;
			}


			if(targetNode_left->level != 0)
			{

				if(verbose)
					std::cout << "Require whole-genome string, but left-hand node only at " << targetNode_left->level << " -- sample missing levels!\n" << std::flush;

				for(unsigned int lI = 0; lI < targetNode_left->level; lI++)
				{
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_nextLevel;
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_endInAffineGap_nextLevel;

					assert(pathsToNodes_left.size() == pathsToNodes_left_endInAffineGap.size());

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_left.begin(); nodeIt != pathsToNodes_left.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
						assert(outgoingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = outgoingEdges.begin(); edgeIt != outgoingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->To;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapOpen + S_gapExtend;
							}
							else
							{
								jumpScore = nodeIt->second.Score + 0;
							}

							if((pathsToNodes_left_nextLevel.count(targetNode) == 0) || (pathsToNodes_left_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level);
								pathsToNodes_left_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_left_endInAffineGap.begin(); nodeIt != pathsToNodes_left_endInAffineGap.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
						assert(outgoingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = outgoingEdges.begin(); edgeIt != outgoingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->To;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapExtend;
							}
							else
							{
								jumpScore = minusInfinity;
							}

							if((pathsToNodes_left_endInAffineGap_nextLevel.count(targetNode) == 0) || (pathsToNodes_left_endInAffineGap_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level);
								pathsToNodes_left_endInAffineGap_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeScoreIt = pathsToNodes_left_endInAffineGap_nextLevel.begin(); nodeScoreIt != pathsToNodes_left_endInAffineGap_nextLevel.end(); nodeScoreIt++)
					{
						Node* n = nodeScoreIt->first;
						assert(pathsToNodes_left_nextLevel.count(n) > 0);
						if(pathsToNodes_left_nextLevel.at(n).Score < pathsToNodes_left_endInAffineGap_nextLevel.at(n).Score)
						{
							pathsToNodes_left_nextLevel.at(n) = pathsToNodes_left_endInAffineGap_nextLevel.at(n);
						}
					}

					pathsToNodes_left = pathsToNodes_left_nextLevel;
					pathsToNodes_left_endInAffineGap = pathsToNodes_left_endInAffineGap_nextLevel;
				}

				assert(pathsToNodes_left.count(targetNode_left) > 0);
				leftMissing_emission = pathsToNodes_left[targetNode_left].emission;
				leftMissing_underlyingGraphLevels = pathsToNodes_left[targetNode_left].underlyingGraphLevels;
				leftMissing_underlyingEdgeLabels = pathsToNodes_left[targetNode_left].edgeLabels;

				assert(leftMissing_emission.size() == leftMissing_underlyingEdgeLabels.size());
				assert(leftMissing_emission.size() == leftMissing_underlyingGraphLevels.size());
			}

			if(targetNode_right->level != (levels - 1))
			{

				if(verbose)
					std::cout << "Require whole-genome string, but right-hand node only at " << targetNode_right->level << " -- sample missing levels!\n" << std::flush;

				for(unsigned int lI = (levels - 1); lI > targetNode_right->level; lI--)
				{
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_nextLevel;
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_endInAffineGap_nextLevel;

					assert(pathsToNodes_right.size() == pathsToNodes_right_endInAffineGap.size());

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_right.begin(); nodeIt != pathsToNodes_right.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> incomingEdges = n->Incoming_Edges;
						assert(incomingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = incomingEdges.begin(); edgeIt != incomingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->From;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapOpen + S_gapExtend;
							}
							else
							{
								jumpScore = nodeIt->second.Score + 0;
							}

							if((pathsToNodes_right_nextLevel.count(targetNode) == 0) || (pathsToNodes_right_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level - 1);
								pathsToNodes_right_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_right_endInAffineGap.begin(); nodeIt != pathsToNodes_right_endInAffineGap.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> incomingEdges = n->Incoming_Edges;
						assert(incomingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = incomingEdges.begin(); edgeIt != incomingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->From;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapExtend;
							}
							else
							{
								jumpScore = minusInfinity;
							}

							if((pathsToNodes_right_endInAffineGap_nextLevel.count(targetNode) == 0) || (pathsToNodes_right_endInAffineGap_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level - 1);
								pathsToNodes_right_endInAffineGap_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeScoreIt = pathsToNodes_right_endInAffineGap_nextLevel.begin(); nodeScoreIt != pathsToNodes_right_endInAffineGap_nextLevel.end(); nodeScoreIt++)
					{
						Node* n = nodeScoreIt->first;
						assert(pathsToNodes_right_nextLevel.count(n) > 0);
						if(pathsToNodes_right_nextLevel.at(n).Score < pathsToNodes_right_endInAffineGap_nextLevel.at(n).Score)
						{
							pathsToNodes_right_nextLevel.at(n) = pathsToNodes_right_endInAffineGap_nextLevel.at(n);
						}
					}

					pathsToNodes_right = pathsToNodes_right_nextLevel;
					pathsToNodes_right_endInAffineGap = pathsToNodes_right_endInAffineGap_nextLevel;
				}

				assert(pathsToNodes_right.count(targetNode_right) > 0);
				rightMissing_emission = pathsToNodes_right[targetNode_right].emission;
				rightMissing_underlyingGraphLevels = pathsToNodes_right[targetNode_right].underlyingGraphLevels;
				rightMissing_underlyingEdgeLabels = pathsToNodes_right[targetNode_right].edgeLabels;
				std::reverse(rightMissing_emission.begin(), rightMissing_emission.end());
				std::reverse(rightMissing_underlyingGraphLevels.begin(), rightMissing_underlyingGraphLevels.end());
				std::reverse(rightMissing_underlyingEdgeLabels.begin(), rightMissing_underlyingEdgeLabels.end());

				assert(rightMissing_emission.size() == rightMissing_underlyingEdgeLabels.size());
				assert(rightMissing_emission.size() == rightMissing_underlyingGraphLevels.size());
			}

			std::string concatenated_emission = leftMissing_emission + total_sequenceLabels + rightMissing_emission ;
			std::vector<int> concatenated_underlyingGraphLevels;
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), leftMissing_underlyingGraphLevels.begin(), leftMissing_underlyingGraphLevels.end());
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), total_underlyingGraphLevels.begin(), total_underlyingGraphLevels.end());
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), rightMissing_underlyingGraphLevels.begin(), rightMissing_underlyingGraphLevels.end());
			std::string concatenated_underlyingEdgeLabels = leftMissing_underlyingEdgeLabels + total_realUnderlyingEdgeLabels + rightMissing_underlyingEdgeLabels;

			assert(concatenated_emission.length() == concatenated_underlyingEdgeLabels.length());
			assert(concatenated_emission.length() == concatenated_underlyingGraphLevels.size());

			sequenceLabels_ret = concatenated_emission;
			string_begin_ret = 0;
			string_end_ret = levels - 1;
			underlyingEdges_levels_ret = concatenated_underlyingGraphLevels;
			underlyingEdgeLabels_ret = concatenated_underlyingEdgeLabels;
		}
		else
		{
			sequenceLabels_ret = total_sequenceLabels;
			string_begin_ret = spannedLevels_begin;
			string_end_ret = spannedLevels_end;
			underlyingEdges_levels_ret = total_underlyingGraphLevels;
			underlyingEdgeLabels_ret = total_realUnderlyingEdgeLabels;
		}
	}
	else
	{
		std::cerr << "recursion...\n" << std::flush;

		return sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
				g,
				kMerSize,
				S_match,
				S_mismatch,
				S_gapOpen,
				S_gapExtend,
				sequenceLabels_ret,
				underlyingEdgeLabels_ret,
				underlyingEdges_levels_ret,
				string_begin_ret,
				string_end_ret,
				minTrunkLength,
				maxTrunkLength,
				minTrunks,
				maxTrunks,
				minSideTrunkExtensions,
				maxSideTrunkExtensions,
				minIntermediatePenalty,
				maxIntermediatePenalty,
				betweenTrunkGapMin,
				betweenTrunkGapMax,
				returnWholeGenomeString
		);
	}

}


void sampleStringFromGraph_forSeedAndExtend(Graph* g, int kMerSize, int S_match, int S_mismatch, int S_gapOpen, int S_gapExtend, std::string& edgeLabels_ret, std::string& underlyingEdges_ret, std::vector<int>& underlyingEdges_levels_ret, int& string_begin_ret, int& string_end_ret, int minTrunkLength, int maxTrunkLength, int minSideTrunkExtensions, int maxSideTrunkExtensions, int minIntermediatePenalty, int maxIntermediatePenalty, int minExtensionStaticRegionLength, int maxExtensionStaticRegionLength, bool guaranteeSmoothExtensibility, bool returnWholeGenomeString)
{
	int levels = g->NodesPerLevel.size();
	assert(levels > 150);
	assert(minTrunkLength >= kMerSize);
	assert(maxTrunkLength >= minTrunkLength);

	if(guaranteeSmoothExtensibility)
		assert((minExtensionStaticRegionLength * S_match) > maxIntermediatePenalty);

	bool error = false;
	bool verbose = false;

	std::vector<Edge*> traversedEgdes_trunk;
	std::string edgeLabels_trunk;

	assert((levels - 30 - maxTrunkLength) > 0);
	int trunk_begin = 30 + Utilities::randomNumber(levels - 60 - maxTrunkLength);
	int want_trunk_characters = minTrunkLength + Utilities::randomNumber(maxTrunkLength - minTrunkLength);
	int have_trunk_characters = 0;
	if(verbose)
	{
		std::cout << "Sample trunk of length " << want_trunk_characters << ", beginning at level " << trunk_begin << " of " << levels << "\n" << std::flush;
	}


	std::set<Node*> nodesFirstLevel = g->NodesPerLevel.at(trunk_begin);
	std::vector<Node*> nodesFirstLevel_vec(nodesFirstLevel.begin(), nodesFirstLevel.end());
	Node* currentNode = nodesFirstLevel_vec.at(Utilities::randomNumber(nodesFirstLevel_vec.size() - 1));

	while(((have_trunk_characters != want_trunk_characters)) && (! error))
	{
		if(!( currentNode->level <= (levels - 2)))
		{
			error = true;
			break;
		}

		assert(((int)currentNode->level >= trunk_begin) && ((int)currentNode->level <= (levels - 2)));

		std::set<Edge*> availableEdges = currentNode->Outgoing_Edges;
		std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

		int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
		Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);

		traversedEgdes_trunk.push_back(selectedEdge);
		string emission = g->CODE.deCode(selectedEdge->locus_id, selectedEdge->emission);
		assert(emission.length() == 1);

		edgeLabels_trunk.append(emission);
		currentNode = selectedEdge->To;

		if(emission != "_")
		{
			have_trunk_characters++;
		}
	}

	auto selectEdge = [&](Node* currentNode, int direction) -> Edge* {
		assert((direction == -1) || (direction == 1));
		std::set<Edge*> availableEdges;
		if(direction == 1)
		{
			availableEdges = currentNode->Outgoing_Edges;
		}
		else
		{
			availableEdges = currentNode->Incoming_Edges;
		}

		std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

		if(availableEdges.size() == 0)
		{
			return 0;
		}

		int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
		Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);
		return selectedEdge;
	};

	std::vector<Edge*> edges_leftExtension;
	std::string leftExtension;
	std::vector<Edge*> edges_rightExtension;
	std::string rightExtension;

	int have_bothsides_extensions = 0;
	int want_bothsides_extensions = minSideTrunkExtensions + Utilities::randomNumber(maxSideTrunkExtensions - minSideTrunkExtensions);

	if(verbose)
		std::cout << "Now extend with " << want_bothsides_extensions << " both sides.\n" << std::flush;

	assert(traversedEgdes_trunk.size() != 0);
	Node* currentNode_left = traversedEgdes_trunk.at(0)->From;
	Node* currentNode_right = traversedEgdes_trunk.back()->To;
	while((! error) && (have_bothsides_extensions != want_bothsides_extensions))
	{
		// left extension
		if(verbose)
			std::cout << "\tExtension round " << have_bothsides_extensions << "\n" << std::flush;

		int proposed_extension_penalty_left = minIntermediatePenalty + Utilities::randomNumber(maxIntermediatePenalty - minIntermediatePenalty);

		bool extension_left_is_gap = false;
		bool gapExtension_left_gapsGraph = false;
		if(proposed_extension_penalty_left >= -1 * (S_gapOpen + S_gapExtend) )
		{
			if(Utilities::randomNumber(1) == 1)
			{
				extension_left_is_gap = true;
				if(Utilities::randomNumber(1) == 1)
				{
					gapExtension_left_gapsGraph = true;
				}
			}
		}

		if(verbose)
		{
			std::cout << "\t\tLeft" << "\n";
			std::cout << "\t\t\tpenalty: " << proposed_extension_penalty_left << "\n";
			std::cout << "\t\t\tgap: " << extension_left_is_gap << "\n" << std::flush;
			std::cout << "\t\t\tif gap, gap is in graph: " << gapExtension_left_gapsGraph << "\n" << std::flush;
		}


		int accumulated_penalty_left = 0;
		bool left_in_affine_gap = false;
		do {

			if(gapExtension_left_gapsGraph)
			{
				assert(extension_left_is_gap);
				left_in_affine_gap = false;
				int number_of_gaps = (proposed_extension_penalty_left - (-1) * S_gapOpen) / ( (-1) * S_gapExtend );
				assert(number_of_gaps > 0);

				for(int i = 0; i < number_of_gaps; i++)
				{
					leftExtension.push_back(Utilities::randomNucleotide());
					edges_leftExtension.push_back(0);
				}

				int realized_penalty = S_gapOpen + S_gapExtend * number_of_gaps;
				accumulated_penalty_left += ( (-1) * realized_penalty);
				break;
			}
			else
			{
				Edge* nextEdge = selectEdge(currentNode_left, -1);
				if(nextEdge == 0)
				{
					error = true;
					break;
				}
				else
				{
					string edgeEmission = g->CODE.deCode(nextEdge->locus_id, nextEdge->emission);
					assert(edgeEmission.length() == 1);

					if(extension_left_is_gap)
					{
						if(edgeEmission == "_")
						{
							left_in_affine_gap = false;
							leftExtension.push_back('_');
							edges_leftExtension.push_back(nextEdge);
							currentNode_left = nextEdge->From;
						}
						else
						{
							int thisStep_penalty;

							if(left_in_affine_gap)
							{
								thisStep_penalty = S_gapExtend;
							}
							else
							{
								thisStep_penalty = S_gapOpen + S_gapExtend;
								left_in_affine_gap = true;
							}

							int penalty_with_thisStep_added = accumulated_penalty_left + (-1)*thisStep_penalty;

							if(penalty_with_thisStep_added <= proposed_extension_penalty_left)
							{

								leftExtension.push_back('_');
								edges_leftExtension.push_back(nextEdge);
								accumulated_penalty_left = penalty_with_thisStep_added;
								currentNode_left = nextEdge->From;
							}
							else
							{
								break;
							}
						}
					}
					else
					{
						int thisStep_penalty = S_mismatch;
						int penalty_with_thisStep_added = accumulated_penalty_left + (-1)*thisStep_penalty;
						if(penalty_with_thisStep_added <= proposed_extension_penalty_left)
						{
							leftExtension.push_back(Utilities::randomNucleotide());
							edges_leftExtension.push_back(nextEdge);
							accumulated_penalty_left = penalty_with_thisStep_added;
							currentNode_left = nextEdge->From;
						}
						else
						{
							break;
						}
					}
				}
			}
		} while((! error) && (accumulated_penalty_left < proposed_extension_penalty_left));

		assert(accumulated_penalty_left <= proposed_extension_penalty_left);

		int proposedStaticExtensionLeft = minExtensionStaticRegionLength + Utilities::randomNumber(maxExtensionStaticRegionLength - minExtensionStaticRegionLength);
		int collectedStaticLeftExtension = 0;
		while((! error) && (collectedStaticLeftExtension < proposedStaticExtensionLeft))
		{
			Edge* nextEdge = selectEdge(currentNode_left, -1);
			if(nextEdge == 0)
			{
				error = true;
				break;
			}

			string edgeEmission = g->CODE.deCode(nextEdge->locus_id, nextEdge->emission);
			assert(edgeEmission.length() == 1);

			if(edgeEmission != "_")
			{
				collectedStaticLeftExtension++;
			}

			leftExtension.append(edgeEmission);
			edges_leftExtension.push_back(nextEdge);
			currentNode_left = nextEdge->From;
		}

		if(verbose)
			std::cout << "\t\tRight\n" << std::flush;

		int proposed_extension_penalty_right = minIntermediatePenalty + Utilities::randomNumber(maxIntermediatePenalty - minIntermediatePenalty);

		bool extension_right_is_gap = false;
		bool gapExtension_right_gapsGraph = false;
		if(proposed_extension_penalty_right >= -1 * (S_gapOpen + S_gapExtend) )
		{
			if(Utilities::randomNumber(1) == 1)
			{
				extension_right_is_gap = true;
				if(Utilities::randomNumber(1) == 1)
				{
					gapExtension_right_gapsGraph = true;
				}
			}
		}

		if(verbose)
		{
			std::cout << "\t\t\tpenalty: " << proposed_extension_penalty_right << "\n";
			std::cout << "\t\t\tgap: " << extension_right_is_gap << "\n" << std::flush;
		}

		int accumulated_penalty_right = 0;
		bool right_in_affine_gap = false;
		do {
			if(gapExtension_right_gapsGraph)
			{
				assert(extension_right_is_gap);
				right_in_affine_gap = false;
				int number_of_gaps = (proposed_extension_penalty_right - (-1) * S_gapOpen) / ( (-1)*S_gapExtend );
				assert(number_of_gaps > 0);

				for(int i = 0; i < number_of_gaps; i++)
				{
					rightExtension.push_back(Utilities::randomNucleotide());
					edges_rightExtension.push_back(0);
				}

				int realized_penalty = S_gapOpen + S_gapExtend * number_of_gaps;
				accumulated_penalty_right += ( -1 * realized_penalty );
				break;
			}
			else
			{
	//			std::cerr << "\t\t\t\t" <<  accumulated_penalty_right << "\n" << std::flush;
	//			std::cerr << "\t\t\t\t currentNode_right: " <<  currentNode_right->level << "\n" << std::flush;

				Edge* nextEdge = selectEdge(currentNode_right, 1);
				if(nextEdge == 0)
				{
					error = true;
					break;
				}
				else
				{
					string edgeEmission = g->CODE.deCode(nextEdge->locus_id, nextEdge->emission);
					assert(edgeEmission.length() == 1);

					if(extension_right_is_gap)
					{
						if(edgeEmission == "_")
						{
							right_in_affine_gap = false;
							rightExtension.push_back('_');
							edges_rightExtension.push_back(nextEdge);
							currentNode_right = nextEdge->To;
						}
						else
						{
							int thisStep_penalty;

							if(right_in_affine_gap)
							{
								thisStep_penalty = S_gapExtend;
							}
							else
							{
								thisStep_penalty = S_gapOpen + S_gapExtend;
								right_in_affine_gap = true;
							}

							int penalty_with_thisStep_added = accumulated_penalty_right + (-1)*thisStep_penalty;

							if(penalty_with_thisStep_added <= proposed_extension_penalty_right)
							{

								rightExtension.push_back('_');
								edges_rightExtension.push_back(nextEdge);
								accumulated_penalty_right = penalty_with_thisStep_added;
								currentNode_right = nextEdge->To;
							}
							else
							{
								break;
							}
						}
					}
					else
					{
						int thisStep_penalty = S_mismatch;
						int penalty_with_thisStep_added = accumulated_penalty_right + (-1) * thisStep_penalty;
						if(penalty_with_thisStep_added <= proposed_extension_penalty_right)
						{

							rightExtension.push_back(Utilities::randomNucleotide());
							edges_rightExtension.push_back(nextEdge);
							accumulated_penalty_right = penalty_with_thisStep_added;
							currentNode_right = nextEdge->To;
						}
						else
						{
							break;
						}
					}
				}
			}
		} while((! error) && (accumulated_penalty_right < proposed_extension_penalty_right));


		assert(accumulated_penalty_right <= proposed_extension_penalty_right);

		int proposedStaticExtensionRight = minExtensionStaticRegionLength + Utilities::randomNumber(maxExtensionStaticRegionLength - minExtensionStaticRegionLength);

		int collectedStaticRightExtension = 0;
		while((! error) && (collectedStaticRightExtension < proposedStaticExtensionRight))
		{
			Edge* nextEdge = selectEdge(currentNode_right, 1);
			if(nextEdge == 0)
			{
				error = true;
				break;
			}

			string edgeEmission = g->CODE.deCode(nextEdge->locus_id, nextEdge->emission);
			assert(edgeEmission.length() == 1);

			if(edgeEmission != "_")
			{
				collectedStaticRightExtension++;
			}

			rightExtension.append(edgeEmission);
			edges_rightExtension.push_back(nextEdge);
			currentNode_right = nextEdge->To;
		}

		have_bothsides_extensions++;
	}

	if((!error) && (have_trunk_characters == want_trunk_characters) && ((edges_leftExtension.size() > 0) || (edges_rightExtension.size() > 0)))
	{
		assert(! error);

		assert(leftExtension.size() == edges_leftExtension.size());
		std::reverse(leftExtension.begin(), leftExtension.end());
		std::reverse(edges_leftExtension.begin(), edges_leftExtension.end());

		assert(edges_rightExtension.size() == rightExtension.size());

		std::vector<Edge*> combinedEdges = edges_leftExtension;
		combinedEdges.insert(combinedEdges.end(), traversedEgdes_trunk.begin(), traversedEgdes_trunk.end());
		combinedEdges.insert(combinedEdges.end(), edges_rightExtension.begin(), edges_rightExtension.end());

		std::string combinedEmission = leftExtension;
		combinedEmission.append(edgeLabels_trunk);
		combinedEmission.append(rightExtension);

		assert(combinedEdges.size() > 0);
		int spannedLevels_begin = combinedEdges.at(0)->From->level;
		int spannedLevels_end = combinedEdges.back()->To->level;
		int spannedLevels = spannedLevels_end - spannedLevels_begin;
		assert(spannedLevels > 0);

		std::string realUnderlyingEdgeLabels;
		std::vector<int> underlyingGraphLevels;
		for(unsigned int eI = 0; eI < combinedEdges.size(); eI++)
		{
			Edge* e = combinedEdges.at(eI);
			if(e == 0)
			{
				realUnderlyingEdgeLabels.append("_");
				underlyingGraphLevels.push_back(-1);
			}
			else
			{
				string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
				assert(edgeEmission.length() == 1);
				realUnderlyingEdgeLabels.append(edgeEmission);

				int level = e->From->level;
				if(eI > 0)
				{
					assert(underlyingGraphLevels.size() > 0);
					assert((underlyingGraphLevels.at(underlyingGraphLevels.size() - 1) == -1) || (level == (underlyingGraphLevels.at(underlyingGraphLevels.size() - 1) + 1)));
				}

				underlyingGraphLevels.push_back(level);
			}
		}

		assert(realUnderlyingEdgeLabels.length() == combinedEmission.length());

		for(int staticCapI = 1; staticCapI <= minExtensionStaticRegionLength; staticCapI++)
		{
			std::string char_l_from_emission = combinedEmission.substr(staticCapI - 1, 1);
			std::string char_l_from_graph = realUnderlyingEdgeLabels.substr(staticCapI - 1, 1);
			assert((char_l_from_graph == "_") || (char_l_from_emission == char_l_from_graph));

			std::string char_r_from_emission = combinedEmission.substr(combinedEmission.size() - staticCapI, 1);
			std::string char_r_from_graph = realUnderlyingEdgeLabels.substr(realUnderlyingEdgeLabels.size() - staticCapI, 1);
			assert((char_r_from_graph == "_") || (char_r_from_emission == char_r_from_graph));
		}

		if(returnWholeGenomeString)
		{
			Node* targetNode_left = currentNode_left;

			std::string leftMissing_emission;
			std::vector<int> leftMissing_underlyingGraphLevels;
			std::string leftMissing_underlyingEdgeLabels;

			Node* targetNode_right = currentNode_right;

			std::string rightMissing_emission;
			std::vector<int> rightMissing_underlyingGraphLevels;
			std::string rightMissing_underlyingEdgeLabels;


			class WGfindExtensionPointer {
			public:
				double Score;
				std::string emission;
				std::string edgeLabels;
				std::vector<int> underlyingGraphLevels;
				WGfindExtensionPointer(){
					Score = 0;
				}
			};

			std::map<Node*, WGfindExtensionPointer> pathsToNodes_left;
			std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_endInAffineGap;

			std::map<Node*, WGfindExtensionPointer> pathsToNodes_right;
			std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_endInAffineGap;

			std::set<Node*> nodesFirstLevel = g->NodesPerLevel.at(0);
			for(std::set<Node*>::iterator nodeIt = nodesFirstLevel.begin(); nodeIt != nodesFirstLevel.end(); nodeIt++)
			{
				Node* n = *nodeIt;
				pathsToNodes_left[n] = WGfindExtensionPointer();
				pathsToNodes_left_endInAffineGap[n] = WGfindExtensionPointer();
				pathsToNodes_left_endInAffineGap[n].Score = minusInfinity;
			}

			std::set<Node*> nodesLastLevel = g->NodesPerLevel.at(levels - 1);
			for(std::set<Node*>::iterator nodeIt = nodesLastLevel.begin(); nodeIt != nodesLastLevel.end(); nodeIt++)
			{
				Node* n = *nodeIt;
				pathsToNodes_right[n] = WGfindExtensionPointer();
				pathsToNodes_right_endInAffineGap[n] = WGfindExtensionPointer();
				pathsToNodes_right_endInAffineGap[n].Score = minusInfinity;
			}

			if(currentNode_left->level != 0)
			{
				for(unsigned int lI = 0; lI < targetNode_left->level; lI++)
				{
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_nextLevel;
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_left_endInAffineGap_nextLevel;

					assert(pathsToNodes_left.size() == pathsToNodes_left_endInAffineGap.size());

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_left.begin(); nodeIt != pathsToNodes_left.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
						assert(outgoingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = outgoingEdges.begin(); edgeIt != outgoingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->To;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapOpen + S_gapExtend;
							}
							else
							{
								jumpScore = nodeIt->second.Score + 0;
							}

							if((pathsToNodes_left_nextLevel.count(targetNode) == 0) || (pathsToNodes_left_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level);
								pathsToNodes_left_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_left_endInAffineGap.begin(); nodeIt != pathsToNodes_left_endInAffineGap.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
						assert(outgoingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = outgoingEdges.begin(); edgeIt != outgoingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->To;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapExtend;
							}
							else
							{
								jumpScore = minusInfinity;
							}

							if((pathsToNodes_left_endInAffineGap_nextLevel.count(targetNode) == 0) || (pathsToNodes_left_endInAffineGap_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level);
								pathsToNodes_left_endInAffineGap_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeScoreIt = pathsToNodes_left_endInAffineGap_nextLevel.begin(); nodeScoreIt != pathsToNodes_left_endInAffineGap_nextLevel.end(); nodeScoreIt++)
					{
						Node* n = nodeScoreIt->first;
						assert(pathsToNodes_left_nextLevel.count(n) > 0);
						if(pathsToNodes_left_nextLevel.at(n).Score < pathsToNodes_left_endInAffineGap_nextLevel.at(n).Score)
						{
							pathsToNodes_left_nextLevel.at(n) = pathsToNodes_left_endInAffineGap_nextLevel.at(n);
						}
					}

					pathsToNodes_left = pathsToNodes_left_nextLevel;
					pathsToNodes_left_endInAffineGap = pathsToNodes_left_endInAffineGap_nextLevel;
				}

				assert(pathsToNodes_left.count(targetNode_left) > 0);
				leftMissing_emission = pathsToNodes_left[targetNode_left].emission;
				leftMissing_underlyingGraphLevels = pathsToNodes_left[targetNode_left].underlyingGraphLevels;
				leftMissing_underlyingEdgeLabels = pathsToNodes_left[targetNode_left].edgeLabels;

				assert(leftMissing_emission.size() == leftMissing_underlyingEdgeLabels.size());
				assert(leftMissing_emission.size() == leftMissing_underlyingGraphLevels.size());
			}

			if(currentNode_right->level != (levels - 1))
			{
				for(unsigned int lI = (levels - 1); lI > targetNode_right->level; lI--)
				{
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_nextLevel;
					std::map<Node*, WGfindExtensionPointer> pathsToNodes_right_endInAffineGap_nextLevel;

					assert(pathsToNodes_right.size() == pathsToNodes_right_endInAffineGap.size());

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_right.begin(); nodeIt != pathsToNodes_right.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> incomingEdges = n->Incoming_Edges;
						assert(incomingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = incomingEdges.begin(); edgeIt != incomingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->From;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapOpen + S_gapExtend;
							}
							else
							{
								jumpScore = nodeIt->second.Score + 0;
							}

							if((pathsToNodes_right_nextLevel.count(targetNode) == 0) || (pathsToNodes_right_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level - 1);
								pathsToNodes_right_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeIt = pathsToNodes_right_endInAffineGap.begin(); nodeIt != pathsToNodes_right_endInAffineGap.end(); nodeIt++)
					{
						Node* n = nodeIt->first;
						std::set<Edge*> incomingEdges = n->Incoming_Edges;
						assert(incomingEdges.size() > 0);
						for(std::set<Edge*>::iterator edgeIt = incomingEdges.begin(); edgeIt != incomingEdges.end(); edgeIt++)
						{
							Edge* e = *edgeIt;
							string emission = g->CODE.deCode(e->locus_id, e->emission);
							assert(emission.length() == 1);
							Node* targetNode = e->From;
							double jumpScore;
							if(emission != "_")
							{
								jumpScore = nodeIt->second.Score + S_gapExtend;
							}
							else
							{
								jumpScore = minusInfinity;
							}

							if((pathsToNodes_right_endInAffineGap_nextLevel.count(targetNode) == 0) || (pathsToNodes_right_endInAffineGap_nextLevel.at(targetNode).Score < jumpScore))
							{
								WGfindExtensionPointer nextExtension = nodeIt->second;
								nextExtension.Score = jumpScore;
								nextExtension.edgeLabels += emission;
								nextExtension.emission += "_";
								nextExtension.underlyingGraphLevels.push_back(n->level - 1);
								pathsToNodes_right_endInAffineGap_nextLevel[targetNode] = nextExtension;
							}
						}
					}

					for(std::map<Node*, WGfindExtensionPointer>::iterator nodeScoreIt = pathsToNodes_right_endInAffineGap_nextLevel.begin(); nodeScoreIt != pathsToNodes_right_endInAffineGap_nextLevel.end(); nodeScoreIt++)
					{
						Node* n = nodeScoreIt->first;
						assert(pathsToNodes_right_nextLevel.count(n) > 0);
						if(pathsToNodes_right_nextLevel.at(n).Score < pathsToNodes_right_endInAffineGap_nextLevel.at(n).Score)
						{
							pathsToNodes_right_nextLevel.at(n) = pathsToNodes_right_endInAffineGap_nextLevel.at(n);
						}
					}

					pathsToNodes_right = pathsToNodes_right_nextLevel;
					pathsToNodes_right_endInAffineGap = pathsToNodes_right_endInAffineGap_nextLevel;
				}

				assert(pathsToNodes_right.count(targetNode_right) > 0);
				rightMissing_emission = pathsToNodes_right[targetNode_right].emission;
				rightMissing_underlyingGraphLevels = pathsToNodes_right[targetNode_right].underlyingGraphLevels;
				rightMissing_underlyingEdgeLabels = pathsToNodes_right[targetNode_right].edgeLabels;
				std::reverse(rightMissing_emission.begin(), rightMissing_emission.end());
				std::reverse(rightMissing_underlyingGraphLevels.begin(), rightMissing_underlyingGraphLevels.end());
				std::reverse(rightMissing_underlyingEdgeLabels.begin(), rightMissing_underlyingEdgeLabels.end());

				assert(rightMissing_emission.size() == rightMissing_underlyingEdgeLabels.size());
				assert(rightMissing_emission.size() == rightMissing_underlyingGraphLevels.size());
			}

			assert(combinedEmission.size() == realUnderlyingEdgeLabels.size());
			assert(underlyingGraphLevels.size() == combinedEmission.size());

			std::string concatenated_emission = leftMissing_emission + combinedEmission + rightMissing_emission ;
			std::vector<int> concatenated_underlyingGraphLevels;
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), leftMissing_underlyingGraphLevels.begin(), leftMissing_underlyingGraphLevels.end());
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), underlyingGraphLevels.begin(), underlyingGraphLevels.end());
				concatenated_underlyingGraphLevels.insert(concatenated_underlyingGraphLevels.end(), rightMissing_underlyingGraphLevels.begin(), rightMissing_underlyingGraphLevels.end());
			std::string concatenated_underlyingEdgeLabels = leftMissing_underlyingEdgeLabels + realUnderlyingEdgeLabels + rightMissing_underlyingEdgeLabels;

			assert(concatenated_emission.length() == concatenated_underlyingEdgeLabels.length());
			assert(concatenated_emission.length() == concatenated_underlyingGraphLevels.size());

			edgeLabels_ret = concatenated_emission;
			string_begin_ret = 0;
			string_end_ret = levels - 1;
			underlyingEdges_levels_ret = concatenated_underlyingGraphLevels;
			underlyingEdges_ret = concatenated_underlyingEdgeLabels;
		}
		else
		{
			edgeLabels_ret = combinedEmission;
			string_begin_ret = spannedLevels_begin;
			string_end_ret = spannedLevels_end;
			underlyingEdges_levels_ret = underlyingGraphLevels;
			underlyingEdges_ret = realUnderlyingEdgeLabels;
		}

		int firstEdgeIndexTrunk_notGap = 0;
		while(g->CODE.deCode(traversedEgdes_trunk.at(firstEdgeIndexTrunk_notGap)->locus_id, traversedEgdes_trunk.at(firstEdgeIndexTrunk_notGap)->emission) == "_")
		{
			firstEdgeIndexTrunk_notGap++;
		}

		lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret = traversedEgdes_trunk.at(firstEdgeIndexTrunk_notGap)->From->level;
		lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret = traversedEgdes_trunk.back()->To->level;
		lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret = leftExtension.length();
		lastCall_sampleString_forSeedAndExtend_trunk_end_sequence_ret = lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret + edgeLabels_trunk.length() - 1;
	}
	else
	{
		std::cerr << "recursion...\n" << std::flush;

		return sampleStringFromGraph_forSeedAndExtend(
				g,
				kMerSize,
				S_match,
				S_mismatch,
				S_gapOpen,
				S_gapExtend,
				edgeLabels_ret,
				underlyingEdges_ret,
				underlyingEdges_levels_ret,
				string_begin_ret,
				string_end_ret,
				minTrunkLength,
				maxTrunkLength,
				minSideTrunkExtensions,
				maxSideTrunkExtensions,
				minIntermediatePenalty,
				maxIntermediatePenalty,
				minExtensionStaticRegionLength,
				maxExtensionStaticRegionLength,
				guaranteeSmoothExtensibility,
				returnWholeGenomeString
		);
	}

}

void sampleStringFromGraph(Graph* g, std::string& edgeLabels_ret, int& string_begin_ret, int& string_end_ret, int min_string_length, int max_string_length)
{
	int levels = g->NodesPerLevel.size();
	assert(levels > 2);
	assert(max_string_length >= min_string_length);
	assert(min_string_length > 0);

	string_begin_ret = Utilities::randomNumber(levels - 2);

	int max_end_level = string_begin_ret + max_string_length;
	if(max_end_level > (levels - 1))
	{
		max_end_level = levels - 1;
	}
	assert(max_end_level > string_begin_ret);

	int min_end_level = string_begin_ret + min_string_length;
	int span_end = max_end_level - min_end_level;
	if(span_end == 0)
	{
		string_end_ret = max_end_level;
	}
	else if(span_end < 0)
	{
		string_end_ret = max_end_level;
	}
	else
	{
		assert(span_end > 0);
		string_end_ret = min_end_level + Utilities::randomNumber(span_end);
	}

	assert(string_begin_ret >= 0);
	assert(string_begin_ret < levels);
	assert(string_end_ret >= 0);
	assert(string_end_ret < levels);

	assert(string_end_ret > string_begin_ret);
	assert((string_end_ret - string_begin_ret) <= max_string_length);

	std::vector<Edge*> traversedEdges;
	std::string edgeLabels;

	std::set<Node*> nodesFirstLevel = g->NodesPerLevel.at(string_begin_ret);
	std::vector<Node*> nodesFirstLevel_vec(nodesFirstLevel.begin(), nodesFirstLevel.end());
	Node* currentNode = nodesFirstLevel_vec.at(Utilities::randomNumber(nodesFirstLevel_vec.size() - 1));
	while((int)currentNode->level != string_end_ret)
	{
		assert(((int)currentNode->level >= string_begin_ret) && ((int)currentNode->level <= string_end_ret));

		std::set<Edge*> availableEdges = currentNode->Outgoing_Edges;
		std::vector<Edge*> availableEdges_vec(availableEdges.begin(), availableEdges.end());

		int selectedEdge_index = Utilities::randomNumber(availableEdges_vec.size() - 1);
		Edge* selectedEdge = availableEdges_vec.at(selectedEdge_index);

		traversedEdges.push_back(selectedEdge);
		string emission = g->CODE.deCode(selectedEdge->locus_id, selectedEdge->emission);
		assert(emission.length() == 1);

		edgeLabels.append(emission);
		currentNode = selectedEdge->To;
	}

	assert((int)edgeLabels.length() == (string_end_ret - string_begin_ret));

	edgeLabels_ret = edgeLabels;
}

diploidGenomeString generateRandomGenome(int minimumLength, bool beginAndEndOneLevel)
{
	int achievedLength = 0;

	int homozygousTractLength_Min = 10;
	int homozygousTractLength_Max = 40;

	int heterozygousTractLength_Min = 1;
	int heterozygousTractLength_Max = 20;

	diploidGenomeString gS;

	while(achievedLength < minimumLength)
	{
		int typeToGenerate = Utilities::randomNumber(2);

		std::vector<std::string> segment;

		if(typeToGenerate == 0)
		{
			int segmentLength = Utilities::randomNumber(homozygousTractLength_Max - homozygousTractLength_Min) + homozygousTractLength_Min;
			segment.push_back(Utilities::generateRandomSequence(segmentLength));
			achievedLength += segmentLength;
		}
		else if(typeToGenerate == 1)
		{
			int segmentLength = Utilities::randomNumber(heterozygousTractLength_Max - heterozygousTractLength_Min) + heterozygousTractLength_Min;
			segment.push_back(Utilities::generateRandomSequenceWithGaps(segmentLength));
			segment.push_back(Utilities::generateRandomSequenceWithGaps(segmentLength));
			achievedLength += segmentLength;
		}
		else if(typeToGenerate == 2)
		{
			int segmentLength = Utilities::randomNumber(heterozygousTractLength_Max - heterozygousTractLength_Min) + heterozygousTractLength_Min;
			segment.push_back(Utilities::generateRandomSequence(segmentLength));
			segment.push_back(Utilities::repeatString("_", segmentLength));
			achievedLength += segmentLength;
		}

		gS.push_back(segment);

	}

	std::vector<std::string> segment;
	segment.push_back(Utilities::generateRandomSequence(2));
	gS.push_back(segment);

	if(beginAndEndOneLevel)
	{
		std::vector<std::string> segment;
		segment.push_back(Utilities::generateRandomSequence(2));
		gS.insert(gS.begin(), segment);
	}

	return gS;
}

void _printDiploidGenomeString(diploidGenomeString& gS)
{
	std::string firstLine_separated;
	std::string middleLine_separated;
	std::string bottomLine_separated;

	std::string firstLine;
	std::string middleLine;
	std::string bottomLine;

	for(int segmentI = 0; segmentI < (int)gS.size(); segmentI++)
	{
		int segmentLength = gS.at(segmentI).at(0).size();
		if(gS.at(segmentI).size() == 1)
		{
			firstLine_separated.append(Utilities::repeatString(" ", segmentLength+1));
			middleLine_separated.append(gS.at(segmentI).at(0)+" ");
			bottomLine_separated.append(Utilities::repeatString(" ", segmentLength+1));

			firstLine.append(Utilities::repeatString(" ", segmentLength));
			middleLine.append(gS.at(segmentI).at(0));
			bottomLine.append(Utilities::repeatString(" ", segmentLength));

		}
		else
		{
			assert(gS.at(segmentI).size() == 2);
			if(!((int)gS.at(segmentI).at(1).size() == (int)segmentLength))
			{
				std::cout << "!((int)gS.at(segmentI).at(1).size() == (int)segmentLength)\n";
				std::cout << "segmentI: " << segmentI << "\n";
				std::cout << "gS.at(segmentI).at(0): " << gS.at(segmentI).at(0) << "\n";
				std::cout << "gS.at(segmentI).at(1): " << gS.at(segmentI).at(1) << "\n" << std::flush;
			}
			assert((int)gS.at(segmentI).at(1).size() == (int)segmentLength);
			firstLine_separated.append(gS.at(segmentI).at(0)+" ");
			middleLine_separated.append(Utilities::repeatString(" ", segmentLength+1));
			bottomLine_separated.append(gS.at(segmentI).at(1)+" ");

			firstLine.append(gS.at(segmentI).at(0));
			middleLine.append(Utilities::repeatString(" ", segmentLength));
			bottomLine.append(gS.at(segmentI).at(1));
		}
	}

	std::cout << firstLine_separated << "\n" << middleLine_separated << "\n" << bottomLine_separated << "\n\n";
	std::cout << firstLine << "\n" << middleLine << "\n" << bottomLine << "\n";

}

void testDiagonalInverseExtendAlignment()
{

	int individualTests = 0;

	for(unsigned int graphIteration = 1; graphIteration <= 100; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		std::cout << "Generate random genome to align to...\n" << std::flush;
		diploidGenomeString gS = generateRandomGenome(80, true);
		_printDiploidGenomeString(gS);

		std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);


		std::cout << "Create GraphAligner...\n" << std::flush;
		GraphAligner_affine gA(gS_graph, 5);

		auto alignDiagonalInverse = [&](std::string sequence_with_gaps, std::vector<int> sequence_characterOrigin) {
			std::string randomString_noGaps;
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < sequence_with_gaps.size(); cI++)
			{
				char string_character = sequence_with_gaps.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_characterOrigin.push_back(sequence_characterOrigin.at(cI));
				}
			}


			// normal affine needleman

			std::string aligned_randomString_noGaps;
			std::string aligned_graph;
			std::vector<int> aligned_graph_levels;
			int score;
			gA.fullNeedleman_affine(randomString_noGaps, score, aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			int expectedScore = gA.score_fullNeedleman_affine(aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			assert(expectedScore == score);

			std::cout << "\t" << aligned_graph << "\n";
			std::cout << "\t" << aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned score: " << score << ", expected: " << expectedScore << "\n\n";


			// parameter bounds
			unsigned int levels = gS_graph->NodesPerLevel.size();
			unsigned int sequenceLength = randomString_noGaps.length();
			int max_levelI = levels - 1;
			int max_seqI = sequenceLength;

			// empty NWt

			VirtualNWTable NWt(&gA, &randomString_noGaps);

			// non-inverse local needleman



			std::vector<localExtension_pathDescription> alignments_fwd = gA.fullNeedleman_affine_diagonal_extension(
					randomString_noGaps, 0, 0, 0, 0, &NWt, true, true);
			assert(alignments_fwd.size() == 1);
			for(unsigned int alignmentI = 0; alignmentI < alignments_fwd.size(); alignmentI++)
			{
				localExtension_pathDescription pD = alignments_fwd.at(alignmentI);
				int expectedScore_fwd = gA.score_fullNeedleman_affine(pD.alignedGraph, pD.alignedGraph_levels, pD.alignedSequence);
				assert(expectedScore_fwd == pD.Score);
				assert(pD.Score == score);
				std::cout << "FWD Needleman with score = " << pD.Score << ":\n";
				std::cout << "\t" << pD.alignedSequence << "\n" << std::flush;

			}

			// non inverse backward needleman!
			std::vector<localExtension_pathDescription> alignments_bwd = gA.fullNeedleman_affine_diagonal_extension(
					randomString_noGaps, max_seqI, max_levelI, 0, 0, &NWt, false, true);

			assert(alignments_bwd.size() == 1);
			for(unsigned int alignmentI = 0; alignmentI < alignments_bwd.size(); alignmentI++)
			{
				localExtension_pathDescription pD = alignments_bwd.at(alignmentI);
				int expectedScore_bwd = gA.score_fullNeedleman_affine(pD.alignedGraph, pD.alignedGraph_levels, pD.alignedSequence);
				assert(expectedScore_bwd == pD.Score);
				std::cout << "BWD Needleman with score = " << pD.Score << " ( expected score: " << expectedScore_bwd << ") :\n";
				std::cout << "\t" << pD.alignedSequence << "\n" << std::flush;
				assert(pD.Score == score);
			}

			/*
			 * for later

			unsigned int covered_noGap_characters = 0;
			unsigned int validatable_noGap_characters = 0;
			unsigned int validatable_noGap_characters_OK = 0;

			for(unsigned int alignedI = 0; alignedI < aligned_randomString_noGaps.size(); alignedI++)
			{
				char alignedC = aligned_randomString_noGaps.at(alignedI);
				int originGraph = aligned_graph_levels.at(alignedI);
				if(alignedC != '_')
				{
					int origin = randomString_characterOrigin.at(covered_noGap_characters);
					if(origin != -1)
					{
						validatable_noGap_characters++;
						if(origin == originGraph)
						{
							validatable_noGap_characters_OK++;
						}
					}
					covered_noGap_characters++;
				}
			}
			assert(covered_noGap_characters == randomString_noGaps.size());

			std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n";
			*/
		};

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			int stringStart;
			int stringStop;

			sampleStringFromGraph(gS_graph, randomString, stringStart, stringStop);

			std::string randomString_aligned = Utilities::repeatString("_", stringStart) +
					randomString + Utilities::repeatString("_", gS_graph->NodesPerLevel.size() - 1 - stringStop);

			// baseline test: align non-modified random string
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(stringStart+cI);
			}

			std::cout << "1) BASELINE" << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment: " << "\n\n\t" << randomString_aligned << "\n\n";

			alignDiagonalInverse(randomString, randomString_characterOrigin);

			// extended test: align
			std::pair<std::string, std::vector<int>> modified_randomString = Utilities::modifySequence(randomString, randomString_characterOrigin);

			std::cout << "\n2) EXTENDED MODIFIED" << "\n\n" << "Original string:" << "\n" << randomString << "\n" << "Modified string:" << "\n" << modified_randomString.first << "\n\n";

			alignDiagonalInverse(modified_randomString.first, modified_randomString.second);

			std::cout << "\n\n";

			individualTests++;

		}

		delete(gS_graph);
	}

	std::cout << "testDiagonalInverseExtendAlignment(): " << individualTests << " tests were succesful!\n" << std::flush;
}



void testExtensionStop()
{
	int individualTests = 0;

	bool verbose = true;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(200);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_forSeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,
					10,
					15,
					3,
					4,
					3,
					7,
					7,
					7
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

//			std::string randomString_aligned = Utilities::repeatString("_", stringStart) +
//					randomString + Utilities::repeatString("_", gS_graph->NodesPerLevel.size() - 1 - stringStop);


			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
				}
			}

			// find chains for un-gapped random string
			std::vector<kMerChain> seq_chains = gA.getGI()->findChains(randomString_noGaps);
			assert(seq_chains.size() > 0);

			VirtualNWTable vNW(&gA, &randomString_noGaps);
			std::vector<NWPath*> found_NWpaths;

			std::cout << "Add found paths to vNW..\n";
			for(unsigned int chainI = 0; chainI < seq_chains.size(); chainI++)
			{
				kMerChain& thisChain = seq_chains.at(chainI);
				std::vector<NWPath*> seed_NW_paths = gA.findVirtualNWExactSeeds(randomString_noGaps, thisChain.sequence_begin, thisChain.sequence_end, thisChain.graph_firstLevel, thisChain.graph_lastLevel);

				if(verbose)
					std::cout << "\tChain " << chainI << ", found seed paths: " << seed_NW_paths.size() << "\n" << std::flush;

				found_NWpaths.insert(found_NWpaths.end(), seed_NW_paths.begin(), seed_NW_paths.end());

				for(unsigned int seedPathI = 0; seedPathI < seed_NW_paths.size(); seedPathI++)
				{
					NWPath* thisSeedPath = seed_NW_paths.at(seedPathI);

					if(verbose)
						std::cout << "\t\tSeed path " << seedPathI << "\n" << std::flush;

						vNW.addPath(thisSeedPath);
				}
			}

			std::cout << "Extend and expect to see hit NW paths...\n" << std::flush;

			for(unsigned int pathI = 0; pathI < found_NWpaths.size(); pathI++)
			{
				std::cout << "\tNWPath " << pathI << "\n";

				NWPath* thisSeedPath = found_NWpaths.at(pathI);

				assert(thisSeedPath->first_edges.size() == 1);
				NWEdge* firstSeedEdge = *(thisSeedPath->first_edges.begin());

				assert(thisSeedPath->last_edges.size() == 1);
				NWEdge* lastSeedEdge = *(thisSeedPath->last_edges.begin());

				std::vector<localExtension_pathDescription> backwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
						randomString_noGaps,
						firstSeedEdge->from_y,
						firstSeedEdge->from_x,
						firstSeedEdge->from_z,
						-20,
						&vNW,
						false);

				std::vector<localExtension_pathDescription> forwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
						randomString_noGaps,
						lastSeedEdge->to_y,
						lastSeedEdge->to_x,
						lastSeedEdge->to_z,
						-20,
						&vNW,
						true);

				if(verbose)
				{
					std::cout << "\tBackward extensions: " << backwardExtensions.size() << "\n" << std::flush;
					std::cout << "\tForward extensions: " << forwardExtensions.size() << "\n";
				}

			}

			individualTests++;
		}

		delete(gS_graph);
	}

	std::cout << "testExtensionStop(): " << individualTests << " tests.\n" << std::flush;
}



void test_Simple_longRange_SeedAndExtend()
{
	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	size_t achievableMatches = 0;
	size_t achievedMatches = 0;

	size_t rightPositionCharacters_counted = 0;
	size_t rightPositionCharacters_correct = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(800);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_nonAffine gA_nonAffine(gS_graph, aligner_kMerSize);
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;


		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,

					10,
					15,

					2,
					4,

					2,
					3,

					0,
					0,

					30,
					40,

					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

			int minimumAchievableScore = gA_nonAffine.score_fullNeedleman(underlyingEdges, underlyingEdges_levels, randomString);
			int minimumAchievableMatches = gA.countMatchesInSequence(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose || true)
				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore <<  " and " << minimumAchievableMatches << " matches of sequence string: " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;

			achievableMatches += minimumAchievableMatches;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

			if(verbose)
				std::cout << "Start full-string alignment...\n" << std::flush;

			seedAndExtend_return wholeString_alignment = gA.seedAndExtend2(randomString_noGaps, true);

			{
				seedAndExtend_return& thisAlignment = wholeString_alignment;

				int thisAlignmentScore = gA_nonAffine.score_fullNeedleman(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);
				int thisAlignmentMatches = gA.countMatchesInSequence(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);

				if(true || verbose)
				{
					std::cout << "\tAlignment  " << " [internal score " << thisAlignment.Score << "]>\n";
					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
					std::cout << "\t\tAffine NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;
				}

				if((thisAlignmentScore >= thisAlignmentScore) && (thisAlignmentMatches >= minimumAchievableMatches))
				{
					individualTests_fullSuccessful++;
				}

//				assert(thisAlignmentMatches >= minimumAchievableMatches);

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

				if(verbose || true)
					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

				achievedMatches += thisAlignmentMatches;

				rightPositionCharacters_counted += validatable_noGap_characters;
				rightPositionCharacters_correct += validatable_noGap_characters_OK;
			}

			individualTests++;
		}

		delete(gS_graph);
	}

	assert(achievableMatches > 0);
	assert(rightPositionCharacters_counted > 0);
	std::cout << "test_Simple_longRange_SeedAndExtend(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
	std::cout << "\t Matches: " << achievedMatches  << " / " << achievableMatches << " => " << ((double)achievedMatches/(double)achievableMatches) << "\n";
	std::cout << "\t Positions: " << rightPositionCharacters_correct  << " / " << rightPositionCharacters_counted << " => " << (double(rightPositionCharacters_correct)/(double)rightPositionCharacters_counted) << "\n";

	std::cout << std::flush;
}


void test_Simple_longRange_SeedAndExtend_3()
{
	// this thing uses "mutation" and extension

	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(800);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;


		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,

					10,
					15,

					2,
					4,

					2,
					3,

					3,
					7,

					30,
					40,

					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

			int minimumAchievableScore = gA.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose || true)
				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\n in alignment with expected score " << minimumAchievableScore <<  " and " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

			if(verbose)
				std::cout << "Start full-string alignment...\n" << std::flush;


			seedAndExtend_return wholeString_alignments_seedAndExtend2 = gA.seedAndExtend2(randomString_noGaps, false, 1);
			seedAndExtend_return wholeString_alignments_seedAndExtend3 = gA.seedAndExtend3(randomString_noGaps, false);

			std::cerr << "wholeString_alignments_seedAndExtend2.Score: " << wholeString_alignments_seedAndExtend2.Score << "\n";
			std::cerr << "wholeString_alignments_seedAndExtend3.Score: " << wholeString_alignments_seedAndExtend3.Score << "\n" << std::flush;

			// std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "    wholeString_alignments_distanceModel_2.Score: " << wholeString_alignments_distanceModel_2.Score << "\n" << std::flush;

			{
				seedAndExtend_return& thisAlignment_sAE2 = wholeString_alignments_seedAndExtend2;
				seedAndExtend_return& thisAlignment_sAE3 = wholeString_alignments_seedAndExtend3;

				// int thisAlignmentScore_sAE2 = gA.score_fullNeedleman_affine(thisAlignment_sAE2.graph_aligned, thisAlignment_sAE2.graph_aligned_levels, thisAlignment_sAE2.sequence_aligned);
				int thisAlignmentScore_sAE3 = gA.score_fullNeedleman_affine(thisAlignment_sAE3.graph_aligned, thisAlignment_sAE3.graph_aligned_levels, thisAlignment_sAE3.sequence_aligned);

				assert(thisAlignment_sAE3.Score == thisAlignmentScore_sAE3);
				assert(thisAlignment_sAE3.Score == thisAlignment_sAE2.Score);
			}

			individualTests_fullSuccessful++;
			individualTests++;
		}

		delete(gS_graph);
	}

	std::cout << "test_Simple_longRange_SeedAndExtend_3(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";

	std::cout << std::flush;
}


void test_Simple_longRange_SeedAndExtend_4()
{
	// this thing uses "mutation" and extension

	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	size_t achievableMatches = 0;
	size_t achievedMatches = 0;

	size_t rightPositionCharacters_counted = 0;
	size_t rightPositionCharacters_correct = 0;

	bool verbose = false;
	bool repeat = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS;
		std::string file_for_gS = "_serializedGenomeString.txt";

		if(repeat)
		{
			std::cerr << "Repeat!\n" << std::flush;
			gS = readGenomeStringFromFile(file_for_gS);
		}
		else
		{
			gS = generateRandomGenome(800);
			storeGenomeStringInFile(gS, file_for_gS);
		}
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_endsFree gA(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;


		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;

			std::string file_for_randomString = "_serializedRandomString.txt";

			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(repeat)
			{
				std::ifstream randomString_fH;
				randomString_fH.open(file_for_randomString.c_str());
				if(! randomString_fH.is_open())
				{
					throw std::runtime_error("Cannot open genome string file for reading: "+ file_for_randomString);
				}
				std::getline(randomString_fH, randomString);
				randomString_fH.close();

				std::string randomString_noGaps;
				for(unsigned int cI = 0; cI < randomString.size(); cI++)
				{
					char string_character = randomString.at(cI);
					if(string_character != '_')
					{
						randomString_noGaps.push_back(string_character);
					}
				}


				seedAndExtend_return wholeString_alignments_seedAndExtend4 = gA.seedAndExtend4(randomString_noGaps, false);
				std::cerr << "\n\nwholeString_alignments_seedAndExtend4.Score: " << wholeString_alignments_seedAndExtend4.Score << "\n\n" << std::flush;
				int thisAlignmentScore_sAE4 = gA.score_endsFree(wholeString_alignments_seedAndExtend4.graph_aligned, wholeString_alignments_seedAndExtend4.graph_aligned_levels, wholeString_alignments_seedAndExtend4.sequence_aligned);
				std::cerr << "\nthisAlignmentScore_sAE4: " << thisAlignmentScore_sAE4 << "\n" << std::flush;
				assert(wholeString_alignments_seedAndExtend4.Score == thisAlignmentScore_sAE4);
			}
			else
			{
				if(verbose)
					std::cout << "Sample possible emission from graph...\n" << std::flush;

				sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
						gS_graph,
						aligner_kMerSize,
						gA.get_S_match(),
						gA.get_S_mismatch(),
						gA.get_S_gapOpen(),
						gA.get_S_gapExtend(),
						randomString,
						underlyingEdges,
						underlyingEdges_levels,
						stringStart,
						stringStop,

						10,
						15,

						2,
						4,

						2,
						3,

						3,
						7,

						30,
						40,

						true
				);

				std::ofstream randomString_fH;
				randomString_fH.open(file_for_randomString.c_str());
				if(! randomString_fH.is_open())
				{
					throw std::runtime_error("Cannot open genome string file for writing: "+ file_for_randomString);
				}
				randomString_fH << randomString;
				randomString_fH.close();

				int minimumAchievableMatches = gA.countMatchesInSequence(underlyingEdges, underlyingEdges_levels, randomString);
				achievableMatches += minimumAchievableMatches;


				if(verbose)
					std::cout << "... sampling done!\n" << std::flush;

				assert(randomString.length() == underlyingEdges.length());
				assert(underlyingEdges_levels.size() == underlyingEdges.length());

				int minimumAchievableScore = gA.score_endsFree(underlyingEdges, underlyingEdges_levels, randomString);

				if(verbose || true)
					std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\n in alignment with expected score " << minimumAchievableScore <<  " and " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;

				// find correct alignment positions for characters in aligned string
				std::vector<int> randomString_characterOrigin;
				std::vector<bool> randomString_unmodifiedCharacterFromGraph;
				for(unsigned int cI = 0; cI < randomString.size(); cI++)
				{
					randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
					if(randomString.at(cI) == underlyingEdges.at(cI))
					{
						randomString_unmodifiedCharacterFromGraph.push_back(true);
					}
					else
					{
						randomString_unmodifiedCharacterFromGraph.push_back(false);
					}
				}

				// remove gaps from generated random str
				std::string randomString_noGaps;
				std::vector<int> randomString_noGaps_characterOrigin;
				std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
				for(unsigned int cI = 0; cI < randomString.size(); cI++)
				{
					char string_character = randomString.at(cI);
					if(string_character != '_')
					{
						randomString_noGaps.push_back(string_character);
						randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
					}
				}

				if(verbose)
					std::cout << "Start full-string alignment...\n" << std::flush;


				// seedAndExtend_return wholeString_alignments_seedAndExtend2 = gA.seedAndExtend2(randomString_noGaps, false, 1);
				seedAndExtend_return wholeString_alignments_seedAndExtend4 = gA.seedAndExtend4(randomString_noGaps, false);

				// std::cerr << "wholeString_alignments_seedAndExtend2.Score: " << wholeString_alignments_seedAndExtend2.Score << "\n";
				std::cerr << "\n\nwholeString_alignments_seedAndExtend4.Score: " << wholeString_alignments_seedAndExtend4.Score << "\n\n" << std::flush;

				// std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "    wholeString_alignments_distanceModel_2.Score: " << wholeString_alignments_distanceModel_2.Score << "\n" << std::flush;


				// seedAndExtend_return& thisAlignment_sAE2 = wholeString_alignments_seedAndExtend2;
				seedAndExtend_return& thisAlignment_sAE4 = wholeString_alignments_seedAndExtend4;

				// int thisAlignmentScore_sAE2 = gA.score_fullNeedleman_affine(thisAlignment_sAE2.graph_aligned, thisAlignment_sAE2.graph_aligned_levels, thisAlignment_sAE2.sequence_aligned);
				int thisAlignmentScore_sAE4 = gA.score_endsFree(thisAlignment_sAE4.graph_aligned, thisAlignment_sAE4.graph_aligned_levels, thisAlignment_sAE4.sequence_aligned);
				std::cerr << "\nthisAlignmentScore_sAE4: " << thisAlignmentScore_sAE4 << "\n" << std::flush;

				assert(thisAlignment_sAE4.Score == thisAlignmentScore_sAE4);
				// assert(thisAlignment_sAE4.Score == thisAlignment_sAE2.Score);

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment_sAE4.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment_sAE4.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment_sAE4.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

				if(verbose || true)
					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

				int thisAlignmentMatches = gA.countMatchesInSequence(thisAlignment_sAE4.graph_aligned, thisAlignment_sAE4.graph_aligned_levels, thisAlignment_sAE4.sequence_aligned);
				achievedMatches += thisAlignmentMatches;

				rightPositionCharacters_counted += validatable_noGap_characters;
				rightPositionCharacters_correct += validatable_noGap_characters_OK;


				individualTests_fullSuccessful++;
				individualTests++;
			}
		}

		delete(gS_graph);
	}

	std::cout << "test_Simple_longRange_SeedAndExtend_4(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
	assert(achievableMatches > 0);
	assert(rightPositionCharacters_counted > 0);
	std::cout << "\t Matches: " << achievedMatches  << " / " << achievableMatches << " => " << ((double)achievedMatches/(double)achievableMatches) << "\n";
	std::cout << "\t Positions: " << rightPositionCharacters_correct  << " / " << rightPositionCharacters_counted << " => " << (double(rightPositionCharacters_correct)/(double)rightPositionCharacters_counted) << "\n";

	std::cout << std::flush;
}


void test_Simple_longRange_SeedAndExtend_2()
{
	// this thing uses "mutation" and extension

	int individualTests = 0;
	int individualTests_fullSuccessful = 0;

	size_t achievableMatches = 0;
	size_t achievedMatches = 0;

	size_t rightPositionCharacters_counted = 0;
	size_t rightPositionCharacters_correct = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(800);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		std::cout << "\t...done!\n" << std::flush;


		int test_iterations = 1;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_for_Simple_longRange_SeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,

					10,
					15,

					2,
					4,

					2,
					3,

					3,
					7,

					30,
					40,

					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

			int minimumAchievableScore = gA.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);
			int minimumAchievableMatches = gA.countMatchesInSequence(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose || true)
				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore <<  " and " << minimumAchievableMatches << " matches of sequence string: " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;


			achievableMatches += minimumAchievableMatches;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

			if(verbose)
				std::cout << "Start full-string alignment...\n" << std::flush;


			seedAndExtend_return wholeString_alignments = gA.seedAndExtend2(randomString_noGaps, false, 1);

			// std::cerr << "wholeString_alignments.Score: " << wholeString_alignments.Score << "    wholeString_alignments_distanceModel_2.Score: " << wholeString_alignments_distanceModel_2.Score << "\n" << std::flush;

			{
				seedAndExtend_return& thisAlignment = wholeString_alignments;

				int thisAlignmentScore = gA.score_fullNeedleman_affine(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);
				int thisAlignmentMatches = gA.countMatchesInSequence(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);

				if(true || verbose)
				{
					std::cout << "\tAlignment [internal score " << thisAlignment.Score << "]>\n";
					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
					std::cout << "\t\tAffine NW score_ " << thisAlignmentScore << ", matches in sequence " << thisAlignmentMatches << "\n" << std::flush;
				}

				if((thisAlignmentScore >= thisAlignmentScore) && (thisAlignmentMatches >= minimumAchievableMatches))
				{
					individualTests_fullSuccessful++;
				}

				assert(thisAlignment.Score == thisAlignmentScore);

//				assert(thisAlignmentMatches >= minimumAchievableMatches);

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

				if(verbose || true)
					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

				achievedMatches += thisAlignmentMatches;

				rightPositionCharacters_counted += validatable_noGap_characters;
				rightPositionCharacters_correct += validatable_noGap_characters_OK;

//				assert(1 == 0);
			}

			individualTests++;
		}

		delete(gS_graph);
	}

	assert(achievableMatches > 0);
	assert(rightPositionCharacters_counted > 0);
	std::cout << "test_Simple_longRange_SeedAndExtend_2(): " << individualTests << " tests, of which " << individualTests_fullSuccessful << " were fully successful.\n";
	std::cout << "\t Matches: " << achievedMatches  << " / " << achievableMatches << " => " << ((double)achievedMatches/(double)achievableMatches) << "\n";
	std::cout << "\t Positions: " << rightPositionCharacters_correct  << " / " << rightPositionCharacters_counted << " => " << (double(rightPositionCharacters_correct)/(double)rightPositionCharacters_counted) << "\n";

	std::cout << std::flush;
}


void testSeedAndExtend()
{
	int individualTests = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(300);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_forSeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,
					10,
					15,
					2,
					2,
					5,
					40,
					4,
					30,
					false,
					true
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

			int minimumAchievableScore = gA.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose || true)
				std::cout << "1) BASELINE " << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore << ": " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
				}
			}

			if(verbose)
				std::cout << "Start full-string alignment...\n" << std::flush;

			seedAndExtend_return wholeString_alignment = gA.seedAndExtend2(randomString_noGaps);

			assert (1 == 0);

			{
				seedAndExtend_return& thisAlignment = wholeString_alignment;

				int thisAlignmentScore = gA.score_fullNeedleman_affine(thisAlignment.graph_aligned, thisAlignment.graph_aligned_levels, thisAlignment.sequence_aligned);

				if(verbose)
				{
					std::cout << "\tAlignment " << " [internal score " << thisAlignment.Score << "]>\n";
					std::cout << "\t\t" << thisAlignment.graph_aligned << "\n";
					std::cout << "\t\t" << thisAlignment.sequence_aligned << "\n";
					std::cout << "\t\tAffine NW score_ " << thisAlignmentScore << "\n" << std::flush;
				}

				unsigned int covered_noGap_characters = 0;
				unsigned int validatable_noGap_characters = 0;
				unsigned int validatable_noGap_characters_OK = 0;

				for(unsigned int alignedI = 0; alignedI < thisAlignment.sequence_aligned.size(); alignedI++)
				{
					char alignedC = thisAlignment.sequence_aligned.at(alignedI);
					int originGraph = thisAlignment.graph_aligned_levels.at(alignedI);
					if(alignedC != '_')
					{
						int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
						if(origin != -1)
						{
							validatable_noGap_characters++;
							if(origin == originGraph)
							{
								validatable_noGap_characters_OK++;
							}
						}
						covered_noGap_characters++;
					}
				}

				if(verbose)
					std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;
			}

			individualTests++;
		}

		delete(gS_graph);
	}

	std::cout << "testSeedAndExtend(): " << individualTests << " tests.\n" << std::flush;
}


void testSeedAndExtend_Algorithm()
{
	int individualTests = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(500);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_forSeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop,
					20,
					30,
					10,
					12,
					5,
					40,
					4,
					30,
					false
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

			int minimumAchievableScore = gA.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose)
				std::cout << "1) BASELINE" << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore << ": " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n" << std::flush;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}


			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
				}
			}

			// find chains for un-gapped random string
			std::vector<kMerChain> seq_chains = gA.getGI()->findChains(randomString_noGaps);
			assert(seq_chains.size() > 0);

			VirtualNWTable vNW(&gA, &randomString_noGaps);
			std::vector<NWPath*> found_NWpaths;

			if(verbose)
				std::cout << "Add found paths to vNW..\n";

			for(unsigned int chainI = 0; chainI < seq_chains.size(); chainI++)
			{
				kMerChain& thisChain = seq_chains.at(chainI);
				std::vector<NWPath*> seed_NW_paths = gA.findVirtualNWExactSeeds(randomString_noGaps, thisChain.sequence_begin, thisChain.sequence_end, thisChain.graph_firstLevel, thisChain.graph_lastLevel);

				if(verbose)
					std::cout << "\tChain " << chainI << ", found seed paths: " << seed_NW_paths.size() << "\n" << std::flush;

				found_NWpaths.insert(found_NWpaths.end(), seed_NW_paths.begin(), seed_NW_paths.end());

				for(unsigned int seedPathI = 0; seedPathI < seed_NW_paths.size(); seedPathI++)
				{
					NWPath* thisSeedPath = seed_NW_paths.at(seedPathI);

					if(verbose)
						std::cout << "\t\tSeed path " << seedPathI << "\n" << std::flush;

						vNW.addPath(thisSeedPath);
				}
			}

			if(verbose)
				std::cout << "Extend and expect to see hit NW paths...\n" << std::flush;

			VirtualNWTable vNW2(&gA, &randomString_noGaps);

			for(unsigned int pathI = 0; pathI < found_NWpaths.size(); pathI++)
			{
				if(verbose)
					std::cout << "\tNWPath " << pathI << "\n";

				NWPath* thisSeedPath = found_NWpaths.at(pathI);

				assert(thisSeedPath->first_edges.size() == 1);
				NWEdge* firstSeedEdge = *(thisSeedPath->first_edges.begin());

				assert(thisSeedPath->last_edges.size() == 1);
				NWEdge* lastSeedEdge = *(thisSeedPath->last_edges.begin());

				std::vector<localExtension_pathDescription> backwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
						randomString_noGaps,
						firstSeedEdge->from_y,
						firstSeedEdge->from_x,
						firstSeedEdge->from_z,
						-20,
						&vNW,
						false);

				std::vector<localExtension_pathDescription> forwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
						randomString_noGaps,
						lastSeedEdge->to_y,
						lastSeedEdge->to_x,
						lastSeedEdge->to_z,
						-20,
						&vNW,
						true);

				NWPath* newPath = thisSeedPath->clonePathWithoutTable();

				if(verbose)
					std::cout << "Created newPath.\n" << std::flush;

				if(backwardExtensions.size() > 0)
				{
					newPath->entry_edges.clear();
					newPath->first_edges.clear();
				}

				if(forwardExtensions.size() > 0)
				{
					newPath->exit_edges.clear();
					newPath->last_edges.clear();
				}

				//newPath->_printPath();

				for(unsigned int bwI = 0; bwI < backwardExtensions.size(); bwI++)
				{
//					std::cout << "\tintegrate bwI " << bwI << " / " << backwardExtensions.size() << "\n" << std::flush;

					localExtension_pathDescription& this_bwE = backwardExtensions.at(bwI);
					newPath->takeInExtensionPath(&this_bwE, -1);
//					newPath->_printPath();
				}

				for(unsigned int fwI = 0; fwI < forwardExtensions.size(); fwI++)
				{
//					std::cout << "\tintegrate fwI " << fwI << " / " << forwardExtensions.size() << "\n" << std::flush;

					localExtension_pathDescription& this_fwE = forwardExtensions.at(fwI);
					newPath->takeInExtensionPath(&this_fwE, 1);
//					newPath->_printPath();
				}

				newPath->recalculateFirstLast();

//				std::cout << "Final print after first/last recalculation!\n" << std::flush;
//				newPath->_printPath();

				vNW2.addPath(newPath);

				if(verbose)
				{
					std::cout << "\tBackward extensions: " << backwardExtensions.size() << "\n" << std::flush;
					std::cout << "\tForward extensions: " << forwardExtensions.size() << "\n";
				}
			}

			vNW2.checkConsistency();

			std::cout << "Added all resulting extensions paths. In total: " << vNW2.getNumPaths() << " paths and " << vNW2.getNumEdges() << " edges.\n" << std::flush;

			individualTests++;

			vNW.freeMemory();
			vNW2.freeMemory();
		}



		delete(gS_graph);
	}

	std::cout << "testSeedAndExtend_Algorithm(): " << individualTests << " tests.\n" << std::flush;
}


void testChainFindingAndExtension()
{
	int individualTests = 0;
	int individualTests_successful = 0;

	bool verbose = false;

	for(unsigned int graphIteration = 1; graphIteration <= 100; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		if(verbose)
			std::cout << "Generate random genome to align to...\n" << std::flush;

		diploidGenomeString gS = generateRandomGenome(200);
		_printDiploidGenomeString(gS);

		if(verbose)
			std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);

		if(verbose)
			std::cout << "Create GraphAligner...\n" << std::flush;

		int aligner_kMerSize = 5;
		GraphAligner_affine gA(gS_graph, aligner_kMerSize);

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			std::string underlyingEdges;
			std::vector<int> underlyingEdges_levels;
			int stringStart;
			int stringStop;

			if(verbose)
				std::cout << "Sample possible emission from graph...\n" << std::flush;

			sampleStringFromGraph_forSeedAndExtend(
					gS_graph,
					aligner_kMerSize,
					gA.get_S_match(),
					gA.get_S_mismatch(),
					gA.get_S_gapOpen(),
					gA.get_S_gapExtend(),
					randomString,
					underlyingEdges,
					underlyingEdges_levels,
					stringStart,
					stringStop
			);

			if(verbose)
				std::cout << "... sampling done!\n" << std::flush;

			assert(randomString.length() == underlyingEdges.length());
			assert(underlyingEdges_levels.size() == underlyingEdges.length());

//			std::string randomString_aligned = Utilities::repeatString("_", stringStart) +
//					randomString + Utilities::repeatString("_", gS_graph->NodesPerLevel.size() - 1 - stringStop);

			int minimumAchievableScore = gA.score_fullNeedleman_affine(underlyingEdges, underlyingEdges_levels, randomString);

			if(verbose)
				std::cout << "1) BASELINE" << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore << ": " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n";

			std::cout << std::flush;

			// find correct alignment positions for characters in aligned string
			std::vector<int> randomString_characterOrigin;
			std::vector<bool> randomString_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(underlyingEdges_levels.at(cI));
				if(randomString.at(cI) == underlyingEdges.at(cI))
				{
					randomString_unmodifiedCharacterFromGraph.push_back(true);
				}
				else
				{
					randomString_unmodifiedCharacterFromGraph.push_back(false);
				}
			}

			// remove gaps from generated random str
			std::string randomString_noGaps;
			std::vector<int> randomString_noGaps_characterOrigin;
			std::vector<int> randomString_noGaps_unmodifiedCharacterFromGraph;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				char string_character = randomString.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_noGaps_characterOrigin.push_back(randomString_characterOrigin.at(cI));
					randomString_noGaps_unmodifiedCharacterFromGraph.push_back(randomString_unmodifiedCharacterFromGraph.at(cI));
				}
			}

			// find chains for un-gapped random string
			std::vector<kMerChain> seq_chains = gA.getGI()->findChains(randomString_noGaps);
			assert(seq_chains.size() > 0);

			bool haveOneGoodChain = false;
			for(int firstIterationLimitChains = 0; firstIterationLimitChains <= 1; firstIterationLimitChains++)
			{
				if(firstIterationLimitChains == 1)
				{
					std::cout << "\n\nALL-CHAIN SEARCH STARTED!\n\n" << std::flush;
					std::cout << "1) BASELINE" << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment with expected score " << minimumAchievableScore << ": " << "\n\n\t" << underlyingEdges << "\n\t" << randomString << "\n\n";

				}
				for(unsigned int chainI = 0; chainI < seq_chains.size(); chainI++)
				{
					kMerChain& thisChain = seq_chains.at(chainI);

					if(firstIterationLimitChains == 1)
					{
						std::cout << "Chain " << chainI << "\n";
						std::cout << "\t" << "lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret" << ": " << lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret << "\n";
						std::cout << "\t" << "lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret" << ": " << lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret << "\n" << std::flush;
						std::cout << "\t" << "lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret" << ": " << lastCall_sampleString_forSeedAndExtend_trunk_begin_sequence_ret << "\n" << std::flush;
						std::cout << "\t" << "lastCall_sampleString_forSeedAndExtend_trunk_end_sequence_ret" << ": " << lastCall_sampleString_forSeedAndExtend_trunk_end_sequence_ret << "\n" << std::flush;

						std::cout << "\t" << "graph_firstLevel" << ": " << thisChain.graph_firstLevel << "\n" << std::flush;
						std::cout << "\t" << "graph_lastLevel" << ": " << thisChain.graph_lastLevel << "\n" << std::flush;
						std::cout << "\t" << "sequence_begin" << ": " << thisChain.sequence_begin << "\n" << std::flush;
						std::cout << "\t" << "sequence_end" << ": " << thisChain.sequence_end << "\n" << std::flush;
					}

					if(!((thisChain.graph_firstLevel <= lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret) && (thisChain.graph_lastLevel >= lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret)))
					{
						if(firstIterationLimitChains == 0)
							continue;
					}

					std::vector<NWPath*> seed_NW_paths = gA.findVirtualNWExactSeeds(randomString_noGaps, thisChain.sequence_begin, thisChain.sequence_end, thisChain.graph_firstLevel, thisChain.graph_lastLevel);

					if(firstIterationLimitChains)
						std::cout << "Chain " << chainI << ", found seed paths: " << seed_NW_paths.size() << "\n" << std::flush;

					for(unsigned int seedPathI = 0; seedPathI < seed_NW_paths.size(); seedPathI++)
					{
						NWPath* thisSeedPath = seed_NW_paths.at(seedPathI);

						if(firstIterationLimitChains)
							std::cout << "\tSeed path " << seedPathI << "\n" << std::flush;

						VirtualNWTable vNW(&gA, &randomString_noGaps);
						vNW.addPath(thisSeedPath);

						std::string seed_reconstructed_sequence;
						std::string seed_reconstructed_graph;
						std::vector<int> seed_reconstructed_graph_levels;
						vNW._testTracePath(thisSeedPath, seed_reconstructed_sequence, seed_reconstructed_graph, seed_reconstructed_graph_levels);
						std::vector<Edge*> seed_alignedGraph_edges = thisSeedPath->graphEdgesPath();
						assert(seed_alignedGraph_edges.size() == seed_reconstructed_sequence.size());

						if(firstIterationLimitChains)
							std::cout << "\t\t" << seed_reconstructed_graph << "\n\t\t" << seed_reconstructed_sequence << "\n" << std::flush;

						std::vector<int> seed_alignedGraph_levels;
	//					for(unsigned int eI = 0; eI < seed_alignedGraph_edges.size(); eI++)
	//					{
	//						Edge* e = seed_alignedGraph_edges.at(eI);
	//						string emission = gS_graph->CODE.deCode(e->locus_id, e->emission);
	//						assert(emission.length() == 1);
	//						if(emission == "_")
	//						{
	//							seed_alignedGraph_levels.push_back(-1);
	//						}
	//						else
	//						{
	//							seed_reconstructed_graph_levels.push_back(e->From->level);
	//						}
	//					}

						assert(thisSeedPath->first_edges.size() == 1);
						NWEdge* firstSeedEdge = *(thisSeedPath->first_edges.begin());

						assert(thisSeedPath->last_edges.size() == 1);
						NWEdge* lastSeedEdge = *(thisSeedPath->last_edges.begin());

	//					std::cerr << "firstSeedEdge:\n";
	//					std::cerr << "\t" << "to_x: " << firstSeedEdge->to_x << "\n";
	//					std::cerr << "\t" << "to_y: " << firstSeedEdge->to_y << "\n";
	//					std::cerr << "\t" << "to_z: " << firstSeedEdge->to_z << "\n";
	//					std::cerr << "chain.sequence_begin: " << thisChain.sequence_begin << "\n";
	//					std::cerr << "chain.graph_firstLevel: " << thisChain.graph_firstLevel << "\n" << std::flush;

	//					std::cerr << "lastSeedEdge:\n";
	//					std::cerr << "\t" << "to_x: " << lastSeedEdge->to_x << "\n";
	//					std::cerr << "\t" << "to_y: " << lastSeedEdge->to_y << "\n";
	//					std::cerr << "\t" << "to_z: " << lastSeedEdge->to_z << "\n";
	//					std::cerr << "thisChain.sequence_end: " << thisChain.sequence_end << "\n";
	//					std::cerr << "thisChain.graph_lastLevel: " << thisChain.graph_lastLevel << "\n" << std::flush;

						std::vector<localExtension_pathDescription> backwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
								randomString_noGaps,
								firstSeedEdge->from_y,
								firstSeedEdge->from_x,
								firstSeedEdge->from_z,
								-20,
								&vNW,
								false);

						std::vector<localExtension_pathDescription> forwardExtensions = gA.fullNeedleman_affine_diagonal_extension(
								randomString_noGaps,
								lastSeedEdge->to_y,
								lastSeedEdge->to_x,
								lastSeedEdge->to_z,
								-20,
								&vNW,
								true);

						if(firstIterationLimitChains)
						{
							std::cout << "\tBackward extensions: " << backwardExtensions.size() << "\n" << std::flush;
							std::cout << "\tForward extensions: " << forwardExtensions.size() << "\n";
						}

						for(int backwardI = -1; backwardI < (int)backwardExtensions.size(); backwardI++)
						{
							for(int forwardI = -1; forwardI < (int)forwardExtensions.size(); forwardI++)
							{
								if(firstIterationLimitChains)
									std::cout << "\t\t" << backwardI << "/" << forwardI << "\n";
								localExtension_pathDescription bwE = (backwardI != -1) ? backwardExtensions.at(backwardI) : localExtension_pathDescription();
								localExtension_pathDescription fwE = (forwardI != -1) ? forwardExtensions.at(forwardI) : localExtension_pathDescription();

								std::string fullAlignedSequence = ((backwardI != -1) ? bwE.alignedSequence : "") + seed_reconstructed_sequence + ((forwardI != -1 ) ? fwE.alignedSequence: "");
								std::string fulllAlignedGraph = ((backwardI != -1) ? bwE.alignedGraph : "") + seed_reconstructed_graph + ((forwardI != -1 ) ? fwE.alignedGraph : "");
								std::vector<int> fullAlignedLevels = ((backwardI != -1) ? bwE.alignedGraph_levels : std::vector<int>());
								fullAlignedLevels.insert(fullAlignedLevels.end(), seed_reconstructed_graph_levels.begin(), seed_reconstructed_graph_levels.end());
								if(forwardI != -1)
									fullAlignedLevels.insert(fullAlignedLevels.end(), fwE.alignedGraph_levels.begin(), fwE.alignedGraph_levels.end());

								if(firstIterationLimitChains)
								{
									std::cout << "\t\t\t" << fulllAlignedGraph << " [graph]" << "\n";
									std::cout << "\t\t\t" << fullAlignedSequence << " [sequence]" << "\n" << std::flush;
								}

								int thisAlignmentScore = gA.score_fullNeedleman_affine(fulllAlignedGraph, fullAlignedLevels, fullAlignedSequence);

								if(firstIterationLimitChains)
									std::cout << "\t\t\t" << thisAlignmentScore << " [score]" << "\n" << std::flush;

								for(unsigned int lI = 0; lI < fullAlignedLevels.size(); lI++)
								{
									if(lI > 0)
									{
										assert((fullAlignedLevels.at(lI) == -1) || (fullAlignedLevels.at(lI-1) == -1) || (fullAlignedLevels.at(lI) == (fullAlignedLevels.at(lI-1)+1)));
									}
								}

								if(thisAlignmentScore >= minimumAchievableScore)
								{
									haveOneGoodChain = true;


									unsigned int covered_noGap_characters = 0;
									unsigned int validatable_noGap_characters = 0;
									unsigned int validatable_noGap_characters_OK = 0;

									for(unsigned int alignedI = 0; alignedI < fullAlignedSequence.size(); alignedI++)
									{
										char alignedC = fullAlignedSequence.at(alignedI);
										int originGraph = fullAlignedLevels.at(alignedI);
										if(alignedC != '_')
										{
											int origin = randomString_noGaps_characterOrigin.at(covered_noGap_characters);
											if(origin != -1)
											{
												validatable_noGap_characters++;
												if(origin == originGraph)
												{
													validatable_noGap_characters_OK++;
												}
	//											std::cerr << alignedI << " " << alignedC << " " << origin << " " << originGraph << "\n" << std::flush;
											}
											covered_noGap_characters++;
										}
									}

									// assert(covered_noGap_characters == randomString_noGaps.size());

									if(firstIterationLimitChains)
										std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n" << std::flush;

								}
							}
						}
					}
				}

				if(haveOneGoodChain)
				{
					firstIterationLimitChains = 2;
				}
			}

			if(haveOneGoodChain)
			{
				individualTests_successful++;
			}

			if(! haveOneGoodChain)
			{
				std::cerr << "! haveOneGoodChain!\n\n" << std::flush;
				_printDiploidGenomeString(gS);
				std::cout << std::flush;
				std::cerr << std::flush;
			}

			// assert(haveOneGoodChain);

			individualTests++;
		}

		delete(gS_graph);
	}

	std::cout << "testChainFindingAndExtension(): " << individualTests << " tests, of which " << individualTests_successful << " were succesful!\n" << std::flush;
}


void testGraphAligner()
{
	int individualTests = 0;

	for(unsigned int graphIteration = 1; graphIteration <= 10; graphIteration++)
	{
		std::cout << "Graph test iteration " << graphIteration << "\n\n=============================================================\n\n=============================================================\n\n";

		std::cout << "Generate random genome to align to...\n" << std::flush;
		diploidGenomeString gS = generateRandomGenome(80);
		_printDiploidGenomeString(gS);

		std::cout << "Generate graph from genome...\n" << std::flush;
		Graph* gS_graph = genomeString2Graph(gS);


		std::cout << "Create GraphAligner...\n" << std::flush;
		GraphAligner_nonAffine gA(gS_graph, 5);
		GraphAligner_affine gA_affine(gS_graph, 5);
		GraphAligner_endsFree gA_endsFree(gS_graph, 5);

		auto alignString = [&](std::string sequence_with_gaps, std::vector<int> sequence_characterOrigin) {
			std::string randomString_noGaps;
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < sequence_with_gaps.size(); cI++)
			{
				char string_character = sequence_with_gaps.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_characterOrigin.push_back(sequence_characterOrigin.at(cI));
				}
			}

			std::string aligned_randomString_noGaps;
			std::string aligned_graph;
			std::vector<int> aligned_graph_levels;
			int score;
			gA.fullNeedleman(randomString_noGaps, score, aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			int expectedScore = gA.score_fullNeedleman(aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			assert(expectedScore == score);

			std::cout << "\t" << aligned_graph << "\n";
			std::cout << "\t" << aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned score: " << score << ", expected: " << expectedScore << "\n\n";

			unsigned int covered_noGap_characters = 0;
			unsigned int validatable_noGap_characters = 0;
			unsigned int validatable_noGap_characters_OK = 0;

			for(unsigned int alignedI = 0; alignedI < aligned_randomString_noGaps.size(); alignedI++)
			{
				char alignedC = aligned_randomString_noGaps.at(alignedI);
				int originGraph = aligned_graph_levels.at(alignedI);
				if(alignedC != '_')
				{
					int origin = randomString_characterOrigin.at(covered_noGap_characters);
					if(origin != -1)
					{
						validatable_noGap_characters++;
						if(origin == originGraph)
						{
							validatable_noGap_characters_OK++;
						}
					}
					covered_noGap_characters++;
				}
			}
			assert(covered_noGap_characters == randomString_noGaps.size());

			std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n";

			std::string diagonally_aligned_graph;
			std::string diagonally_aligned_randomString_noGaps;
			std::vector<int> diagonally_aligned_graph_levels;
			int diagonal_score;
			gA.fullNeedleman_diagonal(randomString_noGaps, diagonal_score, diagonally_aligned_graph, diagonally_aligned_graph_levels, diagonally_aligned_randomString_noGaps);

			std::cout << "\t" << diagonally_aligned_graph << "\n";
			std::cout << "\t" << diagonally_aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned score: " << diagonal_score << "\n\n" << std::flush;

			assert(score == diagonal_score);
		};

		auto alignString_affine = [&](std::string sequence_with_gaps, std::vector<int> sequence_characterOrigin) {
			std::string randomString_noGaps;
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < sequence_with_gaps.size(); cI++)
			{
				char string_character = sequence_with_gaps.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_characterOrigin.push_back(sequence_characterOrigin.at(cI));
				}
			}

			std::string aligned_randomString_noGaps;
			std::string aligned_graph;
			std::vector<int> aligned_graph_levels;
			int score;
			gA_affine.fullNeedleman_affine(randomString_noGaps, score, aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			int expectedScore = gA_affine.score_fullNeedleman_affine(aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);

			std::cout << "\t" << aligned_graph << "\n";
			std::cout << "\t" << aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned (affine normal) score: " << score << ", expected: " << expectedScore << "\n\n" << std::flush;

			assert(expectedScore == score);

			unsigned int covered_noGap_characters = 0;
			unsigned int validatable_noGap_characters = 0;
			unsigned int validatable_noGap_characters_OK = 0;

			for(unsigned int alignedI = 0; alignedI < aligned_randomString_noGaps.size(); alignedI++)
			{
				char alignedC = aligned_randomString_noGaps.at(alignedI);
				int originGraph = aligned_graph_levels.at(alignedI);
				if(alignedC != '_')
				{
					int origin = randomString_characterOrigin.at(covered_noGap_characters);
					if(origin != -1)
					{
						validatable_noGap_characters++;
						if(origin == originGraph)
						{
							validatable_noGap_characters_OK++;
						}
					}
					covered_noGap_characters++;
				}
			}
			assert(covered_noGap_characters == randomString_noGaps.size());

			std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n";

			std::string diagonally_aligned_graph;
			std::string diagonally_aligned_randomString_noGaps;
			std::vector<int> diagonally_aligned_graph_levels;
			int diagonal_score;
			gA_affine.fullNeedleman_affine_diagonal(randomString_noGaps, diagonal_score, diagonally_aligned_graph, diagonally_aligned_graph_levels, diagonally_aligned_randomString_noGaps);

			std::cout << "\t" << diagonally_aligned_graph << "\n";
			std::cout << "\t" << diagonally_aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned (affine diagonal) score: " << diagonal_score << "\n\n" << std::flush;

			assert(score == diagonal_score);
		};

		auto alignString_endsFree = [&](std::string sequence_with_gaps, std::vector<int> sequence_characterOrigin) {
			std::string randomString_noGaps;
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < sequence_with_gaps.size(); cI++)
			{
				char string_character = sequence_with_gaps.at(cI);
				if(string_character != '_')
				{
					randomString_noGaps.push_back(string_character);
					randomString_characterOrigin.push_back(sequence_characterOrigin.at(cI));
				}
			}

			std::string aligned_randomString_noGaps;
			std::string aligned_graph;
			std::vector<int> aligned_graph_levels;
			int score;
			gA_endsFree.fullNeedleman_endsFree(randomString_noGaps, score, aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);
			int expectedScore = gA_endsFree.score_endsFree(aligned_graph, aligned_graph_levels, aligned_randomString_noGaps);

			std::cout << "\t" << aligned_graph << "\n";
			std::cout << "\t" << aligned_randomString_noGaps << "\n";
			std::cout << "\t" << "Returned (ends-free affine) score: " << score << ", expected: " << expectedScore << "\n\n" << std::flush;

			assert(expectedScore == score);

			unsigned int covered_noGap_characters = 0;
			unsigned int validatable_noGap_characters = 0;
			unsigned int validatable_noGap_characters_OK = 0;

			for(unsigned int alignedI = 0; alignedI < aligned_randomString_noGaps.size(); alignedI++)
			{
				char alignedC = aligned_randomString_noGaps.at(alignedI);
				int originGraph = aligned_graph_levels.at(alignedI);
				if(alignedC != '_')
				{
					int origin = randomString_characterOrigin.at(covered_noGap_characters);
					if(origin != -1)
					{
						validatable_noGap_characters++;
						if(origin == originGraph)
						{
							validatable_noGap_characters_OK++;
						}
					}
					covered_noGap_characters++;
				}
			}
			assert(covered_noGap_characters == randomString_noGaps.size());

			std::cout << "Can validate " << validatable_noGap_characters << " characters, " << validatable_noGap_characters_OK << " at right position!\n";
		};

		int test_iterations = 10;
		for(int iteration = 1; iteration <= test_iterations; iteration++)
		{
			std::cout << "Test iteration " << iteration << "\n==============================================\n(individual iteration " << individualTests << ")\n\n" << std::flush;

			// generate random string
			std::string randomString;
			int stringStart;
			int stringStop;

			sampleStringFromGraph(gS_graph, randomString, stringStart, stringStop);

			std::string randomString_aligned = Utilities::repeatString("_", stringStart) +
					randomString + Utilities::repeatString("_", gS_graph->NodesPerLevel.size() - 1 - stringStop);

			// baseline test: align non-modified random string
			std::vector<int> randomString_characterOrigin;
			for(unsigned int cI = 0; cI < randomString.size(); cI++)
			{
				randomString_characterOrigin.push_back(stringStart+cI);
			}

			std::cout << "1) BASELINE" << "\n\n" << "String" << "\n" << randomString << ",\n\nin alignment: " << "\n\n\t" << randomString_aligned << "\n\n";
			alignString_affine(randomString, randomString_characterOrigin);

			// extended test: align
			std::pair<std::string, std::vector<int>> modified_randomString = Utilities::modifySequence(randomString, randomString_characterOrigin);

			std::cout << "\n2) EXTENDED MODIFIED" << "\n\n" << "Original string:" << "\n" << randomString << "\n" << "Modified string:" << "\n" << modified_randomString.first << "\n\n";
			alignString_affine(modified_randomString.first, modified_randomString.second);

			std::cout << "\n3) ENDS-FREE nonmodified" << "\n\n" << "Original string:" << "\n" << randomString << "\n\n";
			alignString_endsFree(randomString, randomString_characterOrigin);

			std::cout << "\n4) ENDS-FREE modified" << "\n\n" << "Original string:" << "\n" << randomString << "\n" << "Modified string:" << "\n" << modified_randomString.first << "\n\n";
			alignString_endsFree(modified_randomString.first, modified_randomString.second);

			std::cout << "\n\n";

			individualTests++;
		}

		delete(gS_graph);
	}

	std::cout << "testGraphAligner(): " << individualTests << " tests were succesful!\n" << std::flush;
}

void testGraphAligner_old()
{
	std::cout << "Define diploidGenomeString...\n" << std::flush;

/*
Graph in text:

                     AT        ATA
ACGTACGATCGTACGACGCGT  ACGTACGT   ACGTACGT
                     GC        ___


Graph in tabs:

 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	A	T	 	 	 	 	 	 	 	 	A	T	A
A	C	G	T	A	C	G	A	T	C	G	T	A	C	G	A	C	G	C	G	T	 	 	A	C	G	T	A	C	G	T	 	 	 	A	C	G	T	A	C	G	T
 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	G	C	 	 	 	 	 	 	 	 	_	_	_


 */
	diploidGenomeString gS;
	std::vector<std::string> firstCompartment;
	firstCompartment.push_back("AGCATCATCGGCAATTCTGCAATACTCGAACCGCTCGT");

	std::vector<std::string> secondCompartment;
	secondCompartment.push_back("AT");
	secondCompartment.push_back("GC");

	std::vector<std::string> thirdCompartment;
	thirdCompartment.push_back("ACGTACGT");

	std::vector<std::string> fourthCompartment;
	fourthCompartment.push_back("ATA");
	fourthCompartment.push_back("___");

	std::vector<std::string> fifthCompartment;
	fifthCompartment.push_back("TTTCACAT");

	gS.push_back(firstCompartment);
	gS.push_back(secondCompartment);
	gS.push_back(thirdCompartment);
	gS.push_back(fourthCompartment);
	gS.push_back(fifthCompartment);

	std::cout << "Build graph from diploidGenomeString...\n" << std::flush;

	Graph* gS_graph = genomeString2Graph(gS);

	std::string seq1 = "ACGTACGTATCTTTCACAT";

	/*
	std::cout << "Index graph...\n" << std::flush;

	GraphAndIndex graph_index(gS_graph, 5);

	std::cout << "\tdone.\n" << std::flush;

	graph_index.printIndex();

	// rC CGTCGTACGATCGTACGT

	std::vector<kMerChain> seq1_chains = graph_index.findChains(seq1);
	std::cout << "Chains for sequence " << seq1 << "\n\n" << std::flush;

	_printChains(seq1_chains);

	*/

	GraphAligner gA(gS_graph, 5);

	// gA.fullNeedleman_affine(seq1);

	// gA.seedAndExtend("ACGTACGTTTTCACAT");


}


