/*
 * GraphAlignernonAffine.cpp
 *
 *  Created on: 27.07.2013
 *      Author: AlexanderDilthey
 */

#include "GraphAlignernonAffine.h"

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

GraphAligner_nonAffine::GraphAligner_nonAffine(Graph* graph, int k) : GraphAligner(graph, k) {
	// TODO Auto-generated constructor stub

}

void GraphAligner_nonAffine::fullNeedleman(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned)
{
	std::vector< std::vector < std::vector <double> > > m;
	std::vector< std::vector < std::vector < backtraceStep > > > m_backtrace;

	bool verbose = false;

	// init matrix
	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	m.resize(levels);
	m_backtrace.resize(levels);
	for(unsigned int levelI = 0; levelI < levels; levelI++)
	{
		m.at(levelI).resize(sequenceLength+1);
		m_backtrace.at(levelI).resize(sequenceLength+1);
		for(unsigned int seqI = 0; seqI <= sequenceLength; seqI++)
		{
			unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();
			m.at(levelI).at(seqI).resize(statesPerLevel, 0);
			m_backtrace.at(levelI).at(seqI).resize(statesPerLevel);
		}
	}


	// matrix needs "gap borders" -- gaps in sequence
	for(unsigned int levelI = 0; levelI < levels; levelI++)
	{
		unsigned int seqI = 0;
		unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();

		for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
		{

			double score;
			backtraceStep bS;
			if(levelI > 0)
			{
				std::vector<double> scores;
				std::vector <backtraceStep> scores_backSteps;

				std::vector<std::pair<int, Edge*> > previousZs = _graph_get_previous_z_values_and_edges(levelI, stateI);
				assert(previousZs.size() > 0);

				for(unsigned int zI = 0; zI < previousZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = previousZs.at(zI);

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					double score_gapInSequence = m.at(levelI-1).at(seqI).at(thisZjump.first) + ((edgeEmission == "_") ? S_graphGap : S_gap);

					backtraceStep backtrack_gapInSequence;
					backtrack_gapInSequence.x = levelI-1;
					backtrack_gapInSequence.y = seqI;
					backtrack_gapInSequence.z = thisZjump.first;
					backtrack_gapInSequence.usedEdge = thisZjump.second;

					scores.push_back(score_gapInSequence);
					scores_backSteps.push_back(backtrack_gapInSequence);
				}

				std::pair<double, unsigned int> maxScore = Utilities::findVectorMax(scores);
				score = maxScore.first;
				bS = scores_backSteps.at(maxScore.second);
			}
			else
			{
				score = 0;
			}

			m.at(levelI).at(seqI).at(stateI) = score;
			m_backtrace.at(levelI).at(seqI).at(stateI) = bS;
		}
	}

	// matrix needs "gap borders" -- gaps in graph
	for(unsigned int seqI = 0; seqI <= sequenceLength; seqI++)
	{
		double score = S_gap * seqI;

		unsigned int levelI = 0;
		unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();

		for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
		{
			m.at(levelI).at(seqI).at(stateI) = score;

			backtraceStep bS;
			if(seqI > 0)
			{
				bS.x = 0;
				bS.y = seqI - 1;
				bS.z = stateI;
				bS.usedEdge = 0;
			}
			m_backtrace.at(levelI).at(seqI).at(stateI) = bS;
		}
	}

	for(unsigned int levelI = 1; levelI < levels; levelI++)
	{
		for(unsigned int seqI = 1; seqI <= sequenceLength; seqI++)
		{
			unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();

			for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
			{
				std::vector<double> scores;
				std::vector <backtraceStep> scores_backSteps;

				std::string sequenceEmission = sequence.substr(seqI - 1, 1);

				// gap in graph
				double score_gapInGraph = m.at(levelI).at(seqI-1).at(stateI) + S_gap;
				backtraceStep backtrack_gapInGraph;
				backtrack_gapInGraph.x = levelI;
				backtrack_gapInGraph.y = seqI-1;
				backtrack_gapInGraph.z = stateI;
				backtrack_gapInGraph.usedEdge = 0;
				scores.push_back(score_gapInGraph);
				scores_backSteps.push_back(backtrack_gapInGraph);

				// match / mismatch
				std::vector<std::pair<int, Edge*> > previousZs = _graph_get_previous_z_values_and_edges(levelI, stateI);
				for(unsigned int zI = 0; zI < previousZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = previousZs.at(zI);

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					double score_MatchMismatch = m.at(levelI-1).at(seqI-1).at(thisZjump.first) + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

					// std::cout << "Sequence: " << (seqI - 1) << " graph: " << (levelI-1) << " sequenceEmission: " << sequenceEmission << " edge emission: " << edgeEmission << " score: " << score_MatchMismatch << "\n";

					backtraceStep backtrack_MatchMismatch;
					backtrack_MatchMismatch.x = levelI-1;
					backtrack_MatchMismatch.y = seqI-1;
					backtrack_MatchMismatch.z = thisZjump.first;
					backtrack_MatchMismatch.usedEdge = thisZjump.second;

					scores.push_back(score_MatchMismatch);
					scores_backSteps.push_back(backtrack_MatchMismatch);
				}

				// gap in sequence
				for(unsigned int zI = 0; zI < previousZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = previousZs.at(zI);

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					if(verbose && (levelI == 19))
					{
						std::cout << "gap in sequence at level " << levelI << " / seqI " << seqI << " / stateI " << stateI << "\n";
						std::cout << "previous score: "  << m.at(levelI-1).at(seqI).at(thisZjump.first) << "\n";
						std::cout << "previous level: " << (levelI - 1) << "\n";
						std::cout << "previous seqI: " << (seqI) << "\n";
						std::cout << "previous z: " << thisZjump.first << "\n";
						std::cout << "jump score: " << ((edgeEmission == "_") ? S_graphGap : S_gap) << "\n";
					}
					double score_gapInSequence = m.at(levelI-1).at(seqI).at(thisZjump.first) + ((edgeEmission == "_") ? S_graphGap : S_gap);

					backtraceStep backtrack_gapInSequence;
					backtrack_gapInSequence.x = levelI-1;
					backtrack_gapInSequence.y = seqI;
					backtrack_gapInSequence.z = thisZjump.first;
					backtrack_gapInSequence.usedEdge = thisZjump.second;

					scores.push_back(score_gapInSequence);
					scores_backSteps.push_back(backtrack_gapInSequence);
				}

				std::pair<double, unsigned int> maxScore = Utilities::findVectorMax(scores);

				if(verbose && (seqI == sequenceLength))
				{
					std::cout << "levelI: " << levelI << "\n";
					std::cout << "seqI: " << seqI << "\n";
					std::cout << "stateI: " << stateI << "\n";

					for(unsigned int i = 0; i < scores.size(); i++)
					{
						std::cout << "Alternative " << i << ":\n";
						std::cout << "\t" << scores.at(i) << "\n";
					}

					std::cout << "SELECT #" << maxScore.second << " with " << maxScore.first << "\n\n";
				}


				m.at(levelI).at(seqI).at(stateI) = maxScore.first;
				m_backtrace.at(levelI).at(seqI).at(stateI) = scores_backSteps.at(maxScore.second);
			}
		}
	}

	if(verbose)
	{
		std::cout << "MATRIX:\n";
		for(unsigned int levelI = 0; levelI < levels; levelI++)
		{
			for(unsigned int seqI = 0; seqI <= sequenceLength; seqI++)
			{
				unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();
				std::vector<std::string> stateScores;
				for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
				{
					int d_x = levelI - m_backtrace.at(levelI).at(seqI).at(stateI).x;
					int d_y = seqI - m_backtrace.at(levelI).at(seqI).at(stateI).y;
					int b_z = m_backtrace.at(levelI).at(seqI).at(stateI).z;

					std::string stateString = Utilities::DtoStr(m.at(levelI).at(seqI).at(stateI)) + "[" +
							Utilities::ItoStr(d_x) + "," +
							Utilities::ItoStr(d_y) + "," +
							Utilities::ItoStr(b_z) + "]";

					stateScores.push_back(stateString);
				}
				std::string statesScores = Utilities::join(stateScores, "; ");
				std::cout << statesScores << "\t";
			}
			std::cout << "\n";
		}
		std::cout << "\n";

		// backtrace
		std::cout << "Levels: " << levels << "\n";
		std::cout << "sequenceLength: " << sequenceLength << "\n";
	}
	std::pair<double, unsigned int> maxFinalState = Utilities::findVectorMax(m.at(levels-1).at(sequenceLength));

	if(verbose)
		std::cout << "Maximum final score " << maxFinalState.first << " in state " << maxFinalState.second << "\n";

	int backtrace_x = levels - 1;
	int backtrace_y = sequenceLength;
	int backtrace_z = maxFinalState.second;

	std::string reconstructed_graph;
	std::string reconstructed_sequence;
	std::vector<int> reconstructed_graph_levels;

	// print matrix

	while((backtrace_x != 0) || (backtrace_y != 0))
	{
		if(verbose)
		{
			std::cout << "backtrace_x: " << backtrace_x << "\n";
			std::cout << "backtrace_y: " << backtrace_y << "\n";
			std::cout << "backtrace_z: " << backtrace_z << "\n";
			std::cout << "Score: " << m.at(backtrace_x).at(backtrace_y).at(backtrace_z) << "\n\n";
		}

		assert((backtrace_x >= 0) && (backtrace_x <= (levels - 1)));
		assert((backtrace_y >= 0) && (backtrace_y <= (int)sequenceLength));
		assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));

		backtraceStep& step = m_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z);

		std::string edgeEmission;
		if(step.usedEdge != 0)
		{
			edgeEmission = g->CODE.deCode(step.usedEdge->locus_id, step.usedEdge->emission);
			assert(edgeEmission.size() == 1);
		}
		std::string sequenceEmission;
		if(backtrace_y >= 1)
		{
			sequenceEmission = sequence.substr(backtrace_y - 1, 1);
		}

		int next_x = step.x;
		int next_y = step.y;
		int next_z = step.z;

		if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
		{
			// match or mismatch
			assert(step.usedEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(backtrace_x - 1);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
		{
			// gap in graph
			reconstructed_graph.append("_");
			reconstructed_graph_levels.push_back(-1);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
		{
			// gap in sequence
			assert(step.usedEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(backtrace_x - 1);
			reconstructed_sequence.append("_");
		}
		else
		{
			assert(1 == 0);
		}

		backtrace_x = next_x;
		backtrace_y = next_y;
		backtrace_z = next_z;
	}

	std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
	std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
	std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());


	if(verbose)
		std::cout << "\n" << reconstructed_graph << "\n" << reconstructed_sequence << "\n\n";

	score = maxFinalState.first;
	graph_aligned = reconstructed_graph;
	graph_aligned_levels = reconstructed_graph_levels;
	sequence_aligned = reconstructed_sequence;
}

int GraphAligner_nonAffine::score_fullNeedleman(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence)
{
	assert(reconstructed_graph.length() == reconstructed_sequence.length());
	int score = 0;
	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::string G = reconstructed_graph.substr(i, 1);
		if(G == "_")
		{
			if(S == "_")
			{
				assert(reconstructed_graph_levels.at(i) != -1);
				score += S_graphGap;
			}
			else
			{
				if(reconstructed_graph_levels.at(i) != -1)
				{
					score += S_mismatch;
				}
				else
				{
					score += S_gap;
				}
			}
		}
		else
		{
			if(S == "_")
			{
				score += S_gap;
			}
			else
			{
				if(S == G)
				{
					score += S_match;
				}
				else
				{
					score += S_mismatch;
				}
			}
		}
	}

	return score;
}


std::map<NWEdge*, std::map<NWEdge*, int> > GraphAligner_nonAffine::naivelyCalculateGraphDistancesForVirtualNW_nonAffine(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, int> >& lastPositionDistances_perZ)
{
	assert(S_graphGap == 0);

	bool verbose = true;
	lastPositionDistances_perZ.clear();

	std::vector<NWEdge*> entryEdges_vector = vNW.getEntryEdges();
	std::vector<NWEdge*> exitEdges_vector = vNW.getExitEdges();

//	std::set<NWEdge*> entryEdges_set(entryEdges_vector.begin(), entryEdges_vector.end());
//	std::set<NWEdge*> exitEdges_set(exitEdges_vector.begin(), exitEdges_vector.end());

	std::map<int, std::map<int, std::set<NWEdge*> > >  entryEdges_from_coordinates;
	std::map<int, std::map<int, std::set<NWEdge*> > >  exitEdges_to_coordinates;

//	std::map<Edge*, std::set<NWEdge*> > graphEdges_2_NW_entry;
	for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
	{
		NWEdge* edge = entryEdges_vector.at(eI);
		entryEdges_from_coordinates[edge->from_x][edge->from_z].insert(edge);
//		Edge* graphEdge = edge->usedGraphEdge;
//		graphEdges_2_NW_entry[graphEdge].insert(edge);
//		assert(exitEdges_set.count(edge) == 0);
	}

//	std::map<Edge*, std::set<NWEdge*> > graphEdges_2_NW_exit;
	for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
	{
		NWEdge* edge = exitEdges_vector.at(eI);
		exitEdges_to_coordinates[edge->to_x][edge->to_z].insert(edge);
//		Edge* graphEdge = edge->usedGraphEdge;
//		graphEdges_2_NW_exit[graphEdge].insert(edge);
	}

	if(verbose)
		std::cout << "naivelyCalculateGraphDistancesForVirtualNW(..) for " << entryEdges_vector.size() <<  " entry NW edges and " << exitEdges_vector.size() << " exit NW edges.\n" << std::flush;

	std::map<Node*, std::map<NWEdge*, int> > runningNodeDistances;
	std::set<Node*> nodes_l0 = g->NodesPerLevel.at(0);

	std::map<NWEdge*, std::map<NWEdge*, int> > NWedge_distances;

	for(unsigned int lI = 0; lI < g->NodesPerLevel.size(); lI++)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<NWEdge*, int> > runningNodeDistances_thisLevel;

		if(lI == (g->NodesPerLevel.size() - 1))
		{
			lastPositionDistances_perZ.resize(nodes_thisLevel.size());
		}

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;
			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == 0)
			{
				runningNodeDistances_thisLevel[n] = std::map<NWEdge*, int>();
				runningNodeDistances_thisLevel[n][0] = 0;
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);
				runningNodeDistances_thisLevel[n] = std::map<NWEdge*, int>();

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;
					assert(runningNodeDistances.count(nodeFrom) > 0);

					for(std::map<NWEdge*, int>::iterator edgeDistIt = runningNodeDistances.at(nodeFrom).begin(); edgeDistIt != runningNodeDistances.at(nodeFrom).end(); edgeDistIt++)
					{
						NWEdge* interestingEdge = edgeDistIt->first;
						int newDistance;
						if(edgeEmission == "_")
						{
							newDistance = edgeDistIt->second;
						}
						else
						{
							newDistance = edgeDistIt->second + 1;
						}
						if((runningNodeDistances_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_thisLevel.at(n).at(interestingEdge) > newDistance))
						{
							runningNodeDistances_thisLevel.at(n)[interestingEdge] = newDistance;
						}
					}

					if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
					{
						std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
						for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
						{
							NWEdge* e = *nwEdgeIt;
							runningNodeDistances_thisLevel.at(n)[e] = 0;
						}
					}
				}
			}

			if(lI != (g->NodesPerLevel.size() -1 ))
			{
				if(entryEdges_from_coordinates.count(lI) && entryEdges_from_coordinates.at(lI).count(z))
				{
					std::set<NWEdge*> nwEdges = entryEdges_from_coordinates.at(lI).at(z);
					for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
					{
						NWEdge* e = *nwEdgeIt;
						for(std::map<NWEdge*, int>::iterator exitEdgeIt = runningNodeDistances_thisLevel.at(n).begin(); exitEdgeIt != runningNodeDistances_thisLevel.at(n).end(); exitEdgeIt++)
						{
							NWEdge* exitEdge = exitEdgeIt->first;
							assert(e != 0);
							if((exitEdge == 0) || (e->from_y >= exitEdge->to_y))
							{
								int distance = exitEdgeIt->second;
								assert(distance >= 0);
								NWedge_distances[e][exitEdge] = distance;
							}
						}
					}
				}
			}
			else
			{
				for(std::map<NWEdge*, int>::iterator edgeDistIt = runningNodeDistances_thisLevel.at(n).begin(); edgeDistIt != runningNodeDistances_thisLevel.at(n).end(); edgeDistIt++)
				{
					NWEdge* interestingEdge = edgeDistIt->first;
					int thisDistance = edgeDistIt->second;

					if((lastPositionDistances_perZ.at(z).count(interestingEdge) == 0) || (lastPositionDistances_perZ.at(z).at(interestingEdge) > thisDistance))
					{
						lastPositionDistances_perZ.at(z)[interestingEdge] = thisDistance;
					}
				}
			}
		}

		runningNodeDistances = runningNodeDistances_thisLevel;
	}

	size_t calculated_distances = 0;
	for(std::map<NWEdge*, std::map<NWEdge*, int> >::iterator edgeIt = NWedge_distances.begin(); edgeIt != NWedge_distances.end(); edgeIt++)
	{
		for(std::map<NWEdge*, int>::iterator edgeIt2 = edgeIt->second.begin(); edgeIt2 != edgeIt->second.end(); edgeIt2++)
		{
			calculated_distances++;
		}
	}

	std::cout << "naivelyCalculateGraphDistancesForVirtualNW_nonAffine(..): Calculated " << calculated_distances << " distances.\n" << std::flush;

	return NWedge_distances;
}

void GraphAligner_nonAffine::seedAndExtend(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned, bool testing)
{
	/*
	 *	Description of the employed seed-and-extend algorithm
	 *
	 *	1)	We find all stretches of kMers in sequence that have a contiguous sequence of matches in the graph
	 *	   	(This does not necessarily imply that there is a path in the graph connecting the kMers)
	 *
	 *	2)	We translate the matched kMer sequences into lists of paths in the virtual NW table. We employ a greedy
	 *	  	algorithm which will follow a kMer sequence through the graph, and abort if it encounters a mismatch (defined
	 *	  	as a character in the kMer sequence which is not present in the graph at the current level). If this does not occur
	 *	  	until the end of the sequence, we have an exact match. If it occurs somewhere in the middle and the matched stretch is
	 *	  	long enough, we retain the matched subpath. After that event, we re-initialize the greedy matcher at the position after
	 *	  	the mismatch and see whether we can identify other subpaths of sufficient length (i.e. going back to the first step, ignoring
	 *	   	the matched bits).
	 *
	 *	3)	We have a data structure to store path information. Each path is defined as a single trunk sequence of (mostly) diagonal edges in the
	 *		virtual NW table, with a (possibly empty) tree of extension paths (which may contain gap edges) either side of the trunk. Each
	 *		path thus consists of a series of edges between nodes in the virtual NW table, which can be ordered by graph and sequence level.
	 *		(We will create an appropriate index data structure).
	 *
	 *	4)	After having translated the exact kMer matches into corresponding paths (the exact matches are path trunks) we carry out local extension
	 *		at both ends of the trunk. We employ the following algorithm:
	 *			- we use a diagonalized affine NW algorithm and initialize the first cell to a value x > 0.
	 *			- we carry out NW until all cells in a diagonal have fallen to a value < 0.
	 *			- if at any point we hit an edge being present in another path, we store the distance betwenen the two paths and retain
	 *			  the corresponding extension leading to the second path as one path in the generated tree of paths. We assign minusInfinity to all edges which
	 *			  are present in other paths.
	 *			- once all values in the last computed diagonal have fallen to below 0, we identify the maximum (maxima, to be more precise) in the
	 *			  covered matrix area. If the maximum(a) is > x, we backtrack and make sure that each path from origin to maximum ends up as a path in
	 *			  the tree.
	 *			- we create the necessary tree data structures for each extended local match.
	 *
	 *	5)	We traverse the leaf nodes of the left-hand-trees of the extended matches from left-to-right, top-to-bottom, and compute a maximum entry value
	 *		for each lead node (these entry values will then of course determine the end point values of the extended matches).
	 *		For each leaf node with coordinates x, y, z:
	 *			We scan for end paths with end coordinates <=x, <=y which are z-compatible
	 *			We scan for paths with start coordinates <x, <y which go through either (x, <=y, z) or (<=x, y, z-compatible)
	 *			Both classes of potential starting coordinates are connected with our extension using (possible two) affine gaps.
	 *			The y-distance for two points is always clear.
	 *			The x-distance has to take graph structure into account.
	 *			We will exactly compute the graph distance between all ending nodes and all starting nodes.
	 *			This, however, does not quite solve the jump problem. We will also compute the distance between all extension starting nodes (from
	 *			which only once will be active at a time) and store the graph distance between selected starting node and extension alongside each edge score.
	 *
	 *
	 *	path data structure:
	 *		bool valuesComputed: whether the scores for the edges of the path have been computed
	 *		set<*edges>
	 *
	 *	edge data structure:
	 *		int x1, int y1, int z1, int x2, int y2, int z2: coordinates of the edge
	 *		int runningScore: running maximum difference score to the beginning.
	 *		int runningGraphDistance: graph distance between beginning of path and end-of-edge
	 *		*path
	 *
	 *
	 *	index-over-table data structure:
	 *		hash int -> edge*: edges starting from x-level int
	 *		hash int-> edge*: edges starting from y-level int
	 *		edgeExists(x1, y1, z1, x2, y2, z2) -> bool: function to test for existence of an edge
	 *
	 *
	 *
	 *
	 *
	 */

	std::vector<kMerChain> seq_chains = gI.findChains(sequence);
	VirtualNWTable vNW(this, &sequence);

	std::cout << "Found " << seq_chains.size() << " chains.\n" << std::flush;

	if(testing)
	{
		std::cerr << "Warning - testing on, will ignore chain because not covering graph from " << lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret << " to " << lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret << "!\n" << std::flush;
	}
	for(unsigned int chainI = 0; chainI < seq_chains.size(); chainI++)
	{
		kMerChain& chain = seq_chains.at(chainI);

		std::string coveredSequence = sequence.substr(chain.sequence_begin, chain.sequence_end - chain.sequence_begin + 1);
//		std::cout << "Chain " << chainI << "\n";
//		std::cout << "\t" << "matchedkMers" << ": " << chain.matchedkMers << "\n";
//		std::cout << "\t" << "coveredSequence" << ": " << coveredSequence << "\n";
//		std::cout << "\t" << "sequence_begin" << ": " << chain.sequence_begin << "\n";
//		std::cout << "\t" << "sequence_end" << ": " << chain.sequence_end << "\n" << std::flush;
//		std::cout << "\t" << "graph_firstLevel" << ": " << chain.graph_firstLevel << "\n" << std::flush;
//		std::cout << "\t" << "graph_lastLevel" << ": " << chain.graph_lastLevel << "\n" << std::flush;


		if(testing)
		{
			if(!((chain.graph_firstLevel <= lastCall_sampleString_forSeedAndExtend_trunk_begin_graph_ret) && (chain.graph_lastLevel >= lastCall_sampleString_forSeedAndExtend_trunk_end_graph_ret)))
			{
				continue;
			}
		}

		// assert(chain.matchedkMers == (coveredSequence.length() - kMerSize + 1));

		//std::vector<std::pair<std::vector<Edge*>, std::string>> edgePaths = findEdgePathForSequence(coveredSequence, chain.graph_firstLevel, chain.graph_lastLevel);

		std::vector<NWPath*> seed_NW_paths = findVirtualNWExactSeeds(sequence, chain.sequence_begin, chain.sequence_end, chain.graph_firstLevel, chain.graph_lastLevel);

		std::cout << "Chain " << chainI << ", covered sequence " << coveredSequence << ", found seed paths: " << seed_NW_paths.size() << "\n" << std::flush;

		if(seed_NW_paths.size() > 0)
		{
			vNW.addPath(seed_NW_paths.at(0));

			std::string reconstructed_sequence;
			std::string reconstructed_graph;
			std::vector<int> reconstructed_graph_levels;
			vNW._testTracePath(seed_NW_paths.at(0), reconstructed_sequence, reconstructed_graph, reconstructed_graph_levels);

			std::cout << "\tChain traceback:\n";
			std::cout << "\t" << reconstructed_graph << " [graph]" << "\n";
			std::cout << "\t" << reconstructed_sequence << " [sequence]" << "\n" << std::flush;

			std::vector<localExtension_pathDescription> backwardExtensions;// = fullNeedleman_affine_diagonal_extension(sequence, chain.sequence_begin, chain.graph_firstLevel, -20, &vNW, false);
			std::vector<localExtension_pathDescription> forwardExtensions;// = fullNeedleman_affine_diagonal_extension(sequence, chain.sequence_end, chain.graph_lastLevel, -20, &vNW, true);

			std::cout << "\tBackward extensions: " << backwardExtensions.size() << "\n";
			std::cout << "\tForward extensions: " << forwardExtensions.size() << "\n";

			for(unsigned int backwardI = 0; backwardI < backwardExtensions.size(); backwardI++)
			{
				for(unsigned int forwardI = 0; forwardI < forwardExtensions.size(); forwardI++)
				{
					std::cout << "\t\t" << backwardI << "/" << forwardI << "\n";
					localExtension_pathDescription& bwE = backwardExtensions.at(backwardI);
					localExtension_pathDescription& fwE = forwardExtensions.at(forwardI);

					std::string fullAlignedSequence = bwE.alignedSequence + reconstructed_sequence + fwE.alignedSequence;
					std::string fulllAlignedGraph = bwE.alignedGraph + reconstructed_graph + fwE.alignedGraph;

					std::cout << "\t\t\t" << fulllAlignedGraph << " [graph]" << "\n";
					std::cout << "\t\t\t" << fullAlignedSequence << " [sequence]" << "\n" << std::flush;

				}
			}

		}

		// extendChain(seq_chains.at(chainI));

	}
}


void GraphAligner_nonAffine::fullNeedleman_diagonal(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned)
{
	std::map<int, std::map<int, std::map<int, double>> > m;
	std::map<int, std::map<int, std::map<int, backtraceStep>> > m_backtrace;

	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(0).size();
	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		m[0][0][stateI] = 0;
		m_backtrace[0][0][stateI] = backtraceStep();
		std::vector<int> existingCoordinates;
		existingCoordinates.push_back(0);
		existingCoordinates.push_back(0);
		existingCoordinates.push_back(stateI);
		m1_diagonal.push_back(existingCoordinates);
	}


	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	int diagonals = sequenceLength + levels - 1;
	int max_levelI = levels - 1;
	int max_seqI = sequenceLength;

	for(int diagonalI = 1; diagonalI <= diagonals; diagonalI++)
	{
		std::map<int, std::map<int, std::map<int, std::vector<double> > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, std::vector<backtraceStep> > > > thisDiagonal_backtrace;

		// extend from m-2 diagonal
		for(int m2I = 0; m2I < (int)m2_diagonal.size(); m2I++)
		{
			std::vector<int>& previous_coordinates = m2_diagonal.at(m2I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			int next_levelI = previous_coordinates.at(0) + 1;
			int next_seqI = previous_coordinates.at(1) + 1;

			if((next_levelI > max_levelI) || (next_seqI > max_seqI))
				continue;

			std::string sequenceEmission = sequence.substr(previous_seqI, 1);

			std::vector<std::pair<int, Edge*> > nextZs = _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI);
			assert(nextZs.size() > 0);

			for(unsigned int zI = 0; zI < nextZs.size(); zI++)
			{
				std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
				int next_stateI = thisZjump.first;

				std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
				assert(edgeEmission.size() == 1);

				double score_MatchMismatch = m.at(previous_levelI).at(previous_seqI).at(previous_stateI) + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;

				thisDiagonal[next_levelI][next_seqI][next_stateI].push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].push_back(backtrack_MatchMismatch);
			}
		}

		// extend from m-1 diagonal
		for(int m1I = 0; m1I < (int)m1_diagonal.size(); m1I++)
		{
			std::vector<int>& previous_coordinates = m1_diagonal.at(m1I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			// gap in graph
			int gapInGraph_next_levelI = previous_coordinates.at(0);
			int gapInGraph_next_seqI = previous_coordinates.at(1) + 1;
			if((gapInGraph_next_levelI <= max_levelI) && (gapInGraph_next_seqI <= max_seqI))
			{
				double score_gapInGraph = m.at(previous_levelI).at(previous_seqI).at(previous_stateI) + S_gap;

				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep backtrack_gapInGraph;
				backtrack_gapInGraph.x = previous_levelI;
				backtrack_gapInGraph.y = previous_seqI;
				backtrack_gapInGraph.z = previous_stateI;
				backtrack_gapInGraph.usedEdge = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].push_back(score_gapInGraph);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].push_back(backtrack_gapInGraph);
			}

			// gap in sequence
			int gapInSequence_next_levelI = previous_coordinates.at(0) + 1;
			int gapInSequence_next_seqI = previous_coordinates.at(1);

			if((gapInSequence_next_levelI <= max_levelI) && (gapInSequence_next_seqI <= max_seqI))
			{
				std::vector<std::pair<int, Edge*> > nextZs = _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI);
				assert(nextZs.size() > 0);

				for(unsigned int zI = 0; zI < nextZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
					int next_stateI = thisZjump.first;

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					double score_gapInSequence = m.at(previous_levelI).at(previous_seqI).at(previous_stateI) + ((edgeEmission == "_") ? S_graphGap : S_gap);

					backtraceStep backtrack_gapInSequence;
					backtrack_gapInSequence.x = previous_levelI;
					backtrack_gapInSequence.y = previous_seqI;
					backtrack_gapInSequence.z = previous_stateI;
					backtrack_gapInSequence.usedEdge = thisZjump.second;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].push_back(score_gapInSequence);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].push_back(backtrack_gapInSequence);
				}
			}
		}

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, std::vector<double> > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, std::vector<double> > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, std::vector<double> >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;
					std::vector<double>& scores = diagIt2->second.at(stateI);
					std::vector<backtraceStep>& scores_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI);

					std::pair<double, unsigned int> maxScore = Utilities::findVectorMax(scores);

					// condition for keeping this maximum
					if(true)
					{
						m[levelI][seqI][stateI] = maxScore.first;
						m_backtrace[levelI][seqI][stateI] = scores_bt.at(maxScore.second);

						std::vector<int> thisCoordinates;
						thisCoordinates.push_back(levelI);
						thisCoordinates.push_back(seqI);
						thisCoordinates.push_back(stateI);

						m_thisDiagonal.push_back(thisCoordinates);
					}
				}
			}
		}

		m2_diagonal = m1_diagonal;
		m1_diagonal = m_thisDiagonal;
	}

	if(verbose)
	{
		std::cout << "MATRIX:\n";
		for(unsigned int levelI = 0; levelI < levels; levelI++)
		{
			for(unsigned int seqI = 0; seqI <= sequenceLength; seqI++)
			{
				unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();
				std::vector<std::string> stateScores;
				for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
				{
					if((m.count(levelI) > 0) && (m.at(levelI).count(seqI) > 0) && (m.at(levelI).at(seqI).count(stateI) > 0))
					{
						int d_x = levelI - m_backtrace.at(levelI).at(seqI).at(stateI).x;
						int d_y = seqI - m_backtrace.at(levelI).at(seqI).at(stateI).y;
						int b_z = m_backtrace.at(levelI).at(seqI).at(stateI).z;

						std::string stateString = Utilities::DtoStr(m.at(levelI).at(seqI).at(stateI)) + "[" +
								Utilities::ItoStr(d_x) + "," +
								Utilities::ItoStr(d_y) + "," +
								Utilities::ItoStr(b_z) + "]";

						stateScores.push_back(stateString);
					}
					else
					{
						stateScores.push_back("-");
					}
				}
				std::string statesScores = Utilities::join(stateScores, "; ");
				std::cout << statesScores << "\t";
			}
			std::cout << "\n";
		}
		std::cout << "\n";

		// backtrace
		std::cout << "Levels: " << levels << "\n";
		std::cout << "sequenceLength: " << sequenceLength << "\n";
	}

	std::pair<double, int> maxFinalState = Utilities::findIntMapMax(m.at(levels-1).at(sequenceLength));

	if(verbose)
		std::cout << "Maximum final score " << maxFinalState.first << " in state " << maxFinalState.second << "\n";

	int backtrace_x = levels - 1;
	int backtrace_y = sequenceLength;
	int backtrace_z = maxFinalState.second;

	std::string reconstructed_graph;
	std::string reconstructed_sequence;
	std::vector<int> reconstructed_graph_levels;

	// print matrix

	while((backtrace_x != 0) || (backtrace_y != 0))
	{
		if(verbose)
		{
			std::cout << "backtrace_x: " << backtrace_x << "\n";
			std::cout << "backtrace_y: " << backtrace_y << "\n";
			std::cout << "backtrace_z: " << backtrace_z << "\n";
			std::cout << "Score: " << m.at(backtrace_x).at(backtrace_y).at(backtrace_z) << "\n\n";
		}

		assert((backtrace_x >= 0) && (backtrace_x <= (levels - 1)));
		assert((backtrace_y >= 0) && (backtrace_y <= (int)sequenceLength));
		assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));

		backtraceStep& step = m_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z);

		std::string edgeEmission;
		if(step.usedEdge != 0)
		{
			edgeEmission = g->CODE.deCode(step.usedEdge->locus_id, step.usedEdge->emission);
			assert(edgeEmission.size() == 1);
		}
		std::string sequenceEmission;
		if(backtrace_y >= 1)
		{
			sequenceEmission = sequence.substr(backtrace_y - 1, 1);
		}

		int next_x = step.x;
		int next_y = step.y;
		int next_z = step.z;

		if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
		{
			// match or mismatch
			assert(step.usedEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(backtrace_x - 1);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
		{
			// gap in graph
			reconstructed_graph.append("_");
			reconstructed_graph_levels.push_back(-1);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
		{
			// gap in sequence
			assert(step.usedEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(backtrace_x - 1);
			reconstructed_sequence.append("_");
		}
		else
		{
			assert(1 == 0);
		}

		backtrace_x = next_x;
		backtrace_y = next_y;
		backtrace_z = next_z;
	}

	std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
	std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
	std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());


	if(verbose)
		std::cout << "\n" << reconstructed_graph << "\n" << reconstructed_sequence << "\n\n";

	score = maxFinalState.first;
	graph_aligned = reconstructed_graph;
	graph_aligned_levels = reconstructed_graph_levels;
	sequence_aligned = reconstructed_sequence;
}




void GraphAligner_nonAffine::findShortGappedGraphConnection_nonAffine(int start_x, int start_z, int stop_x, int stop_z, std::vector<Edge*>& usedEdges, std::string& graphSequence, std::vector<int>& graphSequence_levels)
{
	assert(S_graphGap == 0);
	assert(start_x > stop_x);
	assert(stop_x >= 0);
	assert(stop_x < (int)g->NodesPerLevel.size());

	assert(start_z >= 0);
	assert((int)nodesPerLevel_ordered.size() > start_x);
	assert(start_z < (int)nodesPerLevel_ordered.at(start_x).size());

	class runningDistance_bookkeeper {
	public:
		int Distance;
		std::vector<Edge*> usedEdges;
		std::string graphSequence;
		std::vector<int> graphSequence_levels;
	};
	std::map<Node*, runningDistance_bookkeeper> runningDistance_fromStart;

	std::map<NWEdge*, std::map<NWEdge*, int> > NWedge_distances;

	for( int lI = start_x; lI >= stop_x; lI--)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, runningDistance_bookkeeper> runningNodeDistances_thisLevel;

		if(lI == start_x)
		{
			for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
			{
				Node* n = *nIt;
				int z = nodesPerLevel_ordered_rev.at(lI).at(n);
				if(z == start_z)
				{
					runningNodeDistances_thisLevel[n] = runningDistance_bookkeeper();
					runningNodeDistances_thisLevel[n].Distance = 0;
				}
			}
			assert(runningNodeDistances_thisLevel.size() > 0);
		}
		else
		{
			for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
			{
				Node* n = *nIt;
				assert((int)n->level == lI);
				std::set<Edge*> outgoingEdges = n->Outgoing_Edges;
				assert(outgoingEdges.size() > 0);

				for(std::set<Edge*>::iterator eIt = outgoingEdges.begin(); eIt != outgoingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					if(runningDistance_fromStart.count(e->To))
					{
						runningDistance_bookkeeper& edgeBound_distance = runningDistance_fromStart.at(e->To);

						int newDistance = edgeBound_distance.Distance + ((edgeEmission == "_") ? 0 : 1);

						if((runningNodeDistances_thisLevel.count(n) == 0) || (runningNodeDistances_thisLevel.at(n).Distance > newDistance))
						{
							runningDistance_bookkeeper newDistanceBookkeeper = edgeBound_distance;
							newDistanceBookkeeper.Distance = newDistance;
							newDistanceBookkeeper.graphSequence.append(edgeEmission);
							newDistanceBookkeeper.usedEdges.push_back(e);
							newDistanceBookkeeper.graphSequence_levels.push_back(n->level);

							runningNodeDistances_thisLevel[n] = newDistanceBookkeeper;
						}
					}
				}
			}
		}

		runningDistance_fromStart = runningNodeDistances_thisLevel;
	}

	runningDistance_bookkeeper selectedDistance;
	if(stop_z == -1)
	{
		double minimumDistance;
		Node* selectedMinimumDistanceNode;

		assert(nodesPerLevel_ordered.at(stop_x).size() > 0);
		for(int zI = 0; zI < (int)nodesPerLevel_ordered.at(stop_x).size(); zI++)
		{
			Node* stopNode = nodesPerLevel_ordered.at(stop_x).at(zI);
			assert(runningDistance_fromStart.count(stopNode));
			double D = runningDistance_fromStart.at(stopNode).Distance;
			if((zI == 0) || (D < minimumDistance))
			{
				minimumDistance = D;
				selectedMinimumDistanceNode = stopNode;
			}
		}

		selectedDistance = runningDistance_fromStart.at(selectedMinimumDistanceNode);
	}
	else
	{
		Node* stopNode = nodesPerLevel_ordered.at(stop_x).at(stop_z);
		assert(runningDistance_fromStart.count(stopNode));
		selectedDistance = runningDistance_fromStart.at(stopNode);
	}

	usedEdges = selectedDistance.usedEdges;
	graphSequence = selectedDistance.graphSequence;
	graphSequence_levels = selectedDistance.graphSequence_levels;

}

std::map<NWEdge*, std::map<NWEdge*, int> > GraphAligner_nonAffine::fastCalculateGraphDistancesForVirtualNW_nonAffine(VirtualNWTable& vNW, std::vector<std::map<NWEdge*, int> >& lastPositionDistances_perZ)
{
	assert(S_graphGap == 0);

	bool verbose = false;
	lastPositionDistances_perZ.clear();

	std::vector<NWEdge*> entryEdges_vector = vNW.getEntryEdges();
	std::vector<NWEdge*> exitEdges_vector = vNW.getExitEdges();

//	std::set<NWEdge*> entryEdges_set(entryEdges_vector.begin(), entryEdges_vector.end());
//	std::set<NWEdge*> exitEdges_set(exitEdges_vector.begin(), exitEdges_vector.end());

	std::map<int, std::map<int, std::set<NWEdge*> > >  entryEdges_from_coordinates;
	std::map<int, std::map<int, std::set<NWEdge*> > >  exitEdges_to_coordinates;

//	std::map<Edge*, std::set<NWEdge*> > graphEdges_2_NW_entry;
	for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
	{
		NWEdge* edge = entryEdges_vector.at(eI);
		entryEdges_from_coordinates[edge->from_x][edge->from_z].insert(edge);
//		Edge* graphEdge = edge->usedGraphEdge;
//		graphEdges_2_NW_entry[graphEdge].insert(edge);
//		assert(exitEdges_set.count(edge) == 0);
	}

//	std::map<Edge*, std::set<NWEdge*> > graphEdges_2_NW_exit;
	for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
	{
		NWEdge* edge = exitEdges_vector.at(eI);
		exitEdges_to_coordinates[edge->to_x][edge->to_z].insert(edge);
//		Edge* graphEdge = edge->usedGraphEdge;
//		graphEdges_2_NW_exit[graphEdge].insert(edge);
	}

	if(verbose)
		std::cout << "fastCalculateGraphDistancesForVirtualNW_nonAffine(..) for " << entryEdges_vector.size() <<  " entry NW edges and " << exitEdges_vector.size() << " exit NW edges.\n" << std::flush;

	std::map<Node*, std::map<NWEdge*, int> > runningNodeDistances;
	std::set<Node*> nodes_l0 = g->NodesPerLevel.at(0);

	std::map<NWEdge*, std::map<NWEdge*, int> > NWedge_distances;
	std::map<NWEdge*, std::set<NWEdge*> > deletionPoints_for_edges;

	for(unsigned int lI = 0; lI < g->NodesPerLevel.size(); lI++)
	{
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<NWEdge*, int> > runningNodeDistances_thisLevel;

		if(lI == (g->NodesPerLevel.size() - 1))
		{
			lastPositionDistances_perZ.resize(nodes_thisLevel.size());
		}


		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;
			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			std::map<NWEdge*, std::set<NWEdge*>> thisNode_distancesFromNWEdges;

			if(lI == 0)
			{
				runningNodeDistances_thisLevel[n] = std::map<NWEdge*, int>();
				runningNodeDistances_thisLevel[n][0] = 0;
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);
				runningNodeDistances_thisLevel[n] = std::map<NWEdge*, int>();

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;
					assert(runningNodeDistances.count(nodeFrom) > 0);


					// if the edge we are just traversing to get into node n is one of our NWEdge exit edges,
					// and if we have any NWEdges marked for deletion at this NWEdge exit edge, we will
					// get rid of it now

//					std::set<NWEdge*> nwEdges_to_ignore;
//					if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
//					{
//						// exit edge from this level identified - start logging
//
//						std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
//						for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
//						{
//							NWEdge* nwE = *nwEdgeIt;
//							if(deletionPoints_for_edges.count(nwE))
//							{
//								nwEdges_to_ignore.insert(deletionPoints_for_edges.at(nwE).begin(), deletionPoints_for_edges.at(nwE).end());
//							}
//						}
//					}

					for(std::map<NWEdge*, int>::iterator edgeDistIt = runningNodeDistances.at(nodeFrom).begin(); edgeDistIt != runningNodeDistances.at(nodeFrom).end(); edgeDistIt++)
					{
						// for all logged edges, increase distances as necessary

						NWEdge* interestingEdge = edgeDistIt->first;
//						if(! nwEdges_to_ignore.count(interestingEdge))
//						{
							int newDistance;
							if(edgeEmission == "_")
							{
								newDistance = edgeDistIt->second;
							}
							else
							{
								newDistance = edgeDistIt->second + 1;
							}

							if((runningNodeDistances_thisLevel.at(n).count(interestingEdge) == 0) || (runningNodeDistances_thisLevel.at(n).at(interestingEdge) > newDistance))
							{
								runningNodeDistances_thisLevel.at(n)[interestingEdge] = newDistance;
								thisNode_distancesFromNWEdges[interestingEdge] = exitEdges_to_coordinates.at(lI).at(z);
							}
//						}
					}

					if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
					{
						// exit edge from this level identified - start logging

						std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
						for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
						{
							NWEdge* e = *nwEdgeIt;
							runningNodeDistances_thisLevel.at(n)[e] = 0;
						}
					}
				}
			}

			for(std::map<NWEdge*, std::set<NWEdge*>>::iterator edgeDistIt = thisNode_distancesFromNWEdges.begin(); edgeDistIt != thisNode_distancesFromNWEdges.end(); edgeDistIt++)
			{
				NWEdge* thisEdgeWithDistance = edgeDistIt->first;
				std::set<NWEdge*> potentiallyTraversedExitEdes = edgeDistIt->second;
				bool oneTraversedExitEdgeMarkedAsDeletor = false;
				bool allTraversedExitEdgeMarkedAsDeletor = true;
				assert(potentiallyTraversedExitEdes.size() > 0);
				for(std::set<NWEdge*>::iterator traversedEdgeIt = potentiallyTraversedExitEdes.begin(); traversedEdgeIt != potentiallyTraversedExitEdes.end(); traversedEdgeIt++)
				{
					NWEdge* potentiallyTraversedExitEdge = *traversedEdgeIt;
					if(deletionPoints_for_edges.count(potentiallyTraversedExitEdge) && deletionPoints_for_edges.at(potentiallyTraversedExitEdge).count(thisEdgeWithDistance))
					{
						oneTraversedExitEdgeMarkedAsDeletor = true;
					}
					else
					{
						allTraversedExitEdgeMarkedAsDeletor = false;
					}
				}
				assert((oneTraversedExitEdgeMarkedAsDeletor && allTraversedExitEdgeMarkedAsDeletor) || ((! oneTraversedExitEdgeMarkedAsDeletor) && (! allTraversedExitEdgeMarkedAsDeletor)));
				if(oneTraversedExitEdgeMarkedAsDeletor)
				{
					runningNodeDistances_thisLevel.at(n).erase(thisEdgeWithDistance);
				}
			}


			if(lI != (g->NodesPerLevel.size() -1 ))
			{
				// non-last level

				if(entryEdges_from_coordinates.count(lI) && entryEdges_from_coordinates.at(lI).count(z))
				{
					// new entry edge identified

					std::set<NWEdge*> nwEdges = entryEdges_from_coordinates.at(lI).at(z);
					for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
					{
						NWEdge* e = *nwEdgeIt;

						// look at all existing exit edges...
						for(std::map<NWEdge*, int>::iterator exitEdgeIt = runningNodeDistances_thisLevel.at(n).begin(); exitEdgeIt != runningNodeDistances_thisLevel.at(n).end(); exitEdgeIt++)
						{
							NWEdge* exitEdge = exitEdgeIt->first;
							assert(e != 0);

							// check whether the entry edge is y-compatible with this exit edge
							if((exitEdge == 0) || (e->from_y >= exitEdge->to_y))
							{
								int distance = exitEdgeIt->second;
								assert(distance >= 0);

								// log edge distance
								NWedge_distances[e][exitEdge] = distance;

								// we will now go to the *exit edges* of the path that this *entry* edge belongs to
								// if we get there, we can delete the exit edges that we have just logged
								NWPath* path_for_entryEdge = e->path;
								std::set<NWEdge*> exitEdges_for_path = path_for_entryEdge->exit_edges;
								for(std::set<NWEdge*>::iterator furtherExitEdgeIt = exitEdges_for_path.begin(); furtherExitEdgeIt != exitEdges_for_path.end(); furtherExitEdgeIt++)
								{
									NWEdge* furtherExitEdge = *furtherExitEdgeIt;
									deletionPoints_for_edges[furtherExitEdge].insert(exitEdge);
								}
							}
						}
					}
				}
			}
			else
			{
				// log edge distances to the different nodes in the last graph

				for(std::map<NWEdge*, int>::iterator edgeDistIt = runningNodeDistances_thisLevel.at(n).begin(); edgeDistIt != runningNodeDistances_thisLevel.at(n).end(); edgeDistIt++)
				{
					NWEdge* interestingEdge = edgeDistIt->first;
					int thisDistance = edgeDistIt->second;

					if((lastPositionDistances_perZ.at(z).count(interestingEdge) == 0) || (lastPositionDistances_perZ.at(z).at(interestingEdge) > thisDistance))
					{
						lastPositionDistances_perZ.at(z)[interestingEdge] = thisDistance;
					}
				}
			}
		}

		runningNodeDistances = runningNodeDistances_thisLevel;
	}

	size_t calculated_distances = 0;
	for(std::map<NWEdge*, std::map<NWEdge*, int> >::iterator edgeIt = NWedge_distances.begin(); edgeIt != NWedge_distances.end(); edgeIt++)
	{
		for(std::map<NWEdge*, int>::iterator edgeIt2 = edgeIt->second.begin(); edgeIt2 != edgeIt->second.end(); edgeIt2++)
		{
			calculated_distances++;
		}
	}

	std::cout << "fastCalculateGraphDistancesForVirtualNW_nonAffine(..): Calculated " << calculated_distances << " distances.\n" << std::flush;

	return NWedge_distances;
}
