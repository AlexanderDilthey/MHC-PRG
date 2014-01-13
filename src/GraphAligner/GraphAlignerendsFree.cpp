/*
 * GraphAlignerendsFree.cpp
 *
 *  Created on: 27.07.2013
 *      Author: AlexanderDilthey
 */

#include "GraphAlignerendsFree.h"

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


GraphAligner_endsFree::GraphAligner_endsFree(Graph* graph, int k) : GraphAligner(graph, k) {
	// TODO Auto-generated constructor stub
}


std::vector<localExtension_pathDescription> GraphAligner_endsFree::fullNeedleman_endsFree_diagonal_extension(std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int diagonal_stop_threshold, VirtualNWTable* blockedPathsTable, bool directionPositive,  bool returnGlobalScore)
{
	// this is not really "ends-free"!

	double minusInfinity = -1 * numeric_limits<double>::max();

	if(returnGlobalScore)
	{
		diagonal_stop_threshold = minusInfinity;
	}

	class mScore {
	public:
		double D;
		double GraphGap;
		double SequenceGap;
	};

	class mScore_backtrace {
	public:
		backtraceStep_affine D;
		backtraceStep_affine GraphGap;
		backtraceStep_affine SequenceGap;
	};

	class mScore_alternatives {
	public:
		std::vector<double> D;
		std::vector<double> GraphGap;
		std::vector<double> SequenceGap;
	};

	class mScore_backtrace_alternatives {
	public:
		std::vector<backtraceStep_affine> D;
		std::vector<backtraceStep_affine> GraphGap;
		std::vector<backtraceStep_affine> SequenceGap;
	};

	std::map<int, std::map<int, std::map<int, mScore>> > scores;
	std::map<int, std::map<int, std::map<int, mScore_backtrace>> > scores_backtrace;

	std::vector<std::vector<int> > coordinates_for_backtracking;
	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	if(verbose)
	{
		std::cout << "fullNeedleman_affine_diagonal_extension(..) called.\n" << std::flush;
	}

	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	int diagonals = sequenceLength + levels - 1;
	int max_levelI = levels - 1;
	int max_seqI = sequenceLength;

	assert(startLevel_graph >= 0);
	assert(start_sequence >= 0);
	assert(startLevel_graph <= max_levelI);
	assert(start_sequence <= max_seqI);

	double currentMaximum = 0;
	std::vector<std::vector<int> > currentMaxima_coordinates;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(startLevel_graph).size();
	assert((startZ_graph >= 0) && (startZ_graph < (int)statesPerLevel0));

	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		if((int)stateI == startZ_graph)
		{
			scores[startLevel_graph][start_sequence][stateI].D = 0;
			scores_backtrace[startLevel_graph][start_sequence][stateI] = mScore_backtrace();
		}
		else
		{
			scores[startLevel_graph][start_sequence][stateI].D = minusInfinity;
		}

		scores[startLevel_graph][start_sequence][stateI].GraphGap = minusInfinity;
		scores[startLevel_graph][start_sequence][stateI].SequenceGap = minusInfinity;

		if((int)stateI == startZ_graph)
		{
			std::vector<int> existingCoordinates;
			existingCoordinates.push_back(startLevel_graph);
			existingCoordinates.push_back(start_sequence);
			existingCoordinates.push_back(stateI);
			m1_diagonal.push_back(existingCoordinates);
			currentMaxima_coordinates.push_back(existingCoordinates);
		}
	}

	std::map<NWPath*, std::pair<double, std::vector<int> > > hit_NW_paths;

	for(int diagonalI = 1; diagonalI <= diagonals; diagonalI++)
	{
		std::map<int, std::map<int, std::map<int, mScore_alternatives > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, mScore_backtrace_alternatives > > > thisDiagonal_backtrace;

		// extend from m-2 diagonal
		for(int m2I = 0; m2I < (int)m2_diagonal.size(); m2I++)
		{
			std::vector<int>& previous_coordinates = m2_diagonal.at(m2I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			int next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);

			if((next_levelI > max_levelI) || (next_seqI > max_seqI))
				continue;

			if((next_levelI < 0) || (next_seqI < 0))
				continue;

			std::string sequenceEmission = (directionPositive ? sequence.substr(previous_seqI, 1) : sequence.substr(previous_seqI-1, 1));

			std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
			assert(nextZs.size() > 0);

			for(unsigned int zI = 0; zI < nextZs.size(); zI++)
			{
				std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
				int next_stateI = thisZjump.first;

				std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
				assert(edgeEmission.size() == 1);

				double score_MatchMismatch = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep_affine backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;
				backtrack_MatchMismatch.sourceMatrix = 0;

				thisDiagonal[next_levelI][next_seqI][next_stateI].D.push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].D.push_back(backtrack_MatchMismatch);
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
			int gapInGraph_next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);
			if	(
					(directionPositive && (gapInGraph_next_levelI <= max_levelI) && (gapInGraph_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInGraph_next_levelI >= 0) && (gapInGraph_next_seqI >= 0))
				)
			{
				double score_gapInGraph_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_open = minusInfinity;
				}
				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep_affine backtrack_gapInGraph_open;
				backtrack_gapInGraph_open.x = previous_levelI;
				backtrack_gapInGraph_open.y = previous_seqI;
				backtrack_gapInGraph_open.z = previous_stateI;
				backtrack_gapInGraph_open.usedEdge = 0;
				backtrack_gapInGraph_open.sourceMatrix = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_open);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_open);

				double score_gapInGraph_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).GraphGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_extend = minusInfinity;
				}

				backtraceStep_affine backtrack_gapInGraph_extend;
				backtrack_gapInGraph_extend.x = previous_levelI;
				backtrack_gapInGraph_extend.y = previous_seqI;
				backtrack_gapInGraph_extend.z = previous_stateI;
				backtrack_gapInGraph_extend.usedEdge = 0;
				backtrack_gapInGraph_extend.sourceMatrix = 1;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_extend);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_extend);
			}

			// gap in sequence
			int gapInSequence_next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int gapInSequence_next_seqI = previous_coordinates.at(1);

			if	(
					(directionPositive && (gapInSequence_next_levelI <= max_levelI) && (gapInSequence_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInSequence_next_levelI >= 0) && (gapInSequence_next_seqI >= 0))
				)
			{
				std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
				assert(nextZs.size() > 0);

				for(unsigned int zI = 0; zI < nextZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
					int next_stateI = thisZjump.first;

					std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
					assert(edgeEmission.size() == 1);

					// open sequence gap

					double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
					if(edgeEmission == "_")
					{
						score_gapInSequence_open = minusInfinity;
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_open = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_open;
					backtrack_gapInSequence_open.x = previous_levelI;
					backtrack_gapInSequence_open.y = previous_seqI;
					backtrack_gapInSequence_open.z = previous_stateI;
					backtrack_gapInSequence_open.usedEdge = thisZjump.second;
					backtrack_gapInSequence_open.sourceMatrix = 0;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);

					// extend sequence gap

					double score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap + S_extendGap;
//					if(edgeEmission == "_")
//					{
//						score_gapInSequence_extend = minusInfinity;
//					}
					if(edgeEmission == "_")
					{
						if(scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap == minusInfinity)
						{
							score_gapInSequence_extend = minusInfinity;
						}
						else
						{
							score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap +  S_graphGap;
						}
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_extend = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_extend;
					backtrack_gapInSequence_extend.x = previous_levelI;
					backtrack_gapInSequence_extend.y = previous_seqI;
					backtrack_gapInSequence_extend.z = previous_stateI;
					backtrack_gapInSequence_extend.usedEdge = thisZjump.second;
					backtrack_gapInSequence_extend.sourceMatrix = 2;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_extend);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_extend);

					// non-affine sequence gap
					if(edgeEmission == "_")
					{
						double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_graphGap;

						backtraceStep_affine backtrack_gapInSequence_nonAffine;
						backtrack_gapInSequence_nonAffine.x = previous_levelI;
						backtrack_gapInSequence_nonAffine.y = previous_seqI;
						backtrack_gapInSequence_nonAffine.z = previous_stateI;
						backtrack_gapInSequence_nonAffine.usedEdge = thisZjump.second;
						backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(score_gapInSequence_nonAffine);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(backtrack_gapInSequence_nonAffine);
					}
				}
			}
		}

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, mScore_alternatives > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, mScore_alternatives > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, mScore_alternatives >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;

					// call maxima for GraphGap
					std::vector<double>& scores_GraphGap = diagIt2->second.at(stateI).GraphGap;
					std::vector<backtraceStep_affine>& scores_GraphGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).GraphGap;
					double selectedScore_GraphGap;
					backtraceStep_affine selectedStep_GraphGap;

					if(scores_GraphGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_GraphGap = Utilities::findVectorMax(scores_GraphGap);
						selectedScore_GraphGap = maxScore_GraphGap.first;
						selectedStep_GraphGap = scores_GraphGap_bt.at(maxScore_GraphGap.second);
					}
					else
					{
						selectedScore_GraphGap = minusInfinity;
					}

					// call maxima for SequenceGap
					std::vector<double>& scores_SequenceGap = diagIt2->second.at(stateI).SequenceGap;
					std::vector<backtraceStep_affine>& scores_SequenceGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).SequenceGap;

					double selectedScore_SequenceGap;
					backtraceStep_affine selectedStep_SequenceGap;
					if(scores_SequenceGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_SequenceGap = Utilities::findVectorMax(scores_SequenceGap);
						selectedScore_SequenceGap = maxScore_SequenceGap.first;
						selectedStep_SequenceGap = scores_SequenceGap_bt.at(maxScore_SequenceGap.second);
					}
					else
					{
						selectedScore_SequenceGap = minusInfinity;
					}

					// two additional steps for D, jumping from the two gap matrices

					double score_intoD_fromGraphGap = selectedScore_GraphGap;
					backtraceStep_affine step_intoD_fromGraphGap;
					step_intoD_fromGraphGap.x = levelI;
					step_intoD_fromGraphGap.y = seqI;
					step_intoD_fromGraphGap.z = stateI;
					step_intoD_fromGraphGap.sourceMatrix = 1;
					step_intoD_fromGraphGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromGraphGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromGraphGap);

					double score_intoD_fromSequenceGap = selectedScore_SequenceGap;
					backtraceStep_affine step_intoD_fromSequenceGap;
					step_intoD_fromSequenceGap.x = levelI;
					step_intoD_fromSequenceGap.y = seqI;
					step_intoD_fromSequenceGap.z = stateI;
					step_intoD_fromSequenceGap.sourceMatrix = 2;
					step_intoD_fromSequenceGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromSequenceGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromSequenceGap);

					// find final maximum for D
					std::vector<double>& scores_D = diagIt2->second.at(stateI).D;
					std::vector<backtraceStep_affine>& scores_D_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D;
					assert(scores_D.size() == scores_D_bt.size());

					std::pair<double, unsigned int> maxScore_D = Utilities::findVectorMax(scores_D);

					bool blockOutCell = false;
					if(directionPositive)
					{
						std::set<NWEdge*> emanatingNWEdges = blockedPathsTable->getEdgesEmanatingFrom(levelI, seqI, stateI);
						if((emanatingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = emanatingNWEdges.begin(); edgeIt != emanatingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}
					else
					{
						std::set<NWEdge*> incomingNWEdges = blockedPathsTable->getEdgesGoingInto(levelI, seqI, stateI);
						if((incomingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = incomingNWEdges.begin(); edgeIt != incomingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}

					if(blockOutCell)
					{
						double scoreThatWillBeDeleted = maxScore_D.first;
						std::vector<int> thisCoordinates;
						thisCoordinates.push_back(levelI);
						thisCoordinates.push_back(seqI);
						thisCoordinates.push_back(stateI);
						if(scoreThatWillBeDeleted == currentMaximum)
						{
							currentMaxima_coordinates.push_back(thisCoordinates);
						}
						else if(scoreThatWillBeDeleted > currentMaximum)
						{
							currentMaximum = scoreThatWillBeDeleted;
							currentMaxima_coordinates.clear();
							currentMaxima_coordinates.push_back(thisCoordinates);
						}

						scores[levelI][seqI][stateI].D = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

						assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

						scores[levelI][seqI][stateI].GraphGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

						scores[levelI][seqI][stateI].SequenceGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;
					}
					else
					{
						if(maxScore_D.first >= diagonal_stop_threshold)
						{
							scores[levelI][seqI][stateI].D = maxScore_D.first;
							scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

							assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

							scores[levelI][seqI][stateI].GraphGap = selectedScore_GraphGap;
							scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

							scores[levelI][seqI][stateI].SequenceGap = selectedScore_SequenceGap;
							scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;

							std::vector<int> thisCoordinates;
							thisCoordinates.push_back(levelI);
							thisCoordinates.push_back(seqI);
							thisCoordinates.push_back(stateI);
							m_thisDiagonal.push_back(thisCoordinates);

							backtraceStep_affine oneRealStepBackwards = scores_backtrace[levelI][seqI][stateI].D;
							while((oneRealStepBackwards.x == levelI) && (oneRealStepBackwards.y == seqI))
							{
								assert(oneRealStepBackwards.sourceMatrix != 0);
								if(oneRealStepBackwards.sourceMatrix == 1)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
								}
								else if(oneRealStepBackwards.sourceMatrix == 2)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
								}
								else
								{
									assert(1 == 0);
								}
							}
							int previousScore;
							if(oneRealStepBackwards.sourceMatrix == 0)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].D;
							}
							else if(oneRealStepBackwards.sourceMatrix == 1)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
							}
							else if(oneRealStepBackwards.sourceMatrix == 2)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
							}
							else
							{
								assert( 1 == 0 );
							}
							int scoreDifference = maxScore_D.first - previousScore;

							if(maxScore_D.first == currentMaximum)
							{
								if(scoreDifference != 0)
								{
									currentMaxima_coordinates.push_back(thisCoordinates);
								}
							}
							else if(maxScore_D.first > currentMaximum)
							{
								currentMaximum = maxScore_D.first;
								currentMaxima_coordinates.clear();
								currentMaxima_coordinates.push_back(thisCoordinates);
							}
						}
						else
						{
							// assert(1 == 0);
						}
					}
				}
			}
		}

		m2_diagonal = m1_diagonal;
		m1_diagonal = m_thisDiagonal;
	}

	std::vector<localExtension_pathDescription> forReturn;
	auto backtraceFrom = [&](int start_x, int start_y, int start_z, double StartScore) {
		int backtrace_x = start_x;
		int backtrace_y = start_y;
		int backtrace_z = start_z;
		int backtrace_matrix = 0;

		std::string reconstructed_graph;
		std::string reconstructed_sequence;
		std::vector<int> reconstructed_graph_levels;
		std::vector<std::vector<int>> edge_coordinates;
		std::vector<Edge*> utilizedEdges;

		std::vector<int> startCoordinates;
		startCoordinates.push_back(start_x);
		startCoordinates.push_back(start_y);
		startCoordinates.push_back(start_z);
		edge_coordinates.push_back(startCoordinates);

		if(verbose)
		{
			std::cout << "\tbacktraceFrom() called.\n";

		}
		while((backtrace_x != startLevel_graph) || (backtrace_y != start_sequence))
		{
			if(verbose)
			{
				std::cout << "\t\tbacktrace_x: " << backtrace_x << "\n";
				std::cout << "\t\tbacktrace_y: " << backtrace_y << "\n";
				std::cout << "\t\tbacktrace_z: " << backtrace_z << "\n";
				std::cout << "\t\tselected matrix: " << backtrace_matrix << "\n" << std::flush;
			}

			if(directionPositive)
			{
				assert((backtrace_x >= startLevel_graph) && (backtrace_x <= (levels - 1)));
				assert((backtrace_y >= start_sequence) && (backtrace_y <= (int)sequenceLength));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}
			else
			{
				assert((backtrace_x <= startLevel_graph) && (backtrace_x >= 0));
				assert((backtrace_y <= start_sequence) && (backtrace_y >= 0));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}

			double Score;
			backtraceStep_affine step;
			if(backtrace_matrix == 0)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
			}
			else if(backtrace_matrix == 1)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
			}
			else if(backtrace_matrix == 2)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
			}

			if((backtrace_x == start_x) && (backtrace_y == start_y) && (backtrace_z == start_z))
			{
				Score = StartScore;
			}

			if(verbose)
			{
				std::cout << "\t\tScore: " << Score << "\n\n" << std::flush;
			}

			std::string edgeEmission;
			if(step.usedEdge != 0)
			{
				edgeEmission = g->CODE.deCode(step.usedEdge->locus_id, step.usedEdge->emission);
				assert(edgeEmission.size() == 1);
			}
			std::string sequenceEmission;
			if((backtrace_y >= 1) && (directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y - 1, 1);
			}
			if((backtrace_y < max_seqI) && (! directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y, 1);
			}

			int next_x = step.x;
			int next_y = step.y;
			int next_z = step.z;
			int next_matrix = step.sourceMatrix;

			bool dontAddCoordinates = false; // if we jump from one matrix to the other without changing coordinates, we don't store the coordinates...
			if(directionPositive)
			{
				if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);

				}
				else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}
			else
			{
				if((next_x == (backtrace_x + 1)) && (next_y == (backtrace_y + 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y + 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);
				}
				else if((next_x == (backtrace_x + 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}

			if(verbose)
			{
				std::cout << "\t\t\t" << "next_x: " << next_x << "\n";
				std::cout << "\t\t\t" << "next_y: " << next_y << "\n";
				std::cout << "\t\t\t" << "next_z: " << next_z << "\n";
				std::cout << "\t\t\t" << "next_matrix: " << next_matrix << "\n\n" << std::flush;
			}

			backtrace_x = next_x;
			backtrace_y = next_y;
			backtrace_z = next_z;
			backtrace_matrix = next_matrix;

			if( ! dontAddCoordinates)
			{
				std::vector<int> nextCoordinates;
				nextCoordinates.push_back(next_x);
				nextCoordinates.push_back(next_y);
				nextCoordinates.push_back(next_z);
				edge_coordinates.push_back(nextCoordinates);
			}
		}

		if(directionPositive)
		{
			std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
			std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
			std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());
			std::reverse(utilizedEdges.begin(), utilizedEdges.end());
			std::reverse(edge_coordinates.begin(), edge_coordinates.end());
		}

		assert(reconstructed_graph.size() == reconstructed_sequence.size());
		assert(reconstructed_graph_levels.size() == reconstructed_graph.size());

		localExtension_pathDescription pathReturn;
		pathReturn.Score = StartScore;
		pathReturn.usedEdges = utilizedEdges;
		pathReturn.coordinates = edge_coordinates;
		pathReturn.alignedSequence = reconstructed_sequence;
		pathReturn.alignedGraph = reconstructed_graph;
		pathReturn.alignedGraph_levels = reconstructed_graph_levels;

		forReturn.push_back(pathReturn);
	};

	if(! returnGlobalScore)
	{
		if(currentMaximum > 0)
		{
			if(verbose)
			{
				std::cout << "\tMaximum " << currentMaximum << ", achieved at " << currentMaxima_coordinates.size() << " positions!\n" << std::flush;
			}

			for(int maximumI = 0; maximumI < (int)currentMaxima_coordinates.size(); maximumI++)
			{
				std::vector<int> coordinates = currentMaxima_coordinates.at(maximumI);
				if(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D != minusInfinity)
				{
					if(verbose)
						std::cout << " - Start maximum backtrace from " <<coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

					backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D);
				}
			}
		}

		if(verbose)
			std::cout << "Have hit " << hit_NW_paths.size() << " NW paths, backtrace independent of achieved value!\n" << std::flush;

		for(std::map<NWPath*, std::pair<double, std::vector<int> > >::iterator hitPathsIt = hit_NW_paths.begin(); hitPathsIt != hit_NW_paths.end(); hitPathsIt++)
		{
			// NWPath* path = hitPathsIt->first;
			double score = hitPathsIt->second.first;
			std::vector<int> coordinates = hitPathsIt->second.second;
			assert(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D == minusInfinity);

			if(verbose)
				std::cout << " - Start NW path backtrace from " << coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

			backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), score);
		}
	}
	else
	{

		if(directionPositive)
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(levels-1).at(sequenceLength).begin(); altIt != scores.at(levels-1).at(sequenceLength).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMax(finalState_scores);

			backtraceFrom(levels - 1, sequenceLength, maxFinalState.second, maxFinalState.first);
		}
		else
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(0).at(0).begin(); altIt != scores.at(0).at(0).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;

				if(verbose)
					std::cout << "Final state z-value " << altIt->first << ", D value " << scores.at(0).at(0).at(altIt->first).D << ", GraphGap: " <<  scores.at(0).at(0).at(altIt->first).GraphGap << ", SequenceGap: " <<  scores.at(0).at(0).at(altIt->first).SequenceGap << "\n" << std::flush;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMax(finalState_scores);

			backtraceFrom(0, 0, maxFinalState.second, maxFinalState.first);
		}
	}
	return forReturn;
}



void GraphAligner_endsFree::fullNeedleman_endsFree(std::string sequence, int& score, std::string& graph_aligned, std::vector<int>& graph_aligned_levels, std::string& sequence_aligned)
{

	double minusInfinity = -1 * numeric_limits<double>::max();
	assert(S_graphGap == 0);

	class mScore {
	public:
		double D;
		double GraphGap;
		double SequenceGap;
	};

	class mScore_backtrace {
	public:
		backtraceStep_affine D;
		backtraceStep_affine GraphGap;
		backtraceStep_affine SequenceGap;
	};

	class mScore_alternatives {
	public:
		std::vector<double> D;
		std::vector<double> GraphGap;
		std::vector<double> SequenceGap;
	};

	class mScore_backtrace_alternatives {
	public:
		std::vector<backtraceStep_affine> D;
		std::vector<backtraceStep_affine> GraphGap;
		std::vector<backtraceStep_affine> SequenceGap;
	};

	std::map<int, std::map<int, std::map<int, mScore>> > scores;
	std::map<int, std::map<int, std::map<int, mScore_backtrace>> > scores_backtrace;

	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(0).size();
	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		scores[0][0][stateI].D = 0;
		scores_backtrace[0][0][stateI] = mScore_backtrace();

		scores[0][0][stateI].GraphGap = minusInfinity;
		scores[0][0][stateI].SequenceGap = minusInfinity;

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
		std::map<int, std::map<int, std::map<int, mScore_alternatives > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, mScore_backtrace_alternatives > > > thisDiagonal_backtrace;

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

				double score_MatchMismatch = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep_affine backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;
				backtrack_MatchMismatch.sourceMatrix = 0;

				thisDiagonal[next_levelI][next_seqI][next_stateI].D.push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].D.push_back(backtrack_MatchMismatch);
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
				double score_gapInGraph_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;

				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep_affine backtrack_gapInGraph_open;
				backtrack_gapInGraph_open.x = previous_levelI;
				backtrack_gapInGraph_open.y = previous_seqI;
				backtrack_gapInGraph_open.z = previous_stateI;
				backtrack_gapInGraph_open.usedEdge = 0;
				backtrack_gapInGraph_open.sourceMatrix = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_open);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_open);

				double score_gapInGraph_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).GraphGap + S_extendGap;

				backtraceStep_affine backtrack_gapInGraph_extend;
				backtrack_gapInGraph_extend.x = previous_levelI;
				backtrack_gapInGraph_extend.y = previous_seqI;
				backtrack_gapInGraph_extend.z = previous_stateI;
				backtrack_gapInGraph_extend.usedEdge = 0;
				backtrack_gapInGraph_extend.sourceMatrix = 1;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_extend);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_extend);
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

					if((gapInSequence_next_seqI == 0) || (gapInSequence_next_seqI == max_seqI))
					{
						double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + 0;
						backtraceStep_affine backtrack_gapInSequence_open;
						backtrack_gapInSequence_open.x = previous_levelI;
						backtrack_gapInSequence_open.y = previous_seqI;
						backtrack_gapInSequence_open.z = previous_stateI;
						backtrack_gapInSequence_open.usedEdge = thisZjump.second;
						backtrack_gapInSequence_open.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);
					}
					else
					{
						std::string edgeEmission = g->CODE.deCode(thisZjump.second->locus_id, thisZjump.second->emission);
						assert(edgeEmission.size() == 1);

						// open sequence gap

						double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
						if(edgeEmission == "_")
						{
							score_gapInSequence_open = minusInfinity;
						}

						backtraceStep_affine backtrack_gapInSequence_open;
						backtrack_gapInSequence_open.x = previous_levelI;
						backtrack_gapInSequence_open.y = previous_seqI;
						backtrack_gapInSequence_open.z = previous_stateI;
						backtrack_gapInSequence_open.usedEdge = thisZjump.second;
						backtrack_gapInSequence_open.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);

						// extend sequence gap

						double score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap + S_extendGap;
						if(edgeEmission == "_")
						{
							if(scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap == minusInfinity)
							{
								score_gapInSequence_extend = minusInfinity;
							}
							else
							{
								score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap +  S_graphGap;
							}
						}

						backtraceStep_affine backtrack_gapInSequence_extend;
						backtrack_gapInSequence_extend.x = previous_levelI;
						backtrack_gapInSequence_extend.y = previous_seqI;
						backtrack_gapInSequence_extend.z = previous_stateI;
						backtrack_gapInSequence_extend.usedEdge = thisZjump.second;
						backtrack_gapInSequence_extend.sourceMatrix = 2;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_extend);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_extend);

						// non-affine sequence gap
						if(edgeEmission == "_")
						{
							double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_graphGap;

							backtraceStep_affine backtrack_gapInSequence_nonAffine;
							backtrack_gapInSequence_nonAffine.x = previous_levelI;
							backtrack_gapInSequence_nonAffine.y = previous_seqI;
							backtrack_gapInSequence_nonAffine.z = previous_stateI;
							backtrack_gapInSequence_nonAffine.usedEdge = thisZjump.second;
							backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

							thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(score_gapInSequence_nonAffine);
							thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(backtrack_gapInSequence_nonAffine);
						}
					}
				}
			}
		}

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, mScore_alternatives > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, mScore_alternatives > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, mScore_alternatives >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;

					// call maxima for GraphGap
					std::vector<double>& scores_GraphGap = diagIt2->second.at(stateI).GraphGap;
					std::vector<backtraceStep_affine>& scores_GraphGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).GraphGap;
					double selectedScore_GraphGap;
					backtraceStep_affine selectedStep_GraphGap;

					if(scores_GraphGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_GraphGap = Utilities::findVectorMax(scores_GraphGap);
						selectedScore_GraphGap = maxScore_GraphGap.first;
						selectedStep_GraphGap = scores_GraphGap_bt.at(maxScore_GraphGap.second);
					}
					else
					{
						selectedScore_GraphGap = minusInfinity;
					}

					// call maxima for SequenceGap
					std::vector<double>& scores_SequenceGap = diagIt2->second.at(stateI).SequenceGap;
					std::vector<backtraceStep_affine>& scores_SequenceGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).SequenceGap;

					double selectedScore_SequenceGap;
					backtraceStep_affine selectedStep_SequenceGap;
					if(scores_SequenceGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_SequenceGap = Utilities::findVectorMax(scores_SequenceGap);
						selectedScore_SequenceGap = maxScore_SequenceGap.first;
						selectedStep_SequenceGap = scores_SequenceGap_bt.at(maxScore_SequenceGap.second);
					}
					else
					{
						selectedScore_SequenceGap = minusInfinity;
					}

					// two additional steps for D, jumping from the two gap matrices

					double score_intoD_fromGraphGap = selectedScore_GraphGap;
					backtraceStep_affine step_intoD_fromGraphGap;
					step_intoD_fromGraphGap.x = levelI;
					step_intoD_fromGraphGap.y = seqI;
					step_intoD_fromGraphGap.z = stateI;
					step_intoD_fromGraphGap.sourceMatrix = 1;
					step_intoD_fromGraphGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromGraphGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromGraphGap);

					double score_intoD_fromSequenceGap = selectedScore_SequenceGap;
					backtraceStep_affine step_intoD_fromSequenceGap;
					step_intoD_fromSequenceGap.x = levelI;
					step_intoD_fromSequenceGap.y = seqI;
					step_intoD_fromSequenceGap.z = stateI;
					step_intoD_fromSequenceGap.sourceMatrix = 2;
					step_intoD_fromSequenceGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromSequenceGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromSequenceGap);

					// find final maximum for D
					std::vector<double>& scores_D = diagIt2->second.at(stateI).D;
					std::vector<backtraceStep_affine>& scores_D_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D;
					assert(scores_D.size() == scores_D_bt.size());

					std::pair<double, unsigned int> maxScore_D = Utilities::findVectorMax(scores_D);

					// condition for keeping this maximum
					if(true)
					{
						scores[levelI][seqI][stateI].D = maxScore_D.first;
						scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);
						if(scores_backtrace[levelI][seqI][stateI].D.x == -1)
						{
							std::cout << std::flush;
							std::cerr << "levelI: " << levelI << "\n";
							std::cerr << "seqI: " << seqI << "\n";
							std::cerr << "stateI: " << stateI << "\n";
							std::cerr << "selected score: " << maxScore_D.first << "\n";
							std::cerr << "selected index: " << maxScore_D.second << "\n\n";
							std::cerr << "coming from:\n";
							std::cerr << "\t" << "sourceMatrix: " << scores_backtrace[levelI][seqI][stateI].D.sourceMatrix << "\n";
							std::cerr << "\t" << "x: " << scores_backtrace[levelI][seqI][stateI].D.x << "\n";
							std::cerr << "\t" << "y: " << scores_backtrace[levelI][seqI][stateI].D.y << "\n";
							std::cerr << "\t" << "z: " << scores_backtrace[levelI][seqI][stateI].D.z << "\n";

							std::cerr << std::flush;
						}
						assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

						scores[levelI][seqI][stateI].GraphGap = selectedScore_GraphGap;
						scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

						scores[levelI][seqI][stateI].SequenceGap = selectedScore_SequenceGap;
						scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;

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

//	if(verbose)
//	{
//		std::cout << "MATRIX:\n";
//		for(unsigned int levelI = 0; levelI < levels; levelI++)
//		{
//			for(unsigned int seqI = 0; seqI <= sequenceLength; seqI++)
//			{
//				unsigned int statesPerLevel = g->NodesPerLevel.at(levelI).size();
//				std::vector<std::string> stateScores;
//				for(unsigned int stateI = 0; stateI < statesPerLevel; stateI++)
//				{
//					if((m.count(levelI) > 0) && (m.at(levelI).count(seqI) > 0) && (m.at(levelI).at(seqI).count(stateI) > 0))
//					{
//						int d_x = levelI - m_backtrace.at(levelI).at(seqI).at(stateI).x;
//						int d_y = seqI - m_backtrace.at(levelI).at(seqI).at(stateI).y;
//						int b_z = m_backtrace.at(levelI).at(seqI).at(stateI).z;
//
//						std::string stateString = Utilities::DtoStr(m.at(levelI).at(seqI).at(stateI)) + "[" +
//								Utilities::ItoStr(d_x) + "," +
//								Utilities::ItoStr(d_y) + "," +
//								Utilities::ItoStr(b_z) + "]";
//
//						stateScores.push_back(stateString);
//					}
//					else
//					{
//						stateScores.push_back("-");
//					}
//				}
//				std::string statesScores = Utilities::join(stateScores, "; ");
//				std::cout << statesScores << "\t";
//			}
//			std::cout << "\n";
//		}
//		std::cout << "\n";
//
//		// backtrace
//		std::cout << "Levels: " << levels << "\n";
//		std::cout << "sequenceLength: " << sequenceLength << "\n";
//	}

	std::map<int, double> finalState_scores;
	for(std::map<int, mScore>::iterator altIt = scores.at(levels-1).at(sequenceLength).begin(); altIt != scores.at(levels-1).at(sequenceLength).end(); altIt++)
	{
		finalState_scores[altIt->first] = altIt->second.D;
	}
	std::pair<double, int> maxFinalState = Utilities::findIntMapMax(finalState_scores);

	if(verbose)
		std::cout << "Maximum final score " << maxFinalState.first << " in state " << maxFinalState.second << "\n";

	int backtrace_x = levels - 1;
	int backtrace_y = sequenceLength;
	int backtrace_z = maxFinalState.second;
	int backtrace_matrix = 0;

	std::string reconstructed_graph;
	std::string reconstructed_sequence;
	std::vector<int> reconstructed_graph_levels;

	while((backtrace_x != 0) || (backtrace_y != 0))
	{
		if(verbose)
		{
			std::cout << "backtrace_x: " << backtrace_x << "\n";
			std::cout << "backtrace_y: " << backtrace_y << "\n";
			std::cout << "backtrace_z: " << backtrace_z << "\n";
			std::cout << "selected matrix: " << backtrace_matrix << "\n" << std::flush;
		}

		assert((backtrace_x >= 0) && (backtrace_x <= (levels - 1)));
		assert((backtrace_y >= 0) && (backtrace_y <= (int)sequenceLength));
		assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
		assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));

		double Score;
		backtraceStep_affine step;
		if(backtrace_matrix == 0)
		{
			Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
			step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
		}
		else if(backtrace_matrix == 1)
		{
			Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
			step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
		}
		else if(backtrace_matrix == 2)
		{
			Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
			step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
		}
		if(verbose)
		{
			std::cout << "Score: " << Score << "\n\n" << std::flush;
		}

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
		int next_matrix = step.sourceMatrix;

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
			// assert(1 == 0);
		}

		if(verbose)
		{
			std::cout << "\t" << "next_x: " << next_x << "\n";
			std::cout << "\t" << "next_y: " << next_y << "\n";
			std::cout << "\t" << "next_z: " << next_z << "\n";
			std::cout << "\t" << "next_matrix: " << next_matrix << "\n\n" << std::flush;
		}

		backtrace_x = next_x;
		backtrace_y = next_y;
		backtrace_z = next_z;
		backtrace_matrix = next_matrix;
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


seedAndExtend_return GraphAligner_endsFree::seedAndExtend4(std::string sequence, bool noExtension)
{
	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();

	bool verbose = false;
	bool superquiet = false;

	/*
	 * Step 1: we will identify matching exact chains, and translate those to paths in vNW.
	 */

	if(! superquiet)
		std::cout << Utilities::timestamp() << "GraphAligner_endsFree::seedAndExtend4(..): search chains.\n" << std::flush;
	std::vector<kMerChain> seq_chains = gI.findChains(sequence);

	if(! superquiet)
		std::cout << Utilities::timestamp() << "Done. Have " << seq_chains.size() << " chains.\n" << std::flush;

	assert(seq_chains.size() > 0);

	VirtualNWTable vNW(this, &sequence);
	VirtualNWTable vNW2(this, &sequence);

	std::vector<NWPath*> found_NWpaths;

	if(! superquiet)
		std::cout << Utilities::timestamp() << "Add found paths to vNW..\n" << std::flush;

	for(unsigned int chainI = 0; chainI < seq_chains.size(); chainI++)
	{
		kMerChain& thisChain = seq_chains.at(chainI);

		// std::cerr << "Call findVirtualNWExactSeeds from saE3...\n" << std::flush;
		std::vector<NWPath*> seed_NW_paths = findVirtualNWExactSeeds(sequence, thisChain.sequence_begin, thisChain.sequence_end, thisChain.graph_firstLevel, thisChain.graph_lastLevel);

		if((! superquiet) && ((chainI % 1000) == 0))
			std::cout << "\rChain " << chainI << " / " << seq_chains.size() << ", found seed paths: " << seed_NW_paths.size() << std::flush;

		found_NWpaths.insert(found_NWpaths.end(), seed_NW_paths.begin(), seed_NW_paths.end());

		for(unsigned int seedPathI = 0; seedPathI < seed_NW_paths.size(); seedPathI++)
		{
			NWPath* thisSeedPath = seed_NW_paths.at(seedPathI);

			if(verbose)
				std::cout << "\t\tSeed path " << seedPathI << "\n" << std::flush;

			if(noExtension)
			{
				vNW2.addPath(thisSeedPath);
			}
			else
			{
				vNW.addPath(thisSeedPath);
			}
		}
	}

	if(! superquiet)
			std::cout << "\n" << Utilities::timestamp() << "Extend and expect to see hit NW paths...\n" << std::flush;

	/*
	 *  Step 2: we extend exact matches and add the resulting paths to vNW2
	 */

	if(! noExtension)
	{
		for(unsigned int pathI = 0; pathI < found_NWpaths.size(); pathI++)
		{
			if(verbose)
				std::cout << "\tNWPath " << pathI << "\n";

			if((! superquiet) && ((pathI % 1000) == 0))
				std::cout << "\rPath " << pathI << " / " << found_NWpaths.size() << std::flush;

			NWPath* thisSeedPath = found_NWpaths.at(pathI);

			assert(thisSeedPath->first_edges.size() == 1);
			NWEdge* firstSeedEdge = *(thisSeedPath->first_edges.begin());

			assert(thisSeedPath->last_edges.size() == 1);
			NWEdge* lastSeedEdge = *(thisSeedPath->last_edges.begin());

			// -11 is enough for two mismatches, or three gaps
			std::vector<localExtension_pathDescription> backwardExtensions = fullNeedleman_endsFree_diagonal_extension(
					sequence,
					firstSeedEdge->from_y,
					firstSeedEdge->from_x,
					firstSeedEdge->from_z,
					-11,
					&vNW,
					false);

			std::vector<localExtension_pathDescription> forwardExtensions = fullNeedleman_endsFree_diagonal_extension(
					sequence,
					lastSeedEdge->to_y,
					lastSeedEdge->to_x,
					lastSeedEdge->to_z,
					-11,
					&vNW,
					true);

			NWPath* newPath = thisSeedPath->clonePathWithoutTable();

			if(verbose)
				std::cout << "Created newPath.\n" << std::flush;

//			if(backwardExtensions.size() > 0)
//			{
//				newPath->entry_edges.clear();
//				newPath->first_edges.clear();
//			}
//
//			if(forwardExtensions.size() > 0)
//			{
//				newPath->exit_edges.clear();
//				newPath->last_edges.clear();
//			}

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
	}

	vNW2.checkConsistency();

	if(! superquiet)
	{
			std::cout << "\n" << Utilities::timestamp() << "Extensions done.\n" << std::flush;
			std::cout << "\tAdded all resulting extensions paths. In total: " << vNW2.getNumPaths() << " paths and " << vNW2.getNumEdges() << " edges.\n" << std::flush;
			std::cout << "\tBefore additional exit points added: " << vNW2.getNumEntryEgdes() << " entry edges and " << vNW2.getNumExitEdges() << " exit edges.\n" << std::flush;
			std::cout << "Stats for vNW2:\n";
			vNW.printSequenceCoverageStats();
	}

	/*
	 * Step 2.1: We add additional exit points for edges that might allow for vertical jumps
	 */

	for(std::set<NWPath*>::iterator pathIt = vNW2.paths.begin(); pathIt != vNW2.paths.end(); pathIt++)
	{
		NWPath* path = *pathIt;
		for(std::set<NWEdge*>::iterator entryIt = path->entry_edges.begin(); entryIt != path->entry_edges.end(); entryIt++)
		{
			NWEdge* entryEdge = *entryIt;

			int entry_x = entryEdge->from_x;
			int entry_y = entryEdge->from_y;

			if(vNW2.index_edges_stop_x.count(entry_x))
			{
				std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths;

				std::set<NWEdge*> egdes_into_x = vNW2.index_edges_stop_x.at(entry_x);
				for(std::set<NWEdge*>::iterator sameXEdgeIt = egdes_into_x.begin(); sameXEdgeIt != egdes_into_x.end(); sameXEdgeIt++)
				{
					NWEdge* sameXEdge = *sameXEdgeIt;
					if(sameXEdge->to_y <= entry_y)
					{
						if(sameXEdge->path != entryEdge->path)
						{
							if((hitPaths.count(sameXEdge->path) == 0) || (hitPaths.at(sameXEdge->path).first < sameXEdge->to_y))
							{
								hitPaths[sameXEdge->path] = std::pair<int, NWEdge*>(sameXEdge->to_y, sameXEdge);
							}
						}
					}
				}

				for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths.begin(); makeEdgeExitIt != hitPaths.end(); makeEdgeExitIt++)
				{
					assert(makeEdgeExitIt->first != entryEdge->path);
					assert(makeEdgeExitIt->second.second != entryEdge);
					makeEdgeExitIt->second.second->makeExitEdge();
				}
			}

			if(vNW2.index_edges_stop_y.count(entry_y))
			{
				std::map<NWPath*, std::pair<int, NWEdge*> > hitPaths;

				std::set<NWEdge*> egdes_into_y = vNW2.index_edges_stop_y.at(entry_y);
				for(std::set<NWEdge*>::iterator sameYEdgeIt = egdes_into_y.begin(); sameYEdgeIt != egdes_into_y.end(); sameYEdgeIt++)
				{
					NWEdge* sameYEdge = *sameYEdgeIt;
					if(sameYEdge->to_x <= entry_x)
					{
						if(sameYEdge->path != entryEdge->path)
						{
							if((hitPaths.count(sameYEdge->path) == 0) || (hitPaths.at(sameYEdge->path).first < sameYEdge->to_x))
							{
								hitPaths[sameYEdge->path] = std::pair<int, NWEdge*>(sameYEdge->to_x, sameYEdge);
							}
						}
					}
				}

				for(std::map<NWPath*, std::pair<int, NWEdge*> >::iterator makeEdgeExitIt = hitPaths.begin(); makeEdgeExitIt != hitPaths.end(); makeEdgeExitIt++)
				{
					assert(makeEdgeExitIt->first != entryEdge->path);
					assert(makeEdgeExitIt->second.second != entryEdge);
					makeEdgeExitIt->second.second->makeExitEdge();
				}
			}
		}
	}

	vNW2.checkConsistency();

	/*
	 * Step 3.1: Sort entry edges. We will want to traverse the edges by column.
	 */

	if(! superquiet)
		std::cout << Utilities::timestamp() << "After additional exit points added: " << vNW2.getNumEntryEgdes() << " entry edges and " << vNW2.getNumExitEdges() << " exit edges.\n" << std::flush;

	std::vector<NWEdge*> allEntryEdges = vNW2.getEntryEdges();
	std::list<NWEdge*> allEntryEdges_forSorting(allEntryEdges.begin(), allEntryEdges.end());
	allEntryEdges_forSorting.sort([](NWEdge* i, NWEdge* j){
		if(i->from_y != j->from_y)
		{
			return (i->from_y < j->from_y);
		}
		else
		{
			return (i->from_x < j->from_x);
		}
	});


	std::vector<NWEdge*> allEntryEdges_sorted(allEntryEdges_forSorting.begin(), allEntryEdges_forSorting.end());

	vNW2.checkConsistency();

	if(! superquiet)
		std::cout  << Utilities::timestamp() << "Sorting done." << "\n" << std::flush;


	/*
	 * Step 3.2: Start assinging scores to entry edges
	 */


	assert(S_graphGap == 0);

	bool debug = false;

	std::vector<double> finalScore_alternatives;
	std::vector<NWEdge*> finalBacktrack_alternatives;
	std::vector<int> finalBacktrack_Zs;

	std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startAffineGap;
	std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> > NWedges_graphDist_startNormal;

	std::vector<NWEdge*> entryEdges_vector = vNW2.getEntryEdges();
	std::vector<NWEdge*> exitEdges_vector = vNW2.getExitEdges();

	std::map<int, std::map<int, std::set<NWEdge*> > >  entryEdges_from_coordinates;
	std::map<int, std::map<int, std::set<NWEdge*> > >  exitEdges_to_coordinates;

	for(unsigned int eI = 0; eI < entryEdges_vector.size(); eI++)
	{
		NWEdge* edge = entryEdges_vector.at(eI);
		entryEdges_from_coordinates[edge->from_x][edge->from_z].insert(edge);
	}

	for(unsigned int eI = 0; eI < exitEdges_vector.size(); eI++)
	{
		NWEdge* edge = exitEdges_vector.at(eI);
		exitEdges_to_coordinates[edge->to_x][edge->to_z].insert(edge);
	}

	if(! superquiet)
		std::cout << "seedAndExtend4(..): Now fill NW table for " << entryEdges_vector.size() <<  " entry NW edges and " << exitEdges_vector.size() << " exit NW edges.\n" << std::flush;

	std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_notStartInAffineGap_endInAnything;
	std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_notStartInAffineGap_endInAffineGap;
	std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_startInAffineGap_endInAnything;
	std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_startInAffineGap_endInAffineGap;

	std::set<Node*> nodes_l0 = g->NodesPerLevel.at(0);

	for(unsigned int lI = 0; lI < g->NodesPerLevel.size(); lI++)
	{
		if(debug) std::cout << "lI: " << lI << "\n" << std::flush;

		if(! superquiet)
		{
			if((lI % 1000) == 0)
			{
				std::cout << "\r\t" << "Level: " << lI << "   " << std::flush;
			}
		}
		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;
		std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		std::map<Node*, sAE3_nodePositionStorage > runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

		std::map<Node*, std::set<NWEdge*> > runningNodeDistances_notStartInAffineGap_endInAnything_thisLevelAssigned;
		std::map<Node*, std::set<NWEdge*> > runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevelAssigned;
		std::map<Node*, std::set<NWEdge*> > runningNodeDistances_startInAffineGap_endInAnything_thisLevelAssigned;
		std::map<Node*, std::set<NWEdge*> > runningNodeDistances_startInAffineGap_endInAffineGap_thisLevelAssigned;

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(debug) std::cout << "Node: " << n << "\n" << std::flush;

			// int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == 0)
			{
				if(debug) std::cout << "a1" << "\n" << std::flush;


				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n].addNewEdge(0, 0);
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n].addNewEdge(0, minusInfinity);

				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n].addNewEdge(0, 0);
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n].addNewEdge(0, minusInfinity);

				if(debug) std::cout << "a2" << "\n" << std::flush;
			}
			else
			{
				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevelAssigned[n] = std::set<NWEdge*>();
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevelAssigned[n] = std::set<NWEdge*>();
				runningNodeDistances_startInAffineGap_endInAnything_thisLevelAssigned[n] = std::set<NWEdge*>();
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevelAssigned[n] = std::set<NWEdge*>();

				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				if(incomingEdges.size() == 1)
				{
					Node* fromNode = (*(incomingEdges.begin()))->From;
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = runningNodeDistances_notStartInAffineGap_endInAnything.at(fromNode);
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = runningNodeDistances_notStartInAffineGap_endInAffineGap.at(fromNode);
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = runningNodeDistances_startInAffineGap_endInAnything.at(fromNode);
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = runningNodeDistances_startInAffineGap_endInAffineGap.at(fromNode);
				}
				else
				{

					if(debug) std::cerr << "A" << "\n" << std::flush;

					std::set<sAE3_nodePositionStorage*> existingSources_notStartInAffineGap_endInAnything_set;
					std::vector<sAE3_nodePositionStorage*> existingSources_notStartInAffineGap_endInAffineGap;
					std::vector<sAE3_nodePositionStorage*> existingSources_startInAffineGap_endInAnything;
					std::vector<sAE3_nodePositionStorage*> existingSources_startInAffineGap_endInAffineGap;

					if(debug) std::cerr << "B" << "\n" << std::flush;

					for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
					{
						Edge* e = *eIt;
						Node* eFrom = e->From;
						sAE3_nodePositionStorage* nodeDistances_notStartInAffineGap_endInAnything = &(runningNodeDistances_notStartInAffineGap_endInAnything.at(eFrom));
						if(existingSources_notStartInAffineGap_endInAnything_set.count(nodeDistances_notStartInAffineGap_endInAnything) == 0)
						{
							existingSources_notStartInAffineGap_endInAnything_set.insert(nodeDistances_notStartInAffineGap_endInAnything);

							{
								sAE3_nodePositionStorage* nodeDistances_notStartInAffineGap_endInAffineGap = &(runningNodeDistances_notStartInAffineGap_endInAffineGap.at(eFrom));
								existingSources_notStartInAffineGap_endInAffineGap.push_back(nodeDistances_notStartInAffineGap_endInAffineGap);
							}
							{
								sAE3_nodePositionStorage* nodeDistances_startInAffineGap_endInAnything = &(runningNodeDistances_startInAffineGap_endInAnything.at(eFrom));
								existingSources_startInAffineGap_endInAnything.push_back(nodeDistances_startInAffineGap_endInAnything);
							}
							{
								sAE3_nodePositionStorage* nodeDistances_startInAffineGap_endInAffineGap = &(runningNodeDistances_startInAffineGap_endInAffineGap.at(eFrom));
								existingSources_startInAffineGap_endInAffineGap.push_back(nodeDistances_startInAffineGap_endInAffineGap);
							}
						}
					}

					if(debug) std::cerr << "C" << "\n" << std::flush;

					std::vector<sAE3_nodePositionStorage*> existingSources_notStartInAffineGap_endInAnything(existingSources_notStartInAffineGap_endInAnything_set.begin(), existingSources_notStartInAffineGap_endInAnything_set.end());
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n].initFromExistingSources(existingSources_notStartInAffineGap_endInAnything);
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n].initFromExistingSources(existingSources_notStartInAffineGap_endInAffineGap);
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n].initFromExistingSources(existingSources_startInAffineGap_endInAnything);
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n].initFromExistingSources(existingSources_startInAffineGap_endInAffineGap);
				}

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;
					std::vector<NWEdge*> defined_NW_edges = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).getEdges();

					for(std::vector<NWEdge*>::iterator nwEdgeIt = defined_NW_edges.begin(); nwEdgeIt != defined_NW_edges.end(); nwEdgeIt++)
					{
						NWEdge* interestingEdge = *nwEdgeIt;

						double previousDistance_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).getEdgeValue(interestingEdge);
						double previousDistance_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap.at(nodeFrom).getEdgeValue(interestingEdge);
						double previousDistance_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything.at(nodeFrom).getEdgeValue(interestingEdge);
						double previousDistance_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap.at(nodeFrom).getEdgeValue(interestingEdge);

						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								double newDistance_notStartInAffineGap_affine_extend;
								double newDistance_notStartInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start = minusInfinity;
									// newDistance_notStartInAffineGap_affine_extend = minusInfinity;
									if(previousDistance_notStartInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_notStartInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_graphGap;
									}

								}
								else
								{
									newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap + S_extendGap;
									newDistance_notStartInAffineGap_affine_start = previousDistance_notStartInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								if(interestingEdge == 0)
								{
									newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap;
									newDistance_notStartInAffineGap_affine_start = previousDistance_notStartInAffineGap_endInAnything;
//									assert(newDistance_notStartInAffineGap_affine_extend == 0);
//									assert(newDistance_notStartInAffineGap_affine_start == 0);
								}

								double newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start > newDistance_notStartInAffineGap_affine_extend) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if((runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevelAssigned.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(interestingEdge) < newDistance_notStartInAffineGap_affine_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).updateEdgeValue(interestingEdge, newDistance_notStartInAffineGap_affine_maximum);
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevelAssigned.at(n).insert(interestingEdge);
								}
							}

							// for the ones that did start with an affine gap
							{
								double newDistance_startInAffineGap_affine_extend;
								double newDistance_startInAffineGap_affine_start;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start = minusInfinity;
									//newDistance_startInAffineGap_affine_extend = minusInfinity;
									if(previousDistance_startInAffineGap_endInAffineGap == minusInfinity)
									{
										newDistance_startInAffineGap_affine_extend = minusInfinity;
									}
									else
									{
										newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_graphGap;
									}
								}
								else
								{
									newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap + S_extendGap;
									newDistance_startInAffineGap_affine_start = previousDistance_startInAffineGap_endInAnything + S_openGap + S_extendGap;
								}

								if(interestingEdge == 0)
								{
									newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap;
									newDistance_startInAffineGap_affine_start = previousDistance_startInAffineGap_endInAnything;
//									assert(newDistance_startInAffineGap_affine_extend == 0);
//									assert(newDistance_startInAffineGap_affine_start == 0);
								}

								double newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start > newDistance_startInAffineGap_affine_extend) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if((runningNodeDistances_startInAffineGap_endInAffineGap_thisLevelAssigned.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(interestingEdge) < newDistance_startInAffineGap_affine_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).updateEdgeValue(interestingEdge, newDistance_startInAffineGap_affine_maximum);
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevelAssigned.at(n).insert(interestingEdge);
								}
							}
						}

						if(debug) std::cout << "d3" << "\n" << std::flush;

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								assert(runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevelAssigned.at(n).count(interestingEdge));
								double alternativeScore_endInAffine = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_notStartInAffineGap_endInAnything;
								}

								if(interestingEdge == 0)
								{
									followGraphGap = previousDistance_notStartInAffineGap_endInAnything;
								}

								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if(interestingEdge == 0)
								{
//									assert(endInAnything_maximum == 0);
								}
								if((runningNodeDistances_notStartInAffineGap_endInAnything_thisLevelAssigned.at(n).count(interestingEdge) == 0) || (runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).updateEdgeValue(interestingEdge, endInAnything_maximum);
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevelAssigned.at(n).insert(interestingEdge);
								}
							}

							// for the ones that did start with an affine gap
							{
								assert(runningNodeDistances_startInAffineGap_endInAffineGap_thisLevelAssigned.at(n).count(interestingEdge));
								double alternativeScore_endInAffine = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(interestingEdge);
								double followGraphGap = minusInfinity;
								if(edgeEmission == "_")
								{
									followGraphGap = previousDistance_startInAffineGap_endInAnything;
								}

								if(interestingEdge == 0)
								{
									followGraphGap = previousDistance_startInAffineGap_endInAnything;
								}

								double endInAnything_maximum = (followGraphGap > alternativeScore_endInAffine) ? followGraphGap : alternativeScore_endInAffine;
								if(interestingEdge == 0)
								{
//									assert(endInAnything_maximum == 0);
								}
								if((runningNodeDistances_startInAffineGap_endInAnything_thisLevelAssigned.at(n).count(interestingEdge) == 0) || (runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(interestingEdge) < endInAnything_maximum))
								{
									runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).updateEdgeValue(interestingEdge, endInAnything_maximum);
									runningNodeDistances_startInAffineGap_endInAnything_thisLevelAssigned.at(n).insert(interestingEdge);
								}
							}
						}
					}
				}

				if(runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).edgeScores.count((NWEdge*)0x80234a60))
				{
					if(lI < 480)
					{
						std::cerr << "Running score for 0x80234a60 at node " << n << " (lI " << lI << "): " << runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).edgeScores.at((NWEdge*)0x80234a60) << "\n" << std::flush;
					}
				}
			}
		}


		auto scoreEdgeJump = [&](NWEdge* entryEdge, NWEdge* exitEdge, graphPointDistance graphDistance_startAffineGap, graphPointDistance graphDistance_startNormal, bool& jumpComingFromSequenceGap) -> double {

			if(exitEdge == 0)
			{
				if(verbose)
					std::cout << "\t" << "THE COORDINATE ORIGIN" << "\n\n" << std::flush;
			}
			else
			{
				if(verbose)
				{
					std::cout << "\t" << "exitEdge->to_x: " << exitEdge->to_x << "\n";
					std::cout << "\t" << "exitEdge->to_y: " << exitEdge->to_y << "\n";
					std::cout << "\t" << "exitEdge->to_z: " << exitEdge->to_z << "\n\n" << std::flush;
				}
			}

			bool debug = false;
			if((void*)entryEdge == (void*)0x817fa5d8)
			{
				if(exitEdge != 0)
				{
					if((exitEdge->to_x == 54) && (exitEdge->to_y == 31))
					{
						debug = true;
					}

				}
			}

			bool exitEdge_isAffineSequenceGap = (exitEdge == 0) ? false : exitEdge->isSequenceGap_endsFree();
			bool exitEdge_isGraphGap = (exitEdge == 0) ? false : exitEdge->isGraphGap();
			bool entryEdge_isAffineSequenceGap = (entryEdge == 0) ? false : entryEdge->isSequenceGap_endsFree();
			bool entryEdge_isGraphGap = (entryEdge == 0) ? false : entryEdge->isGraphGap();

			int distance_along_y;
			int raw_distance_along_x;
			int exitEdge_to_y;
			int entryEdge_from_y;

			if((entryEdge == 0) && (exitEdge == 0))
			{
				distance_along_y = sequence.length();
				raw_distance_along_x = g->NodesPerLevel.size() - 1;
				exitEdge_to_y = 0;
				entryEdge_from_y = sequence.length();
			}
			else if(entryEdge == 0)
			{
				assert(exitEdge != 0);
				distance_along_y = sequence.length() - exitEdge->to_y;
				raw_distance_along_x = g->NodesPerLevel.size() - 1 - exitEdge->to_x;
				exitEdge_to_y = exitEdge->to_y;
				entryEdge_from_y = sequence.length();
			}
			else if(exitEdge == 0)
			{
				assert(entryEdge != 0);
				distance_along_y = entryEdge->from_y;
				raw_distance_along_x = entryEdge->from_x;
				exitEdge_to_y = 0;
				entryEdge_from_y = entryEdge->from_y;
			}
			else
			{
				assert(entryEdge != 0);
				assert(exitEdge != 0);
				distance_along_y = entryEdge->from_y - exitEdge->to_y;
				raw_distance_along_x = entryEdge->from_x - exitEdge->to_x;
				exitEdge_to_y = exitEdge->to_y;
				entryEdge_from_y = entryEdge->from_y;
			}

			assert(distance_along_y >= 0);
			assert(graphDistance_startAffineGap.Score_endInAnything >= graphDistance_startAffineGap.Score_endAffinely);
			assert(graphDistance_startAffineGap.Score_endInAnything > minusInfinity);
			assert(graphDistance_startNormal.Score_endInAnything >= graphDistance_startNormal.Score_endAffinely);
			assert(graphDistance_startNormal.Score_endInAnything > minusInfinity);

			double combinedScore;
			if(entryEdge == 0)
			{
				// we jump to the end - we only need to take into account the y-jump, the x-jump becomes irrelevant
				double gapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
				double scoreFromEntryEdge = 0;
				double utilizedDistance = (exitEdge != 0) ? exitEdge->getScoreAfterEdge() : 0;
				double utilizedYDistance = ((distance_along_y > 0) ? gapOpenPenalty : 0) + distance_along_y * S_extendGap;
				combinedScore = utilizedDistance + utilizedYDistance + scoreFromEntryEdge;
				jumpComingFromSequenceGap = false;
			}
			else
			{
				if((raw_distance_along_x != 0) && (distance_along_y != 0))
				{
					// we have two alternatives: either the sequence gap comes first and the graph gap follows, or the other way around
					std::vector<double> scores;
					std::vector<bool> scores_fromAffineSequenceGap;
					{
						// case 1: sequence gap comes first
						{
							// case 1.1: the sequence gap ends in anything

							graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);
							double sequenceGapScore = utilizedDistance.Score_endInAnything;
							double graphGapScore = S_openGap + distance_along_y * S_extendGap;
							double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
							double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;

							scores.push_back(local_combinedScore);
							scores_fromAffineSequenceGap.push_back(false);
						}
						{
							// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
						}
					}
					{
						// case 2: graph gap comes first
						{
							// case 2.1: sequence gap ends in anything
							double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
							double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
							double sequenceGapScore = graphDistance_startNormal.Score_endInAnything;
							double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
							double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
							if(exitEdge == 0)
							{
								local_combinedScore = minusInfinity;
							}
							scores.push_back(local_combinedScore);
							scores_fromAffineSequenceGap.push_back(false);
						}
						{
							// case 2.2: sequence gap ends affinely
							double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
							double graphGapScore = graphGapOpenPenalty + distance_along_y * S_extendGap;
							double sequenceGapScore = graphDistance_startNormal.Score_endAffinely;
							double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
							double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
							if(exitEdge == 0)
							{
								local_combinedScore = minusInfinity;
							}
							scores.push_back(local_combinedScore);
							scores_fromAffineSequenceGap.push_back(true);
						}
					}

					std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMax(scores);
					combinedScore = scoreAlternativesMax.first;
					jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);
				}
				else if((raw_distance_along_x == 0) && (distance_along_y != 0))
				{
					// graph gap
					double gapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);
					assert(graphDistance_startAffineGap.Score_endInAnything == graphDistance_startNormal.Score_endInAnything);
					combinedScore = utilizedDistance.Score_endInAnything + gapOpenPenalty + distance_along_y * S_extendGap + scoreFromEntryEdge;
					jumpComingFromSequenceGap = false;
				}
				else if((raw_distance_along_x != 0) && (distance_along_y == 0))
				{
					// sequence gap
					graphPointDistance& utilizedDistance = (exitEdge_isAffineSequenceGap ? graphDistance_startAffineGap : graphDistance_startNormal);

					std::vector<double> scores;
					std::vector<bool> scores_fromAffineSequenceGap;

					// case 1: we want to end in an affine gap
					{
						double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
						double scoreGraphDistance = utilizedDistance.Score_endAffinely;
						double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
						scores.push_back(local_combinedScore);
						scores_fromAffineSequenceGap.push_back(true);

						if((void*)entryEdge == (void*)0x80276690)
						{
							if((void*)exitEdge == (void*)0x80234a60)
							{
								std::cerr << "\nEntry edge: " << entryEdge << "\n";
								std::cerr << "Exit edge: " << exitEdge << "\n";
								std::cerr << "exitEdge_isAffineSequenceGap: " << exitEdge_isAffineSequenceGap << "\n";

								std::cerr << "scoreFromEntryEdge: " << scoreFromEntryEdge << "\n";
								std::cerr << "utilizedDistance.Score_endAffinely: " << utilizedDistance.Score_endAffinely << "\n";
								std::cerr << "local_combinedScore: " << local_combinedScore << "\n";
								std::cerr << "scores_fromAffineSequenceGap: true" << "\n" << std::flush;
							}
						}
					}

					// case 2: end in anything
					{
						double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
						double scoreGraphDistance = utilizedDistance.Score_endInAnything;
						double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
						scores.push_back(local_combinedScore);
						scores_fromAffineSequenceGap.push_back(false);

						if((void*)entryEdge == (void*)0x80276690)
						{
							if((void*)exitEdge == (void*)0x80234a60)
							{
								std::cerr << "\nEntry edge: " << entryEdge << "\n";
								std::cerr << "Exit edge: " << exitEdge << "\n";
								std::cerr << "exitEdge_isAffineSequenceGap: " << exitEdge_isAffineSequenceGap << "\n";
								std::cerr << "scoreFromEntryEdge: " << scoreFromEntryEdge << "\n";
								std::cerr << "utilizedDistance.Score_endInAnything: " << utilizedDistance.Score_endInAnything << "\n";
								std::cerr << "local_combinedScore: " << local_combinedScore << "\n";
								std::cerr << "scores_fromAffineSequenceGap: false" << "\n" << std::flush;
							}
						}
					}

					std::pair<double, unsigned int> scoreAlternativesMax = Utilities::findVectorMax(scores);
					combinedScore = scoreAlternativesMax.first;
					jumpComingFromSequenceGap = scores_fromAffineSequenceGap.at(scoreAlternativesMax.second);

				}
				else
				{
					assert((raw_distance_along_x == 0) && (distance_along_y == 0));
					assert(graphDistance_startAffineGap.Score_endInAnything == graphDistance_startNormal.Score_endInAnything);
					double scoreFromEntryEdge = (entryEdge != 0) ? entryEdge->calculateScore_endsFree(exitEdge_isAffineSequenceGap, exitEdge_isGraphGap, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
					combinedScore = graphDistance_startAffineGap.Score_endInAnything + scoreFromEntryEdge;
					jumpComingFromSequenceGap = (exitEdge != 0) ? exitEdge->isSequenceGap_endsFree() : false;
				}
			}

			if(verbose)
			{
				std::cout << "\t" << "graph jump scores:\n";
				std::cout << "\t\t" << " affine start, affine end: " << graphDistance_startAffineGap.Score_endAffinely << "\n";
				std::cout << "\t\t" << " non-affine start, affine end: " << graphDistance_startNormal.Score_endAffinely << "\n";
				std::cout << "\t\t" << " affine start, arbitrary end: " << graphDistance_startAffineGap.Score_endInAnything << "\n";
				std::cout << "\t\t" << " non-affine start, arbitrary end: " << graphDistance_startNormal.Score_endInAnything << "\n";
				std::cout << "\t" << "distance_along_y: " << distance_along_y << "\n" << std::flush;
			}

			if((void*)entryEdge == (void*)0x80276690)
			{
				if((void*)exitEdge == (void*)0x80234a60)
				{
					std::cerr << "Entry edge: " << entryEdge << "\n";
					std::cerr << "Exit edge: " << exitEdge << "\n";
					std::cerr << "\tto_x: " << exitEdge->to_x << "\n";
					std::cerr << "\tto_y: " << exitEdge->to_y << "\n";
					std::cerr << "\tto_z: " << exitEdge->to_z << "\n" << std::flush;
					std::cerr << "\texitEdge->getScoreAfterEdge(): " << exitEdge->getScoreAfterEdge() << "\n" << std::flush;
					std::cerr << "Chosen jump score: " << combinedScore << "\n\n" << std::flush;

//					std::cerr << "\t" << "graph jump scores:\n";
//					std::cerr << "\t\t" << " affine start, affine end: " << graphDistance_startAffineGap.Score_endAffinely << "\n";
//					std::cerr << "\t\t" << " non-affine start, affine end: " << graphDistance_startNormal.Score_endAffinely << "\n";
//					std::cerr << "\t\t" << " affine start, arbitrary end: " << graphDistance_startAffineGap.Score_endInAnything << "\n";
//					std::cerr << "\t\t" << " non-affine start, arbitrary end: " << graphDistance_startNormal.Score_endInAnything << "\n";
//					std::cerr << "\t" << "distance_along_y: " << distance_along_y << "\n" << std::flush;
				}
			}

			return combinedScore;
		};

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{
			Node* n = *nIt;

			if(debug) std::cout << "Node: " << n << "\n" << std::flush;

			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			class entryExitEdgeSpecifier {
			public:
				NWEdge* e;
				int type;
				entryExitEdgeSpecifier()
				{
					type = -1;
					e = 0;
				}
			};

			std::list<entryExitEdgeSpecifier> edges2process;

			if(exitEdges_to_coordinates.count(lI) && exitEdges_to_coordinates.at(lI).count(z))
			{
				std::set<NWEdge*> nwEdges = exitEdges_to_coordinates.at(lI).at(z);
				for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
				{
					NWEdge* exitEdge = *nwEdgeIt;
					entryExitEdgeSpecifier thisEdge;
					thisEdge.e = exitEdge;
					thisEdge.type = 1;
					edges2process.push_back(thisEdge);
				}
			}

			if(lI != (g->NodesPerLevel.size() - 1))
			{
				if(entryEdges_from_coordinates.count(lI) && entryEdges_from_coordinates.at(lI).count(z))
				{
					std::set<NWEdge*> nwEdges = entryEdges_from_coordinates.at(lI).at(z);
					for(std::set<NWEdge*>::iterator nwEdgeIt = nwEdges.begin(); nwEdgeIt != nwEdges.end(); nwEdgeIt++)
					{
						NWEdge* entryEdge = *nwEdgeIt;
						assert(entryEdge != 0);
						entryExitEdgeSpecifier thisEdge;
						thisEdge.e = entryEdge;
						thisEdge.type = 2;
						edges2process.push_back(thisEdge);
					}
				}
			}

			edges2process.sort([&](entryExitEdgeSpecifier a, entryExitEdgeSpecifier b) -> bool {
				int y_a;
				int y_b;
				if(a.type == 1)
				{
					if(a.e == 0)
					{
						y_a = 0;
					}
					else
					{
						y_a = a.e->to_y;
					}
				}
				else
				{
					assert(a.type == 2);
					assert(a.e != 0);
					y_a = a.e->from_y;
				}

				if(b.type == 1)
				{
					if(b.e == 0)
					{
						y_b = 0;
					}
					else
					{
						y_b = b.e->to_y;
					}
				}
				else
				{
					assert(b.type == 2);
					assert(b.e != 0);
					y_b = b.e->from_y;
				}

				if(y_a != y_b)
				{
					return (y_a < y_b);
				}
				else
				{
					if(a.type != b.type)
					{
						return (a.type < b.type);
					}
					else
					{
						return (a.e < b.e);
					}
				}
			});


			std::vector<entryExitEdgeSpecifier> edges2process_vector(edges2process.begin(), edges2process.end());

			for(unsigned int eI = 0; eI < edges2process_vector.size(); eI++)
			{
				entryExitEdgeSpecifier& thisES = edges2process_vector.at(eI);
				if(thisES.type == 1)
				{
					NWEdge* exitEdge = thisES.e;

					// std::cerr << "score1\n" << std::flush;
					double edgeScore = exitEdge->getScoreAfterEdge();
					// std::cerr << "score2\n" << std::flush;

					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).addNewEdge(exitEdge, edgeScore);
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).addNewEdge(exitEdge, minusInfinity);

					runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).addNewEdge(exitEdge, edgeScore);
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).addNewEdge(exitEdge, edgeScore);

					if(exitEdge == (NWEdge*)0x80234a60)
					{
						std::cerr << "Running score for 0x80234a60 at node " << exitEdge->usedGraphEdge->From << " (lI " << lI << "): " << runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).edgeScores.at((NWEdge*)0x80234a60) << "\n" << std::flush;
					}
				}
				else
				{
					assert(thisES.type == 2);
					NWEdge* entryEdge = thisES.e;

					std::vector<double> scores_achieved;
					std::vector<NWEdge*> scores_backtrack;
					std::vector<bool> scores_jumpFromAffineSequenceGap;

					std::vector<NWEdge*> potentialExitEdges = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).getEdges();
					omp_set_num_threads(threads);

					long long max_i = potentialExitEdges.size() - 1;
					long long chunk_size = max_i / threads;

					#pragma omp parallel
					{
						std::vector<NWEdge*> local_scores_backtrack;
						std::vector<double> local_scores_achieved;
						std::vector<bool> local_scores_jumpFromAffineSequenceGap;

						assert(omp_get_num_threads() == threads);
						long long thisThread = omp_get_thread_num();
						long long firstPair = thisThread * chunk_size;
						long long lastPair = (thisThread+1) * chunk_size - 1;
						if((thisThread == (threads-1)) && (lastPair < max_i))
						{
							lastPair = max_i;
						}

						for(long long i = firstPair; i <= lastPair; i++)
						{
							NWEdge* exitEdge = potentialExitEdges.at(i);

							assert(entryEdge != 0);
							if((exitEdge == 0) || (entryEdge->from_y >= exitEdge->to_y))
							{
								graphPointDistance distanceSpecifier_startAffineGap;
								graphPointDistance distanceSpecifier_startNormal;
								{

									distanceSpecifier_startAffineGap.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(exitEdge);
									distanceSpecifier_startAffineGap.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(exitEdge);
									distanceSpecifier_startAffineGap.start_in_affine_sequenceGap = true;
								}
								{
									distanceSpecifier_startNormal.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(exitEdge);
									distanceSpecifier_startNormal.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(exitEdge);
									distanceSpecifier_startNormal.start_in_affine_sequenceGap = false;
								}

								if(verbose)
									std::cout << "\n\t" << "-------\n\tConsider source " << exitEdge << "\n" << std::flush;

								bool jumpComingFromSequenceGap;
								double jumpScore = scoreEdgeJump(entryEdge, exitEdge, distanceSpecifier_startAffineGap, distanceSpecifier_startNormal, jumpComingFromSequenceGap);
								local_scores_backtrack.push_back(exitEdge);
								local_scores_achieved.push_back(jumpScore);
								local_scores_jumpFromAffineSequenceGap.push_back(jumpComingFromSequenceGap);
							}
						}

						#pragma omp critical
						{
							scores_achieved.insert(scores_achieved.end(), local_scores_achieved.begin(), local_scores_achieved.end());
							scores_backtrack.insert(scores_backtrack.end(), local_scores_backtrack.begin(), local_scores_backtrack.end());
							scores_jumpFromAffineSequenceGap.insert(scores_jumpFromAffineSequenceGap.end(), local_scores_jumpFromAffineSequenceGap.begin(), local_scores_jumpFromAffineSequenceGap.end());
						}
					}

					std::pair<double, unsigned int> maxOrigin = Utilities::findVectorMax(scores_achieved);
					assert(entryEdge != 0);
					entryEdge->takeScore_endsFree(maxOrigin.first, scores_jumpFromAffineSequenceGap.at(maxOrigin.second), scores_backtrack.at(maxOrigin.second), S_match, S_mismatch, S_openGap, S_extendGap);

//					if((void*)entryEdge == (void*)0x80276690)
//					{
//						std::cerr << "Entry edge: " << entryEdge << "\n";
//						std::cerr << "After jump: " << entryEdge->getScoreAfterEdge() << "\n" << std::flush;
//					}
				}
			}

			if(lI == (g->NodesPerLevel.size() - 1))
			{
				std::vector<NWEdge*> potentialExitEdges = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).getEdges();
				for(std::vector<NWEdge*>::iterator exitEdgeIt = potentialExitEdges.begin(); exitEdgeIt != potentialExitEdges.end(); exitEdgeIt++)
				{
					NWEdge* exitEdge = *exitEdgeIt;
					graphPointDistance distanceSpecifier_startInAffineGap;
					graphPointDistance distanceSpecifier_startNormal;
					{
						distanceSpecifier_startInAffineGap.Score_endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(exitEdge);
						distanceSpecifier_startInAffineGap.Score_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(exitEdge);
						distanceSpecifier_startInAffineGap.start_in_affine_sequenceGap = true;

					}

					{
						distanceSpecifier_startNormal.Score_endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).getEdgeValue(exitEdge);
						distanceSpecifier_startNormal.Score_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).getEdgeValue(exitEdge);
						distanceSpecifier_startNormal.start_in_affine_sequenceGap = false;
					}

					bool jumpComingFromSequenceGap;
					double jumpScore = scoreEdgeJump(0, exitEdge, distanceSpecifier_startInAffineGap, distanceSpecifier_startNormal, jumpComingFromSequenceGap);
					finalScore_alternatives.push_back(jumpScore);
					finalBacktrack_alternatives.push_back(exitEdge);
					finalBacktrack_Zs.push_back(z);
				}
			}
		}

		runningNodeDistances_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		runningNodeDistances_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;
	}

	std::pair<double, unsigned int> maxFinalScore = Utilities::findVectorMax(finalScore_alternatives);
	double finalScore = maxFinalScore.first;
	int finalScore_z = finalBacktrack_Zs.at(maxFinalScore.second);
	NWEdge* finalScore_backtrack = finalBacktrack_alternatives.at(maxFinalScore.second);

	if(! superquiet)
	{
		std::cout << "\n";
	}

	size_t calculated_distances = 0;
	for(std::map<NWEdge*, std::map<NWEdge*, graphPointDistance> >::iterator edgeIt = NWedges_graphDist_startNormal.begin(); edgeIt != NWedges_graphDist_startNormal.end(); edgeIt++)
	{
		for(std::map<NWEdge*, graphPointDistance>::iterator edgeIt2 = edgeIt->second.begin(); edgeIt2 != edgeIt->second.end(); edgeIt2++)
		{
			calculated_distances += 2;
		}
	}

	std::string reconstructedSequence;
	std::string reconstructedGraph;
	std::vector<int> reconstructedGraph_levels;


	if(! superquiet)
		std::cout  << Utilities::timestamp() << "Start backtracking...\n" << std::flush;

	seedAndExtend4_backtrack(
				vNW2,
				sequence,
				finalScore,
				finalScore_z,
				finalScore_backtrack,
				reconstructedSequence,
				reconstructedGraph,
				reconstructedGraph_levels
	);


	if(! superquiet)
		std::cout  << Utilities::timestamp() << "Bactracking done.\n" << std::flush;

	seedAndExtend_return thisBacktrace;
	thisBacktrace.Score = finalScore;
	thisBacktrace.graph_aligned = reconstructedGraph;
	thisBacktrace.graph_aligned_levels = reconstructedGraph_levels;
	thisBacktrace.sequence_aligned = reconstructedSequence;

	return thisBacktrace;
}






void GraphAligner_endsFree::seedAndExtend4_backtrack(VirtualNWTable& vNW2, std::string& sequence, double finalScore, int finalScore_z, NWEdge* finalScore_backtrack, std::string& reconstructedSequence, std::string& reconstructedGraph, std::vector<int>& reconstructedGraph_levels)
{
	double minusInfinity = -1 * numeric_limits<double>::max();

	int completed_x = g->NodesPerLevel.size() - 1;
	int completed_y = sequence.length();
	int completed_z = finalScore_z;
	NWEdge* currentEdge = 0; // = finalScore_backtrack;

	// if(verbose)
		std::cout << "Start backtrace...\n" << std::flush;

	vNW2.checkConsistency();

	std::set<NWEdge*> knownEdges = vNW2.getAllEdges();

	bool verboseBacktrack = false;

	bool firstStep = true;
	do {
		assert(completed_x >= 0);
		assert(completed_y >= 0);

		if(verboseBacktrack)
		{
			std::cout << "Backtrace step...\n";
			std::cout << "\t" << "completed_x: " << completed_x << "\n";
			std::cout << "\t" << "completed_y: " << completed_y << "\n";
			std::cout << "\t" << "completed_z: " << completed_z << "\n";
			std::cout << "\t" << "firstStep: " << firstStep << "\n" << std::flush;
		}

		if(! firstStep)
		{
			// std::cout << "\t" << "currentEdge: " << currentEdge << " with score after edge: " << currentEdge->scoreAfterEdge << "\n" << std::flush;
			assert(currentEdge != 0);
			currentEdge->checkConsistency();
			if(! knownEdges.count(currentEdge))
			{
				std::cerr << "! knownEdges.count(currentEdge): " << currentEdge << "\n" << std::flush;
			}
			assert(knownEdges.count(currentEdge));
		}
		else
		{
			if(verboseBacktrack) std::cout << "\t" << "finalScore: " << finalScore << "\n" << std::flush;
		}

		NWEdge* nextEdge_goingBack;
		if(firstStep)
		{
			nextEdge_goingBack = finalScore_backtrack;
		}
		else
		{
			nextEdge_goingBack = currentEdge->scoreBacktrack;
		}

//		std::cout << "\t" << "nextEdge_goingBack: " << nextEdge_goingBack << "\n" << std::flush;


		if((nextEdge_goingBack != 0) && (completed_x == nextEdge_goingBack->to_x) && (completed_y == nextEdge_goingBack->to_y) && (completed_z == nextEdge_goingBack->to_z))
		{
			// switch from the next edge to the current edge
			if(verboseBacktrack) std::cout << "\t" << "switch from " << currentEdge << " to " <<  nextEdge_goingBack << " as currentEdge.\n" << std::flush;
			currentEdge = nextEdge_goingBack;

		}
		else if((! firstStep) && (completed_x == currentEdge->to_x) && (completed_y == currentEdge->to_y) && (completed_z == currentEdge->to_z))
		{
			// output the current edge
			int diff_current_x = currentEdge->to_x - currentEdge->from_x;
			int diff_current_y = currentEdge->to_y - currentEdge->from_y;
			assert((diff_current_x >= 0) && (diff_current_x <= 1));
			assert((diff_current_y >= 0) && (diff_current_y <= 1));
			assert(diff_current_x || diff_current_y);

			std::string sequence_character;
			std::string graph_character;
			int graph_level;

			if((diff_current_x == 1) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);

				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 1) && (diff_current_y == 0))
			{
				// gap in sequence
				sequence_character = "_";
				Edge* usedGraphEdge = currentEdge->usedGraphEdge;
				assert(usedGraphEdge != 0);
				graph_character = g->CODE.deCode(usedGraphEdge->locus_id, usedGraphEdge->emission);
				assert(graph_character.length() == 1);
				graph_level = usedGraphEdge->From->level;
			}
			else if((diff_current_x == 0) && (diff_current_y == 1))
			{
				sequence_character = sequence.substr(currentEdge->from_y, 1);
				graph_character = "_";
				graph_level = -1;
			}

			reconstructedSequence.append(sequence_character);
			reconstructedGraph.append(graph_character);
			reconstructedGraph_levels.push_back(graph_level);

			if(verboseBacktrack) std::cout << "\t" << "output " << graph_character << " / " << sequence_character << " from current edge.\n" << std::flush;
			if(verboseBacktrack) std::cout << "\t\t" << "score *after* this edge " << currentEdge << " (in normal coordinates): " << currentEdge->scoreAfterEdge <<  " (i.e. after consuming the " << currentEdge->to_x << "-th graph level and the " << currentEdge->to_y << "-th sequence character)\n" << std::flush;

			completed_x = currentEdge->from_x;
			completed_y = currentEdge->from_y;
			completed_z = currentEdge->from_z;
		}
		else if(firstStep || ((completed_x == currentEdge->from_x) && (completed_y == currentEdge->from_y) && (completed_z == currentEdge->from_z)))
		{
			// complete the gap behind the current edge

			int need_gaps_graph;
			int need_gaps_sequence;

			if(nextEdge_goingBack == 0)
			{
				// going back to origin
				need_gaps_sequence = completed_x;
				need_gaps_graph = completed_y;
			}
			else
			{
				// going back to other normal edge
				need_gaps_sequence = completed_x - nextEdge_goingBack->to_x;
				need_gaps_graph = completed_y - nextEdge_goingBack->to_y;
			}

			if(nextEdge_goingBack == 0)
			{
				if(verboseBacktrack) std::cout << "\t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << 0 << " / " << 0  << "\n" << std::flush;
			}
			else
			{
				if(verboseBacktrack) std::cout << "\t" << "make big gap from  " << completed_x << " / " << completed_y <<  " to " << nextEdge_goingBack->to_x << " / " << nextEdge_goingBack->to_y  << "\n" << std::flush;
			}

			assert(need_gaps_sequence >= 0);
			assert(need_gaps_graph >= 0);

			// gap in graph - ie we take sequence, and graph has gaps
			std::string gaps_in_graph_graphCharacters;
			std::string gaps_in_graph_sequenceCharacters;
			gaps_in_graph_graphCharacters.resize(need_gaps_graph, '_');
			std::vector<int> gaps_in_graph_levels;

			if(need_gaps_graph > 0)
			{
				gaps_in_graph_levels.resize(need_gaps_graph, -1);
				int stopYForGapInGraph = (nextEdge_goingBack != 0) ? nextEdge_goingBack->to_y : 0;
				for(int seqI = completed_y - 1; seqI >= stopYForGapInGraph; seqI--)
				{
					std::string sequence_character = sequence.substr(seqI, 1);
					gaps_in_graph_sequenceCharacters.append(sequence_character);
				}
			}

			std::string gaps_in_sequence_graphCharacters;
			std::string gaps_in_sequence_sequenceCharacters;
			std::vector<int> gaps_in_sequence_levels;
			std::vector<Edge*> gaps_in_sequence_usedEdges;

			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startAffinely;
			std::map<int, std::map<int, graphPointDistance_withBacktrack>> distances_startNormally;

			if(need_gaps_sequence > 0)
			{
				int start_x_graph = completed_x;
				int stop_x_graph = (nextEdge_goingBack == 0) ? 0 : nextEdge_goingBack->to_x;

				int start_z_graph = completed_z;
				int stop_z_graph = (nextEdge_goingBack == 0) ? -1 : nextEdge_goingBack->to_z;

				gaps_in_sequence_sequenceCharacters.resize(need_gaps_sequence, '_');

				findShortGappedGraphConnection_endsFree(
						stop_x_graph,
						stop_z_graph,
						start_x_graph,
						start_z_graph,
						distances_startAffinely,
						distances_startNormally
				);
			}

			if(need_gaps_graph && need_gaps_sequence)
			{
				// we need both sequence gap and graph gap. In order to find the optimal sequence, we need
				// to find out which one to place first.
				// this code strongly resembles some code from above, where we compute the actual optimal score.

				bool exitEdge_isAffineSequenceGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isSequenceGap_affine());
				bool exitEdge_isGraphGap = ((nextEdge_goingBack == 0) ? false : nextEdge_goingBack->isGraphGap());

				std::vector<double> scores;
				std::vector<backtrackBookkeeper> backtracks;
				std::vector<bool> sequenceGapFirst;

				assert((nextEdge_goingBack == 0) || (distances_startNormally.size() == 1));
				for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
				{
					int zStart = zStarterIt->first;
					std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
					assert(zStops.size() == 1);

					for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
					{
						int zStop = zStopIt->first;

						{
							// case 1: sequence gap comes first
							{
								// case 1.1: the sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = (exitEdge_isAffineSequenceGap) ? distances_startAffinely.at(zStart).at(zStop) : distances_startNormally.at(zStart).at(zStop);
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double graphGapScore = S_openGap + need_gaps_graph * S_extendGap;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, true, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = sequenceGapScore + graphGapScore + scoreFromEntryEdge;
								if(firstStep)
								{
									local_combinedScore = minusInfinity;
								}
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(true);
							}

							{
								// case 1.2: the sequence gap ends affinely - not relevant, as we know that a graph gap follows
							}
						}
						{
							// case 2: graph gap comes first
							{
								// case 2.1: sequence gap ends in anything
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endInAnything.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								if(nextEdge_goingBack == 0)
								{
									local_combinedScore = minusInfinity;
								}
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
								sequenceGapFirst.push_back(false);
							}
							{
								// case 2.2: sequence gap ends affinely
								graphPointDistance_withBacktrack& utilizedBacktrack = distances_startNormally.at(zStart).at(zStop);

								double graphGapOpenPenalty = (exitEdge_isGraphGap) ? 0 : S_openGap;
								double graphGapScore = graphGapOpenPenalty + need_gaps_graph * S_extendGap;
								double sequenceGapScore = utilizedBacktrack.endAffinely.S;
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double local_combinedScore = graphGapScore + sequenceGapScore + scoreFromEntryEdge;
								if(nextEdge_goingBack == 0)
								{
									local_combinedScore = minusInfinity;
								}
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
								sequenceGapFirst.push_back(false);
							}
						}
					}
				}

				assert(scores.size() > 0);
				std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMax(scores);

				backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);
				bool selectedSequenceGapFirst = sequenceGapFirst.at(maxScore_sequenceGap.second);

				gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
				gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

				std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
				std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				if(verboseBacktrack) std::cout << "\t\t\t" << "selectedSequenceGapFirst: " << selectedSequenceGapFirst << "\n" << std::flush;


				if(selectedSequenceGapFirst)
				{
					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());

					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());
				}
				else
				{
					reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
					reconstructedGraph.append(gaps_in_sequence_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

					reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
					reconstructedGraph.append(gaps_in_graph_graphCharacters);
					reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());
				}

			}
			else
			{
				// we either need the graph gap or the sequence gap. The other fields are going to be empty. The orde
				// of concatentation is thus irrelevant.

				if(need_gaps_sequence > 0)
				{
					// if we go to the origin, we may end up in multiple z values. We traverse all z values to
					// find the best alternative.

					std::vector<double> scores;
					std::vector<backtrackBookkeeper> backtracks;

					assert((nextEdge_goingBack == 0) || (distances_startNormally.size() == 1));

					for(std::map<int, std::map<int, graphPointDistance_withBacktrack>>::iterator zStarterIt = distances_startNormally.begin(); zStarterIt != distances_startNormally.end(); zStarterIt++)
					{
						int zStart = zStarterIt->first;
						std::map<int, graphPointDistance_withBacktrack>& zStops = zStarterIt->second;
						assert(zStops.size() == 1);

						for(std::map<int, graphPointDistance_withBacktrack>::iterator zStopIt = zStops.begin(); zStopIt != zStops.end(); zStopIt++)
						{

							int zStop = zStopIt->first;

							graphPointDistance_withBacktrack utilizedBacktrack = zStopIt->second;
							if((nextEdge_goingBack != 0) && (nextEdge_goingBack->isSequenceGap_affine()))
							{
								utilizedBacktrack = distances_startAffinely.at(zStart).at(zStop);
							}

							// case 1: we want to end in an affine gap
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(true, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endAffinely.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endAffinely);
							}

							// case 2: end in anything
							{
								double scoreFromEntryEdge = (currentEdge != 0) ? currentEdge->calculateScore_endsFree(false, false, S_match, S_mismatch, S_openGap, S_extendGap) : 0;
								double scoreGraphDistance = utilizedBacktrack.endInAnything.S;
								double local_combinedScore = scoreGraphDistance + scoreFromEntryEdge;
								scores.push_back(local_combinedScore);
								backtracks.push_back(utilizedBacktrack.endInAnything);
							}

						}
					}

					assert(scores.size() > 0);
					std::pair<double, unsigned int> maxScore_sequenceGap = Utilities::findVectorMax(scores);
					backtrackBookkeeper selectedBacktrack = backtracks.at(maxScore_sequenceGap.second);

					gaps_in_sequence_graphCharacters = selectedBacktrack.graphSequence;
					gaps_in_sequence_levels = selectedBacktrack.graphSequence_levels;

					std::reverse(gaps_in_sequence_graphCharacters.begin(), gaps_in_sequence_graphCharacters.end());
					std::reverse(gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				}

				reconstructedSequence.append(gaps_in_sequence_sequenceCharacters);
				reconstructedGraph.append(gaps_in_sequence_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_sequence_levels.begin(), gaps_in_sequence_levels.end());

				reconstructedSequence.append(gaps_in_graph_sequenceCharacters);
				reconstructedGraph.append(gaps_in_graph_graphCharacters);
				reconstructedGraph_levels.insert(reconstructedGraph_levels.end(), gaps_in_graph_levels.begin(), gaps_in_graph_levels.end());
			}

			// gap in
			if(nextEdge_goingBack == 0)
			{
				completed_x = 0;
				completed_y = 0;
			}
			else
			{
				completed_x = nextEdge_goingBack->to_x;
				completed_y = nextEdge_goingBack->to_y;
				completed_z = nextEdge_goingBack->to_z;
			}

			if(firstStep)
			{
				currentEdge = nextEdge_goingBack;
			}
		}
		else
		{
			// this case is unforeseen!
			std::cout << "\tUNEXPECTED!\n";
			std::cout << "\t\t" << "currentEdge: " << currentEdge << "\n";
			std::cout << "\t\t" << "currentEdge->from_x: " << currentEdge->from_x << "\n";
			std::cout << "\t\t" << "currentEdge->from_y: " << currentEdge->from_y << "\n";
			std::cout << "\t\t" << "currentEdge->from_z: " << currentEdge->from_z << "\n";
			std::cout << "\t\t" << "currentEdge->to_x: " << currentEdge->to_x << "\n";
			std::cout << "\t\t" << "currentEdge->to_y: " << currentEdge->to_y << "\n";
			std::cout << "\t\t" << "currentEdge->to_z: " << currentEdge->to_z << "\n\n" << std::flush;
			if(nextEdge_goingBack != 0)
			{
				std::cout << "\t\t\t" << "nextEdge: " << nextEdge_goingBack << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_x: " << nextEdge_goingBack->from_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_y: " << nextEdge_goingBack->from_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->from_z: " << nextEdge_goingBack->from_z << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_x: " << nextEdge_goingBack->to_x << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_y: " << nextEdge_goingBack->to_y << "\n";
				std::cout << "\t\t\t" << "nextEdge_goingBack->to_z: " << nextEdge_goingBack->to_z << "\n\n" << std::flush;
			}

			assert( 1 == 0 );
		}

		firstStep = false;

	} while((completed_x != 0) || (completed_y != 0));

	if(verboseBacktrack)	std::cout << "Backtrace done!\n" << std::flush;

	assert(reconstructedGraph.length() == reconstructedSequence.length());
	assert(reconstructedGraph_levels.size() == reconstructedGraph.length());

	std::reverse(reconstructedGraph.begin(), reconstructedGraph.end());
	std::reverse(reconstructedSequence.begin(), reconstructedSequence.end());
	std::reverse(reconstructedGraph_levels.begin(), reconstructedGraph_levels.end());
}


void GraphAligner_endsFree::findShortGappedGraphConnection_endsFree(int start_x, int start_z, int stop_x, int stop_z, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startAffinely, std::map<int, std::map<int, graphPointDistance_withBacktrack>>& distances_startNormally)
{
	// this is not really ends-free - scoring does not take into account that the edges cost 0,
	// but not really relevant either!

	assert(S_graphGap == 0);
	double minusInfinity = -1 * numeric_limits<double>::max();

	assert(stop_x > start_x);
	assert(start_x  >= 0);
	assert(stop_x < (int)g->NodesPerLevel.size());

	assert((start_z == -1) || ((start_z >= 0) && (start_z < (int)g->NodesPerLevel.at(start_x).size())));
	assert((stop_z == -1) || ((stop_z >= 0) && (stop_z < (int)g->NodesPerLevel.at(stop_x).size())));

	// bool verbose = true;
	distances_startAffinely.clear();
	distances_startNormally.clear();

	std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_notStartInAffineGap_endInAnything;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_notStartInAffineGap_endInAffineGap;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_startInAffineGap_endInAnything;
	std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_startInAffineGap_endInAffineGap;

	for(int lI = start_x; lI <= stop_x; lI++)
	{

		std::set<Node*> nodes_thisLevel = g->NodesPerLevel.at(lI);

		std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		std::map<Node*, std::map<Node*, backtrackBookkeeper> > runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

//		if(lI == start_x)
//		{
//			lastPositionDistances_perZ_startInAffineGap.resize(nodes_thisLevel.size());
//			lastPositionDistances_perZ_startNormal.resize(nodes_thisLevel.size());
//		}

		for(std::set<Node*>::iterator nIt = nodes_thisLevel.begin(); nIt != nodes_thisLevel.end(); nIt++)
		{

			Node* n = *nIt;
			int z = nodesPerLevel_ordered_rev.at(lI).at(n);

			if(lI == start_x)
			{
				// jump over node if start_z specified
				bool takeNode = true;
				if(start_z != -1)
				{
					if(z != start_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n][n].S = 0;
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n][n].S = minusInfinity;

					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
					runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n][n].S = 0;
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();;
					runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n][n].S = (lI == 0) ? minusInfinity : 0;
				}
			}
			else
			{
				std::set<Edge*> incomingEdges = n->Incoming_Edges;
				assert(incomingEdges.size() > 0);

				runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_startInAffineGap_endInAnything_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();
				runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel[n] = std::map<Node*, backtrackBookkeeper>();

				for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
				{
					Edge* e = *eIt;
					std::string edgeEmission = g->CODE.deCode(e->locus_id, e->emission);
					assert(edgeEmission.length() == 1);

					Node* nodeFrom = e->From;

					if(! runningNodeDistances_notStartInAffineGap_endInAnything.count(nodeFrom))
					{
						continue;
					}

					for(std::map<Node*, backtrackBookkeeper>::iterator nodeDistIt = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).begin(); nodeDistIt != runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).end(); nodeDistIt++)
					{
						Node* interestingNode = nodeDistIt->first;

						backtrackBookkeeper previousDistance_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything.at(nodeFrom).at(interestingNode);
						backtrackBookkeeper previousDistance_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap.at(nodeFrom).at(interestingNode);
						backtrackBookkeeper previousDistance_startInAffineGap_endInAnything = runningNodeDistances_startInAffineGap_endInAnything.at(nodeFrom).at(interestingNode);
						backtrackBookkeeper previousDistance_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap.at(nodeFrom).at(interestingNode);

						// first part: we end in an affine gap
						{
							// for the ones that did not start with an affine gap
							{
								backtrackBookkeeper newDistance_notStartInAffineGap_affine_extend = previousDistance_notStartInAffineGap_endInAffineGap;
								backtrackBookkeeper newDistance_notStartInAffineGap_affine_start = previousDistance_notStartInAffineGap_endInAnything;

								if(edgeEmission == "_")
								{
									newDistance_notStartInAffineGap_affine_start.S = minusInfinity;

									newDistance_notStartInAffineGap_affine_extend.S += 0;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_notStartInAffineGap_affine_extend.S += S_extendGap;
									newDistance_notStartInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);

									newDistance_notStartInAffineGap_affine_start.S += (S_openGap + S_extendGap);
									newDistance_notStartInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);
								}

								backtrackBookkeeper newDistance_notStartInAffineGap_affine_maximum = (newDistance_notStartInAffineGap_affine_start.S > newDistance_notStartInAffineGap_affine_extend.S) ? newDistance_notStartInAffineGap_affine_start : newDistance_notStartInAffineGap_affine_extend;

								if((runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode).S < newDistance_notStartInAffineGap_affine_maximum.S))
								{
									runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n)[interestingNode] = newDistance_notStartInAffineGap_affine_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								backtrackBookkeeper newDistance_startInAffineGap_affine_extend = previousDistance_startInAffineGap_endInAffineGap;
								backtrackBookkeeper newDistance_startInAffineGap_affine_start = previousDistance_startInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									newDistance_startInAffineGap_affine_start.S = minusInfinity;

									newDistance_startInAffineGap_affine_extend.S += 0;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									newDistance_startInAffineGap_affine_extend.S += S_extendGap;
									newDistance_startInAffineGap_affine_extend.takeNodeAndEdge(nodeFrom, e);

									newDistance_startInAffineGap_affine_start.S += (S_openGap + S_extendGap);
									newDistance_startInAffineGap_affine_start.takeNodeAndEdge(nodeFrom, e);

								}

								backtrackBookkeeper newDistance_startInAffineGap_affine_maximum = (newDistance_startInAffineGap_affine_start.S > newDistance_startInAffineGap_affine_extend.S) ? newDistance_startInAffineGap_affine_start : newDistance_startInAffineGap_affine_extend;

								if((runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode).S < newDistance_startInAffineGap_affine_maximum.S))
								{
									runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n)[interestingNode] = newDistance_startInAffineGap_affine_maximum;
								}
							}
						}

						//second part: we end in anything
						{
							// for the ones that did not start with an affine gap
							{
								backtrackBookkeeper alternativeScore_endInAffine = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode);

								backtrackBookkeeper followGraphGap = previousDistance_notStartInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap.S = previousDistance_notStartInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine.S) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n).at(interestingNode).S < endInAnything_maximum.S))
								{
									runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(n)[interestingNode] = endInAnything_maximum;
								}
							}

							// for the ones that did start with an affine gap
							{
								backtrackBookkeeper alternativeScore_endInAffine = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).at(interestingNode);
								backtrackBookkeeper followGraphGap = previousDistance_startInAffineGap_endInAnything;
								if(edgeEmission == "_")
								{
									// followGraphGap = previousDistance_startInAffineGap_endInAnything;
									followGraphGap.takeNodeAndEdge(nodeFrom, e);
								}
								else
								{
									followGraphGap.S = minusInfinity;
								}

								backtrackBookkeeper endInAnything_maximum = (followGraphGap.S > alternativeScore_endInAffine.S) ? followGraphGap : alternativeScore_endInAffine;
								if((runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).count(interestingNode) == 0) || (runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n).at(interestingNode).S < endInAnything_maximum.S))
								{
									runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(n)[interestingNode] = endInAnything_maximum;
								}
							}
						}
					}
				}
			}

			if(lI == stop_x)
			{
				bool takeNode = true;
				if(stop_z != -1)
				{
					if(z != stop_z)
					{
						takeNode = false;
					}
				}

				if(takeNode)
				{

					for(std::map<Node*, backtrackBookkeeper>::iterator startNodeIt = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).begin(); startNodeIt != runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(n).end(); startNodeIt++)
					{
						Node* startNode = startNodeIt->first;
						Node* stopNode = n;

						int z_startNode = nodesPerLevel_ordered_rev.at(start_x).at(startNode);
						int z_stopNode = nodesPerLevel_ordered_rev.at(stop_x).at(stopNode);

						{
							// start affinely

							graphPointDistance_withBacktrack distanceSpecifier;

							distanceSpecifier.endAffinely = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel.at(stopNode).at(startNode);
							distanceSpecifier.endInAnything = runningNodeDistances_startInAffineGap_endInAnything_thisLevel.at(stopNode).at(startNode);
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startAffinely.count(z_startNode) == 0) || (distances_startAffinely.at(z_startNode).count(z_stopNode) == 0));
							distances_startAffinely[z_startNode][z_stopNode] = distanceSpecifier;
						}

						{
							// start arbitrarily

							graphPointDistance_withBacktrack distanceSpecifier;

							distanceSpecifier.endAffinely = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel.at(stopNode).at(startNode);
							distanceSpecifier.endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel.at(stopNode).at(startNode);
							distanceSpecifier.start_in_affine_sequenceGap = true;
							assert(distanceSpecifier.endAffinely.S <= 0);
							assert(distanceSpecifier.endInAnything.S <= 0);

							assert((distances_startNormally.count(z_startNode) == 0) || (distances_startNormally.at(z_startNode).count(z_stopNode) == 0));
							distances_startNormally[z_startNode][z_stopNode] = distanceSpecifier;
						}
					}
				}
			}
		}

		runningNodeDistances_notStartInAffineGap_endInAnything = runningNodeDistances_notStartInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_notStartInAffineGap_endInAffineGap = runningNodeDistances_notStartInAffineGap_endInAffineGap_thisLevel;

		runningNodeDistances_startInAffineGap_endInAnything	= runningNodeDistances_startInAffineGap_endInAnything_thisLevel;
		runningNodeDistances_startInAffineGap_endInAffineGap = runningNodeDistances_startInAffineGap_endInAffineGap_thisLevel;

	}
}



double GraphAligner_endsFree::score_endsFree(std::string reconstructed_graph, std::vector<int> reconstructed_graph_levels, std::string reconstructed_sequence)
{
	assert(reconstructed_graph.length() == reconstructed_sequence.length());
	assert(reconstructed_graph_levels.size() == reconstructed_graph.length());
	assert(S_graphGap == 0);

	double score = 0;
	bool inAffineGapInGraph = false;
	bool inAffineGapInSequence = false;

	std::vector<double> runningScores;

	std::vector<bool> endsStatus_sequence;
	endsStatus_sequence.resize(reconstructed_graph.length(), false);
	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		if(S == "_")
		{
			endsStatus_sequence.at(i) = true;
		}
		else
		{
			break;
		}
	}
	for(int i = (reconstructed_graph.length() - 1); i >= 0; i--)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		if(S == "_")
		{
			endsStatus_sequence.at(i) = true;
		}
		else
		{
			break;
		}
	}

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::string G = reconstructed_graph.substr(i, 1);

		if((S == "_") && (endsStatus_sequence.at(i) == true))
		{
			score += 0;
			inAffineGapInSequence = false;
			inAffineGapInGraph = false;
		}
		else
		{
			if(S != "_")
			{
				inAffineGapInSequence = false;
			}
			if(G != "_")
			{
				inAffineGapInGraph = false;
			}

			if(G == "_")
			{
				if(reconstructed_graph_levels.at(i) != -1)
				{
					inAffineGapInGraph = false;
					// inAffineGapInSequence = false;

					if(S == "_")
					{
						score += S_graphGap;
					}
					else
					{
						score += S_mismatch;
					}
				}
				else
				{
					assert(S != "_");
					if(inAffineGapInGraph == true)
					{
						score += S_extendGap;
					}
					else
					{
						inAffineGapInGraph = true;
						score += (S_openGap + S_extendGap);
					}
				}
			}
			else
			{
				if(S == "_")
				{
					if(inAffineGapInSequence == true)
					{
						score += S_extendGap;
					}
					else
					{
						inAffineGapInSequence = true;
						score += (S_openGap + S_extendGap);
					}
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
		runningScores.push_back(score);
	}

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		int l = reconstructed_graph_levels.at(i);
		if(l != -1)
			l++;
		std::cout << Utilities::ItoStr(l) << "\t";
	}
	std::cout << "\n";

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string G = reconstructed_graph.substr(i, 1);
		std::cout << G << "\t";
	}
	std::cout << "\n";

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		std::string S = reconstructed_sequence.substr(i, 1);
		std::cout << S << "\t";
	}
	std::cout << "\n";

	for(unsigned int i = 0; i < reconstructed_graph.length(); i++)
	{
		double S = runningScores.at(i);
		std::cout << S << "\t";
	}
	std::cout << "\n" << std::flush;



	return score;
}


