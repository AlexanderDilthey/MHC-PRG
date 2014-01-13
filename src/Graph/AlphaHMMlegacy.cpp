/*
 * AlphaHMMlegacy.cpp
 *
 *  Created on: 22 May 2012
 *      Author: dilthey
 */


#include "AlphaHMM.h"
#include <assert.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <vector>
#include "../Utilities.h"
#include "../MHC-PRG.h"
#include <boost/math/distributions/poisson.hpp>

using namespace boost::math::policies;
using namespace boost::math;

typedef boost::math::poisson_distribution< double, policy < discrete_quantile < integer_round_inwards> > > poisson_up;






double AlphaHMM::fillForwardBackwardTable_5(map<string, long long> globalEmission)
{
	assert(genotypingMode == 5);

	int levels = diploidStatesByLevel.size();

	fw_Viterbi.resize(levels);
	fw_Emission_e1.resize(levels);
	fw_Emission_e2.resize(levels);
	fw_Viterbi_backtrack_int.resize(levels);

	//fw_adaptedEmission.resize(levels);
	//fw_adaptedEmission_e1.resize(levels);
	//fw_adaptedEmission_e2.resize(levels);
	//fw_expectedMissing_e1.resize(levels);
	//fw_expectedMissing_e2.resize(levels);
	//fw_Coverage_e1.resize(levels);
	//fw_Coverage_e2.resize(levels);
	//fw_Viterbi_backtrack.resize(levels);

	fw_underflow_factor = 1;

	//fw_Emission.resize(levels);
	//fw_e1_length.resize(levels);

	// fw_e1_length_valid.resize(levels);
	// fw_e2_length.resize(levels);
	// fw_e2_length_valid.resize(levels);

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		if((level % 1000) == 0)
		{
			int startLevelFill = level;
			int stopLevelFill = level + 1000 - 1;

			cout << "Fill state transitions table from " << startLevelFill << " to " << stopLevelFill << "\n" << flush;

			int previousStepFillStop = startLevelFill - 1;
			int previousStepFillStart = previousStepFillStop - 1000 + 1;

			if(previousStepFillStop >= 0)
			{
				if(previousStepFillStart < 0)
				{
					previousStepFillStart = 0;
				}

				// cout << "\t clear state transitions table from " << previousStepFillStart << " to " << previousStepFillStop << "\n" << flush;

				for(int k = previousStepFillStart; k <= previousStepFillStop; k++)
				{
					assert((int)diploidStateTransitions.size() > k);
					assert((int)diploidStateTransitions_Reverse.size() > k);
					diploidStateTransitions.at(k).clear();
					diploidStateTransitions_Reverse.at(k).clear();
				}
			}

			if (stopLevelFill > (levels - 1))
			{
				stopLevelFill = levels - 1;
			}
			int pieces_to_divide = stopLevelFill - startLevelFill + 1;
			assert(pieces_to_divide > 0);
			int chunk_size = pieces_to_divide / CONFIG.threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == CONFIG.threads);
				int thisThread = omp_get_thread_num();
				int firstLevel = thisThread * chunk_size;
				int lastLevel = (thisThread+1) * chunk_size - 1;

				if((thisThread == (CONFIG.threads-1)) && (lastLevel < (pieces_to_divide-1)))
				{
					lastLevel = pieces_to_divide-1;
				}

				assert(firstLevel >= 0);
				if(lastLevel >= 0)
				{
					#pragma omp critical
					{
						// cout << "\t thread " << thisThread << " from " << firstLevel+startLevelFill << " to " << lastLevel+startLevelFill << "\n" << flush;
					}

					for(int k = firstLevel; k <= lastLevel; k++)
					{
						fillStateTransitions(startLevelFill+k);
					}
				}
			}
		}

		int states = diploidStatesByLevel.at(level);

		//fw.at(level).resize(states);

		fw_Viterbi.at(level).resize(states);
		fw_Emission_e1.at(level).resize(states);
		fw_Emission_e2.at(level).resize(states);
		fw_Viterbi_backtrack_int.at(level).resize(states);

		//fw_Emission.at(level).resize(states);
		//fw_adaptedEmission.at(level).resize(states);
		//fw_adaptedEmission_e1.at(level).resize(states);
		//fw_adaptedEmission_e2.at(level).resize(states);
		//fw_expectedMissing_e1.at(level).resize(states);
		//fw_expectedMissing_e2.at(level).resize(states);

		//fw_Coverage_e1.at(level).resize(states);
		//fw_Coverage_e2.at(level).resize(states);
		//fw_e1_length.at(level).resize(states);
		//fw_e1_length_valid.at(level).resize(states);
		//fw_e2_length.at(level).resize(states);
		//fw_e2_length_valid.at(level).resize(states);
		//fw_Viterbi_backtrack.at(level).resize(states);


		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}

					local_edge_length = local_edge_length / 2;
				}


				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);
				assert(total_edge_length == haploid_edge_lengths.at(level));

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{

							if(observedNumber < (int)observedXunderlying.size())
							{
								for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
								{
									int missing = kMerCount - possibleUnderlyingCount;
									double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
									thisKMer_estimated_missing += p * (double)missing;
								}
							}
							else
							{
								this_kMer_expected_missing = 0;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double expected_missing_coverage = coverage * 1.0;
									double lambda = (double)kMerCount*expected_missing_coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);
					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				if(!(total_kMers_e1 == haploidEdgeValidLength(level, state, 1)))
				{
					cout << "(total_kMers_e1 == haploidEdgeValidLength(level, state, 1))!\n";
					cout << "level: " << level << "\n";
					cout << "diploid state: " << state << "\n";
					cout << "total_kMers_e1: " << total_kMers_e1 << "\n";
					cout << "haploidEdgeValidLength(level, state, 1): " << haploidEdgeValidLength(level, state, 1) << "\n";
				}
				assert(total_kMers_e1 == haploidEdgeValidLength(level, state, 1));
				if(!(total_kMers_e2 == haploidEdgeValidLength(level, state, 2)))
				{
					cout << "(total_kMers_e2 == haploidEdgeValidLength(level, state, 2))!\n";
					cout << "level: " << level << "\n";
					cout << "diploid state: " << state << "\n";
					cout << "total_kMers_e2: " << total_kMers_e2 << "\n";
					cout << "haploidEdgeValidLength(level, state, 2): " << haploidEdgeValidLength(level, state, 2) << "\n";
				}
				assert(total_kMers_e2 == haploidEdgeValidLength(level, state, 2));

				if(level == 0)
				{
					fw_Viterbi.at(level).at(state) = 0;
					fw_Viterbi_backtrack_int.at(level).at(state) = -1;

					//vector<int> states_same_jump_P;
					//fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
						}

					}

					assert(max_jump_state != -1);

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack_int.at(level).at(state) = max_jump_state;

					//fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					//assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
				}

				if(total_kMers_e1 == 0)
				{
					symbol_emission_e1 = 1;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}

				if(total_kMers_e2 == 0)
				{
					symbol_emission_e2 = 1;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				int missing_edges_e1 = total_edge_length - total_kMers_e1;
				int missing_edges_e2 = total_edge_length - total_kMers_e2;

				int informative_back_edges_e1 = 0;
				if (missing_edges_e1 > 0)
				{
					informative_back_edges_e1 = (missing_edges_e1 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e1;
				}
				int informative_back_edges_e2 = 0;
				if (missing_edges_e2 > 0)
				{
					informative_back_edges_e2 = (missing_edges_e2 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e2;
				}

				double adapted_emission_e1 = symbol_emission_e1;
				double adapted_emission_e2 = symbol_emission_e2;

				if(((double)total_kMers_e1/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e1 = total_kMers_e1;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e1;
					weights.push_back(total_kMers_e1);
					values.push_back(symbol_emission_e1);
					while((current_level > 0) && (got_backtracked_edges_e1 < informative_back_edges_e1))
					{
						running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
						current_level--;

						int running_level_edges = haploidEdgeValidLength(current_level, running_state, 1);
						double running_level_emission = fw_Emission_e1.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e1 + running_level_edges) > informative_back_edges_e1)
						{
							int too_many = (got_backtracked_edges_e1 + running_level_edges) - informative_back_edges_e1;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e1 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					if(weights_sum > 0)
					{
						assert(weights_sum > 0);
						adapted_emission_e1 = 0;
						for(int wI = 0; wI < (int)weights.size(); wI++)
						{
							adapted_emission_e1 += (weights.at(wI)/weights_sum) * values.at(wI);
						}
					}
					else
					{
						adapted_emission_e1 = 1;
					}

				}

				if(((double)total_kMers_e2/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e2 = total_kMers_e2;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e2;
					weights.push_back(total_kMers_e2);
					values.push_back(symbol_emission_e2);
					while((current_level > 0) && (got_backtracked_edges_e2 < informative_back_edges_e2))
					{
						running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
						current_level--;

						int running_level_edges = haploidEdgeValidLength(current_level, running_state, 2);
						double running_level_emission = fw_Emission_e2.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e2 + running_level_edges) > informative_back_edges_e2)
						{
							int too_many = (got_backtracked_edges_e2 + running_level_edges) - informative_back_edges_e2;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e2 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					if(weights_sum > 0)
					{
						assert(weights_sum > 0);
						adapted_emission_e2 = 0;
						for(int wI = 0; wI < (int)weights.size(); wI++)
						{
							adapted_emission_e2 += (weights.at(wI)/weights_sum) * values.at(wI);
						}
					}
					else
					{
						adapted_emission_e2 = 1;
					}
				}

				double symbol_emission = 0.5*adapted_emission_e1 + 0.5*adapted_emission_e2;

				double e1_lengthWeight = 1 + 0.10 * ((double)total_kMers_e1/(double)total_edge_length);
				double e2_lengthWeight = 1 + 0.10 * ((double)total_kMers_e2/(double)total_edge_length);

				double symbol_emission_lengthWeighted  =  0.5*adapted_emission_e1*e1_lengthWeight + 0.5*adapted_emission_e2*e2_lengthWeight;

				assert(symbol_emission >= 0);

				//fw_adaptedEmission.at(level).at(state) = symbol_emission;
				//fw_adaptedEmission_e1.at(level).at(state) = adapted_emission_e1*e1_lengthWeight;
				//fw_adaptedEmission_e2.at(level).at(state) = adapted_emission_e2*e2_lengthWeight;

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission_lengthWeighted*(double)total_edge_length);
				//fw_Emission.at(level).at(state) = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				fw_Emission_e1.at(level).at(state) = symbol_emission_e1;
				fw_Emission_e2.at(level).at(state) = symbol_emission_e2;

				//fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				//fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				//assert(fw.at(level).at(state) >= 0);
			}

		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}

double AlphaHMM::fillForwardBackwardTable_6(map<string, long long> globalEmission)
{
	assert(genotypingMode == 6);

	int levels = diploidStatesByLevel.size();

	fw_Viterbi.resize(levels);
	fw_Emission_e1.resize(levels);
	fw_Emission_e2.resize(levels);
	fw_Viterbi_backtrack_int.resize(levels);

	//fw_adaptedEmission.resize(levels);
	//fw_adaptedEmission_e1.resize(levels);
	//fw_adaptedEmission_e2.resize(levels);
	//fw_expectedMissing_e1.resize(levels);
	//fw_expectedMissing_e2.resize(levels);
	//fw_Coverage_e1.resize(levels);
	//fw_Coverage_e2.resize(levels);
	//fw_Viterbi_backtrack.resize(levels);

	fw_underflow_factor = 1;

	//fw_Emission.resize(levels);
	//fw_e1_length.resize(levels);

	// fw_e1_length_valid.resize(levels);
	// fw_e2_length.resize(levels);
	// fw_e2_length_valid.resize(levels);

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	kMerUniquenessInfo uniqueMers = myGraph->kMerUniqueness();

	cout << levels << " Levels\n";
	cout << uniqueMers.levelUniqueKMers.size() << " graph-unique kMers\n";

	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		if((level % 1000) == 0)
		{
			int startLevelFill = level;
			int stopLevelFill = level + 1000 - 1;

			cout << "Fill state transitions table from " << startLevelFill << " to " << stopLevelFill << "\n" << flush;

			int previousStepFillStop = startLevelFill - 1;
			int previousStepFillStart = previousStepFillStop - 1000 + 1;

			if(previousStepFillStop >= 0)
			{
				if(previousStepFillStart < 0)
				{
					previousStepFillStart = 0;
				}

				// cout << "\t clear state transitions table from " << previousStepFillStart << " to " << previousStepFillStop << "\n" << flush;

				for(int k = previousStepFillStart; k <= previousStepFillStop; k++)
				{
					assert((int)diploidStateTransitions.size() > k);
					assert((int)diploidStateTransitions_Reverse.size() > k);
					diploidStateTransitions.at(k).clear();
					diploidStateTransitions_Reverse.at(k).clear();
				}
			}

			if (stopLevelFill > (levels - 1))
			{
				stopLevelFill = levels - 1;
			}
			int pieces_to_divide = stopLevelFill - startLevelFill + 1;
			assert(pieces_to_divide > 0);
			int chunk_size = pieces_to_divide / CONFIG.threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == CONFIG.threads);
				int thisThread = omp_get_thread_num();
				int firstLevel = thisThread * chunk_size;
				int lastLevel = (thisThread+1) * chunk_size - 1;

				if((thisThread == (CONFIG.threads-1)) && (lastLevel < (pieces_to_divide-1)))
				{
					lastLevel = pieces_to_divide-1;
				}

				assert(firstLevel >= 0);
				if(lastLevel >= 0)
				{
					#pragma omp critical
					{
						// cout << "\t thread " << thisThread << " from " << firstLevel+startLevelFill << " to " << lastLevel+startLevelFill << "\n" << flush;
					}

					for(int k = firstLevel; k <= lastLevel; k++)
					{
						fillStateTransitions(startLevelFill+k);
					}
				}
			}
		}

		int states = diploidStatesByLevel.at(level);

		//fw.at(level).resize(states);

		fw_Viterbi.at(level).resize(states);
		fw_Emission_e1.at(level).resize(states);
		fw_Emission_e2.at(level).resize(states);
		fw_Viterbi_backtrack_int.at(level).resize(states);

		//fw_Emission.at(level).resize(states);
		//fw_adaptedEmission.at(level).resize(states);
		//fw_adaptedEmission_e1.at(level).resize(states);
		//fw_adaptedEmission_e2.at(level).resize(states);
		//fw_expectedMissing_e1.at(level).resize(states);
		//fw_expectedMissing_e2.at(level).resize(states);

		//fw_Coverage_e1.at(level).resize(states);
		//fw_Coverage_e2.at(level).resize(states);
		//fw_e1_length.at(level).resize(states);
		//fw_e1_length_valid.at(level).resize(states);
		//fw_e2_length.at(level).resize(states);
		//fw_e2_length_valid.at(level).resize(states);
		//fw_Viterbi_backtrack.at(level).resize(states);


		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;
			map< int, bool > cache_kMerSymbolUnique;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}

					local_edge_length = local_edge_length / 2;
				}


				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);
				assert(total_edge_length == haploid_edge_lengths.at(level));

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;
				double e1_kMer_weight = 0;
				double e2_kMer_weight = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					if(cache_kMerSymbolUnique.count(kMerSymbol) == 0)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						cache_kMerSymbolUnique[kMerSymbol] = (uniqueMers.levelUniqueKMers.count(kMer) > 0) ? true : false;
					}

					bool kMerUnique = cache_kMerSymbolUnique[kMerSymbol];
					double kMerWeight = 1;

					if(kMerUnique)
					{
						kMerWeight = 20;
					}

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{


							if(observedNumber < (int)observedXunderlying.size())
							{
								for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
								{
									int missing = kMerCount - possibleUnderlyingCount;
									double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
									thisKMer_estimated_missing += p * (double)missing;
								}
							}
							else
							{
								this_kMer_expected_missing = 0;
							}


							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double expected_missing_coverage = coverage * 1.0;
									double lambda = (double)kMerCount*expected_missing_coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * kMerWeight * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * kMerWeight * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * kMerWeight * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * kMerWeight * ((double)e2_count / (double)kMerCount);

					e1_kMer_weight += kMerWeight * (double)e1_count;
					e2_kMer_weight += kMerWeight * (double)e2_count;
				}

				if(e1_kMer_weight > 0)
				{
					estimated_missing_e1 = (estimated_missing_e1 / e1_kMer_weight) * total_kMers_e1;
					expected_missing_e1 = (expected_missing_e1 / e1_kMer_weight) * total_kMers_e1;
				}

				if(e2_kMer_weight > 0)
				{
					estimated_missing_e2 = (estimated_missing_e2 / e2_kMer_weight) * total_kMers_e2;
					expected_missing_e2 = (expected_missing_e2 / e2_kMer_weight) * total_kMers_e2;
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				if(!(total_kMers_e1 == haploidEdgeValidLength(level, state, 1)))
				{
					cout << "(total_kMers_e1 == haploidEdgeValidLength(level, state, 1))!\n";
					cout << "level: " << level << "\n";
					cout << "diploid state: " << state << "\n";
					cout << "total_kMers_e1: " << total_kMers_e1 << "\n";
					cout << "haploidEdgeValidLength(level, state, 1): " << haploidEdgeValidLength(level, state, 1) << "\n";
				}
				assert(total_kMers_e1 == haploidEdgeValidLength(level, state, 1));
				if(!(total_kMers_e2 == haploidEdgeValidLength(level, state, 2)))
				{
					cout << "(total_kMers_e2 == haploidEdgeValidLength(level, state, 2))!\n";
					cout << "level: " << level << "\n";
					cout << "diploid state: " << state << "\n";
					cout << "total_kMers_e2: " << total_kMers_e2 << "\n";
					cout << "haploidEdgeValidLength(level, state, 2): " << haploidEdgeValidLength(level, state, 2) << "\n";
				}
				assert(total_kMers_e2 == haploidEdgeValidLength(level, state, 2));

				if(level == 0)
				{
					fw_Viterbi.at(level).at(state) = 0;
					fw_Viterbi_backtrack_int.at(level).at(state) = -1;

					//vector<int> states_same_jump_P;
					//fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
						}

					}

					assert(max_jump_state != -1);

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack_int.at(level).at(state) = max_jump_state;

					//fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					//assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
				}

				if(total_kMers_e1 == 0)
				{
					symbol_emission_e1 = 1;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}

				if(total_kMers_e2 == 0)
				{
					symbol_emission_e2 = 1;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				int missing_edges_e1 = total_edge_length - total_kMers_e1;
				int missing_edges_e2 = total_edge_length - total_kMers_e2;

				int informative_back_edges_e1 = 0;
				if (missing_edges_e1 > 0)
				{
					informative_back_edges_e1 = (missing_edges_e1 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e1;
				}
				int informative_back_edges_e2 = 0;
				if (missing_edges_e2 > 0)
				{
					informative_back_edges_e2 = (missing_edges_e2 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e2;
				}

				double adapted_emission_e1 = symbol_emission_e1;
				double adapted_emission_e2 = symbol_emission_e2;

				if(((double)total_kMers_e1/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e1 = total_kMers_e1;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e1;
					weights.push_back(total_kMers_e1);
					values.push_back(symbol_emission_e1);
					while((current_level > 0) && (got_backtracked_edges_e1 < informative_back_edges_e1))
					{
						running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
						current_level--;

						int running_level_edges = haploidEdgeValidLength(current_level, running_state, 1);
						double running_level_emission = fw_Emission_e1.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e1 + running_level_edges) > informative_back_edges_e1)
						{
							int too_many = (got_backtracked_edges_e1 + running_level_edges) - informative_back_edges_e1;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e1 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					if(weights_sum > 0)
					{
						assert(weights_sum > 0);
						adapted_emission_e1 = 0;
						for(int wI = 0; wI < (int)weights.size(); wI++)
						{
							adapted_emission_e1 += (weights.at(wI)/weights_sum) * values.at(wI);
						}
					}
					else
					{
						adapted_emission_e1 = 1;
					}

				}

				if(((double)total_kMers_e2/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e2 = total_kMers_e2;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e2;
					weights.push_back(total_kMers_e2);
					values.push_back(symbol_emission_e2);
					while((current_level > 0) && (got_backtracked_edges_e2 < informative_back_edges_e2))
					{
						running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
						current_level--;

						int running_level_edges = haploidEdgeValidLength(current_level, running_state, 2);
						double running_level_emission = fw_Emission_e2.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e2 + running_level_edges) > informative_back_edges_e2)
						{
							int too_many = (got_backtracked_edges_e2 + running_level_edges) - informative_back_edges_e2;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e2 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					if(weights_sum > 0)
					{
						assert(weights_sum > 0);
						adapted_emission_e2 = 0;
						for(int wI = 0; wI < (int)weights.size(); wI++)
						{
							adapted_emission_e2 += (weights.at(wI)/weights_sum) * values.at(wI);
						}
					}
					else
					{
						adapted_emission_e2 = 1;
					}
				}

				double symbol_emission = 0.5*adapted_emission_e1 + 0.5*adapted_emission_e2;

				double e1_lengthWeight = 1;
				double e2_lengthWeight = 1;

				double symbol_emission_lengthWeighted  =  0.5*adapted_emission_e1*e1_lengthWeight + 0.5*adapted_emission_e2*e2_lengthWeight;

				assert(symbol_emission >= 0);

				//fw_adaptedEmission.at(level).at(state) = symbol_emission;
				//fw_adaptedEmission_e1.at(level).at(state) = adapted_emission_e1*e1_lengthWeight;
				//fw_adaptedEmission_e2.at(level).at(state) = adapted_emission_e2*e2_lengthWeight;

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission_lengthWeighted*(double)total_edge_length);
				//fw_Emission.at(level).at(state) = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				fw_Emission_e1.at(level).at(state) = symbol_emission_e1;
				fw_Emission_e2.at(level).at(state) = symbol_emission_e2;

				//fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				//fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				//assert(fw.at(level).at(state) >= 0);
			}

		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}


double AlphaHMM::fillForwardBackwardTable_7(map<string, long long> globalEmission)
{
	assert(genotypingMode == 7);

	int levels = diploidStatesByLevel.size();

	fw_Viterbi.resize(levels);
	fw_Emission_e1.resize(levels);
	fw_Emission_e2.resize(levels);
	fw_Viterbi_backtrack_int.resize(levels);

	//fw_adaptedEmission.resize(levels);
	//fw_adaptedEmission_e1.resize(levels);
	//fw_adaptedEmission_e2.resize(levels);
	//fw_expectedMissing_e1.resize(levels);
	//fw_expectedMissing_e2.resize(levels);
	//fw_Coverage_e1.resize(levels);
	//fw_Coverage_e2.resize(levels);
	//fw_Viterbi_backtrack.resize(levels);

	fw_underflow_factor = 1;

	//fw_Emission.resize(levels);
	//fw_e1_length.resize(levels);

	// fw_e1_length_valid.resize(levels);
	// fw_e2_length.resize(levels);
	// fw_e2_length_valid.resize(levels);

	double coverage = estimateCoverage(globalEmission);

	// dMer == duplicate kMer
	kMerNonUniquenessInfo nonUniqueMers = myGraph->kMerNonUniqueness();
	map<int, string> Int2dMer;
	map<string, int> dMer2Int;
	for(set<string>::iterator kMerIt = nonUniqueMers.levelNonUniqueKMers.begin(); kMerIt != nonUniqueMers.levelNonUniqueKMers.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		int number = dMer2Int.size();
		dMer2Int[kMer] = number;
		Int2dMer[number] = kMer;
		assert(nonUniqueMers.kMerMultiplicity.count(kMer) > 0);
		assert(nonUniqueMers.kMerMultiplicity[kMer] < 250);
	}

	vector< map<int, int> > dMerTableLastLevel;
	set<int> kMerConsolidation;
	int lastLeveldMerRecord = -1;

	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map< int, map< int, double > > stateKMerCount_expectedMissing_H0;


	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		if((level % 1000) == 0)
		{
			int startLevelFill = level;
			int stopLevelFill = level + 1000 - 1;

			cout << "Fill state transitions table from " << startLevelFill << " to " << stopLevelFill << "\n" << flush;

			int previousStepFillStop = startLevelFill - 1;
			int previousStepFillStart = previousStepFillStop - 1000 + 1;

			if(previousStepFillStop >= 0)
			{
				if(previousStepFillStart < 0)
				{
					previousStepFillStart = 0;
				}

				// cout << "\t clear state transitions table from " << previousStepFillStart << " to " << previousStepFillStop << "\n" << flush;

				for(int k = previousStepFillStart; k <= previousStepFillStop; k++)
				{
					assert((int)diploidStateTransitions.size() > k);
					assert((int)diploidStateTransitions_Reverse.size() > k);
					diploidStateTransitions.at(k).clear();
					diploidStateTransitions_Reverse.at(k).clear();
				}
			}

			if (stopLevelFill > (levels - 1))
			{
				stopLevelFill = levels - 1;
			}
			int pieces_to_divide = stopLevelFill - startLevelFill + 1;
			assert(pieces_to_divide > 0);
			int chunk_size = pieces_to_divide / CONFIG.threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == CONFIG.threads);
				int thisThread = omp_get_thread_num();
				int firstLevel = thisThread * chunk_size;
				int lastLevel = (thisThread+1) * chunk_size - 1;

				if((thisThread == (CONFIG.threads-1)) && (lastLevel < (pieces_to_divide-1)))
				{
					lastLevel = pieces_to_divide-1;
				}

				assert(firstLevel >= 0);
				if(lastLevel >= 0)
				{
					#pragma omp critical
					{
						// cout << "\t thread " << thisThread << " from " << firstLevel+startLevelFill << " to " << lastLevel+startLevelFill << "\n" << flush;
					}

					for(int k = firstLevel; k <= lastLevel; k++)
					{
						fillStateTransitions(startLevelFill+k);
					}
				}
			}
		}

		int states = diploidStatesByLevel.at(level);

		//fw.at(level).resize(states);

		fw_Viterbi.at(level).resize(states);
		fw_Emission_e1.at(level).resize(states);
		fw_Emission_e2.at(level).resize(states);
		fw_Viterbi_backtrack_int.at(level).resize(states);

		bool thisLevelKeepDMers = false;
		vector< map<int, int> > dMerTableThisLevel;
		if((states < 10000) || (level == 0))
		{
			dMerTableThisLevel.resize(states);
			thisLevelKeepDMers = true;
		}

		//fw_Emission.at(level).resize(states);
		//fw_adaptedEmission.at(level).resize(states);
		//fw_adaptedEmission_e1.at(level).resize(states);
		//fw_adaptedEmission_e2.at(level).resize(states);
		//fw_expectedMissing_e1.at(level).resize(states);
		//fw_expectedMissing_e2.at(level).resize(states);

		//fw_Coverage_e1.at(level).resize(states);
		//fw_Coverage_e2.at(level).resize(states);
		//fw_e1_length.at(level).resize(states);
		//fw_e1_length_valid.at(level).resize(states);
		//fw_e2_length.at(level).resize(states);
		//fw_e2_length_valid.at(level).resize(states);
		//fw_Viterbi_backtrack.at(level).resize(states);


		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map<int, map< int, vector<double> > > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}

					local_edge_length = local_edge_length / 2;
				}


				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);
				assert(total_edge_length == haploid_edge_lengths.at(level));

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				// find out which non-unique kMers we need to calculate this state
				set<string> nonUniqueThisState;
				set<int> backtrackLevels;

				if((nonUniqueMers.level2kMers.count(level) > 0) && (nonUniqueMers.level2kMers[level].size() > 0))
				{
					for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
					{
						int kMerSymbol = kMerIt->first;
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						if(nonUniqueMers.level2kMers[level].count(kMer) > 0)
						{
							nonUniqueThisState.insert(kMer);
							if(lastLeveldMerRecord != (level - 1))
							{
								assert(nonUniqueMers.kMers2Level.count(kMer) > 0);
								set<int> kMerPresentLevels = nonUniqueMers.kMers2Level[kMer];
								for(set<int>::iterator levelIt = kMerPresentLevels.begin(); levelIt != kMerPresentLevels.end(); levelIt++)
								{
									int kMerLevel = *levelIt;
									if(kMerLevel > lastLeveldMerRecord)
									{
										backtrackLevels.insert(kMerLevel);
									}
								}
							}
						}
					}
				}

				// now we find the best previous state to jump from into this state
				double max_jump_P = -1;
				int max_jump_state = -1;
				{
					vector<int> states_same_jump_P;
					set<int> set_jump_P;


					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert((level == 0) || (diploidStateTransitions_Reverse.at(level).at(state).size() > 0));

					int previousStatesMaxIndex = diploidStateTransitions_Reverse.at(level).at(state).size();
					if(level == 0)
					{
						previousStatesMaxIndex = 1;
					}

					for(unsigned int sI2 = 0; sI2 < (unsigned int)previousStatesMaxIndex; sI2++)
					{
						double symbol_emission_e1 = 0;
						double symbol_emission_e2 = 0;
						double estimated_missing_e1 = 0;
						double estimated_missing_e2 = 0;
						double expected_missing_e1 = 0;
						double expected_missing_e2 = 0;
						double total_kMers_e1 = 0;
						double total_kMers_e2 = 0;
						double absolute_coverage_e1 = 0;
						double absolute_coverage_e2 = 0;

						double viterbiJumpP = 0;
						map<string, int> nonUniqueThisStateCount;
						int from = -1;

						if(level > 0)
						{
							stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
							from = jumpState.id;

							/* calculate emission probabilities conditional on previous path */

							// we have a set of non-unique kMers that we are interested in
							// now, either we have a "used kMers" table at the previous level,
							// or we need to carry out some backtracking


							{
								int runningDiploidState = from;
								int runningLevel = level - 1;
								while(runningLevel != lastLeveldMerRecord)
								{
									if(backtrackLevels.count(runningLevel) > 0)
									{
										int haploidStates_runningLevel = haploidStatesByLevel.at(runningLevel).size();
										int s1_runningLevel = runningDiploidState % haploidStates_runningLevel;
										int s2_runningLevel = runningDiploidState / haploidStates_runningLevel;

										Edge* e1_runningLevel = haploidStatesByLevel.at(runningLevel).at(s1_runningLevel).e;
										Edge* e2_runningLevel = haploidStatesByLevel.at(runningLevel).at(s2_runningLevel).e;

										for(set<string>::iterator kMerStringIt = nonUniqueThisState.begin(); kMerStringIt != nonUniqueThisState.end(); kMerStringIt++)
										{
											string kMerString = *kMerStringIt;
											int combined_count = 0;
											if(nonUniqueMers.kMers2Level[kMerString].count(runningLevel) > 0)
											{
												assert(dMer2Int.count(kMerString) > 0);
												int kMerSymbol = myGraph->CODE.doCode(e1_runningLevel->locus_id, kMerString);
												int e1_count = (e1_runningLevel->multiEmission.count(kMerSymbol) > 0) ? e1_runningLevel->multiEmission[kMerSymbol]  : 0;
												int e2_count = (e2_runningLevel->multiEmission.count(kMerSymbol) > 0) ? e2_runningLevel->multiEmission[kMerSymbol]  : 0;
												combined_count = e1_count + e2_count;
											}
											assert(combined_count <= (2*nonUniqueMers.kMerMultiplicity[kMerString]));

											if(nonUniqueThisStateCount.count(kMerString) == 0)
											{
												nonUniqueThisStateCount[kMerString] = 0;
											}
											nonUniqueThisStateCount[kMerString] += combined_count;

											assert(nonUniqueThisStateCount[kMerString] <= (2 * nonUniqueMers.kMerMultiplicity[kMerString]));

										}
									}

									runningDiploidState = fw_Viterbi_backtrack_int.at(runningLevel).at(runningDiploidState);
									runningLevel--;
								}

								for(set<string>::iterator kMerStringIt = nonUniqueThisState.begin(); kMerStringIt != nonUniqueThisState.end(); kMerStringIt++)
								{
									string kMerString = *kMerStringIt;
									assert(dMer2Int.count(kMerString) > 0);
									int number = dMer2Int[kMerString];

									if(nonUniqueThisStateCount.count(kMerString) == 0)
									{
										nonUniqueThisStateCount[kMerString] = 0;
									}

									if(dMerTableLastLevel.at(runningDiploidState).count(number) > 0)
									{
										nonUniqueThisStateCount[kMerString] += dMerTableLastLevel.at(runningDiploidState)[number];
										if(!(nonUniqueThisStateCount[kMerString] <= (2 * nonUniqueMers.kMerMultiplicity[kMerString])))
										{
											cout << "!(nonUniqueThisStateCount[kMerString] <= 2* nonUniqueMers.kMerMultiplicity[kMerString])\n";
											cout << kMerString << "\n";
											cout << "Multiplicity: " <<  nonUniqueMers.kMerMultiplicity[kMerString] << "\n";
											cout << "nonUniqueThisStateCount[kMerString]: " <<  nonUniqueThisStateCount[kMerString] << "\n";
											cout << "dMerTableLastLevel.at(runningDiploidState)[number]: " << dMerTableLastLevel.at(runningDiploidState)[number] << "\n";
										}
										assert(nonUniqueThisStateCount[kMerString] <= (2* nonUniqueMers.kMerMultiplicity[kMerString]));
									}
								}
							}
						}


						for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
						{
							int kMerSymbol = kMerIt->first;
							int kMerCount = kMerIt->second;

							double thisKMer_estimated_missing = 0;
							double this_kMer_expected_missing = 0;

							int e1_count = 0;
							int e2_count = 0;
							if(e1->multiEmission.count(kMerSymbol) > 0)
							{
								e1_count = e1->multiEmission[kMerSymbol];
							}
							if(e2->multiEmission.count(kMerSymbol) > 0)
							{
								e2_count = e2->multiEmission[kMerSymbol];
							}


							assert((e1_count + e2_count) == kMerCount);

							total_kMers_e1 += e1_count;
							total_kMers_e2 += e2_count;

							string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);

							int previousStateCount = 0;
							if(nonUniqueThisState.count(kMer) > 0)
							{
								assert(nonUniqueThisStateCount.count(kMer) > 0);
								previousStateCount = nonUniqueThisStateCount[kMer];
								assert(previousStateCount <= (2 * nonUniqueMers.kMerMultiplicity[kMer]));
								if(thisLevelKeepDMers == true)
								{
									assert(dMer2Int.count(kMer) > 0);
									int dMerNumber = dMer2Int[kMer];
									dMerTableThisLevel.at(state)[dMerNumber] = kMerCount;
								}
							}

							assert(globalEmission.count(kMer) > 0);
							int observedNumber = globalEmission[kMer];
							absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
							absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

							if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0) && (cache_estimated_expected_Missing[kMerSymbol][kMerCount].count(previousStateCount) > 0))
							{
								thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount][previousStateCount].at(0);
								this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount][previousStateCount].at(1);
							}
							else
							{
								string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
								assert(kMer != "*");
								assert(kMer != "_");

								assert(globalEmission.count(kMer) > 0);
								int observedNumber = globalEmission[kMer];

								if(kMerCount > 0)
								{

									double p_evenLower = 0;
									for(int possibleUnderlyingCount = 0; possibleUnderlyingCount <= previousStateCount; possibleUnderlyingCount++)
									{
										if(possibleUnderlyingCount >= (int)observedXunderlying.at(observedNumber).size() )
										{
											cout << "kMer " << kMer << " already observed " << previousStateCount << " times, off-limit! [" << nonUniqueMers.kMerMultiplicity[kMer] << "]\n";
										}
										else
										{
											p_evenLower += observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
										}
									}
									assert(p_evenLower >= 0);
									assert(p_evenLower <= (1 + epsilon));

									int skipFirstState = 0;
									if(previousStateCount != 0)
									{
										int missing = kMerCount;
										thisKMer_estimated_missing += p_evenLower * (double)missing;
										skipFirstState = 1;
									}

									for(int possibleUnderlyingCount = previousStateCount+skipFirstState; possibleUnderlyingCount < previousStateCount+kMerCount; possibleUnderlyingCount++)
									{
										assert(possibleUnderlyingCount < previousStateCount+kMerCount);

										int missing = (kMerCount+previousStateCount) - possibleUnderlyingCount;
										double p;
										if(possibleUnderlyingCount >= (int)observedXunderlying.at(observedNumber).size() )
										{
											p = 0;
										}
										else
										{
											p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
										}

										assert(p >= 0);
										assert(p <= (1 + epsilon));

										thisKMer_estimated_missing += p * (double)missing;

										assert(thisKMer_estimated_missing >= 0);
										assert(thisKMer_estimated_missing <= kMerCount);
									}

									if((stateKMerCount_expectedMissing_H0.count(kMerCount) == 0) && (stateKMerCount_expectedMissing_H0[kMerCount].count(previousStateCount) == 0))
									{
										double expected_missingness = 0;
										for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
										{
											double thisObservation_expectedMissing = 0;

											double expected_missing_coverage = coverage * 1.0;
											double lambda = (double)kMerCount*expected_missing_coverage;
											poisson_up poisson(lambda);
											double thisObservation_weight = pdf(poisson, potentialObserved);

											double p_evenLower = 0;
											for(int possibleUnderlyingCount = 0; possibleUnderlyingCount <= previousStateCount; possibleUnderlyingCount++)
											{
												if(possibleUnderlyingCount <  (int)observedXunderlying.at(potentialObserved).size())
												{
													p_evenLower += observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
												}
											}

											if(!((p_evenLower >= 0) && (p_evenLower <= (1 + epsilon))))
											{
												cout << "p_evenLower: " << p_evenLower << "\n";
												cout << "potentialObserved: " << potentialObserved << "\n";
												cout << "previousStateCount: " << previousStateCount << "\n";

												double t_p_evenLower = 0;
												for(int possibleUnderlyingCount = 0; possibleUnderlyingCount <= previousStateCount; possibleUnderlyingCount++)
												{
													cout << "possibleUnderlyingCount: " << possibleUnderlyingCount << "\n";
													cout << "\t before t_p_evenLower: " << t_p_evenLower << "\n";
													if(possibleUnderlyingCount <  (int)observedXunderlying.at(potentialObserved).size())
													{
														cout << "\t\t add " << observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount) << "\n";
														t_p_evenLower += observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
													}
													cout << "\t old t_p_evenLower: " << t_p_evenLower << "\n";
												}
											}
											assert(p_evenLower >= 0);
											assert(p_evenLower <= (1 + epsilon));

											int skipFirstState = 0;
											if(previousStateCount != 0)
											{
												int missing = kMerCount;
												thisObservation_expectedMissing += p_evenLower * (double)missing;
												skipFirstState = 1;
											}

											for(int possibleUnderlyingCount = previousStateCount + skipFirstState; possibleUnderlyingCount < previousStateCount+kMerCount; possibleUnderlyingCount++)
											{
												int missing = (kMerCount+previousStateCount) - possibleUnderlyingCount;
												double p;
												if(possibleUnderlyingCount >= (int)observedXunderlying.at(potentialObserved).size() )
												{
													p = 0;
												}
												else
												{
													p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
												}
												thisObservation_expectedMissing += p * (double)missing;
											}

											expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
											assert(expected_missingness >= 0);
										}

										#pragma omp critical
										{
											stateKMerCount_expectedMissing_H0[kMerCount][previousStateCount] = 	expected_missingness;
											this_kMer_expected_missing = expected_missingness;
										}
									}
									else
									{
										this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount][previousStateCount];
									}

								}

								vector<double> forCache;
								forCache.push_back(thisKMer_estimated_missing);
								forCache.push_back(this_kMer_expected_missing);
								cache_estimated_expected_Missing[kMerSymbol][kMerCount][previousStateCount] = forCache;
							}

							assert(kMerCount > 0);
							assert(this_kMer_expected_missing >= 0);

							assert(thisKMer_estimated_missing <= kMerCount);

							estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
							expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

							estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
							expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);
						}

						assert(expected_missing_e1 >= 0);
						assert(expected_missing_e2 >= 0);

						estimated_missing_e1 -= expected_missing_e1;
						estimated_missing_e2 -= expected_missing_e2;


						if(!(total_kMers_e1 == haploidEdgeValidLength(level, state, 1)))
						{
							cout << "(total_kMers_e1 == haploidEdgeValidLength(level, state, 1))!\n";
							cout << "level: " << level << "\n";
							cout << "diploid state: " << state << "\n";
							cout << "total_kMers_e1: " << total_kMers_e1 << "\n";
							cout << "haploidEdgeValidLength(level, state, 1): " << haploidEdgeValidLength(level, state, 1) << "\n";
						}
						assert(total_kMers_e1 == haploidEdgeValidLength(level, state, 1));
						if(!(total_kMers_e2 == haploidEdgeValidLength(level, state, 2)))
						{
							cout << "(total_kMers_e2 == haploidEdgeValidLength(level, state, 2))!\n";
							cout << "level: " << level << "\n";
							cout << "diploid state: " << state << "\n";
							cout << "total_kMers_e2: " << total_kMers_e2 << "\n";
							cout << "haploidEdgeValidLength(level, state, 2): " << haploidEdgeValidLength(level, state, 2) << "\n";
						}
						assert(total_kMers_e2 == haploidEdgeValidLength(level, state, 2));

						if(total_kMers_e1 == 0)
						{
							symbol_emission_e1 = 1;
						}
						else
						{
							assert((estimated_missing_e1/total_kMers_e1) <= 1);
							symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
						}
						assert(symbol_emission_e1 >= 0);


						if(total_kMers_e2 == 0)
						{
							symbol_emission_e2 = 1;
						}
						else
						{
							assert((estimated_missing_e2/total_kMers_e2) <= 1);
							symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
						}
						assert(symbol_emission_e2 >= 0);

						int missing_edges_e1 = total_edge_length - total_kMers_e1;
						int missing_edges_e2 = total_edge_length - total_kMers_e2;

						int informative_back_edges_e1 = 0;
						if (missing_edges_e1 > 0)
						{
							informative_back_edges_e1 = (missing_edges_e1 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e1;
						}
						int informative_back_edges_e2 = 0;
						if (missing_edges_e2 > 0)
						{
							informative_back_edges_e2 = (missing_edges_e2 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e2;
						}

						double adapted_emission_e1 = symbol_emission_e1;
						double adapted_emission_e2 = symbol_emission_e2;

						if(((double)total_kMers_e1/(double)total_edge_length) < 0.3)
						{
							int current_level = level;
							int got_backtracked_edges_e1 = total_kMers_e1;
							int running_state = state;

							vector<double> weights;
							vector<double> values;
							double weights_sum = total_kMers_e1;
							weights.push_back(total_kMers_e1);
							values.push_back(symbol_emission_e1);
							while((current_level > 0) && (got_backtracked_edges_e1 < informative_back_edges_e1))
							{
								running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
								current_level--;

								int running_level_edges = haploidEdgeValidLength(current_level, running_state, 1);
								double running_level_emission = fw_Emission_e1.at(current_level).at(running_state);

								int take_edges;
								if((got_backtracked_edges_e1 + running_level_edges) > informative_back_edges_e1)
								{
									int too_many = (got_backtracked_edges_e1 + running_level_edges) - informative_back_edges_e1;
									take_edges = running_level_edges - too_many;
								}
								else
								{
									take_edges = running_level_edges;
								}

								if(take_edges > 0)
								{
									weights.push_back(take_edges);
									weights_sum += take_edges;
									values.push_back(running_level_emission);
									got_backtracked_edges_e1 += take_edges;
								}
							}

							assert(weights.size() == values.size());
							if(weights_sum > 0)
							{
								assert(weights_sum > 0);
								adapted_emission_e1 = 0;
								for(int wI = 0; wI < (int)weights.size(); wI++)
								{
									assert(values.at(wI) >= 0);
									assert((weights.at(wI)/weights_sum) >= 0);

									adapted_emission_e1 += (weights.at(wI)/weights_sum) * values.at(wI);
								}
							}
							else
							{
								adapted_emission_e1 = 1;
							}

						}

						if(((double)total_kMers_e2/(double)total_edge_length) < 0.3)
						{
							int current_level = level;
							int got_backtracked_edges_e2 = total_kMers_e2;
							int running_state = state;

							vector<double> weights;
							vector<double> values;
							double weights_sum = total_kMers_e2;
							weights.push_back(total_kMers_e2);
							values.push_back(symbol_emission_e2);
							while((current_level > 0) && (got_backtracked_edges_e2 < informative_back_edges_e2))
							{
								running_state = fw_Viterbi_backtrack_int.at(current_level).at(running_state);
								current_level--;

								int running_level_edges = haploidEdgeValidLength(current_level, running_state, 2);
								double running_level_emission = fw_Emission_e2.at(current_level).at(running_state);

								int take_edges;
								if((got_backtracked_edges_e2 + running_level_edges) > informative_back_edges_e2)
								{
									int too_many = (got_backtracked_edges_e2 + running_level_edges) - informative_back_edges_e2;
									take_edges = running_level_edges - too_many;
								}
								else
								{
									take_edges = running_level_edges;
								}

								if(take_edges > 0)
								{
									weights.push_back(take_edges);
									weights_sum += take_edges;
									values.push_back(running_level_emission);
									got_backtracked_edges_e2 += take_edges;
								}
							}

							assert(weights.size() == values.size());
							if(weights_sum > 0)
							{
								assert(weights_sum > 0);
								adapted_emission_e2 = 0;
								for(int wI = 0; wI < (int)weights.size(); wI++)
								{
									assert(values.at(wI) >= 0);
									assert((weights.at(wI)/weights_sum) >= 0);

									adapted_emission_e2 += (weights.at(wI)/weights_sum) * values.at(wI);
								}
							}
							else
							{
								adapted_emission_e2 = 1;
							}
						}

						assert(adapted_emission_e1 >= 0);
						assert(adapted_emission_e2 >= 0);

						double symbol_emission = 0.5*adapted_emission_e1 + 0.5*adapted_emission_e2;

						double e1_lengthWeight = 1 + 0.10 * ((double)total_kMers_e1/(double)total_edge_length);
						double e2_lengthWeight = 1 + 0.10 * ((double)total_kMers_e2/(double)total_edge_length);

						assert(e1_lengthWeight >= 1);
						assert(e2_lengthWeight >= 1);

						double symbol_emission_lengthWeighted  =  0.5*adapted_emission_e1*e1_lengthWeight + 0.5*adapted_emission_e2*e2_lengthWeight;

						assert(symbol_emission_lengthWeighted >= 0);
						assert(symbol_emission >= 0);
						assert(total_edge_length >= 0);

						/* end calculation emission probabilties */

						viterbiJumpP = symbol_emission_lengthWeighted*(double)total_edge_length;
						if(level > 0)
						{
							viterbiJumpP = fw_Viterbi.at(level-1).at(from) + viterbiJumpP;
						}

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							assert(viterbiJumpP >= 0);

							fw_Viterbi.at(level).at(state) = viterbiJumpP;
							fw_Emission_e1.at(level).at(state) = symbol_emission_e1;
							fw_Emission_e2.at(level).at(state) = symbol_emission_e2;
							fw_Viterbi_backtrack_int.at(level).at(state) = max_jump_state;
						}
					}
				}

				// we have identified the best state to go from into this state and
				// the optimality value for this state

				if((thisLevelKeepDMers == true) && (level > 0))
				{
					// if we want to build a new "used kMers" table for this state,
					// we now need to start merging with previous states and the last table

					set<string> nonUniquePastStates;

					if(lastLeveldMerRecord < (level - 1))
					{
						for(int lI = lastLeveldMerRecord+1; lI < level; lI++)
						{
							if(nonUniqueMers.level2kMers[lI].size() > 0)
							{
								nonUniquePastStates.insert(nonUniqueMers.level2kMers[lI].begin(), nonUniqueMers.level2kMers[lI].end());
							}
						}
					}


					int runningDiploidState = max_jump_state;
					int runningLevel = level - 1;
					while(runningLevel != lastLeveldMerRecord)
					{
						if((nonUniqueMers.level2kMers.count(runningLevel) > 0) && (nonUniqueMers.level2kMers[runningLevel].size() > 0))
						{
							int haploidStates_runningLevel = haploidStatesByLevel.at(runningLevel).size();
							int s1_runningLevel = runningDiploidState % haploidStates_runningLevel;
							int s2_runningLevel = runningDiploidState / haploidStates_runningLevel;

							Edge* e1_runningLevel = haploidStatesByLevel.at(runningLevel).at(s1_runningLevel).e;
							Edge* e2_runningLevel = haploidStatesByLevel.at(runningLevel).at(s2_runningLevel).e;

							for(set<string>::iterator kMerStringIt = nonUniquePastStates.begin(); kMerStringIt != nonUniquePastStates.end(); kMerStringIt++)
							{
								string kMerString = *kMerStringIt;
								int combined_count = 0;
								if(nonUniqueMers.kMers2Level[kMerString].count(runningLevel) > 0)
								{
									assert(dMer2Int.count(kMerString) > 0);
									int kMerSymbol = myGraph->CODE.doCode(e1_runningLevel->locus_id, kMerString);
									int e1_count = (e1_runningLevel->multiEmission.count(kMerSymbol) > 0) ? e1_runningLevel->multiEmission[kMerSymbol]  : 0;
									int e2_count = (e2_runningLevel->multiEmission.count(kMerSymbol) > 0) ? e2_runningLevel->multiEmission[kMerSymbol]  : 0;
									combined_count = e1_count + e2_count;
								}
								assert(combined_count <= (2 * nonUniqueMers.kMerMultiplicity[kMerString]));

								assert(dMer2Int.count(kMerString) > 0);
								int kMerNum = dMer2Int[kMerString];
								if(dMerTableThisLevel.at(state).count(kMerNum) == 0)
								{
									dMerTableThisLevel.at(state)[kMerNum] = 0;
								}

								assert(dMerTableThisLevel.at(state)[kMerNum] <= (2 *nonUniqueMers.kMerMultiplicity[kMerString]));
								dMerTableThisLevel.at(state)[kMerNum] += combined_count;
								assert(dMerTableThisLevel.at(state)[kMerNum] <= (2 * nonUniqueMers.kMerMultiplicity[kMerString]));

							}
						}

						runningDiploidState = fw_Viterbi_backtrack_int.at(runningLevel).at(runningDiploidState);
						runningLevel--;
					}

					for(map<int, int>::iterator tableNumCountIt = dMerTableLastLevel.at(runningDiploidState).begin(); tableNumCountIt != dMerTableLastLevel.at(runningDiploidState).end(); tableNumCountIt++)
					{
						int kMerNumber = tableNumCountIt->first;
						int count = tableNumCountIt->second;

						if(dMerTableThisLevel.at(state).count(kMerNumber)== 0)
						{
							dMerTableThisLevel.at(state)[kMerNumber] = 0;
						}

						assert(count == dMerTableLastLevel.at(runningDiploidState)[kMerNumber]);
						assert(dMerTableThisLevel.at(state)[kMerNumber] <= (2 * nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]]));
						assert(dMerTableLastLevel.at(runningDiploidState)[kMerNumber] <= (2 * nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]]));

						int before_value = dMerTableThisLevel.at(state)[kMerNumber];
						dMerTableThisLevel.at(state)[kMerNumber] += dMerTableLastLevel.at(runningDiploidState)[kMerNumber];
						assert(Int2dMer.count(kMerNumber) > 0);
						if(!(dMerTableThisLevel.at(state)[kMerNumber] <= (2 * nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]])))
						{
							cout << "!(dMerTableThisLevel.at(state)[kMerNumber] <= (2 * nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]]))\n";
							cout << "kMer: " << Int2dMer[kMerNumber] << "\n";
							cout << "kMerNumber: " << kMerNumber << "\n";
							cout << "dMerTableThisLevel.at(state)[kMerNumber]: " << dMerTableThisLevel.at(state)[kMerNumber] << "\n";
							cout << "nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]]: " << nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]] << "\n";
							cout << "before_value: " << before_value << "\n";
							cout << "dMerTableLastLevel.at(runningDiploidState)[kMerNumber]: " << dMerTableLastLevel.at(runningDiploidState)[kMerNumber] << "\n\n";
							cout << "lastLeveldMerRecord: " << lastLeveldMerRecord << "\n";
							cout << "level: " << level << "\n";
							cout << "Levels at which this kMer occurs:\n";
							vector<int> kMerLevels(nonUniqueMers.kMers2Level[Int2dMer[kMerNumber]].begin(), nonUniqueMers.kMers2Level[Int2dMer[kMerNumber]].end());
							vector<string> kMerLevelStrings = Utilities::ItoStr(kMerLevels);
							cout << "\t" << Utilities::join(kMerLevelStrings, ", ") << "\n\n\n";

						}
						assert(dMerTableThisLevel.at(state)[kMerNumber] <= (2 * nonUniqueMers.kMerMultiplicity[Int2dMer[kMerNumber]]));
					}

					// remove kMers which do not appear in future states or which have count 0

					set<int> removeKMer;
					for(map<int, int>::iterator usedKmerIt = dMerTableThisLevel.at(state).begin(); usedKmerIt != dMerTableThisLevel.at(state).end(); usedKmerIt++)
					{
						int kMerNumber = usedKmerIt->first;
						int count = usedKmerIt->second;

						assert(Int2dMer.count(kMerNumber));
						string kMerString = Int2dMer[kMerNumber];
						assert(nonUniqueMers.kMers2Level.count(kMerString) > 0);
						set<int> presentAtLevels = nonUniqueMers.kMers2Level[kMerString];
						bool oneAhead = false;
						for(set<int>::iterator levelIt = presentAtLevels.begin(); levelIt != presentAtLevels.end(); levelIt++)
						{
							int kMerAtLevel = *levelIt;
							if(kMerAtLevel > level)
							{
								oneAhead = true;
							}
						}
						if((! oneAhead) || (count == 0))
						{
							removeKMer.insert(kMerNumber);
						}
					}

					for(set<int>::iterator kMerNumIt = removeKMer.begin(); kMerNumIt != removeKMer.end(); kMerNumIt++)
					{
						int kMerNumber = *kMerNumIt;
						dMerTableThisLevel.at(state).erase(kMerNumber);
					}
				}
			}
		}

		if(thisLevelKeepDMers == true)
		{
			dMerTableLastLevel = dMerTableThisLevel;
			lastLeveldMerRecord = level;
		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}

double AlphaHMM::fillForwardBackwardTable_2(map<string, long long> globalEmission)
{
	assert(genotypingMode == 2);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);
	fw_Emission.resize(levels);
	fw_Emission_e1.resize(levels);
	fw_Emission_e2.resize(levels);
	fw_expectedMissing_e1.resize(levels);
	fw_expectedMissing_e2.resize(levels);
	fw_Coverage_e1.resize(levels);
	fw_Coverage_e2.resize(levels);
	fw_Viterbi_backtrack.resize(levels);
	fw_underflow_factor = 1;
	fw_e1_length.resize(levels);
	fw_e1_length_valid.resize(levels);
	fw_e2_length.resize(levels);
	fw_e2_length_valid.resize(levels);

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Emission.at(level).resize(states);
		fw_Emission_e1.at(level).resize(states);
		fw_Emission_e2.at(level).resize(states);
		fw_expectedMissing_e1.at(level).resize(states);
		fw_expectedMissing_e2.at(level).resize(states);
		fw_Coverage_e1.at(level).resize(states);
		fw_Coverage_e2.at(level).resize(states);
		fw_e1_length.at(level).resize(states);
		fw_e1_length_valid.at(level).resize(states);
		fw_e2_length.at(level).resize(states);
		fw_e2_length_valid.at(level).resize(states);
		fw_Viterbi_backtrack.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				bool dPrint = false;
				if(((e1->label.substr(0, 10) == "HLAA*02:01") && (e2->label.substr(0, 10) == "HLAA*23:01" )) || ((e2->label.substr(0, 10) == "HLAA*02:01") && (e1->label.substr(0, 10) == "HLAA*23:01" )))
				{
					if(s1 <= s2)
					{
						dPrint = true;
					}
				}

				if(dPrint)
				{
					subst_cout << "LABELS " << e1->label << "/" << e2->label << "\n";
				}

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}

					local_edge_length = local_edge_length / 2;

				}

				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				if(dPrint)
				{
					subst_cout << "\t gaps: " << combinedEdgeEmission[gap_symbol] << "\n";
					subst_cout << "\t stars: " << combinedEdgeEmission[star_symbol] << "\n";

				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{

							for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
							{
								int missing = kMerCount - possibleUnderlyingCount;
								double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
								thisKMer_estimated_missing += p * (double)missing;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double expected_missing_coverage = coverage * 1.0;
									double lambda = (double)kMerCount*expected_missing_coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);

					if(dPrint)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						subst_cout << "\t\t symbol " << kMerSymbol << " " << kMerCount << " [required] "  << globalEmission[kMer] << " [copies in reads] " << thisKMer_estimated_missing  << " [estimated missing] " << this_kMer_expected_missing << " [expected missing] " << "                 " << kMer << "\n";
					}
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				absolute_coverage_e1 = absolute_coverage_e1 / (double)total_kMers_e1;
				absolute_coverage_e2 = absolute_coverage_e2 / (double)total_kMers_e2;
				fw_Coverage_e1.at(level).at(state) = absolute_coverage_e1;
				fw_Coverage_e2.at(level).at(state) = absolute_coverage_e2;

				fw_e1_length.at(level).at(state) = total_edge_length;
				fw_e2_length.at(level).at(state) = total_edge_length;
				fw_e1_length_valid.at(level).at(state) = total_kMers_e1;
				fw_e2_length_valid.at(level).at(state) = total_kMers_e2;

				if(dPrint)
				{
					subst_cout << "\t estimated_missing_e1: " << estimated_missing_e1 << "\t total_kMers_e1: " << total_kMers_e1 << " => symbol_emission_e1: " << symbol_emission_e1 << "\n";
					subst_cout << "\t estimated_missing_e2: " << estimated_missing_e2 << "\t total_kMers_e2: " << total_kMers_e2 << " => symbol_emission_e2: " << symbol_emission_e2 << "\n";
				}


				if(dPrint)
				{
					#pragma omp critical
					{
						cout << subst_cout.str() << flush;
					}
				}

				if(level == 0)
				{
					fw.at(level).at(state) = 0;
					fw_Viterbi.at(level).at(state) = 0;

					vector<int> states_same_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					fw.at(level).at(state) = 0;

					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							states_same_jump_P.clear();
							states_same_jump_P.push_back(max_jump_state);
							assert(set_jump_P.count(max_jump_state) == 0);
							set_jump_P.insert(max_jump_state);
						}
						else if(viterbiJumpP == max_jump_P)
						{
							states_same_jump_P.push_back(from);
							assert(set_jump_P.count(from) == 0);
							set_jump_P.insert(from);
						}
					}

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
					assert(fw_Viterbi_backtrack.at(level).size() > 0);

				}

				if(total_kMers_e1 == 0)
				{
					assert(level > 0);
					symbol_emission_e1 = fw_Emission_e1.at(level-1).at(fw_Viterbi_backtrack.at(level).at(state).at(0));				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}

				if(total_kMers_e2 == 0)
				{
					assert(level > 0);
					symbol_emission_e2 = fw_Emission_e2.at(level-1).at(fw_Viterbi_backtrack.at(level).at(state).at(0));
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				double symbol_emission  = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				assert(symbol_emission >= 0);

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission*(double)total_edge_length);
				fw_Emission.at(level).at(state) = symbol_emission;
				fw_Emission_e1.at(level).at(state) = symbol_emission_e1;
				fw_Emission_e2.at(level).at(state) = symbol_emission_e2;
				fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				assert(fw.at(level).at(state) >= 0);
			}
		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}

double AlphaHMM::fillForwardBackwardTable_3(map<string, long long> globalEmission)
{
	assert(genotypingMode == 3);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);

	fw_Viterbi_e1.resize(levels);
	fw_Viterbi_e2.resize(levels);

	fw_Emission.resize(levels);
	fw_expectedMissing_e1.resize(levels);
	fw_expectedMissing_e2.resize(levels);
	fw_Coverage_e1.resize(levels);
	fw_Coverage_e2.resize(levels);
	fw_Viterbi_backtrack.resize(levels);
	fw_underflow_factor = 1;

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Viterbi_e1.at(level).resize(states);
		fw_Viterbi_e2.at(level).resize(states);

		fw_Emission.at(level).resize(states);
		fw_expectedMissing_e1.at(level).resize(states);
		fw_expectedMissing_e2.at(level).resize(states);
		fw_Coverage_e1.at(level).resize(states);
		fw_Coverage_e2.at(level).resize(states);

		fw_Viterbi_backtrack.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				bool dPrint = false;
				if(((e1->label.substr(0, 10) == "HLAA*02:01") && (e2->label.substr(0, 10) == "HLAA*23:01" )) || ((e2->label.substr(0, 10) == "HLAA*02:01") && (e1->label.substr(0, 10) == "HLAA*23:01" )))
				{
					if(s1 <= s2)
					{
						dPrint = true;
					}
				}

				if(dPrint)
				{
					subst_cout << "LABELS " << e1->label << "/" << e2->label << "\n";
				}

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}
					local_edge_length = local_edge_length / 2;

				}

				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				if(dPrint)
				{
					subst_cout << "\t gaps: " << combinedEdgeEmission[gap_symbol] << "\n";
					subst_cout << "\t stars: " << combinedEdgeEmission[star_symbol] << "\n";

				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{

							for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
							{
								int missing = kMerCount - possibleUnderlyingCount;
								double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
								thisKMer_estimated_missing += p * (double)missing;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double lambda = (double)kMerCount*coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);

					if(dPrint)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						subst_cout << "\t\t symbol " << kMerSymbol << " " << kMerCount << " [required] "  << globalEmission[kMer] << " [copies in reads] " << thisKMer_estimated_missing  << " [estimated missing] " << this_kMer_expected_missing << " [expected missing] " << "                 " << kMer << "\n";
					}
				}


				absolute_coverage_e1 = absolute_coverage_e1 / (double)total_kMers_e1;
				absolute_coverage_e2 = absolute_coverage_e2 / (double)total_kMers_e2;
				fw_Coverage_e1.at(level).at(state) = absolute_coverage_e1;
				fw_Coverage_e2.at(level).at(state) = absolute_coverage_e2;

				if(dPrint)
				{
					subst_cout << "\t estimated_missing_e1: " << estimated_missing_e1 << "\t total_kMers_e1: " << total_kMers_e1 << " => symbol_emission_e1: " << symbol_emission_e1 << "\n";
					subst_cout << "\t estimated_missing_e2: " << estimated_missing_e2 << "\t total_kMers_e2: " << total_kMers_e2 << " => symbol_emission_e2: " << symbol_emission_e2 << "\n";
				}


				if(dPrint)
				{
					#pragma omp critical
					{
						cout << subst_cout.str() << flush;
					}
				}

				if(level == 0)
				{
					fw.at(level).at(state) = 0;
					fw_Viterbi.at(level).at(state) = 0;
					fw_Viterbi_e1.at(level).at(state) = 0;
					fw_Viterbi_e2.at(level).at(state) = 0;

					vector<int> states_same_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					fw.at(level).at(state) = 0;

					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							states_same_jump_P.clear();
							states_same_jump_P.push_back(max_jump_state);
							assert(set_jump_P.count(max_jump_state) == 0);
							set_jump_P.insert(max_jump_state);
						}
						else if(viterbiJumpP == max_jump_P)
						{
							states_same_jump_P.push_back(from);
							assert(set_jump_P.count(from) == 0);
							set_jump_P.insert(from);
						}
					}

					fw_Viterbi_e1.at(level).at(state) = fw_Viterbi_e1.at(level-1).at(states_same_jump_P.at(0));
					fw_Viterbi_e2.at(level).at(state) = fw_Viterbi_e2.at(level-1).at(states_same_jump_P.at(0));

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
					assert(fw_Viterbi_backtrack.at(level).size() > 0);

				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				if(total_kMers_e1 == 0)
				{
					assert(level > 0);
					double current_optimality = fw_Viterbi_e1.at(level).at(state)/((grand_total_edge_length - total_edge_length)/(double)2);
					assert(current_optimality >= 0);
					// assert(current_optimality <= 1);
					symbol_emission_e1 = current_optimality;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}

				if(total_kMers_e2 == 0)
				{
					assert(level > 0);
					double current_optimality = fw_Viterbi_e2.at(level).at(state)/((grand_total_edge_length - total_edge_length)/(double)2);
					assert(current_optimality >= 0);
					// assert(current_optimality <= 1);
					symbol_emission_e2 = current_optimality;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				double symbol_emission = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				assert(symbol_emission >= 0);

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission*(double)total_edge_length);

				fw_Viterbi_e1.at(level).at(state) = fw_Viterbi_e1.at(level).at(state) + (symbol_emission_e1*(double)(total_edge_length/(double)2));
				fw_Viterbi_e2.at(level).at(state) = fw_Viterbi_e2.at(level).at(state) + (symbol_emission_e2*(double)(total_edge_length/(double)2));
				assert((abs(fw_Viterbi_e1.at(level).at(state) + fw_Viterbi_e2.at(level).at(state)) - fw_Viterbi.at(level).at(state)) < epsilon);

				fw_Emission.at(level).at(state) = (symbol_emission*(double)total_edge_length);
				fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				assert(fw.at(level).at(state) >= 0);
			}
		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}


double AlphaHMM::fillForwardBackwardTable_4(map<string, long long> globalEmission)
{
	assert(genotypingMode == 4);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);
	fw_Emission.resize(levels);
	fw_Emission_e1.resize(levels);
	fw_Emission_e2.resize(levels);
	fw_adaptedEmission.resize(levels);
	fw_adaptedEmission_e1.resize(levels);
	fw_adaptedEmission_e2.resize(levels);
	fw_expectedMissing_e1.resize(levels);
	fw_expectedMissing_e2.resize(levels);
	fw_Coverage_e1.resize(levels);
	fw_Coverage_e2.resize(levels);
	fw_Viterbi_backtrack.resize(levels);
	fw_underflow_factor = 1;
	fw_e1_length.resize(levels);
	fw_e1_length_valid.resize(levels);
	fw_e2_length.resize(levels);
	fw_e2_length_valid.resize(levels);

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Emission.at(level).resize(states);
		fw_Emission_e1.at(level).resize(states);
		fw_Emission_e2.at(level).resize(states);
		fw_adaptedEmission.at(level).resize(states);
		fw_adaptedEmission_e1.at(level).resize(states);
		fw_adaptedEmission_e2.at(level).resize(states);
		fw_expectedMissing_e1.at(level).resize(states);
		fw_expectedMissing_e2.at(level).resize(states);
		fw_Coverage_e1.at(level).resize(states);
		fw_Coverage_e2.at(level).resize(states);
		fw_e1_length.at(level).resize(states);
		fw_e1_length_valid.at(level).resize(states);
		fw_e2_length.at(level).resize(states);
		fw_e2_length_valid.at(level).resize(states);
		fw_Viterbi_backtrack.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				bool dPrint = false;
				if(((e1->label.substr(0, 10) == "HLAA*02:01") && (e2->label.substr(0, 10) == "HLAA*23:01" )) || ((e2->label.substr(0, 10) == "HLAA*02:01") && (e1->label.substr(0, 10) == "HLAA*23:01" )))
				{
					if(s1 <= s2)
					{
						dPrint = true;
					}
				}

				if(dPrint)
				{
					subst_cout << "LABELS " << e1->label << "/" << e2->label << "\n";
				}

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}

					local_edge_length = local_edge_length / 2;

				}

				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				if(dPrint)
				{
					subst_cout << "\t gaps: " << combinedEdgeEmission[gap_symbol] << "\n";
					subst_cout << "\t stars: " << combinedEdgeEmission[star_symbol] << "\n";

				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{

							for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
							{
								int missing = kMerCount - possibleUnderlyingCount;
								double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
								thisKMer_estimated_missing += p * (double)missing;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double expected_missing_coverage = coverage * 1.0;
									double lambda = (double)kMerCount*expected_missing_coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);

					if(dPrint)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						subst_cout << "\t\t symbol " << kMerSymbol << " " << kMerCount << " [required] "  << globalEmission[kMer] << " [copies in reads] " << thisKMer_estimated_missing  << " [estimated missing] " << this_kMer_expected_missing << " [expected missing] " << "                 " << kMer << "\n";
					}
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				absolute_coverage_e1 = absolute_coverage_e1 / (double)total_kMers_e1;
				absolute_coverage_e2 = absolute_coverage_e2 / (double)total_kMers_e2;
				fw_Coverage_e1.at(level).at(state) = absolute_coverage_e1;
				fw_Coverage_e2.at(level).at(state) = absolute_coverage_e2;

				fw_e1_length.at(level).at(state) = total_edge_length;
				fw_e2_length.at(level).at(state) = total_edge_length;
				fw_e1_length_valid.at(level).at(state) = total_kMers_e1;
				fw_e2_length_valid.at(level).at(state) = total_kMers_e2;

				if(dPrint)
				{
					subst_cout << "\t estimated_missing_e1: " << estimated_missing_e1 << "\t total_kMers_e1: " << total_kMers_e1 << " => symbol_emission_e1: " << symbol_emission_e1 << "\n";
					subst_cout << "\t estimated_missing_e2: " << estimated_missing_e2 << "\t total_kMers_e2: " << total_kMers_e2 << " => symbol_emission_e2: " << symbol_emission_e2 << "\n";
				}


				if(dPrint)
				{
					#pragma omp critical
					{
						cout << subst_cout.str() << flush;
					}
				}

				if(level == 0)
				{
					fw.at(level).at(state) = 0;
					fw_Viterbi.at(level).at(state) = 0;

					vector<int> states_same_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					fw.at(level).at(state) = 0;

					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							states_same_jump_P.clear();
							states_same_jump_P.push_back(max_jump_state);
							assert(set_jump_P.count(max_jump_state) == 0);
							set_jump_P.insert(max_jump_state);
						}
						else if(viterbiJumpP == max_jump_P)
						{
							states_same_jump_P.push_back(from);
							assert(set_jump_P.count(from) == 0);
							set_jump_P.insert(from);
						}
					}

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
					assert(fw_Viterbi_backtrack.at(level).size() > 0);

				}

				if(total_kMers_e1 == 0)
				{
					symbol_emission_e1 = 1;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}

				if(total_kMers_e2 == 0)
				{
					symbol_emission_e2 = 1;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				int missing_edges_e1 = total_edge_length - total_kMers_e1;
				int missing_edges_e2 = total_edge_length - total_kMers_e2;

				int informative_back_edges_e1 = 0;
				if (missing_edges_e1 > 0)
				{
					informative_back_edges_e1 = (missing_edges_e1 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e1;
				}
				int informative_back_edges_e2 = 0;
				if (missing_edges_e2 > 0)
				{
					informative_back_edges_e2 = (missing_edges_e2 > myGraph->kMerSize) ?  myGraph->kMerSize : missing_edges_e2;
				}

				double adapted_emission_e1 = symbol_emission_e1;
				double adapted_emission_e2 = symbol_emission_e2;

				if(((double)total_kMers_e1/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e1 = total_kMers_e1;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e1;
					weights.push_back(total_kMers_e1);
					values.push_back(symbol_emission_e1);
					while((current_level > 0) && (got_backtracked_edges_e1 < informative_back_edges_e1))
					{
						running_state = fw_Viterbi_backtrack.at(current_level).at(running_state).at(0);
						current_level--;

						int running_level_edges = (int)fw_e1_length_valid.at(current_level).at(running_state);
						double running_level_emission = fw_Emission_e1.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e1 + running_level_edges) > informative_back_edges_e1)
						{
							int too_many = (got_backtracked_edges_e1 + running_level_edges) - informative_back_edges_e1;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e1 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					assert(weights_sum > 0);
					adapted_emission_e1 = 0;
					for(int wI = 0; wI < (int)weights.size(); wI++)
					{
						adapted_emission_e1 += (weights.at(wI)/weights_sum) * values.at(wI);
					}
				}

				if(((double)total_kMers_e2/(double)total_edge_length) < 0.3)
				{
					int current_level = level;
					int got_backtracked_edges_e2 = total_kMers_e2;
					int running_state = state;

					vector<double> weights;
					vector<double> values;
					double weights_sum = total_kMers_e2;
					weights.push_back(total_kMers_e2);
					values.push_back(symbol_emission_e2);
					while((current_level > 0) && (got_backtracked_edges_e2 < informative_back_edges_e2))
					{
						running_state = fw_Viterbi_backtrack.at(current_level).at(running_state).at(0);
						current_level--;

						int running_level_edges = (int)fw_e2_length_valid.at(current_level).at(running_state);
						double running_level_emission = fw_Emission_e2.at(current_level).at(running_state);

						int take_edges;
						if((got_backtracked_edges_e2 + running_level_edges) > informative_back_edges_e2)
						{
							int too_many = (got_backtracked_edges_e2 + running_level_edges) - informative_back_edges_e2;
							take_edges = running_level_edges - too_many;
						}
						else
						{
							take_edges = running_level_edges;
						}

						if(take_edges > 0)
						{
							weights.push_back(take_edges);
							weights_sum += take_edges;
							values.push_back(running_level_emission);
							got_backtracked_edges_e2 += take_edges;
						}
					}

					assert(weights.size() == values.size());
					assert(weights_sum > 0);
					adapted_emission_e2 = 0;
					for(int wI = 0; wI < (int)weights.size(); wI++)
					{
						adapted_emission_e2 += (weights.at(wI)/weights_sum) * values.at(wI);
					}
				}

				double symbol_emission = 0.5*adapted_emission_e1 + 0.5*adapted_emission_e2;

				double e1_lengthWeight = 1 + 0.10 * ((double)total_kMers_e1/(double)total_edge_length);
				double e2_lengthWeight = 1 + 0.10 * ((double)total_kMers_e2/(double)total_edge_length);

				double symbol_emission_lengthWeighted  =  0.5*adapted_emission_e1*e1_lengthWeight + 0.5*adapted_emission_e2*e2_lengthWeight;

				assert(symbol_emission >= 0);

				fw_adaptedEmission.at(level).at(state) = symbol_emission;
				fw_adaptedEmission_e1.at(level).at(state) = adapted_emission_e1*e1_lengthWeight;
				fw_adaptedEmission_e2.at(level).at(state) = adapted_emission_e2*e2_lengthWeight;

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission_lengthWeighted*(double)total_edge_length);
				fw_Emission.at(level).at(state) = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				fw_Emission_e1.at(level).at(state) = symbol_emission_e1;
				fw_Emission_e2.at(level).at(state) = symbol_emission_e2;
				fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				assert(fw.at(level).at(state) >= 0);
			}

		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}


double AlphaHMM::fillForwardBackwardTable_1(map<string, long long> globalEmission)
{

	assert(genotypingMode == 1);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);
	fw_Emission.resize(levels);
	fw_expectedMissing_e1.resize(levels);
	fw_expectedMissing_e2.resize(levels);
	fw_Coverage_e1.resize(levels);
	fw_Coverage_e2.resize(levels);
	fw_Viterbi_backtrack.resize(levels);

	fw_e1_present.resize(levels);
	fw_e2_present.resize(levels);
	fw_e1_length.resize(levels);
	fw_e2_length.resize(levels);

	double coverage = estimateCoverage(globalEmission);
	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Emission.at(level).resize(states);
		fw_expectedMissing_e1.at(level).resize(states);
		fw_expectedMissing_e2.at(level).resize(states);
		fw_Coverage_e1.at(level).resize(states);
		fw_Coverage_e2.at(level).resize(states);

		fw_e1_present.at(level).resize(states);
		fw_e2_present.at(level).resize(states);
		fw_e1_length.at(level).resize(states);
		fw_e2_length.at(level).resize(states);

		fw_Viterbi_backtrack.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				bool dPrint = false;
				if(((e1->label.substr(0, 10) == "HLAA*02:01") && (e2->label.substr(0, 10) == "HLAA*23:01" )) || ((e2->label.substr(0, 10) == "HLAA*02:01") && (e1->label.substr(0, 10) == "HLAA*23:01" )))
				{
					if(s1 <= s2)
					{
						dPrint = true;
					}
				}

				if(dPrint)
				{
					subst_cout << "LABELS " << e1->label << "/" << e2->label << "\n";
				}

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}
					assert((local_edge_length % 2) == 0);
					local_edge_length = local_edge_length / 2;

				}

				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				if(dPrint)
				{
					subst_cout << "\t gaps: " << combinedEdgeEmission[gap_symbol] << "\n";
					subst_cout << "\t stars: " << combinedEdgeEmission[star_symbol] << "\n";
				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];


						if(kMerCount > 0)
						{

							for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
							{
								int missing = kMerCount - possibleUnderlyingCount;
								double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
								thisKMer_estimated_missing += p * (double)missing;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double lambda = (double)kMerCount*coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				double estimated_present_e1 = total_kMers_e1 - estimated_missing_e1;
				double estimated_present_e2 = total_kMers_e2 - estimated_missing_e2;

				assert(estimated_present_e1 >= 0);
				assert(estimated_present_e2 >= 0);

				if(dPrint)
				{
					subst_cout << "\t estimated_missing_e1: " << estimated_missing_e1 << "\t total_kMers_e1: " << total_kMers_e1 << " => symbol_emission_e1: " << symbol_emission_e1 << "\n";
					subst_cout << "\t estimated_missing_e2: " << estimated_missing_e2 << "\t total_kMers_e2: " << total_kMers_e2 << " => symbol_emission_e2: " << symbol_emission_e2 << "\n";
				}


				if(dPrint)
				{
					#pragma omp critical
					{
						cout << subst_cout.str() << flush;
					}
				}

				if(level == 0)
				{


					vector<int> states_same_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;

					double e1_optimality = 0;
					double e2_optimality = 0;
					if(total_kMers_e1 != 0)
					{
						e1_optimality = estimated_present_e1/total_kMers_e1;
					}
					if(total_kMers_e2 != 0)
					{
						e2_optimality = estimated_present_e2/total_kMers_e2;
					}

					fw_Viterbi.at(level).at(state) = 0.5 * e1_optimality + 0.5 * e2_optimality;

					fw_e1_present.at(level).at(state) = estimated_present_e1;
					fw_e2_present.at(level).at(state) = estimated_present_e2;
					fw_e1_length.at(level).at(state) = total_kMers_e1;
					fw_e2_length.at(level).at(state) = total_kMers_e2;

				}
				else
				{
					fw.at(level).at(state) = 0;

					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP =
								0.5 * (fw_e1_present.at(level-1).at(from)+estimated_present_e1)/(fw_e1_length.at(level-1).at(from)+total_kMers_e1) +
								0.5 * (fw_e2_present.at(level-1).at(from)+estimated_present_e2)/(fw_e2_length.at(level-1).at(from)+total_kMers_e2);

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							states_same_jump_P.clear();
							states_same_jump_P.push_back(max_jump_state);
							assert(set_jump_P.count(max_jump_state) == 0);
							set_jump_P.insert(max_jump_state);
						}
					}

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;

					fw_e1_present.at(level).at(state) = fw_e1_present.at(level-1).at(states_same_jump_P.at(0))+estimated_present_e1;
					fw_e2_present.at(level).at(state) = fw_e2_present.at(level-1).at(states_same_jump_P.at(0))+estimated_present_e2;
					fw_e1_length.at(level).at(state) = fw_e1_length.at(level-1).at(states_same_jump_P.at(0))+total_kMers_e1;
					fw_e2_length.at(level).at(state) = fw_e2_length.at(level-1).at(states_same_jump_P.at(0))+total_kMers_e2;

					assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);

				}
				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
					assert(fw_Viterbi_backtrack.at(level).size() > 0);

				}

				// some debugging information
				if(total_kMers_e1 == 0)
				{
					absolute_coverage_e1 = -1;
				}
				else
				{
					absolute_coverage_e1 = absolute_coverage_e1 / (double)total_kMers_e1;
				}
				if(total_kMers_e2 == 0)
				{
					absolute_coverage_e2 = -1;
				}
				else
				{
					absolute_coverage_e2 = absolute_coverage_e2 / (double)total_kMers_e2;
				}
				fw_Coverage_e1.at(level).at(state) = absolute_coverage_e1;
				fw_Coverage_e2.at(level).at(state) = absolute_coverage_e2;
				if(total_kMers_e1 == 0)
				{
					symbol_emission_e1 = 1;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}
				if(total_kMers_e2 == 0)
				{
					symbol_emission_e2 = 1;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}
				double symbol_emission =  0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;

				assert(symbol_emission >= 0);

				fw_Emission.at(level).at(state) = symbol_emission;
				fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				assert(fw.at(level).at(state) >= 0);
			}
		}

		double local_max = 0;
		int local_max_state = -1;
		for(int state2 = 0; state2 < states; state2++)
		{
			if((state2 == 0) || (fw_Viterbi.at(level).at(state2) > local_max))
			{
				local_max = fw_Viterbi.at(level).at(state2);
				local_max_state = state2;
			}
		}

		// cout << "\n\nLevel " << level << ", viterbi  max state: " << local_max_state << "\n";

		//cout << "\tlabel1: " << e1->label << "\n";
		//cout << "\tlabel2: " << e2->label << "\n";

		assert(grand_total_edge_length != 0);
		//cout << "\tRunning emission value " << fw_Viterbi.at(level).at(local_max_state)/grand_total_edge_length << " (" << fw_Viterbi.at(level).at(local_max_state) <<  ")\n";

		for(int state2 = 0; state2 < states; state2++)
		{
			//fw.at(level).at(state2) = fw.at(level).at(state2)/local_max;
			//assert(fw.at(level).at(state2) >= 0);
		}

		//fw_underflow_factor *= local_max;
	}

	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}


double AlphaHMM::fillForwardBackwardTable_0(map<string, long long> globalEmission)
{
	assert(genotypingMode == 0);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);
	fw_Emission.resize(levels);
	fw_expectedMissing_e1.resize(levels);
	fw_expectedMissing_e2.resize(levels);
	fw_Coverage_e1.resize(levels);
	fw_Coverage_e2.resize(levels);
	fw_Viterbi_backtrack.resize(levels);
	fw_underflow_factor = 1;

	double coverage = estimateCoverage(globalEmission);


	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	cout << levels << " Levels\n";
	int grand_total_edge_length = 0;

	for(int level = 0; level < levels; level++)
	{
		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Emission.at(level).resize(states);
		fw_expectedMissing_e1.at(level).resize(states);
		fw_expectedMissing_e2.at(level).resize(states);
		fw_Coverage_e1.at(level).resize(states);
		fw_Coverage_e2.at(level).resize(states);

		fw_Viterbi_backtrack.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();
			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;
			if((thisThread == (CONFIG.threads-1)) && (lastPair < (states-1)))
			{
				lastPair = states - 1;
			}
			assert(lastPair <= (states - 1));

			map< int, map< int, vector<double> > > cache_estimated_expected_Missing;

			int total_edge_length = 0;

			for(int state = firstPair; state <= lastPair; state++)
			{
				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				bool dPrint = false;
				if(((e1->label.substr(0, 10) == "HLAA*02:01") && (e2->label.substr(0, 10) == "HLAA*23:01" )) || ((e2->label.substr(0, 10) == "HLAA*02:01") && (e1->label.substr(0, 10) == "HLAA*23:01" )))
				{
					if(s1 <= s2)
					{
						dPrint = true;
					}
				}

				if(dPrint)
				{
					subst_cout << "LABELS " << e1->label << "/" << e2->label << "\n";
				}

				map<int, int> combinedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedEdgeEmission.count(emissionSymbol) == 0)
						combinedEdgeEmission[emissionSymbol] = 0;

					combinedEdgeEmission[emissionSymbol] += number;
				}

				int local_edge_length = 0;
				if((state == firstPair) || (state == (firstPair + 1)))
				{
					for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
					{
						int number = emissionIt->second;
						local_edge_length += number;
					}
					local_edge_length = local_edge_length / 2;

				}

				if(state == firstPair)
				{
					assert(local_edge_length != 0);
					total_edge_length = local_edge_length;
					if(state == 0)
					{
						grand_total_edge_length += total_edge_length;
					}
				}

				assert(total_edge_length != 0);

				if(state == (firstPair+1))
				{
					assert(total_edge_length == local_edge_length);
				}

				if(dPrint)
				{
					subst_cout << "\t gaps: " << combinedEdgeEmission[gap_symbol] << "\n";
					subst_cout << "\t stars: " << combinedEdgeEmission[star_symbol] << "\n";

				}

				combinedEdgeEmission.erase(gap_symbol);
				combinedEdgeEmission.erase(star_symbol);

				double symbol_emission_e1 = 0;
				double symbol_emission_e2 = 0;
				double estimated_missing_e1 = 0;
				double estimated_missing_e2 = 0;
				double expected_missing_e1 = 0;
				double expected_missing_e2 = 0;
				double total_kMers_e1 = 0;
				double total_kMers_e2 = 0;
				double absolute_coverage_e1 = 0;
				double absolute_coverage_e2 = 0;

				for(map<int, int>::iterator kMerIt = combinedEdgeEmission.begin(); kMerIt != combinedEdgeEmission.end(); kMerIt++)
				{
					int kMerSymbol = kMerIt->first;
					int kMerCount = kMerIt->second;

					double thisKMer_estimated_missing = 0;
					double this_kMer_expected_missing = 0;

					int e1_count = 0;
					int e2_count = 0;
					if(e1->multiEmission.count(kMerSymbol) > 0)
					{
						e1_count = e1->multiEmission[kMerSymbol];
					}
					if(e2->multiEmission.count(kMerSymbol) > 0)
					{
						e2_count = e2->multiEmission[kMerSymbol];
					}

					assert((e1_count + e2_count) == kMerCount);

					total_kMers_e1 += e1_count;
					total_kMers_e2 += e2_count;

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
					assert(globalEmission.count(kMer) > 0);
					int observedNumber = globalEmission[kMer];
					absolute_coverage_e1 += observedNumber*((double)e1_count/(double)kMerCount);
					absolute_coverage_e2 += observedNumber*((double)e2_count/(double)kMerCount);

					if((cache_estimated_expected_Missing.count(kMerSymbol) > 0) && (cache_estimated_expected_Missing[kMerSymbol].count(kMerCount) > 0))
					{
						thisKMer_estimated_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(0);
						this_kMer_expected_missing = cache_estimated_expected_Missing[kMerSymbol][kMerCount].at(1);
					}
					else
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						assert(kMer != "*");
						assert(kMer != "_");

						assert(globalEmission.count(kMer) > 0);
						int observedNumber = globalEmission[kMer];

						if(kMerCount > 0)
						{

							for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
							{
								int missing = kMerCount - possibleUnderlyingCount;
								double p = observedXunderlying.at(observedNumber).at(possibleUnderlyingCount);
								thisKMer_estimated_missing += p * (double)missing;
							}

							if(stateKMerCount_expectedMissing_H0.count(kMerCount) == 0)
							{
								double expected_missingness = 0;
								for(int potentialObserved = 0; potentialObserved < (int)observedXunderlying.size(); potentialObserved++)
								{
									double thisObservation_expectedMissing = 0;

									double lambda = (double)kMerCount*coverage;
									poisson_up poisson(lambda);
									double thisObservation_weight = pdf(poisson, potentialObserved);
									if(kMerCount > 0)
									{
										for(int possibleUnderlyingCount = 0; possibleUnderlyingCount < kMerCount; possibleUnderlyingCount++)
										{
											int missing = kMerCount - possibleUnderlyingCount;
											double p = observedXunderlying.at(potentialObserved).at(possibleUnderlyingCount);
											thisObservation_expectedMissing += p * (double)missing;
										}
									}
									expected_missingness += thisObservation_expectedMissing * thisObservation_weight;
								}

								#pragma omp critical
								{
									stateKMerCount_expectedMissing_H0[kMerCount] = 	expected_missingness;
									this_kMer_expected_missing = expected_missingness;
								}
							}
							else
							{
								this_kMer_expected_missing = stateKMerCount_expectedMissing_H0[kMerCount];
							}

						}

						vector<double> forCache;
						forCache.push_back(thisKMer_estimated_missing);
						forCache.push_back(this_kMer_expected_missing);
						cache_estimated_expected_Missing[kMerSymbol][kMerCount] = forCache;
					}

					assert(kMerCount > 0);

					estimated_missing_e1 += thisKMer_estimated_missing * ((double)e1_count / (double)kMerCount);
					expected_missing_e1 += this_kMer_expected_missing * ((double)e1_count / (double)kMerCount);

					estimated_missing_e2 += thisKMer_estimated_missing * ((double)e2_count / (double)kMerCount);
					expected_missing_e2 += this_kMer_expected_missing * ((double)e2_count / (double)kMerCount);

					if(dPrint)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerSymbol);
						subst_cout << "\t\t symbol " << kMerSymbol << " " << kMerCount << " [required] "  << globalEmission[kMer] << " [copies in reads] " << thisKMer_estimated_missing  << " [estimated missing] " << this_kMer_expected_missing << " [expected missing] " << "                 " << kMer << "\n";
					}
				}

				estimated_missing_e1 -= expected_missing_e1;
				estimated_missing_e2 -= expected_missing_e2;

				if(total_kMers_e1 == 0)
				{
					symbol_emission_e1 = 1;
				}
				else
				{
					symbol_emission_e1 = 1 - (estimated_missing_e1/total_kMers_e1);
				}
				if(total_kMers_e2 == 0)
				{
					symbol_emission_e2 = 1;
				}
				else
				{
					symbol_emission_e2 = 1 - (estimated_missing_e2/total_kMers_e2);
				}

				double symbol_emission = 0.5*symbol_emission_e1 + 0.5*symbol_emission_e2;
				assert(symbol_emission >= 0);

				absolute_coverage_e1 = absolute_coverage_e1 / (double)total_kMers_e1;
				absolute_coverage_e2 = absolute_coverage_e2 / (double)total_kMers_e2;
				fw_Coverage_e1.at(level).at(state) = absolute_coverage_e1;
				fw_Coverage_e2.at(level).at(state) = absolute_coverage_e2;

				if(dPrint)
				{
					subst_cout << "\t estimated_missing_e1: " << estimated_missing_e1 << "\t total_kMers_e1: " << total_kMers_e1 << " => symbol_emission_e1: " << symbol_emission_e1 << "\n";
					subst_cout << "\t estimated_missing_e2: " << estimated_missing_e2 << "\t total_kMers_e2: " << total_kMers_e2 << " => symbol_emission_e2: " << symbol_emission_e2 << "\n";
					subst_cout << "\t\t symbol_emission " << symbol_emission << "\n";
				}


				if(dPrint)
				{
					#pragma omp critical
					{
						cout << subst_cout.str() << flush;
					}
				}

				if(level == 0)
				{
					fw.at(level).at(state) = 0;
					fw_Viterbi.at(level).at(state) = 0;

					vector<int> states_same_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
				}
				else
				{
					fw.at(level).at(state) = 0;

					double max_jump_P = -1;
					int max_jump_state = -1;
					vector<int> states_same_jump_P;
					set<int> set_jump_P;

					set<int> possibleOrigins;
					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						possibleOrigins.insert(from);
					}
					assert(possibleOrigins.size() == diploidStateTransitions_Reverse.at(level).at(state).size());

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) * 1.0;

						if(viterbiJumpP > max_jump_P)
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
							states_same_jump_P.clear();
							states_same_jump_P.push_back(max_jump_state);
							assert(set_jump_P.count(max_jump_state) == 0);
							set_jump_P.insert(max_jump_state);
						}
						else if(viterbiJumpP == max_jump_P)
						{
							states_same_jump_P.push_back(from);
							assert(set_jump_P.count(from) == 0);
							set_jump_P.insert(from);
						}
					}

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack.at(level).at(state) = states_same_jump_P;
					assert(fw_Viterbi_backtrack.at(level).at(state).size() > 0);
				}

				assert(fw_Viterbi.at(level).at(state) != -1);
				if(level != 0)
				{
					assert(fw_Viterbi.at(level).at(state) > 0);
					assert(fw_Viterbi_backtrack.at(level).size() > 0);

				}

				fw_Viterbi.at(level).at(state) = fw_Viterbi.at(level).at(state) + (symbol_emission*(double)total_edge_length);
				fw_Emission.at(level).at(state) = (symbol_emission*(double)total_edge_length);
				fw_expectedMissing_e1.at(level).at(state) = expected_missing_e1;
				fw_expectedMissing_e2.at(level).at(state) = expected_missing_e2;

				assert(fw.at(level).at(state) >= 0);
			}
		}

	}

	assert(grand_total_edge_length != 0);
	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		fw_Viterbi.at(levels-1).at(state) = fw_Viterbi.at(levels-1).at(state) / (double)grand_total_edge_length;
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure " << max_viterbi << " (total length " << grand_total_edge_length << ", unnormalized optimality " << max_viterbi*(double)grand_total_edge_length << ")\n";

	return max_viterbi;
}


