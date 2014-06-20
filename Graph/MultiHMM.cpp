#include "MultiHMM.h"
#include <assert.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "../Utilities.h"
#include "../MHC-PRG.h"
#include <boost/math/distributions/poisson.hpp>

int multiUnderflowProtection = 1;

map<string, long long> MultiHMM::estimateEmissions(string kMerCountsSamplePath, string kMerCountsGenomePath)
{
	// this function is more or less legacy code
	// it used to be used for computing an estimated number of kMers over
	// the area covered by the haplotype graph

	// this is not necessary anymore due to the new emission model
	// we only generate a nice and suitable vector of emissions

	// read in kMer count in sample

	map<string, long long> kMerCountsInSample;
	ifstream kMerCountsSampleFile;
	kMerCountsSampleFile.open (kMerCountsSamplePath.c_str(), ios::in);
	int lineCounter = 0;
	if(kMerCountsSampleFile.is_open())
	{
		string line;
		while(kMerCountsSampleFile.good())
		{
			lineCounter++;
			getline (kMerCountsSampleFile, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			vector<string> fields = Utilities::split(line, ' ');
			if(fields.size() != 2)
				errEx("Strange format in kMerCountsInSample file. Expect two fields (kMer and count). Problem occured in line: "+Utilities::ItoStr(lineCounter));

			if((fields.at(0) != "MeanReadLen") && (fields.at(0) != "TotalKMerCoverage") && (fields.at(0) != "TotalSeq"))
			{
				assert((int)fields.at(0).length() == myGraph->kMerSize);
			}

			kMerCountsInSample[fields.at(0)] = Utilities::StrtoLongLong(fields.at(1));
		}
		kMerCountsSampleFile.close();
	}
	else
	{
		errEx("Cannot open kMer counts file: "+kMerCountsSamplePath);
	}

	return kMerCountsInSample;
}

MultiHMM::MultiHMM(MultiGraph* g, bool verbose)
{
	myGraph = g;

	Edge2Int.clear();
	Edge2Level.clear();
	
	int levels = g->NodesPerLevel.size();
	haploidStatesByLevel.resize(levels-1);
	Int2Edge.resize(levels);
	
	int stats_multiNomial = 0;
	int stats_mereCompatible = 0;
	set<string> stats_kMersInGraph;

	for(int level = 0; level < levels; level++)
	{
		vector<Node*> nodes(g->NodesPerLevel.at(level).begin(), g->NodesPerLevel.at(level).end());
		for(unsigned int nI = 0; nI < g->NodesPerLevel.at(level).size(); nI++)
		{
			Node* n = nodes.at(nI);
			for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
			{
				Edge* e = (*outgoingIt);
				assert(Edge2Int.count(e) == 0);
				if(Edge2Int.count(e) == 0)
				{
					int size_now = Int2Edge.at(level).size();
					Edge2Int[e] = size_now;
					Edge2Level[e] = level;
					Int2Edge.at(level)[size_now] = e;
				}
				
				haploidState thisState;
				thisState.level = level;
				thisState.id = Edge2Int[e];
				thisState.e = e;

				haploidStatesByLevel.at(level).push_back(thisState);

				if(level == 0)
				{
					assert(e->From->Sum_Outgoing() > 0);
					double init_p = e->count/e->From->Sum_Outgoing();
					haploid_initialProbabilities.push_back(init_p);
				}

				// compatibility kMers

				int number_of_multiKMers = 0;
				for(map<int, int>::iterator emissionIt = e->multiEmission.begin(); emissionIt != e->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					assert(emissionIt->second >= 0);

					if(emissionIt->second > 0)
					{
						string kMer = g->CODE.deCode(e->locus_id, emissionSymbol);
						if(kMer != "_")
						{
							stats_kMersInGraph.insert(kMer);
							number_of_multiKMers++;
						}
					}
				}
				if(number_of_multiKMers > 0)
				{
					stats_mereCompatible++;
				}

				// unique kMers

				if(SubLevelToUniqueKMers.count(level) == 0)
				{
					map< int, set<string> > oneLevel;
					SubLevelToUniqueKMers[level] = oneLevel;
				}
				for(int l_i = 0; l_i < (int)e->multiEmission_full.size(); l_i++)
				{
					int symbol = e->multiEmission_full.at(l_i);
					string kMer = g->CODE.deCode(e->locus_id, symbol);
					if(SubLevelToUniqueKMers[level].count(l_i) == 0)
					{
						set<string> newSet;
						SubLevelToUniqueKMers[level][l_i] = newSet;
					}
					if(kMer != "_")
					{
						stats_kMersInGraph.insert(kMer);
						number_of_multiKMers++;
						SubLevelToUniqueKMers[level][l_i].insert(kMer);
					}
					stats_multiNomial++;
				}

			}
		}
	}

	// connect haploid states
	for(int level = 0; level < (levels-1); level++)
	{
		for(unsigned int haploidStateI = 0; haploidStateI < haploidStatesByLevel.at(level).size(); haploidStateI++)
		{
			Node* To = haploidStatesByLevel.at(level).at(haploidStateI).e->To;
			
			double To_Outgoing_Sum = To->Sum_Outgoing();
			double realized_leaveState = 0.0;
			if(To_Outgoing_Sum > 0)
			{
				for(set<Edge*>::iterator outgoingIt = To->Outgoing_Edges.begin(); outgoingIt != To->Outgoing_Edges.end(); outgoingIt++)
				{
					Edge* e = (*outgoingIt);	

					assert(Edge2Int.count(e) > 0);
					assert(Edge2Level[e] == level+1);
					
					int indexTo = Edge2Int[e];

					stateAlternative forward;
					forward.id = indexTo;
					forward.p = e->count / To_Outgoing_Sum;
					realized_leaveState += forward.p;
					
					stateAlternative backward;
					backward.id = haploidStateI;
					backward.p = e->count / To_Outgoing_Sum;
					
					haploidStatesByLevel.at(level+1).at(indexTo).previous.push_back(backward);					
					haploidStatesByLevel.at(level).at(haploidStateI).next.push_back(forward);
				}
				
				assert(abs(realized_leaveState - 1) < epsilon);
			}
		}
	}
	
	// create diploid states

	diploidStateTransitions.resize(levels-1);
	diploidStateTransitions_Reverse.resize(levels-1);
	
	if (verbose)
		cout << "Diploid Graph HMM\n";
	diploidStatesByLevel.resize(levels-1);
	
	int chunk_size = (levels-1) / CONFIG.threads;
	
	#pragma omp parallel
	{
		assert(omp_get_num_threads() == CONFIG.threads);
		int thisThread = omp_get_thread_num();
		int firstPair = thisThread * chunk_size;
		int lastPair = (thisThread+1) * chunk_size - 1;
		if((thisThread == (CONFIG.threads-1)) && (lastPair < (levels-2)))
		{ 
			lastPair = levels-2;				
		}
		
		for(int level = firstPair; level <= lastPair; level++)
		{
			int haploidStates = haploidStatesByLevel.at(level).size();
			int diploidStates = int(pow((double)haploidStates, 2)+.005);
			
			// int diploidStates = pow(haploidStates, 2);
			// cerr << "Level " << level << " " << diploidStates << " diploid states\n";
			
			diploidStatesByLevel.at(level) = diploidStates;
			diploidStateTransitions.at(level).resize(diploidStates);
		
			if (verbose)
				cout << "Level " << level << ": " << diploidStates << " diploid states.\n" << flush;
			
			if(level == 0)
			{
				diploidStateTransitions_Reverse.at(level).resize(diploidStates);
			}
			if(level < (levels-2))
			{
				int diploidStates_n1 = int(pow((double)haploidStatesByLevel.at(level+1).size(), 2)+0.005);
				diploidStateTransitions_Reverse.at(level+1).resize(diploidStates_n1);
			}
				
			for(int stateI = 0; stateI < diploidStates; stateI++)
			{
				int s1 = stateI % haploidStates;
				int s2 = stateI / haploidStates;
				
				if(level == 0)
				{
					diploid_initialProbabilities.push_back(haploid_initialProbabilities.at(s1) * haploid_initialProbabilities.at(s2));
				}

				if(level < (levels-2))
				{
					double summed_jumps = 0.0;
					int haploidStates_nextLevel = haploidStatesByLevel.at(level+1).size();
					for(unsigned int i = 0; i < haploidStatesByLevel.at(level).at(s1).next.size(); i++)
					{
						stateAlternative jump1 = haploidStatesByLevel.at(level).at(s1).next.at(i);
						for(unsigned int j = 0; j < haploidStatesByLevel.at(level).at(s2).next.size(); j++)
						{
							stateAlternative jump2 = haploidStatesByLevel.at(level).at(s2).next.at(j);
							int diploidTarget = jump1.id + (jump2.id * haploidStates_nextLevel);
							double p = jump1.p * jump2.p;
							summed_jumps += p;
							
							stateAlternative jumpStateCombined;
							jumpStateCombined.id = diploidTarget;
							jumpStateCombined.p = p;
							
							stateAlternative jumpStateCombinedBackwards;
							jumpStateCombinedBackwards.id = stateI;
							jumpStateCombinedBackwards.p = p;						
							
							diploidStateTransitions.at(level).at(stateI).push_back(jumpStateCombined);
							diploidStateTransitions_Reverse.at(level+1).at(diploidTarget).push_back(jumpStateCombinedBackwards);
						}
					}
					if(! ((abs(summed_jumps - 1) < epsilon) || (abs(summed_jumps) < epsilon)))
					{
						cerr << "\n\n" << flush;
						cerr << "State: " << stateI << "\n";
						cerr << "\tLevel: " << level << "\n";
						cerr << "Summed jumps: " << summed_jumps << "\n" << flush;
						assert((abs(summed_jumps - 1) < epsilon) || (abs(summed_jumps) < epsilon)) ;				
					}
				}	
			}			
		}
	}

	if(verbose)
	{
		cout << "Levels in HMM: " << levels-1 << "\n";
		cout << "Total number of implemented k-mers: " << stats_kMersInGraph.size() << "\n";
		cout << "Total haploid number of compatibility states: " << stats_mereCompatible << "\n";
		cout << "Total haploid number of multinomial states: " << stats_multiNomial << "\n";
	}

}
multiHaploLabelPair MultiHMM::retrieveProbabilisticSample(bool labelOnly)
{
	vector<int> backtrack_states;
	int levels = diploidStatesByLevel.size();
	int runningState = -1;
	for(int level = levels-1; level >= 0; level--)
	{		
		//int states = diploidStatesByLevel.at(level);
		if(level == (levels-1))
		{
			vector<double> stateWeights;		
			stateWeights = fw.at(level);
			runningState = Utilities::chooseFromVector(stateWeights);
			
		}
		else
		{
			vector<double> stateWeights;		
			vector<int> trueIndices;
			int possibleStatesThisLevel = diploidStateTransitions_Reverse.at(level+1).at(runningState).size();
			for(int sI = 0; sI < possibleStatesThisLevel; sI++)
			{
				stateAlternative thisState = diploidStateTransitions_Reverse.at(level+1).at(runningState).at(sI);
				double weight = thisState.p * fw.at(level).at(thisState.id);
				stateWeights.push_back(weight);
				trueIndices.push_back(thisState.id);
			}
			int selection = Utilities::chooseFromVector(stateWeights);
			runningState = trueIndices.at(selection);
		}
		backtrack_states.push_back(runningState);
	}
	
	reverse(backtrack_states.begin(), backtrack_states.end());
	
	multiHaploLabelPair forReturn;
	for(unsigned int level = 0; level < backtrack_states.size(); level++)
	{
		int diploidState = backtrack_states.at(level);
		int haploidStates = haploidStatesByLevel.at(level).size();
		int s1 = diploidState % haploidStates;
		int s2 = diploidState / haploidStates;
		
		Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
		Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
		
		if((e1->label != "") || (e2->label != ""))
		{
			forReturn.h1.push_back(e1->label);
			forReturn.h2.push_back(e2->label);
		}
	}
	
	return forReturn;	
}

double MultiHMM::fillForwardBackwardTable(map<string, long long> globalEmission)
{
	//assert((globalEmission.size()) == diploidStatesByLevel.size());

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_underflow_factor = 1;

	double emission_probability = 0.0;
	double expectedCoverage = 5;
	assert((globalEmission.count("MeanReadLen") > 0) && (globalEmission.count("TotalSeq") > 0));

	double normal_coverage = globalEmission["TotalSeq"];
	double correction_factor = 0.5 * (double)myGraph->kMerSize/(double)globalEmission["MeanReadLen"];
	cout << "Expected kMer coverage: " << normal_coverage*correction_factor << "\n";

	// TODO test the following
	// expectedCoverage = normal_coverage*correction_factor;

	for(int level = 0; level < levels; level++)
	{
		emission_probability = 0.0;

		int states = diploidStatesByLevel.at(level);
		fw.at(level).resize(states);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;

		if(level == 0)
		{
			double initial_sum = 0.0;
			for(int i = 0; i < states; i++)
			{
				assert((diploid_initialProbabilities.at(i) <= 1) && (diploid_initialProbabilities.at(i) >= 0));
				initial_sum += diploid_initialProbabilities.at(i);
				fw.at(0).at(i) = diploid_initialProbabilities.at(i);
				assert(fw.at(0).at(i) >= 0);
			}
			assert(abs(initial_sum - 1) < epsilon);
		}

		int chunk_size = states / CONFIG.threads;

		vector<long long> unique_count_sublevel;
		assert(SubLevelToUniqueKMers.count(level) > 0);
		unique_count_sublevel.resize(SubLevelToUniqueKMers[level].size());
		for(int sublevel_i = 0; sublevel_i < (int)SubLevelToUniqueKMers[level].size(); sublevel_i++)
		{
			assert(SubLevelToUniqueKMers[level].count(sublevel_i) > 0);
			int total_coverage_this_sublevel = 0;
			for(set<string>::iterator kMerIt = SubLevelToUniqueKMers[level][sublevel_i].begin(); kMerIt != SubLevelToUniqueKMers[level][sublevel_i].end(); kMerIt++)
			{
				string kMer = *kMerIt;
				assert(globalEmission.count(kMer) > 0);
				total_coverage_this_sublevel += globalEmission[kMer];
			}
			unique_count_sublevel.at(sublevel_i) = total_coverage_this_sublevel;
		}

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

			for(int state = firstPair; state <= lastPair; state++)
			{
				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;

				assert(e1->locus_id == e2->locus_id);

				double symbol_emission_prob = 0;

				map<int, int> combinedKMerExpectedEmission_raw;
				int expected_sum = 0;

				// Contributions from compatibility kMers

				for(map<int, int>::iterator emissionIt = e1->multiEmission.begin(); emissionIt != e1->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;
					assert(number >= 0);
					if(combinedKMerExpectedEmission_raw.count(emissionSymbol) == 0)
						combinedKMerExpectedEmission_raw[emissionSymbol] = 0;

					combinedKMerExpectedEmission_raw[emissionSymbol] += number;
					expected_sum += number;
				}

				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;
					assert(number >= 0);
					if(combinedKMerExpectedEmission_raw.count(emissionSymbol) == 0)
						combinedKMerExpectedEmission_raw[emissionSymbol] = 0;

					combinedKMerExpectedEmission_raw[emissionSymbol] += number;
					expected_sum += number;
				}

				// sum up and find out whether above expected threshold
				double P_Compatibility = 0;
				int kMerItCounter = 0;
				for(map<int, int>::iterator emissionIt = combinedKMerExpectedEmission_raw.begin(); emissionIt != combinedKMerExpectedEmission_raw.end(); emissionIt++)
				{
					kMerItCounter++;
					string kMer = myGraph->CODE.deCode(e1->locus_id, emissionIt->first);

					long long count = globalEmission[kMer];
					int expected_underlying_instances_atLeast = emissionIt->second;

					assert(expected_underlying_instances_atLeast > 0);

					if((kMer != "_") && (kMer != "*"))
					{
						assert(globalEmission.count(kMer) > 0);

						double expected_coverage_thisKMer_atLeast = expectedCoverage*expected_underlying_instances_atLeast;
						int count_threshold_thisKMer = expected_coverage_thisKMer_atLeast/3;
						assert(count_threshold_thisKMer > 0);

						boost::math::poisson_distribution<> myPoisson(expected_coverage_thisKMer_atLeast);
						double not_met_threshold_P = cdf(myPoisson, count_threshold_thisKMer);
						assert(not_met_threshold_P >= 0);
						assert(not_met_threshold_P <= 1);

						// cout << "level " << level << ", kMer " << kMer << ", count " << count <<  ", expected_coverage_thisKMer_atLeast " << expected_coverage_thisKMer_atLeast << ", count_threshold_thisKMer " << count_threshold_thisKMer << ", not_met_threshold_P " << not_met_threshold_P << "\n";

						if (count > count_threshold_thisKMer)
						{
							//P_Compatibility *= (1-not_met_threshold_P);
							P_Compatibility += log(1-not_met_threshold_P);
							//assert(P_Compatibility >= 0);
						}
						else
						{
							//P_Compatibility *= not_met_threshold_P;
							P_Compatibility += log(not_met_threshold_P);
							//assert(P_Compatibility >= 0);
						}

						// cout << "P_Compatibility: " << P_Compatibility << "     " << kMerItCounter << "/" << combinedKMerExpectedEmission_raw.size() << "\n";

					}
					else
					{
						if((kMer == "_") || (kMer == "*"))
						{
							// try to reduce the probability for gaps a bit
							// 1 is probably too much, as a concurrent edge with no gaps
							// would also have a probabilty < 1 of emission (random coverage fluctuations)
							int num_of_gaps = expected_underlying_instances_atLeast;
							int count_threshold_uniqueMer = expectedCoverage/3;
							assert(count_threshold_uniqueMer > 0);
							boost::math::poisson_distribution<> myPoisson(expectedCoverage);
							double not_met_threshold_P = cdf(myPoisson, count_threshold_uniqueMer);
							assert(not_met_threshold_P >= 0);
							assert(not_met_threshold_P <= 1);

							double expected_number_notMetP = floor(not_met_threshold_P * num_of_gaps + 0.5);
							assert(num_of_gaps >= expected_number_notMetP);

							double gap_probability = log(not_met_threshold_P) * expected_number_notMetP + log(1 - not_met_threshold_P) * (num_of_gaps - expected_number_notMetP);

							// TODO test the following
							// double gap_probability = log(1 - not_met_threshold_P) * num_of_gaps
							P_Compatibility += gap_probability;
						}
						else
						{
							assert(1 == 0);
						}
					}

					// TODO test the following
					// probability which does not sum over multiple kMers in the same section
					/*
					if((kMer != "_") && (kMer != "*"))
					{
						assert(globalEmission.count(kMer) > 0);

						double expected_coverage_thisKMer_atLeast = expectedCoverage;
						int count_threshold_thisKMer = expected_coverage_thisKMer_atLeast/3;
						assert(count_threshold_thisKMer > 0);

						boost::math::poisson_distribution<> myPoisson(expected_coverage_thisKMer_atLeast);
						double not_met_threshold_P = cdf(myPoisson, count_threshold_thisKMer);
						assert(not_met_threshold_P >= 0);
						assert(not_met_threshold_P <= 1);

						if (count > count_threshold_thisKMer)
						{
							P_Compatibility += log(1-not_met_threshold_P)*expected_underlying_instances_atLeast;
						}
						else
						{
							P_Compatibility += log(not_met_threshold_P)*expected_underlying_instances_atLeast;
						}
					}
					else
					{
						if((kMer == "_") || (kMer == "*"))
						{
							int num_of_gaps = expected_underlying_instances_atLeast;
							int count_threshold_uniqueMer = expectedCoverage/3;
							assert(count_threshold_uniqueMer > 0);
							boost::math::poisson_distribution<> myPoisson(expectedCoverage);
							double not_met_threshold_P = cdf(myPoisson, count_threshold_uniqueMer);
							assert(not_met_threshold_P >= 0);
							assert(not_met_threshold_P <= 1);

							double expected_number_notMetP = floor(not_met_threshold_P * num_of_gaps + 0.5);
							assert(num_of_gaps >= expected_number_notMetP);

							double gap_probability = log(not_met_threshold_P) * expected_number_notMetP + log(1 - not_met_threshold_P) * (num_of_gaps - expected_number_notMetP);
						}
						else
						{
							assert(1 == 0);
						}
					}
					*/

				}

				symbol_emission_prob += P_Compatibility;

				// Contributions from multinomial levels
				// TODO think about the right likelihood
				double P_multinom = 0;
				assert(e1->multiEmission_full.size() == e2->multiEmission_full.size());
				for(int multin_l = 0; multin_l < (int)e1->multiEmission_full.size(); multin_l++)
				{
					int e1_s = e1->multiEmission_full.at(multin_l);
					int e2_s = e2->multiEmission_full.at(multin_l);
					string kMer1 = myGraph->CODE.deCode(e1->locus_id, e1_s);
					string kMer2 = myGraph->CODE.deCode(e1->locus_id, e2_s);

					long long total_kmers_covered = 0;
					double p_emit_correct = 0.995;

					double local_p = 0;

					assert(kMer1 != "*");
					assert(kMer2 != "*");

					if(kMer1 == kMer2)
					{
						if(kMer1 != "_")
						{
							assert(globalEmission.count(kMer1) > 0);
							local_p += (double)globalEmission[kMer1] * log(p_emit_correct);
							//local_p = pow(p_emit_correct, (double)globalEmission[kMer1]);
							total_kmers_covered += globalEmission[kMer1];
						}
					}
					else
					{
						if((kMer1 != "_") && (kMer2 != "_"))
						{
							assert(globalEmission.count(kMer1) > 0);
							// cout << "level " << level << ", multin_l " << multin_l << ", kMer1 " << kMer1 << ", globalEmission[kMer1] " << globalEmission[kMer1] << ", pow(p_emit_correct, (double)globalEmission[kMer1]) " << pow(p_emit_correct, (double)globalEmission[kMer1]) << "\n";
							//local_p = pow(p_emit_correct/2.0, (double)globalEmission[kMer1]);
							local_p += (double)globalEmission[kMer1] * log(p_emit_correct/2);
							total_kmers_covered += globalEmission[kMer1];

							assert(globalEmission.count(kMer2) > 0);
							// cout << "level " << level << ", multin_l " << multin_l << ", kMer2 " << kMer2 << ", globalEmission[kMer2] " << globalEmission[kMer2] << ", pow(p_emit_correct, (double)globalEmission[kMer2]) " << pow(p_emit_correct, (double)globalEmission[kMer2]) << "\n";
							//local_p = pow(p_emit_correct/2.0, (double)globalEmission[kMer2]);
							local_p += (double)globalEmission[kMer2] * log(p_emit_correct/2);
							total_kmers_covered += globalEmission[kMer2];
						}
						else
						{
							if(kMer1 != "_")
							{
								assert(globalEmission.count(kMer1) > 0);
								// cout << "level " << level << ", multin_l " << multin_l << ", kMer1 " << kMer1 << ", globalEmission[kMer1] " << globalEmission[kMer1] << ", pow(p_emit_correct, (double)globalEmission[kMer1]) " << pow(p_emit_correct, (double)globalEmission[kMer1]) << "\n";
								//local_p = pow(p_emit_correct, (double)globalEmission[kMer1]);
								local_p += (double)globalEmission[kMer1] * log(p_emit_correct);
								total_kmers_covered += globalEmission[kMer1];
							}
							else if(kMer2 != "_")
							{
								assert(globalEmission.count(kMer2) > 0);
								// cout << "level " << level << ", multin_l " << multin_l << ", kMer2 " << kMer2 << ", globalEmission[kMer2] " << globalEmission[kMer2] << ", pow(p_emit_correct, (double)globalEmission[kMer2]) " << pow(p_emit_correct, (double)globalEmission[kMer2]) << "\n";
								//local_p = pow(p_emit_correct, (double)globalEmission[kMer2]);
								local_p += (double)globalEmission[kMer2] * log(p_emit_correct);
								total_kmers_covered += globalEmission[kMer2];
							}
							else
							{
								assert(1 == 0);
							}
						}
					}
					assert(total_kmers_covered <= unique_count_sublevel.at(multin_l));

					if(total_kmers_covered != unique_count_sublevel.at(multin_l))
					{
						local_p += (double)(unique_count_sublevel.at(multin_l) - total_kmers_covered) * log((double)1 - p_emit_correct);

						// P_multinom *= pow(1 - p_emit_correct, unique_count_sublevel.at(multin_l) - total_kmers_covered);
						// cout << "level " << level << ", total_kmers_covered " << total_kmers_covered << ", unique_count_sublevel.at(multin_l) " << unique_count_sublevel.at(multin_l) << ", pow(1 - p_emit_correct, unique_count_sublevel.at(multin_l) - total_kmers_covered) " << pow(1 - p_emit_correct, unique_count_sublevel.at(multin_l) - total_kmers_covered) << "\n";
					}

					P_multinom += local_p;

					// cout << "P_multinom: " << P_multinom << "\n";
				}

				symbol_emission_prob += P_multinom;

				//assert(symbol_emission_prob >= 0);

				if(level == 0)
				{
					fw.at(level).at(state) = diploid_initialProbabilities.at(state);
				}
				else
				{
					fw.at(level).at(state) = 0;

					for(unsigned int s2 = 0; s2 < diploidStateTransitions_Reverse.at(level).at(state).size(); s2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(s2);
						int from = jumpState.id;

						if(!(fw.at(level-1).at(from) >= 0))
						{
							cout << "fw.at(level-1).at(from): " << fw.at(level-1).at(from) << "\n";
							cout << "level-1: " << level-1 << "\n";
							cout << "from: " << from << "\n";

						}
						if(!(jumpState.p >= 0))
						{
							cout << "jumpState.p: " << jumpState.p << "\n";
						}

						assert(fw.at(level-1).at(from) >= 0);
						assert(jumpState.p >= 0);

						fw.at(level).at(state) += fw.at(level-1).at(from) * jumpState.p;
						assert(fw.at(level).at(state) >= 0);
					}
				}
				assert(fw.at(level).at(state) >= 0);

				fw.at(level).at(state) = log(fw.at(level).at(state)) + symbol_emission_prob;

				//assert(fw.at(level).at(state) >= 0);
			}
		}

		/*
		if((level > 0) && (multiUnderflowProtection > 0) && ((level % multiUnderflowProtection) == 0))
		{
			double local_sum = 0.0;
			for(int state2 = 0; state2 < states; state2++)
			{
				local_sum += fw.at(level).at(state2);
			}

			if(!(local_sum > 0))
			{
				cout << "local_sum: " << local_sum << "\n";
				cout << "level: " << level << "\n";
			}
			assert(local_sum > 0);

			double factor = 1.0/local_sum;

			for(int state2 = 0; state2 < states; state2++)
			{
				fw.at(level).at(state2) = fw.at(level).at(state2)*factor;
				assert(fw.at(level).at(state2) >= 0);
			}

			fw_underflow_factor *= factor;
		}

		if(level == (levels-1))
		{
			for(int state = 0; state < states; state++)
			{
				emission_probability += fw.at(level).at(state);
			}
		}
		*/

			double local_max = 0;
			for(int state2 = 0; state2 < states; state2++)
			{
				if((state2 == 0) || (fw.at(level).at(state2) > local_max))
				{
					local_max = fw.at(level).at(state2);
				}
			}

			for(int state2 = 0; state2 < states; state2++)
			{
				fw.at(level).at(state2) = exp(fw.at(level).at(state2)-local_max);
				assert(fw.at(level).at(state2) >= 0);
			}

			fw_underflow_factor *= exp(local_max);

			if(level == (levels-1))
			{
				for(int state = 0; state < states; state++)
				{
					emission_probability += fw.at(level).at(state);
				}
			}
	}

	// emission_probability = emission_probability * (1.0/fw_underflow_factor);

	return emission_probability;


}
