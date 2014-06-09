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


void AlphaHMM::init(MultiGraph* g, int pGenotypingMode, bool verbose)
{

	myGraph = g;
	genotypingMode = pGenotypingMode;

	Edge2Int.clear();
	Edge2Level.clear();

	int levels = g->NodesPerLevel.size();
	haploidStatesByLevel.resize(levels-1);
	Int2Edge.resize(levels);

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

				//assert(e->multiEmission_full.size() == 0);

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

			if(!(diploidStates > 0))
			{
				cout << "Error!\ndiploidStates: " << diploidStates << "\nhaploidStates: " << haploidStates << "\n\n" << flush;
			}

			assert(diploidStates > 0);

			diploidStatesByLevel.at(level) = diploidStates;

			if(genotypingMode < 5)
			{
				fillStateTransitions(level);
			}
		}
	}

	maxStateKMerCount = 0;

	int min_states = 0;
	int max_states = 0;
	int sum_states = 0;

	for(int level = 0; level <= levels-2; level++)
	{
		int diploidStates = diploidStatesByLevel.at(level);

		assert(diploidStates > 0);
		int chunk_size = diploidStates / CONFIG.threads;
		assert(chunk_size >= 0);

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		assert(diploidStates > 0);

		if(level == 0)
		{
			min_states = diploidStates;
			max_states = diploidStates;
		}
		if(diploidStates > max_states)
		{
			max_states = diploidStates;
		}
		if(diploidStates < min_states)
		{
			min_states = diploidStates;
		}
		sum_states += diploidStates;

		assert(diploidStates > 0);


		#pragma omp parallel
		{
			assert(omp_get_num_threads() == CONFIG.threads);
			int thisThread = omp_get_thread_num();

			assert(thisThread >= 0);
			assert(chunk_size >= 0);

			int firstPair = thisThread * chunk_size;
			int lastPair = (thisThread+1) * chunk_size - 1;

			if((thisThread == (CONFIG.threads-1)) && (lastPair < (diploidStates-1)))
			{
				assert(diploidStates > 0);
				lastPair = diploidStates - 1;
			}

			assert(lastPair <= (diploidStates - 1));
			assert(diploidStates > 0);
			//assert(lastPair >= 0);

			int local_maxStateKMerCount = 0;
			set<int> local_interestingKMerSymbols;

			for(int state = firstPair; state <= lastPair; state++)
			{
				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				assert(diploidStates > 0);
				assert(lastPair >= 0);


				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
				assert(e1->locus_id == e2->locus_id);

				assert(diploidStates > 0);
				assert(lastPair >= 0);


				for(map<int, int>::iterator emissionIt = e1->multiEmission.begin(); emissionIt != e1->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(diploidStates > 0);
					assert(lastPair >= 0);


					if((emissionSymbol != gap_symbol) && (emissionSymbol != star_symbol))
					{
						assert(number >= 0);
						if(local_interestingKMerSymbols.count(emissionSymbol) == 0)
						{
							local_interestingKMerSymbols.insert(emissionSymbol);
						}

						assert(diploidStates > 0);
						assert(lastPair >= 0);

						if(e2->multiEmission.count(emissionSymbol) > 0)
						{
							number += e2->multiEmission[emissionSymbol];
						}

						assert(diploidStates > 0);
						assert(lastPair >= 0);

						if(local_maxStateKMerCount < number)
						{
							local_maxStateKMerCount = number;
						}

						assert(diploidStates > 0);
						assert(lastPair >= 0);
					}
				}

				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(diploidStates > 0);
					assert(lastPair >= 0);

					if((emissionSymbol != gap_symbol) && (emissionSymbol != star_symbol))
					{
						assert(number >= 0);
						if(local_interestingKMerSymbols.count(emissionSymbol) == 0)
						{
							local_interestingKMerSymbols.insert(emissionSymbol);
						}

						assert(diploidStates > 0);
						assert(lastPair >= 0);

						if(local_maxStateKMerCount < number)
						{
							local_maxStateKMerCount = number;
						}

						assert(diploidStates > 0);
						assert(lastPair >= 0);
					}
				}

			}

			#pragma omp critical
			{
				assert(diploidStates > 0);
				if(lastPair >= 0)
				{
					if(maxStateKMerCount < local_maxStateKMerCount)
					{
						maxStateKMerCount = local_maxStateKMerCount;
					}

					assert(diploidStates > 0);
					assert(lastPair >= 0);

					for(set<int>::iterator kMerSymbolIt = local_interestingKMerSymbols.begin(); kMerSymbolIt != local_interestingKMerSymbols.end(); kMerSymbolIt++)
					{
						assert(diploidStates > 0);
						assert(lastPair >= 0);

						int kMerSymbol = *kMerSymbolIt;

						assert(diploidStates > 0);
						assert(lastPair >= 0);

						string kMer = myGraph->CODE.deCode(locusID, kMerSymbol);

						assert(diploidStates > 0);
						assert(lastPair >= 0);

						interestingKMerStrings.insert(kMer);

						assert(diploidStates > 0);
						assert(lastPair >= 0);

					}
				}
			}
		}
	}

	haploid_edge_lengths.resize(haploidStatesByLevel.size());
	haploid_edge_lengths_valid.resize(haploidStatesByLevel.size());

	for(int level = 0; level < (int)haploidStatesByLevel.size(); level++)
	{
		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		int haploidStates = haploidStatesByLevel.at(level).size();

		haploid_edge_lengths_valid.at(level).resize(haploidStates);

		int expected_length_total = -1;

		for(int hI = 0; hI < haploidStates; hI++)
		{
			Edge* e = haploidStatesByLevel.at(level).at(hI).e;
			assert(e->locus_id == locusID);

			int length_total = 0;
			int length_valid = 0;

			for(map<int, int>::iterator emissionIt = e->multiEmission.begin(); emissionIt != e->multiEmission.end(); emissionIt++)
			{
				int emissionSymbol = emissionIt->first;
				int number = emissionIt->second;

				length_total += number;

				if((emissionSymbol != gap_symbol) && (emissionSymbol != star_symbol))
				{
					length_valid += number;
				}
			}

			if(expected_length_total == -1)
			{
				expected_length_total = length_total;
				haploid_edge_lengths.at(level) = length_total;
			}

			assert(expected_length_total == length_total);

			haploid_edge_lengths_valid.at(level).at(hI) = length_valid;
		}
	}

	if(verbose)
	{
		double average_states_per_level = (double) sum_states / (double) (levels-1);

		cout << "Levels in HMM: " << levels-1 << "\n";
		cout << "\t min_states (per level): " << min_states << "\n";
		cout << "\t max_states (per level): " << max_states << "\n";
		cout << "\t avg_states (per level): " << average_states_per_level << "\n";
	}

}

AlphaHMM::AlphaHMM(MultiGraph* g, int pGenotypingMode, bool verbose)
{
	init(g, pGenotypingMode, verbose);
}


void AlphaHMM::fillStateTransitions(int level)
{
	int levels = myGraph->NodesPerLevel.size();

	assert((int)haploidStatesByLevel.size() > level);
	assert((int)diploidStateTransitions.size() > level);

	int haploidStates = haploidStatesByLevel.at(level).size();
	int diploidStates = int(pow((double)haploidStates, 2)+.005);

	diploidStateTransitions.at(level).clear();
	diploidStateTransitions.at(level).resize(diploidStates);

	if(level == 0)
	{
		diploidStateTransitions_Reverse.at(level).clear();
		diploidStateTransitions_Reverse.at(level).resize(diploidStates);
	}
	if(level < (levels-2))
	{
		int diploidStates_n1 = int(pow((double)haploidStatesByLevel.at(level+1).size(), 2)+0.005);
		diploidStateTransitions_Reverse.at(level+1).clear();
		diploidStateTransitions_Reverse.at(level+1).resize(diploidStates_n1);
	}

	for(int stateI = 0; stateI < diploidStates; stateI++)
	{
		int s1 = stateI % haploidStates;
		int s2 = stateI / haploidStates;

		Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
		Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;
		assert(e1->locus_id == e2->locus_id);

		assert(diploidStateTransitions.at(level).at(stateI).size() == 0);

		if(level == 0)
		{
			diploid_initialProbabilities.push_back(haploid_initialProbabilities.at(s1) * haploid_initialProbabilities.at(s2));
		}

		if(level < (levels-2))
		{
			double summed_jumps = 0.0;
			int haploidStates_nextLevel = haploidStatesByLevel.at(level+1).size();

			assert((int)haploidStatesByLevel.size() > level);
			assert((int)haploidStatesByLevel.at(level).size() > s1);

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

map<string, long long> AlphaHMM::estimateEmissions(string kMerCountsSamplePath, string kMerCountsGenomePath, int kMerSize)
{
	// this function reads in kMer counts from kMerCountsSamplePath

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
				assert((int)fields.at(0).length() == kMerSize);
			}

			kMerCountsInSample[fields.at(0)] = Utilities::StrtoLongLong(fields.at(1));
		}
		kMerCountsSampleFile.close();
	}
	else
	{
		errEx("Cannot open kMer counts file: "+kMerCountsSamplePath);
	}

	ifstream kMerCountsGenomeFile;
	kMerCountsGenomeFile.open (kMerCountsGenomePath.c_str(), ios::in);
	lineCounter = 0;
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
				errEx("Strange format in kMerGenomeCountsFile file. Expect two fields (kMer and count). Problem occured in line: "+Utilities::ItoStr(lineCounter));

			if((fields.at(0) != "MeanReadLen") && (fields.at(0) != "TotalKMerCoverage") && (fields.at(0) != "TotalSeq"))
			{
				assert((int)fields.at(0).length() == kMerSize);

				if(Utilities::StrtoLongLong(fields.at(1)) != 0)
				{
					kMerCountsInSample[fields.at(0)] = -1;
				}
			}

		}
		kMerCountsGenomeFile.close();
	}
	else
	{
		// cerr << "Cannot open genome counts file: "+kMerCountsGenomePath+" -- this is not a problem if you are not using a Cortex-like likelihood!\n\n" << flush;
	}


	return kMerCountsInSample;
}


AlphaHMM::AlphaHMM(MultiGraph* g, bool verbose)
{
	init(g, 4, verbose);
}

int AlphaHMM::haploidEdgeLength(int level)
{
	return haploid_edge_lengths.at(level);
}

int AlphaHMM::haploidEdgeValidLength(int level, int diploidIndex, int e12)
{
	int haploidStates = haploidStatesByLevel.at(level).size();
	int s1 = diploidIndex % haploidStates;
	int s2 = diploidIndex / haploidStates;
	if(e12 == 1)
	{
		return haploid_edge_lengths_valid.at(level).at(s1);
	}
	else if(e12 == 2)
	{
		return haploid_edge_lengths_valid.at(level).at(s2);
	}
	else
	{
		assert(1 == 0);
		return 0;
	}
}


vector<diploidEdgePointerPath> AlphaHMM::retrieveAllBestPaths(map<string, long long>& observedEmissions)
{
	return retrieveAllBestPaths(observedEmissions, false);
}

vector<diploidEdgePointerPath> AlphaHMM::retrieveAllBestPaths(map<string, long long>& observedEmissions, bool getAllPaths)
{
	int levels = diploidStatesByLevel.size();

	vector<int> equally_good_states;
	int max_state = -1;
	double max_state_P = 0;

	for(int s = 0; s < (int)fw_Viterbi.at(levels-1).size(); s++)
	{
		double local_p = fw_Viterbi.at(levels-1).at(s);

		if((local_p > max_state_P) || (s == 0))
		{
			max_state = s;
			max_state_P = local_p;
			equally_good_states.clear();
			equally_good_states.push_back(s);
		}

		if((local_p == max_state_P) && (s != 0))
		{
			equally_good_states.push_back(s);
		}
	}

	vector< vector<int> > diploidStates;
	//for(int i = 0; i < equally_good_states.size(); i++)
	for(int i = 0; i < 1; i++)
	{
		assert((int)equally_good_states.size() > i);
		vector< vector<int> > recursionReturn = retrieveAllBestPaths_rec(levels-1, equally_good_states.at(i), observedEmissions, getAllPaths);
		diploidStates.insert(diploidStates.end(), recursionReturn.begin(), recursionReturn.end());
	}

	vector<diploidEdgePointerPath> forReturn;
	for(unsigned int i = 0; i < diploidStates.size(); i++)
	{
		diploidEdgePointerPath path;
		assert((int)diploidStates.at(i).size() == levels);
		for(int level = 0; level < levels; level++)
		{
			int diploidState = diploidStates.at(i).at(level);

			int haploidStates = haploidStatesByLevel.at(level).size();
			int s1 = diploidState % haploidStates;
			int s2 = diploidState / haploidStates;

			Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
			Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;

			path.h1.push_back(e1);
			path.h2.push_back(e2);
		}

		forReturn.push_back(path);
	}

	return forReturn;
}

vector<diploidEdgePath> AlphaHMM::retrieveAllBestPathsLargeGraph(map<string, long long>& observedEmissions)
{
	vector<diploidEdgePointerPath> paths_multiGraph = retrieveAllBestPaths(observedEmissions);
	vector<diploidEdgePath> forReturn;
	for(unsigned int i = 0; i < paths_multiGraph.size(); i++)
	{
		diploidEdgePointerPath path_largeGraph = myGraph->MultiGraphEdgesToLargeGraphEdges(paths_multiGraph.at(i));
		diploidEdgePath path_largeGraph_indices = myGraph->underlyingLargeGraph->edgePointerPathToEdgeIndexPath(path_largeGraph);
		forReturn.push_back(path_largeGraph_indices);
	}
	return forReturn;
}

vector<diploidNucleotidePath> AlphaHMM::retrieveAllBestNucleotidePaths(map<string, long long>& observedEmissions)
{
	vector<diploidEdgePath> allbestEdgePathsLargeGraph = retrieveAllBestPathsLargeGraph(observedEmissions);
	vector<diploidNucleotidePath> forReturn;
	for(int pI = 0; pI < (int)allbestEdgePathsLargeGraph.size(); pI++)
	{
		diploidEdgePath path_lG = allbestEdgePathsLargeGraph.at(pI);
		diploidNucleotidePath path_nucleotides;

		assert(path_lG.h1.size() == path_lG.h2.size());

		path_nucleotides.h1 = myGraph->underlyingLargeGraph->haploidPathToNucleotides(path_lG.h1);
		path_nucleotides.h2 = myGraph->underlyingLargeGraph->haploidPathToNucleotides(path_lG.h2);

		forReturn.push_back(path_nucleotides);
	}

	return forReturn;
}

int AlphaHMM::getGenotypingMode()
{
	return genotypingMode;
}

vector< vector<int> > AlphaHMM::retrieveAllBestPaths_rec(int level, int runningState, map<string, long long>& observedEmissions, bool getAllPaths)
{
	vector< vector<int> > forReturn;

	vector<int> followStates;

	if(genotypingMode >= 5)
	{
		assert((int)fw_Viterbi_backtrack_int.size() > level);
		assert((int)fw_Viterbi_backtrack_int.at(level).size() > runningState);

		followStates.resize(1);

		assert(level >= 0);
		assert(runningState >= 0);
		assert(followStates.size() == 1);
		followStates.at(0) = fw_Viterbi_backtrack_int.at(level).at(runningState);
	}
	else
	{
		followStates = fw_Viterbi_backtrack.at(level).at(runningState);
	}

	if(level > 0)
	{

		//for(int i = 0; i < followStates.size(); i++)

		if(getAllPaths && (followStates.size() > 1))
		{
			cout << "Level " << level << ", would like to descend into " << followStates.size() << " levels of recursion.\n";

			for(int i = 0; i < (int)followStates.size(); i++)
			{
				int haploidStates_previousLevel = haploidStatesByLevel.at(level-1).size();
				int previousLevel = followStates.at(i);
				int s1_previous = previousLevel % haploidStates_previousLevel;
				int s2_previous = previousLevel / haploidStates_previousLevel;

				Edge* e1 = haploidStatesByLevel.at(level-1).at(s1_previous).e;
				Edge* e2 = haploidStatesByLevel.at(level-1).at(s2_previous).e;

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

				cout << "\t alternative " << i << ", is underlying state " << previousLevel << "\n";
				cout << "\t\t" << s1_previous << "/" << s2_previous << "\n";
				cout << "\t\t" << e1->label << "/" << e2->label << "\n";


				for(map<int, int>::iterator emissionIt = combinedEdgeEmission.begin(); emissionIt != combinedEdgeEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;
					string kMer = myGraph->CODE.deCode(e1->locus_id, emissionSymbol);
					cout << "\t\t\t" << emissionSymbol << " combinedEdge: " << number << ", observed in genome: " << observedEmissions[kMer] << "   [" << kMer <<  "]\n";
				}
			}
		}

		assert(followStates.size() > 0);
		if(followStates.size() > 0)
		{
			for(int i = 0; i < 1; i++)
			{
				vector< vector<int> > recursionReturn = retrieveAllBestPaths_rec(level-1, followStates.at(i), observedEmissions, getAllPaths);
				forReturn.insert(forReturn.end(), recursionReturn.begin(), recursionReturn.end());
			}
		}
	}
	else
	{
		forReturn.resize(1);
	}

	for(int retI = 0; retI < (int)forReturn.size(); retI++)
	{
		forReturn.at(retI).push_back(runningState);
	}

	return forReturn;

}

multiHaplotypePair_and_P AlphaHMM::retrieveViterbiSample(bool labelOnly)
{
	vector<int> backtrack_states;

	// determine Viterbi P (we do not get this from other functions)
	int lastLevel = diploidStatesByLevel.size() - 1;
	double viterbi_P = 0;

	int max_state = -1;
	double max_state_P = 0;
	for(int s = 0; s < (int)fw_Viterbi.at(lastLevel).size(); s++)
	{
		double local_p = fw_Viterbi.at(lastLevel).at(s);
		if((local_p > max_state_P) || (s == 0))
		{
			max_state = s;
			max_state_P = local_p;
		}
	}

	viterbi_P = max_state_P;

	map<string, long long> observedEmissions_null;
	vector<diploidEdgePointerPath> bestPaths = retrieveAllBestPaths(observedEmissions_null, false);
	assert(bestPaths.size() == 1);
	diploidEdgePointerPath bestPath = bestPaths.at(0);


	multiHaplotypePair_and_P forReturn;
	forReturn.P = viterbi_P;
	forReturn.diploidEdgePointers = bestPath;
	forReturn.haploPair = myGraph->edgePointerPathToLabels(bestPath);

	return forReturn;
}

double AlphaHMM::fillForwardBackwardTable(map<string, long long> globalEmission)
{
	if(genotypingMode == 0)
	{
		return fillForwardBackwardTable_0(globalEmission);
	}
	else if(genotypingMode == 1)
	{
		return fillForwardBackwardTable_1(globalEmission);
	}
	else if(genotypingMode == 2)
	{
		return fillForwardBackwardTable_2(globalEmission);
	}
	else if(genotypingMode == 3)
	{
		return fillForwardBackwardTable_3(globalEmission);
	}
	else if(genotypingMode == 4)
	{
		return fillForwardBackwardTable_4(globalEmission);
	}
	else if(genotypingMode == 5)
	{
		return fillForwardBackwardTable_5(globalEmission);
	}
	else if(genotypingMode == 6)
	{
		return fillForwardBackwardTable_6(globalEmission);
	}
	else if(genotypingMode == 7)
	{
		return fillForwardBackwardTable_7(globalEmission);
	}
	else if(genotypingMode == 8)
	{
		return fillForwardBackwardTable_8(globalEmission);
	}
	else
	{
		assert("Unknown genotypingMode for alphaHMM" == "");
		return 0;
	}
}

vector< vector<double> > AlphaHMM::populatekMerMatrix(double coverage, map<string, long long>& globalEmission)
{
	cout << "Expected haploid kMer coverage: " << coverage << "\n";
	cout << "Number of interesting kMers: " << interestingKMerStrings.size() << "\n";

	assert(interestingKMerStrings.size() > 0);

	int maxObservedKMer = 0;
	for(set<string>::iterator kMerIt = interestingKMerStrings.begin(); kMerIt != interestingKMerStrings.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		if(globalEmission.count(kMer) == 0)
		{
			cerr << "kMer " << kMer << " is not defined in coverage file!\n";
		}
		assert(globalEmission.count(kMer) > 0);
		int kMerCount = globalEmission[kMer];
		if(kMerCount > maxObservedKMer)
		{
			maxObservedKMer = kMerCount;
		}
	}
	interestingKMerStrings.clear();

	cout << "Maximum observed kMer count from sample: " << maxObservedKMer << "\n";
	cout << "Maximum count of single kMer in single state: " << maxStateKMerCount << "\n";

	//if((maxObservedKMer < maxStateKMerCount*coverage*1.5) || (maxObservedKMer > 100))
	//{
		maxObservedKMer = maxStateKMerCount*coverage*1.5;
		cout << "We manually fix the maxObservedkMer count to maxStateKMerCount*coverage*1.5: " << maxObservedKMer << "\n";
	//}
	// implication maxObservedKMer >= maxStateKMerCount*coverage*1.5

	vector< vector<double> > observedXunderlying;
	int maxUnderlyingKMer = (double)maxObservedKMer/coverage * 1.5+1;
	//if(maxUnderlyingKMer < (maxStateKMerCount*2))
	//{
	//	maxUnderlyingKMer = maxStateKMerCount*2;
	//}
	//if(maxObserve_forMatrix < maxUnderlyingKMer * coverage * 2)
	//{
	//	maxObserve_forMatrix = maxUnderlyingKMer * coverage * 2 + 1;
	//}
	//assert(maxUnderlyingKMer >= (maxStateKMerCount*2));
	//assert(maxUnderlyingKMer >= ((double)maxObservedKMer/coverage * 2));
	//assert(maxObserve_forMatrix >= (maxUnderlyingKMer * coverage * 2));

	cout << "Assume maximum number of possible underlying kMer count: " << maxUnderlyingKMer << "\n";
	cout << "Assume maximum observable kMer count: " << maxObservedKMer << "\n";

	int maxObserve_forMatrix = maxObservedKMer;

	cout << "Populate matrix " << maxObserve_forMatrix+1 << " x " << maxUnderlyingKMer+1 << " (X x Y)\n";
	observedXunderlying.resize(maxObserve_forMatrix+1);
	int chunk_size = (maxObserve_forMatrix+1) / CONFIG.threads;

	#pragma omp parallel
	{
		assert(omp_get_num_threads() == CONFIG.threads);
		int thisThread = omp_get_thread_num();
		int firstObserved = thisThread * chunk_size;
		int lastObserved = (thisThread+1) * chunk_size - 1;
		if((thisThread == (CONFIG.threads-1)) && (lastObserved < maxObserve_forMatrix))
		{
			lastObserved = maxObserve_forMatrix;
		}

		for(int observedNumber = firstObserved; observedNumber <= lastObserved; observedNumber++)
		{
			observedXunderlying.at(observedNumber).resize(maxUnderlyingKMer+1);
			double col_sum = 0;
			for(int underlyingCount = 0; underlyingCount <= maxUnderlyingKMer; underlyingCount++)
			{
				double lambda;
				if(underlyingCount == 0)
				{
					lambda = 0.01;
				}
				else
				{
					lambda = underlyingCount * coverage;
				}
				poisson_up poisson(lambda);
				double likelihood = pdf(poisson, observedNumber);
				col_sum += likelihood;
				observedXunderlying.at(observedNumber).at(underlyingCount) = likelihood;
			}
			assert(col_sum != 0);
			for(int underlyingCount = 0; underlyingCount <= maxUnderlyingKMer; underlyingCount++)
			{
				observedXunderlying.at(observedNumber).at(underlyingCount) = observedXunderlying.at(observedNumber).at(underlyingCount)/col_sum;
				assert(observedXunderlying.at(observedNumber).at(underlyingCount) >= 0);
				assert(observedXunderlying.at(observedNumber).at(underlyingCount) <= 1);
			}
		}
	}

	return observedXunderlying;

}

void AlphaHMM::debug(map<string, long long> globalEmission)
{
	vector<string> interestingSNPs;
	interestingSNPs.push_back("rs114471122");
	interestingSNPs.push_back("rs114407469");
	interestingSNPs.push_back("rs116159376");
	interestingSNPs.push_back("rs114979331");
	interestingSNPs.push_back("rs115815844");
	interestingSNPs.push_back("rs111631325");
	interestingSNPs.push_back("rs115538520");
	interestingSNPs.push_back("rs116340620");
	interestingSNPs.push_back("rs115593407");
	interestingSNPs.push_back("rs115760093");
	interestingSNPs.push_back("rs115919721");
	interestingSNPs.push_back("rs114310176");
	interestingSNPs.push_back("rs115689245");
	interestingSNPs.push_back("rs114631018");
	interestingSNPs.push_back("rs116685512");
	interestingSNPs.push_back("rs114171774");
	interestingSNPs.push_back("rs116629007");
	interestingSNPs.push_back("rs116429229");
	interestingSNPs.push_back("rs116667638");
	interestingSNPs.push_back("rs114607474");
	interestingSNPs.push_back("rs116114006");
	interestingSNPs.push_back("rs116348417");
	interestingSNPs.push_back("rs115156025");
	interestingSNPs.push_back("rs115314263");
	interestingSNPs.push_back("rs116793085");
	interestingSNPs.push_back("rs116302129");
	interestingSNPs.push_back("rs116377187");
	interestingSNPs.push_back("rs116340898");
	interestingSNPs.push_back("rs114618774");
	interestingSNPs.push_back("rs114673216");
	interestingSNPs.push_back("rs114525617");
	interestingSNPs.push_back("rs115495371");
	interestingSNPs.push_back("rs116839689");
	interestingSNPs.push_back("rs116182315");
	interestingSNPs.push_back("rs116348849");
	interestingSNPs.push_back("rs114524182");
	interestingSNPs.push_back("rs116189721");
	interestingSNPs.push_back("rs116102014");
	interestingSNPs.push_back("rs114931354");
	interestingSNPs.push_back("rs114716686");
	interestingSNPs.push_back("rs115033670");


	map<string, long long> noEmissions;
	vector<diploidEdgePointerPath> bestPaths_Edges = retrieveAllBestPaths(noEmissions);
	diploidEdgePointerPath firstEstimatedEdgePath = bestPaths_Edges.at(0);

	for(int level = firstEstimatedEdgePath.h1.size()-1; level >= 0; level--)
	{
		Edge* selected_e1 = firstEstimatedEdgePath.h1.at(level);
		Edge* selected_e2 = firstEstimatedEdgePath.h2.at(level);

		bool interestingState = false;
		string foundLocus;
		for(unsigned int sI = 0; sI < interestingSNPs.size(); sI++)
		{
			string locusID = interestingSNPs.at(sI);
			if(selected_e1->label.find(locusID) != string::npos)
			{
				interestingState = true;
				foundLocus = locusID;
			}
			if(selected_e2->label.find(locusID) != string::npos)
			{
				interestingState = true;
				foundLocus = locusID;
			}
		}

		if(interestingState)
		{
			cout << "Level " << level << "is interesting (for SNP " << foundLocus << " and possibly other!\n";

			set<string> utilizeKMersatLevel;
			set<int> utilizeCodedKMersatLevel;

			double likelihood_all_kMers_0 = 0;
			map<string, double> logL_count0;
			string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;

			assert(kMersInGraphInfo.Level2Kmers.count(level) > 0);
			for(set<string>::iterator kMerIt = kMersInGraphInfo.Level2Kmers.at(level).begin(); kMerIt != kMersInGraphInfo.Level2Kmers.at(level).end(); kMerIt++)
			{
				string kMer = *kMerIt;
				if(utilizekMers.count(kMer) > 0)
				{
					utilizeKMersatLevel.insert(kMer);
					utilizeCodedKMersatLevel.insert(myGraph->CODE.doCode(locusID, kMer));
					assert(globalEmission.count(kMer) > 0);
					poisson_up poisson(error_rate);

					double likelihood_thisMer_if_underlying_0 = pdf(poisson, globalEmission[kMer]);
					if(likelihood_thisMer_if_underlying_0 == 0)
					{
						likelihood_thisMer_if_underlying_0 = 1e-200;
					}
					double log_likelihood_thisMer_if_underlying_0 = log(likelihood_thisMer_if_underlying_0);
					assert(log_likelihood_thisMer_if_underlying_0 <= 0);
					likelihood_all_kMers_0 += log_likelihood_thisMer_if_underlying_0;
					logL_count0[kMer] = log_likelihood_thisMer_if_underlying_0;
				}
				else
				{
					assert(ignorekMers.count(kMer) > 0);
				}
			}

			assert(likelihood_all_kMers_0 <= 0);

			int states = diploidStatesByLevel.at(level);

			int gap_symbol = myGraph->CODE.doCode(locusID, "_");
			int star_symbol = myGraph->CODE.doCode(locusID, "*");

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

				map<int, map<int, double> > _cache_poisson_PDF;

				for(int state = firstPair; state <= lastPair; state++)
				{
					std::stringstream subst_cout;

					int haploidStates = haploidStatesByLevel.at(level).size();
					int s1 = state % haploidStates;
					int s2 = state / haploidStates;

					Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
					Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;

					subst_cout << "State " << s1 << "/" << s2 << "\n";
					subst_cout << "\tViterbi optimality: " << fw_Viterbi.at(level).at(state) << "\n";
					if((e1 == selected_e1) && (e2 == selected_e2))
					{
						subst_cout << "\tTHIS IS THE SELECTED STATE******************\n";
					}

					assert(e1->locus_id == e2->locus_id);

					map<int, int> combinedCodedEdgeEmission = e1->multiEmission;
					for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
					{
						int emissionSymbol = emissionIt->first;
						int number = emissionIt->second;

						assert(number >= 0);
						if(combinedCodedEdgeEmission.count(emissionSymbol) == 0)
							combinedCodedEdgeEmission[emissionSymbol] = 0;

						combinedCodedEdgeEmission[emissionSymbol] += number;
					}

					combinedCodedEdgeEmission.erase(gap_symbol);
					combinedCodedEdgeEmission.erase(star_symbol);

					double log_emission_p = likelihood_all_kMers_0;
					assert(log_emission_p <= 0);

					for(map<int, int>::iterator kMerIt = combinedCodedEdgeEmission.begin(); kMerIt != combinedCodedEdgeEmission.end(); kMerIt++)
					{
						string kMer = myGraph->CODE.deCode(e1->locus_id, kMerIt->first);

						if(utilizeCodedKMersatLevel.count(kMerIt->first) == 0)
						{
							subst_cout << "\t\t" << kMer << ": ignore.\n";
							continue;
						}

						int underlyingCopyCount = kMerIt->second;
						assert(underlyingCopyCount != 0);

						subst_cout << "\t\t" << kMer << ": " << underlyingCopyCount << " on underlying copies " << underlyingCopyCount << ".\n";


						assert(logL_count0.count(kMer) > 0);
						assert(logL_count0[kMer] <= 0);

						log_emission_p -= logL_count0[kMer];
						assert(log_emission_p <= 0);

						double rate = (double)underlyingCopyCount * coverage;
						assert(globalEmission.count(kMer) > 0);

						double thiskMer_P;
						if((_cache_poisson_PDF.count(underlyingCopyCount) > 0) && (_cache_poisson_PDF.at(underlyingCopyCount).count(globalEmission[kMer]) > 0))
						{
							thiskMer_P = _cache_poisson_PDF[underlyingCopyCount][globalEmission[kMer]];
						}
						else
						{
							poisson_up poisson(rate);
							thiskMer_P = pdf(poisson, globalEmission[kMer]);
							_cache_poisson_PDF[underlyingCopyCount][globalEmission[kMer]] = thiskMer_P;
						}

						if(thiskMer_P == 0)
						{
							thiskMer_P = 1e-200;
						}

						assert(thiskMer_P > 0);
						assert(thiskMer_P <= 1);

						log_emission_p += log(thiskMer_P);
						assert(log_emission_p <= 0);
					}

					subst_cout << "\tlog likelihood: " << log_emission_p << "\n\n";

					#pragma omp critical
					{
						cout << subst_cout.str();
					}
				}
			}

		}
	}
}

double AlphaHMM::fillForwardBackwardTable_8(map<string, long long> globalEmission)
{
	assert(genotypingMode == 8);

	int levels = diploidStatesByLevel.size();

	fw.resize(levels);
	fw_Viterbi.resize(levels);
	fw_Viterbi_backtrack_int.resize(levels);
	fw_Emission.resize(levels);

	fw_underflow_factor = 0;

	double coverage_wG = estimateCoverage(globalEmission);

	vector< vector<double> > observedXunderlying = populatekMerMatrix(coverage_wG, globalEmission);
	map <int, double> stateKMerCount_expectedMissing_H0;

	kMersInGraphInfo = myGraph->getkMerPositions();
	kMersInGraph = kMersInGraphInfo.kMers;
	uniqueMers = myGraph->kMerUniqueness();

	// find out which kMers we want to keep

	int kickedOutBecauseGraphDuplicate = 0;
	int kickedOutBecauseWholeGenomeDuplicate = 0;
	int kickedOutBecauseAssumeOtherDuplication = 0;
	int kickedOutBecauseWayTooManyReads = 0;

	for(set<string>::iterator kMerIt = kMersInGraph.begin(); kMerIt != kMersInGraph.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		bool usekMer = true;

		// non-unique within graph
		if(uniqueMers.levelUniqueKMers.count(kMer) == 0)
		{
			usekMer = false;
			kickedOutBecauseGraphDuplicate++;
		}

		// non-unique within genome
		assert(globalEmission.count(kMer) > 0);
		if(globalEmission[kMer] == -1)
		{
			kickedOutBecauseWholeGenomeDuplicate++;
			usekMer = false;
		}

		// estimate number of underlying copies - consistent with
		// our expectation from graph?

		assert(kMersInGraphInfo.kMerMaxPerLevel.count(kMer) > 0);
		int maxLevelCount = kMersInGraphInfo.kMerMaxPerLevel[kMer];
		assert(maxLevelCount > 0);
		assert(maxLevelCount < maxStateKMerCount);

		long long observedReadsOnKMer = globalEmission[kMer];
		if(globalEmission[kMer] != -1)
		{
			if(observedReadsOnKMer < (int)observedXunderlying.size())
			{
				double p_within_range = 0.0;
				for(int possibleUnderlyingCount = 0; possibleUnderlyingCount <= maxLevelCount+2; possibleUnderlyingCount++)
				{
					double p = observedXunderlying.at(observedReadsOnKMer).at(possibleUnderlyingCount);
					p_within_range += p;
				}

				if(p_within_range < 0.5)
				{
					kickedOutBecauseAssumeOtherDuplication++;
					usekMer = false;
				}
			}
			else
			{
				kickedOutBecauseWayTooManyReads++;
				usekMer = false;
			}
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

	// summary stats of kMer usage
	cout << "AlphaHMM in mode 8!\n";
	cout << levels << " Levels\n";
	cout << "Utilize " << utilizekMers.size() << " of " << kMersInGraph.size() << "kMers\n";
	cout << "\t Kicked out because (multiple reasons possible):\n";
	cout << "\t\tDuplication in graph: " << kickedOutBecauseGraphDuplicate << "\n";
	cout << "\t\tDuplication in genome outside xMHC: " << kickedOutBecauseWholeGenomeDuplicate << "\n";
	cout << "\t\tPosterior copy count estimate exceeds what we expect: " << kickedOutBecauseAssumeOtherDuplication << "\n";
	cout << "\t\tNumber of reads too high ( read number > " << observedXunderlying.size() << ", maxStateKMerCount " << maxStateKMerCount << ") : " << kickedOutBecauseWayTooManyReads << "\n";

	// assert(kickedOutBecauseWholeGenomeDuplicate > 0);

	cout << uniqueMers.levelUniqueKMers.size() << " graph-unique kMers\n";

	// use graph to estimate coverage

	int kMers_used_for_coverage_estimation = 0;
	int levels_used_for_coverage_estimation = 0;
	double sum_kMer_coverage = 0.00;

	for(int level = 0; level < levels; level++)
	{
		if(haploidStatesByLevel.at(level).size() == 1)
		{
			levels_used_for_coverage_estimation++;

			string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
			int gap_symbol = myGraph->CODE.doCode(locusID, "_");
			int star_symbol = myGraph->CODE.doCode(locusID, "*");

			Edge* e1 = haploidStatesByLevel.at(level).at(0).e;
			map<int, int> edgeEmissions = e1->multiEmission;
			edgeEmissions.erase(gap_symbol);
			edgeEmissions.erase(star_symbol);

			for(map<int, int>::iterator emissionIt = edgeEmissions.begin(); emissionIt != edgeEmissions.end(); emissionIt++)
			{
				string kMer = myGraph->CODE.deCode(locusID, emissionIt->first);
				if(ignorekMers.count(kMer) > 0)
				{
					continue;
				}
				else
				{
					kMers_used_for_coverage_estimation += emissionIt->second;
					assert(globalEmission.count(kMer) > 0);
					assert(utilizekMers.count(kMer) > 0);
					sum_kMer_coverage += globalEmission[kMer];
				}
			}
		}
	}

	if(kMers_used_for_coverage_estimation > 0)
	{
		coverage = 0.5 * sum_kMer_coverage / (double)kMers_used_for_coverage_estimation ;
	}
	else
	{
		cout << "\nCannot estimate coverage from graph -- use whole-genome estimate";
		coverage = coverage_wG;
	}


	cout << "\n\nRate estimation\n";
	cout << "\tHaploid kMer coverage estimated from graph: " << coverage << ", based on " << kMers_used_for_coverage_estimation << " kMer copies and " << levels_used_for_coverage_estimation << "/" << levels << " levels \n";
	cout << "\t\t\t (whole-genome estimate: " << coverage_wG << ")\n";

	// use graph to estimate error rates

	error_rate = 0.1;
	bool update_error_rate = true;

	double previous_likelihood = 1;
	set<int> considerLevelsForErrorEstimation;

	while(update_error_rate == true)
	{
		double current_likelihood = 0;

		int kMers_used_for_error_estimation = 0;
		int levels_used_for_error_estimation = 0;
		double sum_error_kMer_rates = 0.00;

		for(int level = 0; level < levels; level++)
		{
			if((previous_likelihood != 1) && (considerLevelsForErrorEstimation.count(level) == 0))
			{
				continue;
			}

			if(haploidStatesByLevel.at(level).size() == 2)
			{
				string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
				int gap_symbol = myGraph->CODE.doCode(locusID, "_");
				int star_symbol = myGraph->CODE.doCode(locusID, "*");

				vector<int> use_kMers;
				use_kMers.resize(haploidStatesByLevel.at(level).size());
				for(int k = 0; k < (int)use_kMers.size(); k++)
				{
					use_kMers.at(k) = 0;
				}

				for(int eI = 0; eI < (int)haploidStatesByLevel.at(level).size(); eI++)
				{
					Edge* e = haploidStatesByLevel.at(level).at(eI).e;
					map<int, int> edgeEmissions = e->multiEmission;
					edgeEmissions.erase(gap_symbol);
					edgeEmissions.erase(star_symbol);

					map<int, int> exclusive_kMers_local = edgeEmissions;

					for(map<int, int>::iterator emissionIt = edgeEmissions.begin(); emissionIt != edgeEmissions.end(); emissionIt++)
					{
						string kMer = myGraph->CODE.deCode(locusID, emissionIt->first);
						if(utilizekMers.count(kMer) > 0)
						{
							use_kMers.at(eI) += emissionIt->second;
						}
					}
				}

				int use_kMers_number = -1;
				bool all_edges_equal_use = true;
				for(int eI = 0; eI < (int)haploidStatesByLevel.at(level).size(); eI++)
				{
					if(use_kMers_number == -1)
					{
						use_kMers_number = use_kMers.at(eI);
					}
					else
					{
						if(use_kMers.at(eI) != use_kMers_number)
						{
							all_edges_equal_use = false;
						}
					}
				}

				if(all_edges_equal_use)
				{
					considerLevelsForErrorEstimation.insert(level);

					//cout << "Level " << level << " error rate estimation\n";

					double max_emission = 1;
					int max_p1 = -1;
					int max_p2 = -1;
					double max_error_rate = -1;
					int max_error_rate_weight = -1;

					for(int p1 = 0; p1 < (int)haploidStatesByLevel.at(level).size(); p1++)
					{
						for(int p2 = p1; p2 < (int)haploidStatesByLevel.at(level).size(); p2++)
						{
							Edge* e1 = haploidStatesByLevel.at(level).at(p1).e;
							map<int, int> edgeEmissions = e1->multiEmission;
							Edge* e2 = haploidStatesByLevel.at(level).at(p2).e;

							for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
							{
								int emissionSymbol = emissionIt->first;
								int emissionCount = emissionIt->second;
								if(edgeEmissions.count(emissionSymbol) == 0)
								{
									edgeEmissions[emissionSymbol] = 0;
								}
								edgeEmissions[emissionSymbol] += emissionCount;
							}

							edgeEmissions.erase(gap_symbol);
							edgeEmissions.erase(star_symbol);

							// get likelihood for this pair of edges

							// ... kMers implied...

							double log_emission_p = 0;
							int _used_kMers = 0;
							for(map<int, int>::iterator kMerIt = edgeEmissions.begin(); kMerIt != edgeEmissions.end(); kMerIt++)
							{
								string kMer = myGraph->CODE.deCode(locusID, kMerIt->first);

								if(utilizekMers.count(kMer) > 0)
								{
									int underlyingCopyCount = kMerIt->second;
									double rate = (double)underlyingCopyCount * coverage;

									poisson_up poisson(rate);
									assert(globalEmission.count(kMer) > 0);
									double thiskMer_P = pdf(poisson, globalEmission[kMer]);

									assert(thiskMer_P >= 0);
									assert(thiskMer_P <= 1);

									log_emission_p += log(thiskMer_P);

									_used_kMers++;
								}
							}

							// .. other kMers from this level not covered by pair of edges

							double other_coverage = 0;
							int other_kMers = 0;

							for(set<string>::iterator kMerIt = kMersInGraphInfo.Level2Kmers.at(level).begin(); kMerIt != kMersInGraphInfo.Level2Kmers.at(level).end(); kMerIt++)
							{
								string kMer = *kMerIt;

								if(utilizekMers.count(kMer) > 0)
								{
									int symbol = myGraph->CODE.doCode(locusID, kMer);

									if(edgeEmissions.count(symbol) == 0)
									{
										poisson_up poisson(error_rate);
										assert(globalEmission.count(kMer) > 0);
										double thiskMer_P = pdf(poisson, globalEmission[kMer]);

										assert(thiskMer_P >= 0);
										assert(thiskMer_P <= 1);

										log_emission_p += log(thiskMer_P);

										other_coverage += globalEmission[kMer];
										other_kMers++;

									}
								}
							}

							double other_kMers_rate;
							if(other_kMers == 0)
							{
								other_kMers_rate = 0;
							}
							else
							{
								if(other_coverage == 0)
								{
									other_kMers_rate = 0.00000000001;
								}
								else
								{
									other_kMers_rate = other_coverage / (double)other_kMers;
								}
							}

							// cout << "\tPair " << p1 << " and " << p2 << " log likelihood " << log_emission_p << " based on " << _used_kMers << " kMers (assumed error rate " << other_kMers_rate << "\n";

							if((max_emission == 1) || (log_emission_p > max_emission))
							{
								max_emission = log_emission_p;
								max_p1 = p1;
								max_p2 = p2;
								max_error_rate = other_kMers_rate;
								max_error_rate_weight = other_kMers;
							}
						}
					}

					//cout << "\tSELECTED ERROR RATE: " << max_error_rate << "\n";

					current_likelihood += max_emission;


					//cout << "\thaploidStatesByLevel.at(level).size(): " << haploidStatesByLevel.at(level).size() << "\n";
					//cout << "\tSelected max pair p1 " << max_p1 << " and p2 " << max_p2 << "\n";

					if(max_error_rate_weight != 0)
					{
						sum_error_kMer_rates += max_error_rate;
						levels_used_for_error_estimation++;
						kMers_used_for_error_estimation += max_error_rate_weight;

						if(levels_used_for_error_estimation > 200)
						{
							break;
						}
					}
				}
			}
		}

		if(kMers_used_for_error_estimation == 0)
		{
			update_error_rate = false;
			cout << "\nCannot estimate error rates from the graph - use fixed estimate!\n";
			error_rate = 0.0001;
		}
		else
		{
			// cout << "Current likelihood with error rate " << error_rate << ": " << current_likelihood << ", based on " << levels_used_for_error_estimation << "/" << levels << " levels and " << kMers_used_for_error_estimation << " kMers \n";

			if(((current_likelihood - previous_likelihood) < 2) && (previous_likelihood != 1))
			{
				update_error_rate = false;
			}
			else
			{
				error_rate = sum_error_kMer_rates / (double)kMers_used_for_error_estimation;
				// cout << "\terror rate updated to " << error_rate << "\n";
				update_error_rate = true;
			}

			previous_likelihood = current_likelihood;
		}
	}

	cout << "Error rate estimated from graph: " << error_rate << "\n\n";

	for(int level = 0; level < levels; level++)
	{
		set<string> utilizeKMersatLevel;
		set<int> utilizeCodedKMersatLevel;

		double likelihood_all_kMers_0 = 0;
		map<string, double> logL_count0;

		string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;

		assert(kMersInGraphInfo.Level2Kmers.count(level) > 0);
		for(set<string>::iterator kMerIt = kMersInGraphInfo.Level2Kmers.at(level).begin(); kMerIt != kMersInGraphInfo.Level2Kmers.at(level).end(); kMerIt++)
		{
			string kMer = *kMerIt;
			if(utilizekMers.count(kMer) > 0)
			{
				utilizeKMersatLevel.insert(kMer);
				utilizeCodedKMersatLevel.insert(myGraph->CODE.doCode(locusID, kMer));
				assert(globalEmission.count(kMer) > 0);
				poisson_up poisson(error_rate);



				double likelihood_thisMer_if_underlying_0 = pdf(poisson, globalEmission[kMer]);
				if(likelihood_thisMer_if_underlying_0 == 0)
				{
					likelihood_thisMer_if_underlying_0 = 1e-200;
				}
				double log_likelihood_thisMer_if_underlying_0 = log(likelihood_thisMer_if_underlying_0);
				assert(log_likelihood_thisMer_if_underlying_0 <= 0);
				likelihood_all_kMers_0 += log_likelihood_thisMer_if_underlying_0;
				logL_count0[kMer] = log_likelihood_thisMer_if_underlying_0;
			}
			else
			{
				assert(ignorekMers.count(kMer) > 0);
			}
		}

		assert(likelihood_all_kMers_0 <= 0);


		if((level % 1000) == 0)
		{
			int startLevelFill = level;
			int stopLevelFill = level + 1000 - 1;

			cout << "Fill state transitions table from " << startLevelFill << " to " << stopLevelFill << "\n" << flush;

			int previousStepFillStop = startLevelFill - 1;
			int previousStepFillStart = previousStepFillStop - 1000 + 1;

			fillStateTransitions(startLevelFill, stopLevelFill, previousStepFillStart, previousStepFillStop);
		}

		int states = diploidStatesByLevel.at(level);

		fw.at(level).resize(states);
		fw_Viterbi.at(level).resize(states);
		fw_Emission.at(level).resize(states);
		fw_Viterbi_backtrack_int.at(level).resize(states);

		//string locusID = haploidStatesByLevel.at(level).at(0).e->locus_id;
		int gap_symbol = myGraph->CODE.doCode(locusID, "_");
		int star_symbol = myGraph->CODE.doCode(locusID, "*");

		if((level % 1000) == 0)
			cout << "fillForwardBackwardTable level " << level << "/" << (levels-1) << " (first thread) \n" << flush;

		int chunk_size = states / CONFIG.threads;

		double max_forward = 0;
		bool set_max_forward_first_value = false;

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

			map<int, map<int, double> > _cache_poisson_PDF;
			double local_max_forward = 0;
			bool have_local_max_forward_value = false;

			for(int state = firstPair; state <= lastPair; state++)
			{

				std::stringstream subst_cout;

				int haploidStates = haploidStatesByLevel.at(level).size();
				int s1 = state % haploidStates;
				int s2 = state / haploidStates;

				Edge* e1 = haploidStatesByLevel.at(level).at(s1).e;
				Edge* e2 = haploidStatesByLevel.at(level).at(s2).e;

				assert(e1->locus_id == e2->locus_id);

				map<int, int> combinedCodedEdgeEmission = e1->multiEmission;
				for(map<int, int>::iterator emissionIt = e2->multiEmission.begin(); emissionIt != e2->multiEmission.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;

					assert(number >= 0);
					if(combinedCodedEdgeEmission.count(emissionSymbol) == 0)
						combinedCodedEdgeEmission[emissionSymbol] = 0;

					combinedCodedEdgeEmission[emissionSymbol] += number;
				}

				combinedCodedEdgeEmission.erase(gap_symbol);
				combinedCodedEdgeEmission.erase(star_symbol);

				/*
				map<string, int> combinedEdgeEmission;
				for(set<string>::iterator kMerIt = utilizeKMersatLevel.begin(); kMerIt != utilizeKMersatLevel.end(); kMerIt++)
				{
					string kMerToUse = *kMerIt;
					combinedEdgeEmission[kMerToUse] = 0;
					int kMerSymbol = myGraph->CODE.doCode(e1->locus_id, kMerToUse);
					if(combinedCodedEdgeEmission.count(kMerSymbol) > 0)
					{
						combinedEdgeEmission[kMerToUse] = combinedCodedEdgeEmission[kMerSymbol];
					}
				}
				*/

				double log_emission_p = likelihood_all_kMers_0;
				assert(log_emission_p <= 0);

				for(map<int, int>::iterator kMerIt = combinedCodedEdgeEmission.begin(); kMerIt != combinedCodedEdgeEmission.end(); kMerIt++)
				{
					if(utilizeCodedKMersatLevel.count(kMerIt->first) == 0)
					{
						continue;
					}

					string kMer = myGraph->CODE.deCode(e1->locus_id, kMerIt->first);

					int underlyingCopyCount = kMerIt->second;
					assert(underlyingCopyCount != 0);

					assert(logL_count0.count(kMer) > 0);
					assert(logL_count0[kMer] <= 0);

					log_emission_p -= logL_count0[kMer];
					assert(log_emission_p <= 0);

					double rate = (double)underlyingCopyCount * coverage;
					assert(globalEmission.count(kMer) > 0);

					double thiskMer_P;
					if((_cache_poisson_PDF.count(underlyingCopyCount) > 0) && (_cache_poisson_PDF.at(underlyingCopyCount).count(globalEmission[kMer]) > 0))
					{
						thiskMer_P = _cache_poisson_PDF[underlyingCopyCount][globalEmission[kMer]];
					}
					else
					{
						poisson_up poisson(rate);
						thiskMer_P = pdf(poisson, globalEmission[kMer]);
						_cache_poisson_PDF[underlyingCopyCount][globalEmission[kMer]] = thiskMer_P;
					}

					if(thiskMer_P == 0)
					{
						thiskMer_P = 1e-200;
					}

					assert(thiskMer_P > 0);
					assert(thiskMer_P <= 1);

					log_emission_p += log(thiskMer_P);
					assert(log_emission_p <= 0);
				}

				// std::stringstream s_stdout;


				if(level == 0)
				{
					fw.at(level).at(state) = log(diploid_initialProbabilities.at(state));
					assert(fw.at(level).at(state) <= 0);

					fw_Viterbi.at(level).at(state) = log(diploid_initialProbabilities.at(state));
					fw_Viterbi_backtrack_int.at(level).at(state) = -1;
				}
				else
				{
					double max_jump_P = -1;
					int max_jump_state = -1;
					double combined_jump_P = 0;

					assert(diploidStateTransitions_Reverse.at(level).at(state).size() > 0);

					// s_stdout << "Compute combined_jump_P for state " << state << " at level " << level << "\n";

					for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(level).at(state).size(); sI2++)
					{
						stateAlternative jumpState = diploidStateTransitions_Reverse.at(level).at(state).at(sI2);
						int from = jumpState.id;
						double jumpState_P = jumpState.p;
						jumpState_P  = 1.0; // if we change this, we need to change it in the backwards sampling procedure as well!

						double viterbiJumpP = fw_Viterbi.at(level-1).at(from) + log(jumpState_P);

						if((max_jump_state == -1) || (viterbiJumpP > max_jump_P))
						{
							max_jump_state = from;
							max_jump_P = viterbiJumpP;
						}

						if(!(fw.at(level-1).at(from) >= 0))
						{
							// cout << "Error: !(fw.at(level-1).at(from) >= 0). Level: " << level << " from: " << from << " value: " << fw_Viterbi.at(level-1).at(from) << "\n" << flush;
						}

						assert(fw.at(level-1).at(from) >= 0);
						assert( jumpState_P >= 0);

						// s_stdout << "\t\t from " << from << ": " << fw.at(level-1).at(from) << ", jumpState_P " << jumpState_P << "\n";

						combined_jump_P += (fw.at(level-1).at(from) * jumpState_P);
						assert(combined_jump_P >= 0);
					}

					assert(max_jump_state != -1);

					assert(combined_jump_P >= 0);

					// s_stdout << "\tcombined_jump_P: " << combined_jump_P << "\n";

					if(combined_jump_P == 0)
					{
						fw.at(level).at(state) = -1 * 1e100;
					}
					else
					{
						fw.at(level).at(state) = log(combined_jump_P);
					}

					// s_stdout << "\tfw.at(level).at(state) : " << fw.at(level).at(state)  << "\n\n";


					//assert(fw.at(level).at(state) >= 0);

					fw_Viterbi.at(level).at(state) = max_jump_P;
					fw_Viterbi_backtrack_int.at(level).at(state) = max_jump_state;

				}

				assert(fw_Viterbi.at(level).at(state) <= 0);
				assert(log_emission_p <= 0);

				fw_Viterbi.at(level).at(state) += log_emission_p;
				fw.at(level).at(state) += log_emission_p;

				// assert(fw.at(level).at(state) <= 0);

				/*
					#pragma omp critical
					{
						subst_cout << "log_emission_p: " << log_emission_p << "\n";
						subst_cout << "fw.at(level).at(state): " << log_emission_p << "\n";
						cout << s_stdout.str() << flush;
					}
				*/


				if((state == firstPair) || (fw.at(level).at(state) > local_max_forward))
				{
					local_max_forward = fw.at(level).at(state);
					have_local_max_forward_value = true;
				}
			}

			//assert(local_max_forward <= 0);


			#pragma omp critical
			{
				if(have_local_max_forward_value)
				{
					if((set_max_forward_first_value == false) || (local_max_forward > max_forward))
					{
						max_forward = local_max_forward;
						set_max_forward_first_value = true;
					}
				}
			}
		}

		// cout << "\n\nmax_forward: " << max_forward << "\n";


		/*
		if(!(max_forward <= 0))
		{
			cout << "Level: " << level << "\n";
			cout << "max_forward: " << max_forward << "\n";
			cout << "states: " << states << "\n" << flush;
		}
		*/

		// assert(max_forward <= 0);

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

			// std::stringstream subst_cout;

			for(int state = firstPair; state <= lastPair; state++)
			{
				// subst_cout << "Re-compute forward for state " << state << " at level " << level << "\n";

				double fw_before = fw.at(level).at(state);

				/*
				subst_cout << "\tfw_before " << fw_before <<  "\n";
				subst_cout << "\tmax_forward " << max_forward <<  "\n";
				subst_cout << "\texp(fw.at(level).at(state) - max_forward): " <<  exp(fw.at(level).at(state) - max_forward) << "\n";
				*/

				fw.at(level).at(state) = exp(fw.at(level).at(state) - max_forward);

				// subst_cout << "\fw.at(level).at(state) " << fw.at(level).at(state) <<  "\n";

				if(!(fw.at(level).at(state) >= 0))
				{
					/*
					#pragma omp critical
					{
					cout << "!(fw.at(level).at(state) >= 0)\n";
					cout << "level: " << level << "\n";
					cout << "state: " << state << "\n";
					cout << "fw_before: " << fw_before << "\n";
					cout << "max_forward: " << max_forward << "\n";
					cout << "fw.at(level).at(state): " << fw.at(level).at(state) << "\n";
					cout << "\n" << flush;
					}
					*/
				}
				assert(fw.at(level).at(state) >= 0);
			}

			#pragma omp critical
			{
				//cout << subst_cout.str() << "\n" << flush;
			}
		}

		fw_underflow_factor += max_forward;
	}

	double max_viterbi = -1;
	int last_states = diploidStatesByLevel.at(levels-1);
	for(int state = 0; state < last_states; state++)
	{
		if((max_viterbi == -1) || (fw_Viterbi.at(levels-1).at(state) > max_viterbi))
		{
			max_viterbi = fw_Viterbi.at(levels-1).at(state);
		}
	}

	cout << "Best path optimality measure (log): " << max_viterbi << "\n";

	return max_viterbi;
}

void AlphaHMM::fillStateTransitions(int startLevelFill, int stopLevelFill, int previousStepFillStop, int previousStepFillStart)
{
	int levels = diploidStatesByLevel.size();

	if(previousStepFillStop >= 0)
	{
		if(previousStepFillStart < 0)
		{
			previousStepFillStart = 0;
		}

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


diploidEdgePointerPath AlphaHMM::sampleFromPosterior(map<string, long long>& observedEmissions)
{
	diploidEdgePointerPath forReturn;
	int levels = diploidStatesByLevel.size();
	int currentState = -1;
	int currentLevel = levels - 1;

	int currentFillStart = -1;
	int currentFillStop = -1;

	while(currentLevel >= 0)
	{
		int states = diploidStatesByLevel.at(currentLevel);

		if((currentFillStart == -1) || (currentLevel == (currentFillStart + 1)))
		{
			int newFillStart = currentLevel - 1000;
			int newFillStop = currentLevel;
			if(newFillStart < 0)
			{
				newFillStart = 0;
			}
			if(newFillStop < 0)
			{
				newFillStop = 0;
			}
			if(newFillStop >= levels)
			{
				newFillStop = levels - 1;
			}
			fillStateTransitions(newFillStart, newFillStop, currentFillStart, currentFillStop);
			currentFillStart = newFillStart;
			currentFillStop = newFillStop;
		}

		if(currentLevel == (levels - 1))
		{
			currentState = Utilities::chooseFromVector(fw.at(currentLevel));
		}
		else
		{
			int previousState = currentState;
			vector<double> alternatives_P;
			vector<int> alternatives_stateI;
			for(unsigned int sI2 = 0; sI2 < diploidStateTransitions_Reverse.at(currentLevel+1).at(previousState).size(); sI2++)
			{
				stateAlternative jumpState = diploidStateTransitions_Reverse.at(currentLevel+1).at(previousState).at(sI2);
				int from = jumpState.id;
				double jumpState_P = jumpState.p;
				jumpState_P  = 1.0; // if we change this, we need to change it in the backwards sampling procedure as well!

				alternatives_P.push_back(jumpState_P * fw.at(currentLevel).at(from));
				alternatives_stateI.push_back(from);
			}

			int currentState_vectorI = Utilities::chooseFromVector(alternatives_P);
			currentState = alternatives_stateI.at(currentState_vectorI);
		}

		int haploidStates = haploidStatesByLevel.at(currentLevel).size();
		int s1 = currentState % haploidStates;
		int s2 = currentState / haploidStates;
		Edge* e1 = haploidStatesByLevel.at(currentLevel).at(s1).e;
		Edge* e2 = haploidStatesByLevel.at(currentLevel).at(s2).e;

		forReturn.h1.push_back(e1);
		forReturn.h2.push_back(e2);

		currentLevel--;
	}

	assert(forReturn.h1.size() == levels);
	assert(forReturn.h1.size() == forReturn.h2.size());

	reverse(forReturn.h1.begin(), forReturn.h1.end());
	reverse(forReturn.h2.begin(), forReturn.h2.end());

	return forReturn;
}

double AlphaHMM::estimateCoverage(map<string, long long>& globalEmission)
{
	assert((globalEmission.count("MeanReadLen") > 0) && (globalEmission.count("TotalSeq") > 0) && (globalEmission.count("TotalKMerCoverage") > 0));

	double totalSeq = globalEmission["TotalSeq"];
	double correction_factor = ((double)globalEmission["MeanReadLen"]-(double)myGraph->kMerSize+1)/(double)globalEmission["MeanReadLen"];
	double total_kmers = globalEmission["TotalKMerCoverage"];
	double coverage = total_kmers / 3200000000.0;
	coverage = totalSeq / (3200000000.0 * 2.0) * correction_factor;
	return coverage;
}

