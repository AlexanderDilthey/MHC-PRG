/*
 * simulationSuite.cpp
 *
 *  Created on: 17 Apr 2012
 *      Author: dilthey
 */

#include "simulationSuite.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include "assert.h"
#include "NextGen.h"
#include <stdlib.h>
#include "../Data/HaplotypePanel.h"
#include "../Data/GenotypePanel.h"
#include "../Graph/Graph.h"
#include "../Graph/Node.h"
#include "../Graph/AlphaHMM.h"
#include "../MHC-PRG.h"
#include "../Utilities.h"
#include "../Graph/LargeGraph.h"
#include "../Graph/MultiGraph.h"
#include <iomanip>

void describeNucleotideGraph(string graph_file)
{
	Graph* g = new Graph();
	g->readFromFile(graph_file);
	int levels = g->NodesPerLevel.size();
	std::cout << "Graph " << graph_file << "\n";
	std::cout << "Levels: " << levels << "\n";
	int nodes_sum = 0;
	int nodes_max = 0;
	for(int l = 0; l < levels; l++)
	{
		int nodes = g->NodesPerLevel.at(l).size();
		nodes_sum += nodes;
		if(nodes > nodes_max)
		{
			nodes_max = nodes;
		}
	}
	double avg_nodes = (double)nodes_sum/(double)levels;
	
	std::cout << "Average number nodes / level: " << avg_nodes << "\n";
	std::cout << "Max number of nodes / level: " << nodes_max << "\n";
}
   
void describeGraph(string graph_file, string temp_dir, string temp_label)
{
	LargeGraph kMerG;
	kMerG.readFromFile(graph_file);
	MultiGraph* multiG = multiBeautifyForAlpha2(&kMerG, "");

	vector< vector<string> > level_categories = kMerG.categorizeEdgeLevels();
	set<string> possible_categories;
	for(int i = 0; i < (int)level_categories.size(); i++)
	{
		for(int j = 0; j < (int)level_categories.at(i).size(); j++)
		{
			possible_categories.insert(level_categories.at(i).at(j));
		}
	}

	multiG->kMerDiagnostics();

	vector<string> ordered_categories(possible_categories.begin(), possible_categories.end());

	vector<levelInfo> levelInformation = kMerG.getLevelInfo();

	string fn_globalstats = temp_dir + "/graph_globalstats_" + temp_label + ".txt";
	ofstream output_stats;
	output_stats.open (fn_globalstats.c_str(), ios::out | ios::trunc);
	if (output_stats.is_open())
	{
		output_stats << "Level\tNodes\tEdges\tSymbols\t" << Utilities::join(ordered_categories, "\t") << "\n";
		for(int i = 0; i < (int)level_categories.size(); i++)
		{
			set<string> categories(level_categories.at(i).begin(), level_categories.at(i).end());
			vector<string> fields;
			fields.push_back(Utilities::ItoStr(i));
			fields.push_back(Utilities::ItoStr(levelInformation.at(i).nodes));
			fields.push_back(Utilities::ItoStr(levelInformation.at(i).edges));
			fields.push_back(Utilities::ItoStr(levelInformation.at(i).symbols));
			for(int k = 0; k < (int)ordered_categories.size(); k++)
			{
				string cat = ordered_categories.at(k);
				if(categories.count(cat) > 0)
				{
					fields.push_back(Utilities::ItoStr(1));
				}
				else
				{
					fields.push_back(Utilities::ItoStr(0));
				}
			}
			output_stats << Utilities::join(fields, "\t") << "\n";
		}
	}
	output_stats.close();
}

void simulationSuite(string graph_file, string temp_dir, string temp_label, int genotypingMode)
{
	assert(agreementStats("A", "C", "A", "C") == 2);
	assert(agreementStats("A", "C", "C", "A") == 2);
	assert(agreementStats("A", "C", "C", "T") == 1);
	assert(agreementStats("A", "C", "G", "A") == 1);
	assert(agreementStats("A", "C", "G", "C") == 1);
	assert(agreementStats("A", "C", "A", "H") == 1);
	assert(agreementStats("A", "A", "A", "H") == 1);
	assert(agreementStats("A", "C", "A", "A") == 1);
	assert(agreementStats("A", "B", "C", "D") == 0);

	assert(agreementStarStats("A", "C", "A", "C").at(0) == 2);
	assert(agreementStarStats("A", "C", "C", "A").at(0) == 2);
	assert(agreementStarStats("A", "C", "C", "T").at(0) == 2);
	assert(agreementStarStats("A", "C", "G", "A").at(0) == 2);
	assert(agreementStarStats("A", "C", "G", "C").at(0) == 2);
	assert(agreementStarStats("A", "C", "A", "H").at(0) == 2);
	assert(agreementStarStats("A", "A", "A", "H").at(0) == 2);
	assert(agreementStarStats("A", "C", "A", "A").at(0) == 2);
	assert(agreementStarStats("A", "B", "C", "D").at(0) == 2);

	assert(agreementStarStats("A", "C", "A", "C").at(1) == 2);
	assert(agreementStarStats("A", "C", "C", "A").at(1) == 2);
	assert(agreementStarStats("A", "C", "C", "T").at(1) == 1);
	assert(agreementStarStats("A", "C", "G", "A").at(1) == 1);
	assert(agreementStarStats("A", "C", "G", "C").at(1) == 1);
	assert(agreementStarStats("A", "C", "A", "H").at(1) == 1);
	assert(agreementStarStats("A", "A", "A", "H").at(1) == 1);
	assert(agreementStarStats("A", "C", "A", "A").at(1) == 1);
	assert(agreementStarStats("A", "B", "C", "D").at(1) == 0);

	assert(agreementStarStats("A", "B", "*", "D").at(0) == 1);
	assert(agreementStarStats("A", "B", "*", "D").at(1) == 0);
	assert(agreementStarStats("A", "B", "*", "A").at(0) == 1);
	assert(agreementStarStats("A", "B", "*", "A").at(1) == 1);
	assert(agreementStarStats("A", "B", "*", "B").at(0) == 1);
	assert(agreementStarStats("A", "B", "*", "B").at(1) == 1);
	assert(agreementStarStats("A", "B", "*", "*").at(0) == 0);
	assert(agreementStarStats("A", "B", "*", "*").at(1) == 0);
	assert(agreementStarStats("*", "B", "A", "B").at(0) == 1);
	assert(agreementStarStats("*", "B", "A", "B").at(1) == 1);

	assert(agreementStarStats("A", "B", "*", "*").at(0) == 0);
	assert(agreementStarStats("A", "B", "*", "*").at(1) == 0);

	assert(agreementStarStats("*", "B", "*", "*").at(0) == 0);
	assert(agreementStarStats("*", "B", "*", "*").at(1) == 0);

	assert(agreementStarStats("*", "*", "B", "*").at(0) == 0);
	assert(agreementStarStats("*", "*", "B", "*").at(1) == 0);

	assert(agreementStarStats("*", "B", "*", "C").at(0) == 1);
	assert(agreementStarStats("*", "B", "*", "C").at(1) == 0);

	assert(agreementStarStats("A", "*", "B", "*").at(0) == 1);
	assert(agreementStarStats("A", "*", "B", "*").at(1) == 0);

	assert(agreementStarStats("A", "*", "A", "*").at(0) == 1);
	assert(agreementStarStats("A", "*", "A", "*").at(1) == 1);


	assert(agreementStarStats("A", "*", "*", "A").at(0) == 1);
	assert(agreementStarStats("A", "*", "*", "A").at(1) == 1);


	LargeGraph kMerG;
	kMerG.readFromFile(graph_file);

	MultiGraph* multiG = multiBeautifyForAlpha2(&kMerG, "");

	//multiG->stats();

	int iterations = 20;
	map<string, map< string, int> > total_comparison_statistics;
	set<string> keys;

	for(int iter = 1; iter <= iterations; iter++)
	{
		cout << "Simulation iteration " << iter << "\n";

		diploidEdgePath simulatedPath = kMerG.simulateRandomPath();
		diploidEdgePointerPath simulatedPath_edgePointers;
		for(int eI = 0; eI < (int)simulatedPath.h1.size(); eI++)
		{
			assert(kMerG.idx2Edge.count(simulatedPath.h1.at(eI)) > 0);
			assert(kMerG.idx2Edge.count(simulatedPath.h2.at(eI)) > 0);
			simulatedPath_edgePointers.h1.push_back(kMerG.idx2Edge[simulatedPath.h1.at(eI)]);
			simulatedPath_edgePointers.h2.push_back(kMerG.idx2Edge[simulatedPath.h2.at(eI)]);
		}

		diploidNucleotidePath realNucleotidePath = kMerG.diploidPathToNucleotides(simulatedPath);

		readSimulationResults readSimulation = kMerG.simulateReadsForPath(simulatedPath);
		map<string, long long> simulatedReads = readSimulation.simulatedReads;

		if(genotypingMode == 8)
		{
			int modifiedkMers = 0;
			for(map<string, long long>::iterator kMerIt = simulatedReads.begin(); kMerIt != simulatedReads.end(); kMerIt++)
			{
				if(modifiedkMers > 2)
				{
					break;
				}
				modifiedkMers++;
				kMerIt->second = -1;
			}
		}

		set<int> emptySet;

		cout << "Original multiG minEdges " << multiG->minEdgeNumber() << "\n";
		cout << "Original multiG maxEdges " << multiG->maxEdgeNumber() << "\n";

		MultiGraph* simpleMultiG = simplifyAccordingToCoverage(multiG, simulatedReads, emptySet);

		cout << "Simplified multiG minEdges " << simpleMultiG->minEdgeNumber() << "\n";
		cout << "Simplified multiG maxEdges " << simpleMultiG->maxEdgeNumber() << "\n" << flush;

		assert(readSimulation.haplotypes.at(0).size() == readSimulation.haplotypeStartPositions.at(0).size());
		assert(readSimulation.haplotypes.at(1).size() == readSimulation.haplotypeStartPositions.at(1).size());


		string temp_fn_haplos = temp_dir + "/" + temp_label + "_iteration_haplotypes_" + Utilities::ItoStr(iter);
		ofstream output_haplos;
		output_haplos.open (temp_fn_haplos.c_str(), ios::out | ios::trunc);
		if (output_haplos.is_open())
		{
			output_haplos << readSimulation.haplotypes.at(0) << "\n";
			for(int haploPos = 0; haploPos < (int)readSimulation.haplotypes.at(0).size(); haploPos++)
			{
				if((haploPos % 10) == 0)
				{
					output_haplos << left << setw(10) << readSimulation.haplotypeStartPositions.at(0).at(haploPos);
				}
			}
			output_haplos << "\n";

			output_haplos << readSimulation.haplotypes.at(1) << "\n";
			for(int haploPos = 0; haploPos < (int)readSimulation.haplotypes.at(1).size(); haploPos++)
			{
				if((haploPos % 10) == 0)
				{
					output_haplos << left << setw(10) << readSimulation.haplotypeStartPositions.at(1).at(haploPos);
				}
			}

			output_haplos << "\n";

		}
		output_haplos.close();

		string temp_fn_haplos_untrimmed = temp_dir + "/" + temp_label + "_iteration_haplotypes_untrimmed_" + Utilities::ItoStr(iter);
		ofstream output_haplos_untrimmed;
		output_haplos_untrimmed.open (temp_fn_haplos_untrimmed.c_str(), ios::out | ios::trunc);
		if (output_haplos_untrimmed.is_open())
		{
			output_haplos_untrimmed << readSimulation.haplotypes_untrimmed.at(0) << "\n";
			output_haplos_untrimmed << readSimulation.haplotypes_untrimmed.at(1) << "\n";
		}
		output_haplos_untrimmed.close();

		string temp_fn = temp_dir + "/" + temp_label + "_iteration_" + Utilities::ItoStr(iter);
		ofstream output;
		output.open (temp_fn.c_str(), ios::out | ios::trunc);
		if (output.is_open())
		{
			for(map<string, long long>::iterator mapIt = simulatedReads.begin(); mapIt != simulatedReads.end(); mapIt++)
			{
				output << mapIt->first << " " << mapIt->second << "\n";
			}
		}
		output.close();

		AlphaHMM aHMM(simpleMultiG, genotypingMode);
		aHMM.fillForwardBackwardTable(simulatedReads);

		//vector<diploidEdgePath> bestPaths_Edges = aHMM.retrieveAllBestPathsLargeGraph(simulatedReads);
		//vector<diploidNucleotidePath> bestPaths = aHMM.retrieveAllBestNucleotidePaths(simulatedReads);


		diploidEdgePointerPath bestPath_Edges_mG = aHMM.sampleFromPosterior(simulatedReads);
		diploidEdgePointerPath bestPath_Edges = simpleMultiG->MultiGraphEdgesToLargeGraphEdges(bestPath_Edges_mG);
		vector<diploidEdgePath> bestPaths_Edges;
		bestPaths_Edges.push_back(kMerG.edgePointerPathToEdgeIndexPath(bestPath_Edges));

		vector<diploidNucleotidePath> bestPaths;
		diploidNucleotidePath bestPath;
		assert(bestPaths_Edges.at(0).h1.size() == bestPaths_Edges.at(0).h2.size());
		bestPath.h1 = kMerG.haploidPathToNucleotides(bestPaths_Edges.at(0).h1);
		bestPath.h2 = kMerG.haploidPathToNucleotides(bestPaths_Edges.at(0).h2);
		bestPaths.push_back(bestPath);


		cout << "\t " << bestPaths.size() << " best paths\n";

		diploidEdgePath firstEstimatedEdgePath = bestPaths_Edges.at(0);
		diploidNucleotidePath firstEstimatedPath = bestPaths.at(0);
		assert(firstEstimatedPath.h1.size() == realNucleotidePath.h1.size());
		assert(firstEstimatedPath.h2.size() == realNucleotidePath.h2.size());

		diploidEdgePath firstEstimatedPath_Edges = bestPaths_Edges.at(0);
		assert(firstEstimatedPath_Edges.h1.size() == realNucleotidePath.h1.size());
		assert(firstEstimatedPath_Edges.h2.size() == realNucleotidePath.h2.size());

		diploidEdgePointerPath firstEstimatedPath_edgePointers;
		for(int eI = 0; eI < (int)firstEstimatedPath_Edges.h1.size(); eI++)
		{
			assert(kMerG.idx2Edge.count(firstEstimatedPath_Edges.h1.at(eI)) > 0);
			assert(kMerG.idx2Edge.count(firstEstimatedPath_Edges.h2.at(eI)) > 0);
			firstEstimatedPath_edgePointers.h1.push_back(kMerG.idx2Edge[firstEstimatedPath_Edges.h1.at(eI)]);
			firstEstimatedPath_edgePointers.h2.push_back(kMerG.idx2Edge[firstEstimatedPath_Edges.h2.at(eI)]);
		}

		/*
		diploidEdgePointerPath firstEstimatedPath_multiEdges = simpleMultiG->LargeEdgePointerPathToMultiEdgePointerPath(firstEstimatedPath_edgePointers, true);
		diploidEdgePointerPath simulatedPath_multiEdges = simpleMultiG->LargeEdgePointerPathToMultiEdgePointerPath(simulatedPath_edgePointers, true);
		assert(firstEstimatedPath_multiEdges.h1.size() == realNucleotidePath.h1.size());
		assert(firstEstimatedPath_multiEdges.h2.size() == realNucleotidePath.h2.size());
		assert(simulatedPath_multiEdges.h1.size() == realNucleotidePath.h1.size());
		assert(simulatedPath_multiEdges.h2.size() == realNucleotidePath.h2.size());

		diploidEdgePointerPath simulatedPath_multiEdges_short = multiG->LargeEdgePointerPathToMultiEdgePointerPath(simulatedPath_edgePointers, false);
		assert(simulatedPath_multiEdges_short.h1.size() == simulatedPath_multiEdges_short.h2.size());
		*/

		map<string, map< string, int> > comparison_statistics;

		vector< vector<string> > level_categories = kMerG.categorizeEdgeLevels();
		assert(firstEstimatedPath.h1.size() == level_categories.size());

		bool inPrintingCoverage = false;
		vector<bool> printPosCoverage;
		vector<int> agreementPerPos;

		vector< set<string> > indelCategories;
		vector<int> indelLengths;

		int thisIndelStart = -1;
		bool lastPositionIndel = false;
		indelLengths.resize(simulatedPath.h1.size());
 		for(int pos = 0; pos < (int)simulatedPath.h1.size(); pos++)
		{
 			indelLengths.at(pos) = 0;

 			bool inIndel = false;

			string real_1 = realNucleotidePath.h1.at(pos);
			string real_2 = realNucleotidePath.h2.at(pos);

			if((real_1 == "_") || (real_2 == "_"))
			{
				inIndel = true;
			}

			if(inIndel && (lastPositionIndel == false))
			{
				thisIndelStart = pos;
			}

			if((! inIndel) && lastPositionIndel)
			{
				// INDEL ends here
				int thisIndelStop = pos - 1;
				int indelLength = thisIndelStop - thisIndelStart + 1;
				for(int k = thisIndelStop; k >= thisIndelStart; k--)
				{
					indelLengths.at(k) = indelLength;
				}
			}

			if((inIndel) && (pos == (simulatedPath.h1.size() - 1)))
			{
				int thisIndelStop = pos;
				int indelLength = thisIndelStop - thisIndelStart + 1;
				for(int k = thisIndelStop; k >= thisIndelStart; k--)
				{
					indelLengths.at(k) = indelLength;
				}
			}

			lastPositionIndel = inIndel;
		}

 		for(int pos = 0; pos < (int)simulatedPath.h1.size(); pos++)
		{
			set<string> thisPos_indelCategories;

			if(indelLengths.at(pos) > 0)
			{
				string real_1 = realNucleotidePath.h1.at(pos);
				string real_2 = realNucleotidePath.h2.at(pos);

				int indelCount = 0;

				if(real_1 == "_")
				{
					indelCount++;
				}
				if(real_2 == "_")
				{
					indelCount++;
				}

				assert(indelCount > 0);

				thisPos_indelCategories.insert("realIndel_"+Utilities::ItoStr(indelCount));
				string strat_Index;
				if(indelLengths.at(pos) < 10)
				{
					strat_Index = Utilities::ItoStr(indelLengths.at(pos));
				}
				else if(indelLengths.at(pos) < 50)
				{
					strat_Index =  Utilities::ItoStr((indelLengths.at(pos) / 10) * 10)+"-"+Utilities::ItoStr(((indelLengths.at(pos) / 10)+1) * 10 - 1);
				}
				else if(indelLengths.at(pos) < 200)
				{
					strat_Index =  Utilities::ItoStr((indelLengths.at(pos) / 50) * 50)+"-"+Utilities::ItoStr(((indelLengths.at(pos) / 50)+1) * 50 - 1);
				}
				else
				{
					strat_Index =  Utilities::ItoStr((indelLengths.at(pos) / 1000) * 1000)+"-"+Utilities::ItoStr(((indelLengths.at(pos) / 1000)+1) * 1000 - 1);
				}

				thisPos_indelCategories.insert("realIndel_L"+strat_Index);
				thisPos_indelCategories.insert("realIndel_"+Utilities::ItoStr(indelCount)+"_L"+strat_Index);
			}

			indelCategories.push_back(thisPos_indelCategories);
		}

		for(int pos = 0; pos < (int)simulatedPath.h1.size(); pos++)
		{
			string estimated_1 = firstEstimatedPath.h1.at(pos);
			string estimated_2 = firstEstimatedPath.h2.at(pos);

			string real_1 = realNucleotidePath.h1.at(pos);
			string real_2 = realNucleotidePath.h2.at(pos);

			int agreement = agreementStats(estimated_1, estimated_2, real_1, real_2);
			agreementPerPos.push_back(agreement);

			bool printThisPositionCoverage = false;
			if(agreement < 2)
			{
				printThisPositionCoverage = true;
				if(inPrintingCoverage == false)
				{
					for(int bI = 1; bI < 60; bI++)
					{
						int bPos = pos - bI;
						if(bPos >= 0)
						{
							printPosCoverage.at(bPos) = true;
						}
					}
				}
				inPrintingCoverage = true;
			}
			else
			{
				inPrintingCoverage = false;
			}
			assert(kMerG.idx2Edge.count(firstEstimatedPath_Edges.h1.at(pos)) > 0);
			assert(kMerG.idx2Edge.count(firstEstimatedPath_Edges.h2.at(pos)) > 0);

			if(agreement < 2)
			{
				//cout << "AGREEMENT " << agreement << " position " << pos << ": " << estimated_1 << "/" << estimated_2 << " vs. real path " << real_1 << "/" << real_2 << "\n";
				//cout << "\t alternative edges: " << estimated_edge_1->To->Outgoing_Edges.size() << "/" << estimated_edge_2->To->Outgoing_Edges.size() << "\n";
			}

			vector<string> categories;

			vector<int> agreementNoStars = agreementStarStats(estimated_1, estimated_2, real_1, real_2);

			categories.push_back("ALL");
			categories.insert(categories.end(), level_categories.at(pos).begin(), level_categories.at(pos).end());
			categories.insert(categories.end(), indelCategories.at(pos).begin(), indelCategories.at(pos).end());


			if(real_1 != real_2)
			{
				vector<string> categories_het = categories;
				for(int cI = 0; cI < (int)categories_het.size(); cI++)
				{
					categories_het.at(cI) = categories_het.at(cI) + "_HET";
				}
				categories.insert(categories.end(), categories_het.begin(), categories_het.end());
			}

			for(int cI = 0; cI < (int)categories.size(); cI++)
			{
				string category = categories.at(cI);

				string key = "agreement_"+Utilities::ItoStr(agreement);
				string key2 = "nS_"+Utilities::ItoStr(agreementNoStars.at(1))+"_of_"+Utilities::ItoStr(agreementNoStars.at(0));

				if(comparison_statistics[category].count(key) == 0)
				{
					comparison_statistics[category][key] = 0;
				}
				comparison_statistics[category][key]++;
				if(total_comparison_statistics[category].count(key) == 0)
				{
					total_comparison_statistics[category][key] = 0;
				}
				total_comparison_statistics[category][key]++;

				keys.insert(key);


				if(comparison_statistics[category].count(key2) == 0)
				{
					comparison_statistics[category][key2] = 0;
				}
				comparison_statistics[category][key2]++;
				if(total_comparison_statistics[category].count(key2) == 0)
				{
					total_comparison_statistics[category][key2] = 0;
				}
				total_comparison_statistics[category][key2]++;

				keys.insert(key2);
			}

			printPosCoverage.push_back(printThisPositionCoverage);
		}
		assert(printPosCoverage.size() == simulatedPath.h1.size());

		/*
		Edge* lastSimulatedMultiEdgeH1 = 0;
		for(int pos = 0; pos < (int)simulatedPath.h1.size(); pos++)
		{
			if(aHMM.getGenotypingMode() < 5)
			{
				if(! (printPosCoverage.at(pos) == true))
				{
					if((pos > 0) && (printPosCoverage.at(pos-1) == true))
					{
						cout << "\n\n===================================================================================\n\n";
					}

					lastSimulatedMultiEdgeH1 = 0;

					continue;
				}

				if(lastSimulatedMultiEdgeH1 != simulatedPath_multiEdges.h1.at(pos))
				{
					Edge* multiE1_simulated = simulatedPath_multiEdges.h1.at(pos);
					Edge* multiE2_simulated = simulatedPath_multiEdges.h2.at(pos);

					Edge* multiE1_estimated = firstEstimatedPath_multiEdges.h1.at(pos);
					Edge* multiE2_estimated = firstEstimatedPath_multiEdges.h2.at(pos);

					int level = multiE1_estimated->From->level;

					assert(aHMM.Edge2Int.count(multiE1_estimated) > 0);
					int multiE1_simulated_haploidIndex = aHMM.Edge2Int[multiE1_simulated];
					assert(aHMM.Edge2Int.count(multiE2_estimated) > 0);
					int multiE2_simulated_haploidIndex = aHMM.Edge2Int[multiE2_simulated];
					assert(aHMM.Edge2Int.count(multiE1_estimated) > 0);
					int multiE1_estimated_haploidIndex = aHMM.Edge2Int[multiE1_estimated];
					assert(aHMM.Edge2Int.count(multiE2_estimated) > 0);
					int multiE2_estimated_haploidIndex = aHMM.Edge2Int[multiE2_estimated];

					int haploidStatesLevel = aHMM.haploidStatesByLevel.at(level).size();
					int diploidSimulatedIndex = multiE2_simulated_haploidIndex*haploidStatesLevel+multiE1_simulated_haploidIndex;
					int multiE1_simulated_haploidIndex_expect = diploidSimulatedIndex % haploidStatesLevel;
					int multiE2_simulated_haploidIndex_expect = diploidSimulatedIndex / haploidStatesLevel;
					assert(multiE1_simulated_haploidIndex_expect == multiE1_simulated_haploidIndex);
					assert(multiE2_simulated_haploidIndex_expect == multiE2_simulated_haploidIndex);

					int diploidEstimatedIndex = multiE2_estimated_haploidIndex*haploidStatesLevel+multiE1_estimated_haploidIndex;
					int multiE1_estimated_haploidIndex_expect = diploidEstimatedIndex % haploidStatesLevel;
					int multiE2_estimated_haploidIndex_expect = diploidEstimatedIndex / haploidStatesLevel;
					assert(multiE1_estimated_haploidIndex_expect == multiE1_estimated_haploidIndex);
					assert(multiE2_estimated_haploidIndex_expect == multiE2_estimated_haploidIndex);

					double optimalitySimulated = aHMM.fw_Viterbi.at(level).at(diploidSimulatedIndex);
					double optimalityEstimated = aHMM.fw_Viterbi.at(level).at(diploidEstimatedIndex);
					double emissionSimulated = aHMM.fw_Emission.at(level).at(diploidSimulatedIndex);
					double emissionEstimated = aHMM.fw_Emission.at(level).at(diploidEstimatedIndex);
					double expectedMissingSimulated_e1 = aHMM.fw_expectedMissing_e1.at(level).at(diploidSimulatedIndex);
					double expectedMissingEstimated_e1 = aHMM.fw_expectedMissing_e1.at(level).at(diploidEstimatedIndex);
					double expectedMissingSimulated_e2 = aHMM.fw_expectedMissing_e2.at(level).at(diploidSimulatedIndex);
					double expectedMissingEstimated_e2 = aHMM.fw_expectedMissing_e2.at(level).at(diploidEstimatedIndex);

					double estimatedEdgeCoverage_e1 = aHMM.fw_Coverage_e1.at(level).at(diploidEstimatedIndex);
					double estimatedEdgeCoverage_e2 = aHMM.fw_Coverage_e2.at(level).at(diploidEstimatedIndex);
					double simulatedEdgeCoverage_e1 = aHMM.fw_Coverage_e1.at(level).at(diploidSimulatedIndex);
					double simulatedEdgeCoverage_e2 = aHMM.fw_Coverage_e2.at(level).at(diploidSimulatedIndex);

					cout << "alphaHMM optimality stats level " << level << "\n";
					cout << "\tsimulated state " << diploidSimulatedIndex << " optimality: " << optimalitySimulated << "\n";
					cout << "\t\t (classical) emission " << emissionSimulated << "  (real under gt4) emission: " << aHMM.fw_adaptedEmission.at(level).at(diploidSimulatedIndex) << "\n";
					cout << "\t\t e1 emission: " << aHMM.fw_Emission_e1.at(level).at(diploidSimulatedIndex) << " adapted emission: " << aHMM.fw_adaptedEmission_e1.at(level).at(diploidSimulatedIndex) << " coverage: " << simulatedEdgeCoverage_e1 << ", expected missing e1: " << expectedMissingSimulated_e1 << ", total length: " << aHMM.fw_e1_length.at(level).at(diploidSimulatedIndex) << ", of which valid: " << aHMM.fw_e1_length_valid.at(level).at(diploidSimulatedIndex) << "\n";
					cout << "\t\t e2 emission: " << aHMM.fw_Emission_e2.at(level).at(diploidSimulatedIndex) << " adapted emission: " << aHMM.fw_adaptedEmission_e2.at(level).at(diploidSimulatedIndex) << " coverage: " << simulatedEdgeCoverage_e2 << ", expected missing e2: " << expectedMissingSimulated_e2 << ", total length: " << aHMM.fw_e2_length.at(level).at(diploidSimulatedIndex) << ", of which valid: " << aHMM.fw_e2_length_valid.at(level).at(diploidSimulatedIndex) << "\n";

					cout << "\testimated state " << diploidEstimatedIndex << " optimality: " << optimalityEstimated << "\n";
					cout << "\t\t (classical) emission " << emissionEstimated << "  (real under gt4) emission: " << aHMM.fw_adaptedEmission.at(level).at(diploidEstimatedIndex) << "\n";
					cout << "\t\t e1 emission: " << aHMM.fw_Emission_e1.at(level).at(diploidEstimatedIndex) << " adapted emission: " << aHMM.fw_adaptedEmission_e1.at(level).at(diploidEstimatedIndex) << " coverage: " << estimatedEdgeCoverage_e1 << ", expected missing e1: " << expectedMissingEstimated_e1 << ", total length: " << aHMM.fw_e1_length.at(level).at(diploidEstimatedIndex) << ", of which valid: " << aHMM.fw_e1_length_valid.at(level).at(diploidEstimatedIndex) << "\n";
					cout << "\t\t e2 emission: " << aHMM.fw_Emission_e2.at(level).at(diploidEstimatedIndex) << " adapted emission: " << aHMM.fw_adaptedEmission_e2.at(level).at(diploidEstimatedIndex) << " coverage: " << estimatedEdgeCoverage_e2 << ", expected missing e2: " << expectedMissingEstimated_e2 << ", total length: " << aHMM.fw_e2_length.at(level).at(diploidEstimatedIndex) << ", of which valid: " << aHMM.fw_e2_length_valid.at(level).at(diploidEstimatedIndex) << "\n";

					lastSimulatedMultiEdgeH1 = simulatedPath_multiEdges.h1.at(pos);
				}
				assert(kMerG.idx2Edge.count(firstEstimatedEdgePath.h1.at(pos)) > 0);
				assert(kMerG.idx2Edge.count(firstEstimatedEdgePath.h2.at(pos)) > 0);
				Edge* estimated_1 = kMerG.idx2Edge[firstEstimatedEdgePath.h1.at(pos)];
				Edge* estimated_2 = kMerG.idx2Edge[firstEstimatedEdgePath.h2.at(pos)];

				assert(kMerG.idx2Edge.count(simulatedPath.h1.at(pos)) > 0);
				assert(kMerG.idx2Edge.count(simulatedPath.h2.at(pos)) > 0);
				Edge* real_1 = kMerG.idx2Edge[simulatedPath.h1.at(pos)];
				Edge* real_2 = kMerG.idx2Edge[simulatedPath.h2.at(pos)];

				map<int, int> combinedEdgeEmission_estimated;
				if(combinedEdgeEmission_estimated.count(estimated_1->largeEmission) == 0)
				{
					combinedEdgeEmission_estimated[estimated_1->largeEmission] = 0;
				}
				if(combinedEdgeEmission_estimated.count(estimated_2->largeEmission) == 0)
				{
					combinedEdgeEmission_estimated[estimated_2->largeEmission] = 0;
				}
				combinedEdgeEmission_estimated[estimated_1->largeEmission]++;
				combinedEdgeEmission_estimated[estimated_2->largeEmission]++;

				map<int, int> combinedEdgeEmission_real;
				if(combinedEdgeEmission_real.count(real_1->largeEmission) == 0)
				{
					combinedEdgeEmission_real[real_1->largeEmission] = 0;
				}
				if(combinedEdgeEmission_real.count(real_2->largeEmission) == 0)
				{
					combinedEdgeEmission_real[real_2->largeEmission] = 0;
				}
				combinedEdgeEmission_real[real_1->largeEmission]++;
				combinedEdgeEmission_real[real_2->largeEmission]++;

				string n_estimated_1 = firstEstimatedPath.h1.at(pos);
				string n_estimated_2 = firstEstimatedPath.h2.at(pos);

				string n_real_1 = realNucleotidePath.h1.at(pos);
				string n_real_2 = realNucleotidePath.h2.at(pos);

				cout << "\t\tDetail-position " << pos << " agreement: " << agreementPerPos.at(pos) << "\n";
				cout << "\t\t\tReal underlying edge " << n_real_1 << "/" << n_real_2 << "\n";

				for(map<int, int>::iterator emissionIt = combinedEdgeEmission_real.begin(); emissionIt != combinedEdgeEmission_real.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;
					string kMer = kMerG.CODE.deCode(real_1->locus_id, emissionSymbol);
					int observedNumber = simulatedReads[kMer];
					string kMerLocalization = Utilities::join(Utilities::ItoStr(vector<int>(readSimulation.kMerOccurenceLocalization[kMer].begin(), readSimulation.kMerOccurenceLocalization[kMer].end())), ", ");

					cout << "\t\t\t\texpected: " << number << "   observed: " << observedNumber << "    in underlying haplotypes: " << readSimulation.kMerOccurence[kMer] << "(" << kMerLocalization << ")  [" << kMer << "]\n";
				}

				cout << "\t\t\tEstimated edge " << n_estimated_1 << "/" << n_estimated_2 << "\n";

				for(map<int, int>::iterator emissionIt = combinedEdgeEmission_estimated.begin(); emissionIt != combinedEdgeEmission_estimated.end(); emissionIt++)
				{
					int emissionSymbol = emissionIt->first;
					int number = emissionIt->second;
					string kMer = kMerG.CODE.deCode(real_1->locus_id, emissionSymbol);
					int observedNumber = simulatedReads[kMer];
					string kMerLocalization = Utilities::join(Utilities::ItoStr(vector<int>(readSimulation.kMerOccurenceLocalization[kMer].begin(), readSimulation.kMerOccurenceLocalization[kMer].end())), ", ");


					cout << "\t\t\t\texpected: " << number << "   observed: " << observedNumber << "    in underlying haplotypes: " << readSimulation.kMerOccurence[kMer] << "(" << kMerLocalization << ")  [" << kMer << "]\n";
				}
			}
		}

		if(aHMM.getGenotypingMode() < 5)
		{
			for(int pos = 1; pos < (int)simulatedPath_multiEdges_short.h1.size(); pos++)
			{
				Edge* multiE1_simulated = simulatedPath_multiEdges_short.h1.at(pos);
				Edge* multiE2_simulated = simulatedPath_multiEdges_short.h2.at(pos);

				Edge* multiE1_simulated_minus1 = simulatedPath_multiEdges_short.h1.at(pos-1);
				Edge* multiE2_simulated_minus1 = simulatedPath_multiEdges_short.h2.at(pos-1);

				int level = multiE1_simulated->From->level;
				int level_minus1 = level - 1;

				assert(aHMM.Edge2Int.count(multiE1_simulated) > 0);
				int multiE1_simulated_haploidIndex = aHMM.Edge2Int[multiE1_simulated];
				assert(aHMM.Edge2Int.count(multiE2_simulated) > 0);
				int multiE2_simulated_haploidIndex = aHMM.Edge2Int[multiE2_simulated];

				assert(aHMM.Edge2Int.count(multiE1_simulated_minus1) > 0);
				int multiE1_simulated_minus1_haploidIndex = aHMM.Edge2Int[multiE1_simulated_minus1];
				assert(aHMM.Edge2Int.count(multiE2_simulated_minus1) > 0);
				int multiE2_simulated_minus1_haploidIndex = aHMM.Edge2Int[multiE2_simulated_minus1];

				int haploidStatesLevel = aHMM.haploidStatesByLevel.at(level).size();
				int haploidStatesLevel_minus1 = aHMM.haploidStatesByLevel.at(level_minus1).size();

				int diploidSimulatedIndex = multiE2_simulated_haploidIndex*haploidStatesLevel+multiE1_simulated_haploidIndex;
				int multiE1_simulated_haploidIndex_expect = diploidSimulatedIndex % haploidStatesLevel;
				int multiE2_simulated_haploidIndex_expect = diploidSimulatedIndex / haploidStatesLevel;
				assert(multiE1_simulated_haploidIndex_expect == multiE1_simulated_haploidIndex);
				assert(multiE2_simulated_haploidIndex_expect == multiE2_simulated_haploidIndex);

				int diploidSimulatedIndex_minus1 = multiE2_simulated_minus1_haploidIndex*haploidStatesLevel_minus1+multiE1_simulated_minus1_haploidIndex;
				int multiE1_simulated_minus1_haploidIndex_expect = diploidSimulatedIndex_minus1 % haploidStatesLevel_minus1;
				int multiE2_simulated_minus1_haploidIndex_expect = diploidSimulatedIndex_minus1 / haploidStatesLevel_minus1;
				assert(multiE1_simulated_minus1_haploidIndex_expect == multiE1_simulated_minus1_haploidIndex);
				assert(multiE2_simulated_minus1_haploidIndex_expect == multiE2_simulated_minus1_haploidIndex);

				int previousPossibleStates = aHMM.diploidStateTransitions_Reverse.at(level).at(diploidSimulatedIndex).size();
				int bestAgreement = 0;
				for(int vI = 0; vI < (int)aHMM.fw_Viterbi_backtrack.at(level).at(diploidSimulatedIndex).size(); vI++)
				{
					int previousBestState = aHMM.fw_Viterbi_backtrack.at(level).at(diploidSimulatedIndex).at(vI);
					int previousBestState_e1 = previousBestState % haploidStatesLevel_minus1;
					int previousBestState_e2 = previousBestState / haploidStatesLevel_minus1;
					assert(previousBestState_e1 < (int)aHMM.haploidStatesByLevel.at(level_minus1).size());
					assert(previousBestState_e2 < (int)aHMM.haploidStatesByLevel.at(level_minus1).size());
					int edge_agreement = agreementStats(previousBestState_e1, previousBestState_e2, multiE1_simulated_minus1_haploidIndex, multiE2_simulated_minus1_haploidIndex);
					if(edge_agreement > bestAgreement)
					{
						bestAgreement = edge_agreement;
					}
				}

				vector<string> categories;
				categories.push_back("localEdge_ALL");
				if(previousPossibleStates > 1)
				{
					categories.push_back("localEdge_MultipleIncoming");
				}

				for(int cI = 0; cI < (int)categories.size(); cI++)
				{
					string category = categories.at(cI);

					string key = "agreement_"+Utilities::ItoStr(bestAgreement);
					if(comparison_statistics[category].count(key) == 0)
					{
						comparison_statistics[category][key] = 0;
					}
					comparison_statistics[category][key]++;
					if(total_comparison_statistics[category].count(key) == 0)
					{
						total_comparison_statistics[category][key] = 0;
					}
					total_comparison_statistics[category][key]++;

					keys.insert(key);
				}

			}

			for(map<string, map< string, int> >::iterator mapIt = comparison_statistics.begin(); mapIt != comparison_statistics.end(); mapIt++)
			{
				cout << "\t\t\t " << mapIt->first << "\n";
				int category_total = 0;
				for(set<string>::iterator keyIt = keys.begin(); keyIt != keys.end(); keyIt++)
				{
					if(mapIt->second.count(*keyIt) == 0)
					{
						mapIt->second[*keyIt] = 0;
					}
					category_total += mapIt->second[*keyIt];
				}
				for(set<string>::iterator keyIt = keys.begin(); keyIt != keys.end(); keyIt++)
				{
					cout << "\t\t\t\t " << *keyIt << " " <<  mapIt->second[*keyIt] << "/" << category_total << "\n";
				}
			}
		}
		*/

		delete simpleMultiG;
	}

	cout << "\nTotal summary\n";
	for(map<string, map< string, int> >::iterator mapIt = total_comparison_statistics.begin(); mapIt != total_comparison_statistics.end(); mapIt++)
	{
		cout << "\t " << mapIt->first << "\n";
		int category_total = 0;
		for(set<string>::iterator keyIt = keys.begin(); keyIt != keys.end(); keyIt++)
		{
			if(mapIt->second.count(*keyIt) == 0)
			{
				mapIt->second[*keyIt] = 0;
			}
			category_total += mapIt->second[*keyIt];
		}
		for(set<string>::iterator keyIt = keys.begin(); keyIt != keys.end(); keyIt++)
		{
			cout << "\t\t " << *keyIt << " " <<  mapIt->second[*keyIt] << "/" << category_total << "\n";
		}
	}

	delete multiG;
}

vector<int> agreementStarStats(string estimated_1, string estimated_2, string real_1, string real_2)
{

	vector<string> estimated;
	if(estimated_1 != "*")
	{
		estimated.push_back(estimated_1);
	}
	if(estimated_2 != "*")
	{
		estimated.push_back(estimated_2);
	}
	vector<string> real;
	if(real_1 != "*")
	{
		real.push_back(real_1);
	}
	if(real_2 != "*")
	{
		real.push_back(real_2);
	}

	int compared_alleles = 0;
	int agreement_alleles = 0;

	if((real.size() == 2) && (estimated.size() == 2))
	{
		compared_alleles = 2;
		if(((estimated_1 == real_1) && (estimated_2 == real_2)) || ((estimated_1 == real_2) && (estimated_2 == real_1)))
		{
			agreement_alleles = 2;
		}
		else if (((estimated_1 == real_1) && (estimated_2 != real_2)) || ((estimated_1 == real_2) && (estimated_2 != real_1)) ||
				((estimated_1 != real_1) && (estimated_2 == real_2)) || ((estimated_1 != real_2) && (estimated_2 == real_1))
				)
		{
			agreement_alleles = 1;

		}
		else
		{
			agreement_alleles = 0;
		}
	}
	else if((real.size() == 0) || (estimated.size() == 0))
	{
		agreement_alleles = 0;
		compared_alleles = 0;
	}
	else
	{
		compared_alleles = 1;
		if((real.size() == 1) && (estimated.size() == 1))
		{
			if(real.at(0) == estimated.at(0))
			{
				agreement_alleles = 1;
			}
			else
			{
				agreement_alleles = 0;
			}
		}
		else if(real.size() == 1)
		{
			assert(estimated.size() == 2);
			if((real.at(0) == estimated.at(0)) || (real.at(0) == estimated.at(1)))
			{
				agreement_alleles = 1;
			}
			else
			{
				agreement_alleles = 0;
			}
		}
		else if(estimated.size() == 1)
		{
			assert(real.size() == 2);
			if((real.at(0) == estimated.at(0)) || (real.at(1) == estimated.at(0)))
			{
				agreement_alleles = 1;
			}
			else
			{
				agreement_alleles = 0;
			}
		}
		else
		{
			assert(1 == 0);
		}
	}

	vector<int> agreement;
	assert(agreement_alleles <= compared_alleles);
	agreement.push_back(compared_alleles);
	agreement.push_back(agreement_alleles);
	return agreement;
}

int agreementStats(string estimated_1, string estimated_2, string real_1, string real_2)
{
	int agreement;

	if(((estimated_1 == real_1) && (estimated_2 == real_2)) || ((estimated_1 == real_2) && (estimated_2 == real_1)))
	{
		agreement = 2;
	}
	else if (((estimated_1 == real_1) && (estimated_2 != real_2)) || ((estimated_1 == real_2) && (estimated_2 != real_1)) ||
			((estimated_1 != real_1) && (estimated_2 == real_2)) || ((estimated_1 != real_2) && (estimated_2 == real_1))
			)
	{
		agreement = 1;

	}
	else
	{
		agreement = 0;
	}

	return agreement;
}


int agreementStats(int estimated_1, int estimated_2, int real_1, int real_2)
{
	int agreement;

	if(((estimated_1 == real_1) && (estimated_2 == real_2)) || ((estimated_1 == real_2) && (estimated_2 == real_1)))
	{
		agreement = 2;
	}
	else if (((estimated_1 == real_1) && (estimated_2 != real_2)) || ((estimated_1 == real_2) && (estimated_2 != real_1)) ||
			((estimated_1 != real_1) && (estimated_2 == real_2)) || ((estimated_1 != real_2) && (estimated_2 == real_1))
			)
	{
		agreement = 1;

	}
	else
	{
		agreement = 0;
	}

	return agreement;
}

