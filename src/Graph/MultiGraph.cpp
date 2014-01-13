/*
 * Graph.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "MultiGraph.h"
#include "../Utilities.h"

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef std::basic_string <unsigned char> ustring;

MultiGraph::MultiGraph() {

}

MultiGraph::~MultiGraph() {
	freeMemory();
}

int MultiGraph::maxEdgeNumber()
{
	int maxEdgeNum = 0;

	for(int level = 0; level < NodesPerLevel.size(); level++)
	{
		set<Node*> nodesAtLevel = NodesPerLevel.at(level);

		int edgeNum = 0;

		for(set<Node*>::iterator nIt = nodesAtLevel.begin(); nIt != nodesAtLevel.end(); nIt++)
		{
			Node* n = *nIt;


			for(set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				edgeNum++;
			}
		}

		if(edgeNum > maxEdgeNum)
		{
			maxEdgeNum = edgeNum;
		}
	}

	return maxEdgeNum;
}


int MultiGraph::minEdgeNumber()
{
	int minEdgeNum = 10000000;

	for(int level = 0; level < NodesPerLevel.size() - 1; level++)
	{
		set<Node*> nodesAtLevel = NodesPerLevel.at(level);

		int edgeNum = 0;

		for(set<Node*>::iterator nIt = nodesAtLevel.begin(); nIt != nodesAtLevel.end(); nIt++)
		{
			Node* n = *nIt;


			for(set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				edgeNum++;
			}
		}

		if(edgeNum < minEdgeNum)
		{
			minEdgeNum = edgeNum;
		}
	}

	return minEdgeNum;
}

kMerPositionInfo MultiGraph::getkMerPositions()
{
	kMerPositionInfo forReturn;

	for(int l = 0; l < ((int)NodesPerLevel.size()-1); l++)
	{
		string locusID = (*(*(NodesPerLevel.at(l).begin()))->Outgoing_Edges.begin())->locus_id;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;

			if(l != (NodesPerLevel.size()-1))
			{
				string local_locus = (*(n->Outgoing_Edges.begin()))->locus_id;
				assert(local_locus == locusID);

				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					Edge* e = *outgoingIt;
					assert(e->locus_id == local_locus);

					for(map<int, int>::iterator emissionIt = e->multiEmission.begin(); emissionIt != e->multiEmission.end(); emissionIt++)
					{
						int emissionSymbol = emissionIt->first;
						int number = emissionIt->second;
						string kMer = CODE.deCode(e->locus_id, emissionSymbol);
						if((kMer != "*") && (kMer != "_"))
						{
							forReturn.kMers.insert(kMer);
							forReturn.kMers2Level[kMer].insert(l);
							if(forReturn.kMerMultiplicity.count(kMer) == 0)
							{
								forReturn.kMerMultiplicity[kMer] = 0;
							}
							forReturn.kMerMultiplicity[kMer] += number;

							if(forReturn.kMerEdgeNum.count(kMer) == 0)
							{
								forReturn.kMerEdgeNum[kMer] = 0;
							}
							forReturn.kMerEdgeNum[kMer]++;

							if(forReturn.kMerMaxPerLevel.count(kMer) == 0)
							{
								forReturn.kMerMaxPerLevel[kMer] = 0;
							}
							if(number > forReturn.kMerMaxPerLevel[kMer])
							{
								forReturn.kMerMaxPerLevel[kMer] = number;
							}

							forReturn.Level2Kmers[l].insert(kMer);
						}
					}
				}
			}
		}

		if(forReturn.Level2Kmers.count(l) == 0)
		{
			set<string> empty;
			forReturn.Level2Kmers[l] = empty;
		}
	}

	return forReturn;
}

kMerUniquenessInfo MultiGraph::kMerUniqueness()
{
	kMerUniquenessInfo forReturn;

	kMerPositionInfo kMerPos = getkMerPositions();

	set<string> levelUniqueKMers;
	for(set<string>::iterator kMerIt = kMerPos.kMers.begin(); kMerIt != kMerPos.kMers.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		assert(kMerPos.kMers2Level.count(kMer) > 0);
		int levels = kMerPos.kMers2Level[kMer].size();

		assert(levels >= 1);

		if(levels == 1)
		{
			forReturn.levelUniqueKMers.insert(kMer);
		}
	}

	for(set<string>::iterator kMerIt = forReturn.levelUniqueKMers.begin(); kMerIt != forReturn.levelUniqueKMers.end(); kMerIt++)
	{
		forReturn.kMers2Level[*kMerIt] = kMerPos.kMers2Level[*kMerIt];
	}

	return forReturn;
}

kMerNonUniquenessInfo MultiGraph::kMerNonUniqueness()
{
	kMerNonUniquenessInfo forReturn;

	kMerPositionInfo kMerPos = getkMerPositions();

	set<string> levelNonUniqueKMers;
	for(set<string>::iterator kMerIt = kMerPos.kMers.begin(); kMerIt != kMerPos.kMers.end(); kMerIt++)
	{
		string kMer = *kMerIt;
		assert(kMerPos.kMers2Level.count(kMer) > 0);
		int levels = kMerPos.kMers2Level[kMer].size();

		assert(levels >= 1);

		if(levels != 1)
		{
			forReturn.levelNonUniqueKMers.insert(kMer);
		}
	}

	for(set<string>::iterator kMerIt = forReturn.levelNonUniqueKMers.begin(); kMerIt != forReturn.levelNonUniqueKMers.end(); kMerIt++)
	{
		forReturn.kMers2Level[*kMerIt] = kMerPos.kMers2Level[*kMerIt];
		for(set<int>::iterator levelIt = forReturn.kMers2Level[*kMerIt] .begin(); levelIt != forReturn.kMers2Level[*kMerIt].end(); levelIt++)
		{
			int level = *levelIt;
			forReturn.level2kMers[level].insert(*kMerIt);
		}
	}

	forReturn.kMerMultiplicity = kMerPos.kMerMultiplicity;

	return forReturn;
}


void MultiGraph::kMerDiagnostics()
{
	kMerPositionInfo kMerPos = getkMerPositions();

	set<string> kMers = kMerPos.kMers;
	map<string, set<int> > kMers2Level = kMerPos.kMers2Level;
	set<string> levelUniqueKMers = kMerUniqueness().levelUniqueKMers;

	cout << "GRAPH DESCRIPTION\n";
	cout << "Total # kMers: " << kMers.size() << "\n";
	cout << "\tof which " << levelUniqueKMers.size() << " (";
	printf("%.2f", (double)levelUniqueKMers.size()/(double)kMers.size() );
	cout << ") are unique\n\n";

	map<int, avg_struct> uniquenessByLength;

	for(int l = 0; l < ((int)NodesPerLevel.size()-1); l++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;

			if(l != (NodesPerLevel.size()-1))
			{
				string local_locus = (*(n->Outgoing_Edges.begin()))->locus_id;

				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					Edge* e = *outgoingIt;
					assert(e->locus_id == local_locus);
					int edge_length = 0;
					int edge_levelUnique = 0;

					for(map<int, int>::iterator emissionIt = e->multiEmission.begin(); emissionIt != e->multiEmission.end(); emissionIt++)
					{
						int emissionSymbol = emissionIt->first;
						int number = emissionIt->second;
						string kMer = CODE.deCode(e->locus_id, emissionSymbol);
						if((kMer != "_") && (kMer != "*"))
						{
							edge_length += number;
							if(levelUniqueKMers.count(kMer) > 0)
							{
								edge_levelUnique += number;
							}
						}
					}

					if(edge_length > 0)
					{
						int key = (int)log(edge_length);
						if(uniquenessByLength.count(key) == 0)
						{
							uniquenessByLength[key].count = 0;
							uniquenessByLength[key].sum = 0;
						}

						uniquenessByLength[key].count += edge_levelUnique;
						uniquenessByLength[key].sum += edge_length;
					}
				}
			}
		}
	}

	cout << "Average unique kMer proportion by edge length:\n";
	for(map<int, avg_struct>::iterator uIt = uniquenessByLength.begin(); uIt != uniquenessByLength.end(); uIt++)
	{
		int index = uIt->first;
		avg_struct data = uIt->second;
		assert(data.sum > 0);
		double avg = (double)data.count/(double)data.sum;
		cout << index << "\t" << "below " << exp(index+1) << ": " << avg << "\n";
	}

	cout << "\n\n";

}


void MultiGraph::unRegisterNode(Node* n)
{
	assert(Nodes.count(n) > 0);
	int l = n->level;
	Nodes.erase(n);
	NodesPerLevel.at(l).erase(n);
	delete(n);
}

void MultiGraph::unRegisterEdge(Edge* e)
{
	assert(Edges.count(e) > 0);
	Edges.erase(e);
	delete(e);
}

void MultiGraph::registerNode(Node* n, unsigned int level)
{
	assert(Nodes.count(n) == 0);
	assert(n->level == level);
	if((level+1) > NodesPerLevel.size())
	{
		NodesPerLevel.resize(level+1);
	}
	NodesPerLevel.at(level).insert(n);
	Nodes.insert(n);
	n->mG = this;
}


multiHaploLabelPair MultiGraph::edgePointerPathToLabels(diploidEdgePointerPath& multiEdgesPath)
{
	multiHaploLabelPair forReturn;

	for(unsigned int level = 0; level < multiEdgesPath.h1.size(); level++)
	{
		Edge* e1 = multiEdgesPath.h1.at(level);
		Edge* e2 = multiEdgesPath.h2.at(level);

		if((e1->label != "") || (e2->label != ""))
		{
			forReturn.h1.push_back(e1->label);
			forReturn.h2.push_back(e2->label);
		}
	}

	return forReturn;

}

diploidEdgePointerPath MultiGraph::LargeEdgePointerPathToMultiEdgePointerPath(diploidEdgePointerPath largeEdgePointerPath, bool sameLength)
{
	diploidEdgePointerPath forReturn;
	int processedLargeEdges = 0;
	int positionInMultiGraph = 0;

	while(processedLargeEdges != (int)largeEdgePointerPath.h1.size())
	{
		assert(largeEdgePointerPath.h1.at(processedLargeEdges)->From->lG == underlyingLargeGraph);
		assert(largeEdgePointerPath.h2.at(processedLargeEdges)->From->lG == underlyingLargeGraph);

		Edge* multiH1 = 0;
		Edge* multiH2 = 0;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(positionInMultiGraph).begin(); nodeIt != NodesPerLevel.at(positionInMultiGraph).end(); nodeIt++)
		{
			for(set<Edge*>::iterator edgeIt = (*nodeIt)->Outgoing_Edges.begin(); edgeIt != (*nodeIt)->Outgoing_Edges.end(); edgeIt++)
			{
				Edge* multiGraphEdge = *edgeIt;
				if(multiGraphEdge->_largeGraphEdge == largeEdgePointerPath.h1.at(processedLargeEdges))
				{
					multiH1 = multiGraphEdge;
				}
				if(multiGraphEdge->_largeGraphEdge == largeEdgePointerPath.h2.at(processedLargeEdges))
				{
					multiH2 = multiGraphEdge;
				}
			}
		}

		assert(multiH1 != 0);
		assert(multiH2 != 0);
		assert(multiH1->multiEmission_full.size() == multiH2->multiEmission_full.size());

		processedLargeEdges += multiH1->multiEmission_full.size();
		positionInMultiGraph++;

		int edgeMultiplicity;
		if(sameLength == true)
		{
			edgeMultiplicity = multiH1->multiEmission_full.size();
		}
		else
		{
			edgeMultiplicity = 1;
		}

		for(int eIm = 1; eIm <= edgeMultiplicity; eIm++)
		{
			forReturn.h1.push_back(multiH1);
			forReturn.h2.push_back(multiH2);
		}

	}

	return forReturn;
}
vector<string> MultiGraph::getAssignedLoci()
{
	vector<string> loci(NodesPerLevel.size()-1, "");
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		//Node* n2 = e->To;
		
		int level = n1->level;
		string locusID = e->locus_id;
		
		if(loci.at(level) != "")
		{
			assert(loci.at(level) == locusID);
		}
		else
		{
			loci.at(level) = locusID;
		}
	}
	return loci;
}

void MultiGraph::checkConsistency(bool terminalCheck)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;
    
	// all nodes implied by the Edges set are in the Nodes set
	// all incoming/outgoing relations implied by the Edges are reflected in the nodes
	// all edges have incoming and outgoing nodes
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		Node* n2 = e->To;
		
		assert(e->count >= 0);
		for(map<int, int>::iterator emissionIt = e->multiEmission.begin(); emissionIt !=  e->multiEmission.end(); emissionIt++)
		{
			if(!(CODE.knowCode(e->locus_id, emissionIt->first)))
			{
				cerr << "Coding problem: at " << e->locus_id << ", do not know what the code " << emissionIt->first << " stands for\n";
			}
			assert(CODE.knowCode(e->locus_id, emissionIt->first));
		}


		assert(n1 != NULL);
		if(n1 != NULL)
		{
			assert(count(n1->Outgoing_Edges.begin(), n1->Outgoing_Edges.end(), e) > 0);
		}

		assert(n2 != NULL);
		if(n2 != NULL)
		{
			assert(count(n2->Incoming_Edges.begin(), n2->Incoming_Edges.end(), e) > 0);
		}

		assert(n1->level == (n2->level-1));

		assert(Nodes.count(n1) > 0);
		assert(Nodes.count(n2) > 0);
	}

	// The NodesPerLevel set and the Nodes set contain the same nodes
	// The edges implied by the Nodes are in the Edges set
	set<Node*> saw_node;
	for(unsigned int level = 0; level < NodesPerLevel.size(); level++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* node = *nodeIt;
			saw_node.insert(node);

			assert(node->level == level);
			assert(Nodes.count(node) > 0);


			for(set<Edge*>::iterator eIt = node->Incoming_Edges.begin(); eIt != node->Incoming_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				assert(Edges.count(e) > 0);
			}

			for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				assert(Edges.count(e) > 0);
			}
		}
	}


	// there are no loneley nodes
	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* node = *nodeIt;
		assert(saw_node.count(node) > 0);
		assert((! terminalCheck) || (node->terminal) || (node->Outgoing_Edges.size() > 0));
		assert((node->level == 0) || (node->Incoming_Edges.size() > 0));
	}
}

void MultiGraph::freeMemory()
{
	vector<Edge*> edges_for_free(Edges.begin(), Edges.end());
	vector<Node*> nodes_for_free(Nodes.begin(), Nodes.end());

	for(unsigned int i = 0; i < edges_for_free.size(); i++)
	{
		Edge* e = edges_for_free.at(i);
		delete(e);
	}

	for(unsigned int i = 0; i < nodes_for_free.size(); i++)
	{
		Node* n = nodes_for_free.at(i);
		delete(n);
	}
}

void MultiGraph::registerEdge(Edge* e)
{
	assert(Edges.count(e) == 0);
	Edges.insert(e);
}

diploidEdgePointerPath MultiGraph::MultiGraphEdgesToLargeGraphEdges(diploidEdgePointerPath& multiEdgesPath)
{
	diploidEdgePointerPath forReturn;

	diploidEdgePointerPath path_mG = multiEdgesPath;

	assert(path_mG.h1.size() == path_mG.h2.size());
	int path_length = path_mG.h1.size();

	int allEmissions_expected = 0;

	for(int path_pos = 0; path_pos < path_length; path_pos++)
	{
		Edge* e1 = path_mG.h1.at(path_pos);

		int found_lG_edges = 1;
		Edge* runningEdge = e1->_largeGraphEdge;
		forReturn.h1.push_back(runningEdge);

		while((int)e1->multiEmission_full.size() != found_lG_edges)
		{
			assert(runningEdge->To->Outgoing_Edges.size() == 1);
			runningEdge = *(runningEdge->To->Outgoing_Edges.begin());
			forReturn.h1.push_back(runningEdge);
			found_lG_edges++;
		}

		assert(found_lG_edges == (int)e1->multiEmission_full.size());

		Edge* e2 = path_mG.h2.at(path_pos);

		found_lG_edges = 1;
		runningEdge = e2->_largeGraphEdge;
		forReturn.h2.push_back(runningEdge);

		while((int)e2->multiEmission_full.size() != found_lG_edges)
		{
			assert(runningEdge->To->Outgoing_Edges.size() == 1);
			runningEdge = *(runningEdge->To->Outgoing_Edges.begin());
			forReturn.h2.push_back(runningEdge);
			found_lG_edges++;
		}

		assert(found_lG_edges == (int)e2->multiEmission_full.size());
		assert(e1->multiEmission_full.size() == (int)e2->multiEmission_full.size());
		allEmissions_expected += e1->multiEmission_full.size();
	}

	assert((int)forReturn.h1.size() == allEmissions_expected);
	assert((int)forReturn.h2.size() == allEmissions_expected);

	return forReturn;
}


void MultiGraph::stats()
{
	cout << "Graph statistics\n" << flush;
	for(int l = 0; l < (int)NodesPerLevel.size(); l++)
	{
		int nodes = NodesPerLevel.at(l).size();
		int edges = 0;
		int symbols = 0;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			edges = edges + n->Outgoing_Edges.size();

			if(l != (NodesPerLevel.size()-1))
			{
				string locus = (*(n->Outgoing_Edges.begin()))->locus_id;
				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					assert((*outgoingIt)->locus_id == locus);
				}

				if(nodeIt == NodesPerLevel.at(l).begin())
				{
					symbols = CODE.getAlleles(locus).size();
				}
			}
		}

		cout << "Level " << l << ", " << nodes << " nodes, " << edges << " edges, " << symbols << "kMers\n";
	}
}

