/*
 * Graph.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "Graph.h"
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

typedef std::basic_string <unsigned char> ustring;

Graph::Graph() {

}

void Graph::unRegisterNode(Node* n)
{
	assert(Nodes.count(n) > 0);
	int l = n->level;
	Nodes.erase(n);
	NodesPerLevel.at(l).erase(n);
	delete(n);
}

void Graph::makeEdgesGaps(double proportion)
{
	srand ( time(NULL) );

	assert(proportion >= 0);
	assert(proportion <= 1);

	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = (*edgeIt);
		double f = (double)rand() / RAND_MAX;
		if(f <= proportion)
		{
			unsigned char newEmission = CODE.doCode(e->locus_id, "_");
			e->emission = newEmission;
		}
	}
}

void Graph::simulateHaplotypes(int number)
{
	srand ( time(NULL) );
	assert(NodesPerLevel.at(0).size()==1);

	for(int iteration = 0; iteration < number; iteration++)
	{
		vector<string> symbols;
		int levels = NodesPerLevel.size();
		Node* currentNode = *(NodesPerLevel.at(0).begin());
		while(currentNode->level != (levels-1))
		{
			double To_Outgoing_Sum = currentNode->Sum_Outgoing();
			assert(To_Outgoing_Sum > 0);
			vector<double> stateP;
			vector<Edge*> outEdges( currentNode->Outgoing_Edges.begin(), currentNode->Outgoing_Edges.end());
			assert(outEdges.size() > 0);
			for(int eI = 0; eI < (int)outEdges.size(); eI++)
			{
				stateP.push_back(outEdges.at(eI)->count/To_Outgoing_Sum);
			}
			int chosenEdgeI = Utilities::chooseFromVector(stateP);
			Edge* chosenEdge = outEdges.at(chosenEdgeI);
			string symbolToEmit = CODE.deCode(chosenEdge->locus_id, chosenEdge->emission);
			if(symbolToEmit != "_")
				symbols.push_back(symbolToEmit);

			currentNode = chosenEdge->To;
		}
		cout << iteration << " " << Utilities::join(symbols, "") << "\n";
	}
}

int Graph::trimGraph(bool remove2DHLA)
{
	int removeNodes = 0;
	int removeEdges = 0;

	set<Node*> leadToEnd;
	set<Node*> leadToStart;

	set<Edge*> zeroEdges;

	for(int l = NodesPerLevel.size()-1; l >= 0; l--)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			if(l == (int)NodesPerLevel.size()-1)
			{
				leadToEnd.insert(n);
			}
			else
			{
				double To_Outgoing_Sum = n->Sum_Outgoing();
				if(To_Outgoing_Sum > 0)
				{
					for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
					{
						Edge* e = (*outgoingIt);
						double p = e->count / To_Outgoing_Sum;
						if(p > 0)
						{
							if(leadToEnd.count(e->To) > 0)
							{
								if(! remove2DHLA)
								{
									leadToEnd.insert(n);
								}
								else
								{
									if((e->locus_id.length()>3) && ((e->locus_id.substr(0, 3) == "HLA") || (e->locus_id.substr(0, 3) == "KIR")))
									{
										if(CODE.allele4D(e->locus_id, e->emission))
										{
											leadToEnd.insert(n);
										}
										else
										{
											zeroEdges.insert(e);
										}
									}
									else
									{
										leadToEnd.insert(n);
									}
								}
							}
						}
						else
						{
							zeroEdges.insert(e);
						}
					}
				}
			}
		}
	}


	for(int l = 0; l < (int)NodesPerLevel.size(); l++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			if(l == 0)
			{
				leadToStart.insert(n);
			}

			double To_Outgoing_Sum = n->Sum_Outgoing();
			if(To_Outgoing_Sum > 0)
			{
				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					Edge* e = (*outgoingIt);
					double p = e->count / To_Outgoing_Sum;
					if(p > 0)
					{
						if(leadToStart.count(n) > 0)
						{
							if(! remove2DHLA)
							{
								leadToStart.insert(e->To);
							}
							else
							{
								if((e->locus_id.length()>3) && ((e->locus_id.substr(0, 3) == "HLA") || (e->locus_id.substr(0, 3) == "KIR")))
								{
									if(CODE.allele4D(e->locus_id, e->emission))
									{
										leadToStart.insert(e->To);
									}
									else
									{
										zeroEdges.insert(e);
									}
								}
								else
								{
									leadToStart.insert(e->To);
								}
							}

						}
					}
					else
					{
						zeroEdges.insert(e);
					}
				}
			}
		}
	}


	for(set<Edge*>::iterator edgeRemoveIt = zeroEdges.begin(); edgeRemoveIt != zeroEdges.end(); edgeRemoveIt++)
	{
		Edge* e = *edgeRemoveIt;
		e->From->Outgoing_Edges.erase(e);
		e->To->Incoming_Edges.erase(e);
		unRegisterEdge(e);
		removeEdges++;
	}

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		if((leadToEnd.count(n) == 0) || (leadToStart.count(n) == 0))
		{
			for(set<Edge*>::iterator edgeRemoveIt = n->Outgoing_Edges.begin(); edgeRemoveIt != n->Outgoing_Edges.end(); edgeRemoveIt++)
			{
				Edge* e = *edgeRemoveIt;
				if(Edges.count(e) > 0)
				{
					e->From->Outgoing_Edges.erase(e);
					e->To->Incoming_Edges.erase(e);
					unRegisterEdge(e);
					removeEdges++;
				}
			}
			for(set<Edge*>::iterator edgeRemoveIt = n->Incoming_Edges.begin(); edgeRemoveIt != n->Incoming_Edges.end(); edgeRemoveIt++)
			{
				Edge* e = *edgeRemoveIt;
				if(Edges.count(e) > 0)
				{
					e->From->Outgoing_Edges.erase(e);
					e->To->Incoming_Edges.erase(e);
					unRegisterEdge(e);
					removeEdges++;
				}
			}
			unRegisterNode(n);
			removeNodes++;
		}
	}

	if(remove2DHLA)
	{
		int edges_2D = 0;
		int hla_Edges = 0;
		for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
		{
			Edge* e = (*edgeIt);
			if((e->locus_id.length()>3) && ((e->locus_id.substr(0, 3) == "HLA") || (e->locus_id.substr(0, 3) == "KIR")))
			{
				hla_Edges++;
				if(! CODE.allele4D(e->locus_id, e->emission))
				{
					edges_2D++;
				}
			}
		}
		cout << "Remaining HLA edges: " << hla_Edges << " (of which " << edges_2D << " are 2-digit)\n";
		assert(edges_2D == 0);
	}
	checkConsistency(false);

	cerr << "Trimming: removed " << removeEdges << " edges and " << removeNodes << " nodes.\n";
	return removeEdges+removeNodes;
}

void Graph::unRegisterEdge(Edge* e)
{
	assert(Edges.count(e) > 0);
	Edges.erase(e);
	delete(e);
}

void Graph::registerNode(Node* n, unsigned int level)
{
	assert(Nodes.count(n) == 0);
	assert(n->level == level);
	if((level+1) > NodesPerLevel.size())
	{
		NodesPerLevel.resize(level+1);
	}
	NodesPerLevel.at(level).insert(n);
	Nodes.insert(n);
	n->g = this;
}

vector<string> Graph::getAssignedLoci()
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

void Graph::checkLocusOrderConsistency(vector<string> loci)
{
	assert((NodesPerLevel.size()-1) == loci.size());
	bool all_ok = true;
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		//Node* n2 = e->To;
		
		int level = n1->level;
		
		if(loci.at(level) != e->locus_id)
		{
			all_ok = false;
			cerr << "Level " << level << " edge says locus " << e->locus_id << ", but locus list says " << loci.at(level) << "\n"; 
		}
	}
	
	if(! all_ok)
	{
		cerr << "Locus list:\n";
		for(int j = 0; j < (int)loci.size(); j++)
		{
			cerr << j << " " << loci.at(j) << "\n";
		}
	}

	assert(all_ok == true);
}


void Graph::printComplexity (string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;
    
	ofstream output;
	output.open (filename.c_str(), ios::out | ios::trunc);
	if (output.is_open())
	{
		output << "Level\tNodes\tEdges\n";

		for(unsigned int level = 0; level < NodesPerLevel.size(); level++)
		{
			int nodes = NodesPerLevel.at(level).size();
			int edges = 0;
			for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;
				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					edges++;
				}
			}

			output << level << "\t total nodes: " << nodes << "\t total edges: " << edges << "\n";

			for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;
				output << "\t\t Outgoing edges: " << node->Outgoing_Edges.size() << "\t Incoming edges: " << node->Incoming_Edges.size() << "\n";
				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					output << "\t\t\t\t" << CODE.deCode(e->locus_id, e->emission) << "          from " << e->From << "to " << e->To << "    [edge " << e << "]\n";
				}
			}
		}
	}
	else
	{
		errEx("Cannot open output file for graph serialization: "+filename);
	}



}



void Graph::checkConsistency(bool terminalCheck)
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

		string emission = lexical_cast<string>(e->emission);
		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->emission << "(end)\n";
			errEx("Unsigned Char Conversion problem!");
		}

		assert(e->count >= 0);
		if(!(CODE.knowCode(e->locus_id, e->emission)))
		{
			cerr << "Coding problem: at " << e->locus_id << ", do not know the code (or locus id) for " << e->emission << "\n";
		}
		assert(CODE.knowCode(e->locus_id, e->emission));

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

void Graph::freeMemory()
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

void Graph::registerEdge(Edge* e)
{
	assert(Edges.count(e) == 0);
	Edges.insert(e);
}


void Graph::removeStarPaths()
{

	// Find edges that are the start of paths that we potentially may want to remove
	// Such paths are defined as: there is an edge with a star symbol, and there are
	// other edges with non-star symbols emanating from the same node

	map<Edge*, bool> edgesPotentialPathStarts;
	vector<int> starSymbols;
	for(int level = 0; level < (int)NodesPerLevel.size(); level++)
	{
		int star_symbol;
		string locusID = "";
		Node* firstNode = *(NodesPerLevel.at(level).begin());

		if(firstNode->Outgoing_Edges.size() > 0)
		{
			Edge* firstEdge = *(firstNode->Outgoing_Edges.begin());
			locusID = firstEdge->locus_id;
			star_symbol = CODE.doCode(locusID, "*");
			starSymbols.push_back(star_symbol);
		}
		else
		{
			assert(level == (NodesPerLevel.size() - 1));
		}

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			set<Node*> has_non_star_edge_incoming;
			set<Edge*> star_edges_outgoing;
			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				assert(locusID != "");
				Edge* e = *edgeIt;
				if(e->emission == star_symbol)
				{
					star_edges_outgoing.insert(e);
				}
				else
				{
					has_non_star_edge_incoming.insert(e->To);
				}
			}

			if(n->Outgoing_Edges.size() > star_edges_outgoing.size())
			{
				for(set<Edge*>::iterator edgeIt = star_edges_outgoing.begin(); edgeIt != star_edges_outgoing.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					if(! e->pgf_protect)
					{
						edgesPotentialPathStarts[e] = true;
					}
				}
			}
		}
	}

	// Find out whether these potential start nodes actually
	// define a path

	set<Node*> delete_Nodes;
	set<Edge*> delete_Edges;
	set<Edge*> edgeFromOtherPath;

	for(map<Edge*, bool>::iterator edgeIt = edgesPotentialPathStarts.begin(); edgeIt != edgesPotentialPathStarts.end(); edgeIt++)
	{
		Edge* e = edgeIt->first;

		if(edgeIt->second != true)
		{
			continue;
		}
		if(delete_Edges.count(e) > 0)
		{
			continue;
		}

		// openPaths contains the paths that we may still want to elongate
		// closedPaths contains the paths that we are still processing
		// the path criterion is: on the way only starred edges + ending in a
		// node that has other incoming edges than the one from the path

		vector< vector<Edge*> > openPaths;
		vector< vector<Edge*> > closedPaths;

		vector<Edge*> firstPath;

		firstPath.push_back(e);
		openPaths.push_back(firstPath);
		set<Edge*> edgeFromPath;
		edgeFromPath.insert(e);

		while(openPaths.size() > 0)
		{
			int extendedPathLength = -1;

			int openPathsSize = openPaths.size();
			for(int openPathI = 0; openPathI < openPathsSize; openPathI++)
			{
				int lastEdgeIndex = openPaths.at(openPathI).size()-1;
				Edge* lastEdge = openPaths.at(openPathI).at(lastEdgeIndex);
				Node* lastEdgeTo = lastEdge->To;

				// legitimate stop nodes are only those which have another
				// incoming edge which does not come this path

				// if we find a legitimate stop node,
				// we end this path and break
				bool stopNode = false;
				for(set<Edge*>::iterator edgeIt = lastEdgeTo->Incoming_Edges.begin(); edgeIt != lastEdgeTo->Incoming_Edges.end(); edgeIt++)
				{
					if((edgeFromPath.count(*edgeIt) == 0) && (edgeFromOtherPath.count(*edgeIt) == 0))
					{
						stopNode = true;
						closedPaths.push_back(openPaths.at(openPathI));
						openPaths.at(openPathI).clear();
						break;
					}
				}

				// if the last node of this path is not legitimate stop node,
				// we extend
				if(! stopNode)
				{
					// find edges to iterate over... if any of the edges is no star symbol, we
					// need to abort ...
					bool pursueEdges = true;
					for(set<Edge*>::iterator edgeIt = lastEdgeTo->Outgoing_Edges.begin(); edgeIt != lastEdgeTo->Outgoing_Edges.end(); edgeIt++)
					{
						Edge* e = *edgeIt;
						int star_symbol = starSymbols.at(lastEdgeTo->level);
						if(e->emission != star_symbol)
						{
							assert(CODE.doCode(e->locus_id, "*") == star_symbol);
							pursueEdges = false;
						}
						if(e->pgf_protect)
						{
							pursueEdges = false;
						}
					}
					// ... same thing if we have no further edges to follow.
					if(lastEdgeTo->Outgoing_Edges.size() == 0)
					{
						pursueEdges = false;
					}

					// if we abort, we clear all found paths ... none is valid; and we break.
					if(! pursueEdges)
					{
						openPaths.clear();
						closedPaths.clear();
						break;
					}

					// on the other hand, if we do not have to abort, we follow the edges
					vector<Edge*> pursueEdgesVector(lastEdgeTo->Outgoing_Edges.begin(), lastEdgeTo->Outgoing_Edges.end());
					vector<Edge*> originalPath = openPaths.at(openPathI);
					for(int pursueEdgesI = 0; pursueEdgesI < (int)pursueEdgesVector.size(); pursueEdgesI++)
					{
						if(pursueEdgesI == 0)
						{
							openPaths.at(openPathI).push_back(pursueEdgesVector.at(pursueEdgesI));
							if(extendedPathLength == -1)
							{
								extendedPathLength = openPaths.at(openPathI).size();
							}
							else
							{
								assert(extendedPathLength == (int)openPaths.at(openPathI).size());
							}
						}
						else
						{
							openPaths.push_back(originalPath);
							openPaths.at(openPaths.size()-1).push_back(pursueEdgesVector.at(pursueEdgesI));
							assert(extendedPathLength == (int)openPaths.at(openPaths.size()-1).size());
						}

						edgeFromPath.insert(pursueEdgesVector.at(pursueEdgesI));
					}
				}
			}

			// if there are any open paths left, we will have to follow them
			for(int openPathI = openPaths.size(); openPathI > 0; openPathI--)
			{
				if(openPaths.at(openPathI-1).size() == 0)
				{
					openPaths.erase(openPaths.begin()+openPathI-1);
				}
				else
				{
					assert(extendedPathLength != -1);
					assert((int)openPaths.at(openPathI-1).size() == extendedPathLength);
				}
			}
		}

		assert(openPaths.size() == 0);

		for(int pI = 0; pI < (int)closedPaths.size(); pI++)
		{
			for(int eI = 0; eI < (int)closedPaths.at(pI).size(); eI++)
			{
				Edge* e = closedPaths.at(pI).at(eI);
				delete_Edges.insert(e);
				edgeFromOtherPath.insert(e);
				if(eI == 0)
				{
					assert(e->From->Outgoing_Edges.size() > 1);
				}
				if(eI != 0)
				{
					delete_Nodes.insert(e->From);
				}
				if(eI != (closedPaths.at(pI).size() - 1))
				{
					delete_Nodes.insert(e->To);
				}
				if(eI == (closedPaths.at(pI).size() - 1))
				{
					assert(e->To->Incoming_Edges.size() > 1);
				}
			}
		}
	}

	int removedEdges = 0;
	int removedNodes = 0;
	map<Node*, vector<string> > deletedKmers_cache;
	for(set<Edge*>::iterator edgeRemoveIt = delete_Edges.begin(); edgeRemoveIt != delete_Edges.end(); edgeRemoveIt++)
	{
		Edge* e = *edgeRemoveIt;
		assert(CODE.deCode(e->locus_id, e->emission) == "*");
		deletedKmers_cache[e->To].push_back(CODE.deCode(e->locus_id, e->emission));
		e->From->Outgoing_Edges.erase(e);
		e->To->Incoming_Edges.erase(e);
		unRegisterEdge(e);
		removedEdges++;
	}
	for(set<Node*>::iterator nodeRemoveIt = delete_Nodes.begin(); nodeRemoveIt != delete_Nodes.end(); nodeRemoveIt++)
	{
		Node* n = *nodeRemoveIt;
		assert(n->Outgoing_Edges.size() == 0);
		assert(n->Incoming_Edges.size() == 0);
		unRegisterNode(n);
		removedNodes++;
	}

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* node = *nodeIt;
		if(!(((node->level == 0) || (node->Incoming_Edges.size() > 0))))
		{
			cout << node << "\n";
			assert(deletedKmers_cache.count(node) > 0);
			for(int i = 0; i < (int)deletedKmers_cache[node].size(); i++)
			{
				cout << "\t" << deletedKmers_cache[node].at(i) << "\n";
			}
		}

		assert((node->level == 0) || (node->Incoming_Edges.size() > 0));

	}

	cout << "Star removal: removed " << removedNodes << " nodes and " << removedEdges << " edges.\n";

	checkConsistency(false);
}



void Graph::writeToFile(string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	checkConsistency(false);

	vector<string> linesFromCode = CODE.serializeIntoVector();
	vector<string> linesForNodes;
	vector<string> linesForEdges;

	map<Edge*, int> EdgeCounter;
	map<Node*, int> NodeCounter;

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;

		if(NodeCounter.count(n) == 0)
		{
			int s = NodeCounter.size();
			NodeCounter[n] = s + 1;
		}
		int n_index = NodeCounter[n];

		vector<string> line_fields;
		line_fields.push_back(lexical_cast<string>(n_index));
		line_fields.push_back(lexical_cast<string>(n->level));
		line_fields.push_back(lexical_cast<string>(n->terminal));

		linesForNodes.push_back(boost::join(line_fields, separatorForSerialization));
	}

	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;

		if(EdgeCounter.count(e) == 0)
		{
			int s = EdgeCounter.size();
			EdgeCounter[e] = s + 1;
		}
		int e_index = EdgeCounter[e];

		vector<string> line_fields;
		line_fields.push_back(lexical_cast<string>(e_index));
		line_fields.push_back(lexical_cast<string>(e->locus_id));
		line_fields.push_back(lexical_cast<string>(e->count));
		string emission = lexical_cast<string>(e->emission);

		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->emission << "(end)\n"; 
			errEx("Conversion problem!");
		}
		line_fields.push_back(emission);


		Node* n1 = e->From;
		Node* n2 = e->To;
		assert(NodeCounter.count(n1) > 0);
		assert(NodeCounter.count(n2) > 0);

		int n1_i = NodeCounter[n1];
		int n2_i = NodeCounter[n2];

		assert(n1_i < NodeCounter.size()+1);
		assert(n2_i < NodeCounter.size()+1);


		line_fields.push_back(lexical_cast<string>(n1_i));
		line_fields.push_back(lexical_cast<string>(n2_i));

		//if(e->label != "")
		//{
			line_fields.push_back(e->label);

		//}

		line_fields.push_back(lexical_cast<string>(e->pgf_protect));

			
		linesForEdges.push_back(boost::join(line_fields, separatorForSerialization));
	}

	ofstream output;
	output.open (filename.c_str(), ios::out | ios::trunc);
	if (output.is_open())
	{
		output << "CODE:" << "\n";
		output << boost::join(linesFromCode, "\n") << "\n";
		output << "NODES:" << "\n";
		output << boost::join(linesForNodes, "\n") << "\n";
		output << "EDGES:" << "\n";
		output << boost::join(linesForEdges, "\n");
	}
	else
	{
		errEx("Cannot open output file for graph serialization: "+filename);
	}
}

void Graph::readFromFile(string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	vector<string> linesForCode;
	vector<string> linesForNodes;
	vector<string> linesForEdges;

	string problematic_part = "|||||||";
	string subsitute_problem = "|||SLASH|||";
	string substitue_indicator = "SLASH";

	int mode = -1;
	ifstream graphStr;
	graphStr.open (filename.c_str(), ios::in);
	if(graphStr.is_open())
	{
		string line;
		int lineCounter = 0;
		while(graphStr.good())
		{
			lineCounter++;

			if((lineCounter % 100000) == 0)
				cout << "\r" << lineCounter << flush;
			
			getline (graphStr, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			if(line.find(problematic_part) != string::npos)
			{
				line.replace(line.find(problematic_part),problematic_part.length(),subsitute_problem);
			}
			
			if(line == "CODE:")
			{
				mode = 1;
			}
			else if(line == "NODES:")
			{
				mode = 2;
			}
			else if(line == "EDGES:")
			{
				mode = 3;
			}
			else
			{
				assert(mode > 0);
				if(mode == 1)
				{
					linesForCode.push_back(line);
				}
				else if(mode == 2)
				{
					linesForNodes.push_back(line);
				}
				else if(mode == 3)
				{
					linesForEdges.push_back(line);
				}
				else
				{
					errEx("Undefined mode state.");
				}
			}
		}
		graphStr.close();
	}
	else
	{
		errEx("Cannot open graph file: "+filename);
	}
		
	CODE.readFromVector(linesForCode);

	map<int, Node*> idx2Node;
	map<int, Edge*> idx2Edge;
	for(unsigned int i = 0; i < linesForNodes.size(); i++)
	{
		if((i % 100000) == 0)
			cout << "\r" << i << "/" << linesForNodes.size() << "[N]" << flush;
		
		string line = linesForNodes.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if(fields.size() != 3)
		{
			errEx("Cannot node-parse this line (expext 6 fields): "+line);
		}
		int n_idx = lexical_cast<int>(fields.at(0));
		int level = lexical_cast<int>(fields.at(1));
		bool terminal = lexical_cast<bool>(fields.at(2));

		Node* n = new Node();
		n->level = level;
		n->terminal = terminal;
		registerNode(n, level);

		idx2Node[n_idx] = n;
	}
	
	bool encounteredErr = false;
	for(unsigned int i = 0; i < linesForEdges.size(); i++)
	{
		if((i % 100000) == 0)
			cout << "\r" << i << "/" << linesForEdges.size() << "[E]" << flush;
		
		string line = linesForEdges.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if((fields.size() != 6) &&  (fields.size() != 8))
		{
			errEx("Cannot edge-parse this line (expext 6/8 fields): "+line);
		}
		int e_idx = lexical_cast<int>(fields.at(0));
		string locusID = fields.at(1);
		double count = lexical_cast<double>(fields.at(2));
		if(fields.at(3) == substitue_indicator)
		{
			fields.at(3) = "|";
		}
		
		unsigned char emission;
        try
        {
        	emission = lexical_cast<unsigned char>(fields.at(3));
        }
        catch(bad_lexical_cast &)
        {
            cerr << "Cannot cast to unsigned char: " << fields.at(3) << "--" << line << "\n";
            for(int i = 0; i < (int)line.length(); i++)
            {
            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
            }
            encounteredErr = true;
            
            exit(1);
        }
        

		int from_idx;
		int to_idx;

		try {
		from_idx = lexical_cast<int>(fields.at(4));
		to_idx = lexical_cast<int>(fields.at(5));
		}
		catch(...)
		{
			cout << "Exception at line\n" << line << "\n\n";
			exit(1);
		}

		string label = "";
		bool protected_pgf = false;
		
		if(fields.size()>6)
		{
			label = fields.at(6);
	        try
	        {
	        	protected_pgf = lexical_cast<bool>(fields.at(7));
	        }
	        catch(bad_lexical_cast &)
	        {
	            cerr << "Cannot cast to bool: " << fields.at(7) << "--" << line << "\n";
	            for(int i = 0; i < (int)line.length(); i++)
	            {
	            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
	            }
	            
				//encounteredErr = true;

				protected_pgf = false;
				
	            //exit(1);
	        }
		
		}

		Edge* e = new Edge();
		e->count = count;
		e->emission = emission;
		e->locus_id = locusID;
		e->label = label;
		e->pgf_protect = protected_pgf;
		idx2Edge[e_idx] = e;
		registerEdge(e);

		assert(idx2Node.count(from_idx) > 0);
		if(!(idx2Node.count(to_idx) > 0))
		{
			cout << "Node index " << to_idx << " is unknown\n";
			cout << "Edge line: " << line << "\n\n" << flush;
		}
		assert(idx2Node.count(to_idx) > 0);

		Node* n1 = idx2Node[from_idx];
		Node* n2 = idx2Node[to_idx];

		e->From = n1;
		e->To = n2;

		n1->Outgoing_Edges.insert(e);
		n2->Incoming_Edges.insert(e);
	}

	if(encounteredErr)
	{
		cerr << "Error!";
		exit(1);
	}
	
	//TODO reactivate
	/*
	cout << "\rCONSISTENCY...             " << flush;
	checkConsistency(false);
	*/

	cout << "\n\n" << flush;
	
	
}


