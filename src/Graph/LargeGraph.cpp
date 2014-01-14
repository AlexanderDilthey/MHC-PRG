/*
 * Graph.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "LargeGraph.h"

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

#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <vector>
#include <utility>

using namespace boost;
using namespace boost::algorithm;

typedef std::basic_string <unsigned char> ustring;
string newickSubTree(Node* n);

LargeGraph::LargeGraph() {
	kMerSize = -1;
}

void LargeGraph::unRegisterNode(Node* n)
{
	assert(Nodes.count(n) > 0);
	int l = n->level;
	Nodes.erase(n);
	NodesPerLevel.at(l).erase(n);
	delete(n);
}

void LargeGraph::reverseGraph()
{
	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;

		Node* oldFrom = e->From;
		e->From = e->To;
		e->To = oldFrom;
	}

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		set<Edge*> oldIncoming = n->Incoming_Edges;
		n->Incoming_Edges = n->Outgoing_Edges;
		n->Outgoing_Edges = oldIncoming;
	}

	int frontIndex = 0;
	int backIndex = NodesPerLevel.size()-1;

	while(frontIndex < backIndex)
	{
		// we will now swap the node levels between frontIndex and backIndex

		set<Node*> oldFrontSet = NodesPerLevel.at(frontIndex);
		NodesPerLevel.at(frontIndex) = NodesPerLevel.at(backIndex);
		NodesPerLevel.at(backIndex) = oldFrontSet;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(frontIndex).begin(); nodeIt != NodesPerLevel.at(frontIndex).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			assert((int)n->level == backIndex);
			n->level = frontIndex;
			n->terminal = false;
		}

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(backIndex).begin(); nodeIt != NodesPerLevel.at(backIndex).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			assert((int)n->level == frontIndex);
			n->level = backIndex;
			n->terminal = (backIndex == (NodesPerLevel.size()-1)) ? true : false;
		}

		// increment / decrement
		frontIndex++;
		backIndex--;
	}

}

LargeGraph::LargeGraph(LargeGraph& otherGraph)
{
	CODE = otherGraph.CODE;
	kMerSize = otherGraph.kMerSize;
	NodesPerLevel.resize(otherGraph.NodesPerLevel.size());

	map<Node*, Node*> otherNodeToThisNode;
	map<Edge*, Edge*> otherEdgeToThisEdge;

	for(set<Node*>::iterator nodeIt = otherGraph.Nodes.begin(); nodeIt != otherGraph.Nodes.end(); nodeIt++)
	{
		Node* otherNode = *nodeIt;

		// allocate new memory
		Node* myNode = new Node();
		otherNodeToThisNode[otherNode] = myNode;

		// call assignment constructor -- will NOT work for
		// pointers in node object
		*myNode = *otherNode;

		// delete graph and edge assignment
		myNode->g = 0;
		myNode->mG = 0;
		myNode->lG = this;

		myNode->Incoming_Edges.clear();
		myNode->Outgoing_Edges.clear();

		assert(myNode->level == otherNode->level);

		// reconstruct the node-based indices for the node
		// into this graph
		Nodes.insert(myNode);
		NodesPerLevel.at(myNode->level).insert(myNode);

		assert(otherGraph.NodeCounter.count(otherNode) > 0);
		NodeCounter[myNode] = otherGraph.NodeCounter[otherNode];
		assert(otherGraph.idx2Node[otherGraph.NodeCounter[otherNode]] == otherNode);

		idx2Node[otherGraph.NodeCounter[otherNode]] = myNode;
	}

	for(set<Edge*>::iterator edgeIt = otherGraph.Edges.begin(); edgeIt != otherGraph.Edges.end(); edgeIt++)
	{
		Edge* otherEdge = *edgeIt;

		Edge* myEdge = new Edge();
		otherEdgeToThisEdge[otherEdge] = myEdge;

		*myEdge = *otherEdge;

		Edges.insert(myEdge);

		assert(otherGraph.EdgeCounter.count(otherEdge) > 0);
		EdgeCounter[myEdge] = otherGraph.EdgeCounter[otherEdge];
		assert(otherGraph.idx2Edge[otherGraph.EdgeCounter[otherEdge]] == otherEdge);

		idx2Edge[otherGraph.EdgeCounter[otherEdge]] = myEdge;

		myEdge->_largeGraphEdge = 0;
		assert(otherNodeToThisNode.count(otherEdge->From) > 0);
		assert(otherNodeToThisNode.count(otherEdge->To) > 0);
		myEdge->From = otherNodeToThisNode[otherEdge->From];
		myEdge->To = otherNodeToThisNode[otherEdge->To];
	}

	for(set<Node*>::iterator nodeIt = otherGraph.Nodes.begin(); nodeIt != otherGraph.Nodes.end(); nodeIt++)
	{
		Node* otherNode = *nodeIt;
		assert(otherNodeToThisNode.count(otherNode) > 0);
		Node* myNode = otherNodeToThisNode[otherNode];

		for(set<Edge*>::iterator edgeIt = otherNode->Incoming_Edges.begin(); edgeIt != otherNode->Incoming_Edges.end(); edgeIt++)
		{
			Edge* otherEdge = *edgeIt;
			assert(otherEdgeToThisEdge.count(otherEdge) > 0);
			Edge* myEdge = otherEdgeToThisEdge[otherEdge];

			myNode->Incoming_Edges.insert(myEdge);
		}

		for(set<Edge*>::iterator edgeIt = otherNode->Outgoing_Edges.begin(); edgeIt != otherNode->Outgoing_Edges.end(); edgeIt++)
		{
			Edge* otherEdge = *edgeIt;
			assert(otherEdgeToThisEdge.count(otherEdge) > 0);
			Edge* myEdge = otherEdgeToThisEdge[otherEdge];

			myNode->Outgoing_Edges.insert(myEdge);
		}
	}

	checkConsistency(false);
}

int LargeGraph::trimGraph()
{
	int removeNodes = 0;
	int removeEdges = 0;

	set<Node*> leadToEnd;
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
								leadToEnd.insert(n);
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
		if(leadToEnd.count(n) == 0)
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

	checkConsistency(false);

	cerr << "Trimming: removed " << removeEdges << " edges and " << removeNodes << " nodes.\n";
	return removeEdges+removeNodes;
}

void LargeGraph::unRegisterEdge(Edge* e)
{
	assert(Edges.count(e) > 0);
	Edges.erase(e);
	delete(e);
}

void LargeGraph::registerNode(Node* n, unsigned int level)
{
	assert(Nodes.count(n) == 0);
	assert(n->level == level);
	if((level+1) > NodesPerLevel.size())
	{
		NodesPerLevel.resize(level+1);
	}
	NodesPerLevel.at(level).insert(n);
	Nodes.insert(n);
	n->lG = this;
}

string LargeGraph::newick()
{
	assert(NodesPerLevel.size()>0);
	assert(NodesPerLevel.at(0).size() == 1);
	Node* n0 = *(NodesPerLevel.at(0).begin());
	return newickSubTree(n0)+";";
}

string newickSubTree(Node* n)
{
	assert(n->Incoming_Edges.size() <= 1);
	vector<string> singleNodesOutput;
	string subNodesOutput;
	if(n->Outgoing_Edges.size() != 0)
	{
		cout << "newick descend into " << n->Outgoing_Edges.size() << "\n";

		for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
		{
			Edge* outEdge = *edgeIt;
			Node* targetNode = outEdge->To;
			singleNodesOutput.push_back(newickSubTree(targetNode));
		}
		subNodesOutput = "("+Utilities::join(singleNodesOutput, ",")+")";
	}

	string symbolGtLt;
	string thisNodeLabel;
	if(n->Incoming_Edges.size()>0)
	{
		assert(n->Incoming_Edges.size()==1);

		Edge* incomingE = *(n->Incoming_Edges.begin());
		symbolGtLt = ">=";
		if(incomingE->largeEmission < 0)
		{
			incomingE->largeEmission *= -1;
			symbolGtLt = "<";
		}
		assert(incomingE->locus_id != "");
		string deCodedLargeEmission = n->lG->CODE.deCode(incomingE->locus_id, incomingE->largeEmission);
		thisNodeLabel = deCodedLargeEmission+">="+Utilities::ItoStr((int)incomingE->emission)+" " +incomingE->label;
	}

	return subNodesOutput+thisNodeLabel;
}
void LargeGraph::simulateHaplotypes(int number)
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
			string symbol_emission = CODE.deCode(chosenEdge->locus_id, chosenEdge->largeEmission);

			if(symbol_emission != "_")
			{
				if(currentNode->level != (levels-2))
				{
					if((symbol_emission.substr(0, 1) == "A") || (symbol_emission.substr(0, 1) == "C") || (symbol_emission.substr(0, 1) == "G") || (symbol_emission.substr(0, 1) == "T"))
					{
						symbol_emission = symbol_emission.substr(0, 1);
					}
					else
					{
						symbol_emission = symbol_emission.substr(0, 4);
					}
				}
				symbols.push_back(symbol_emission);
			}
			currentNode = chosenEdge->To;
		}
		cout << iteration << " " << Utilities::join(symbols, "") << "\n";
	}
}

vector<string> LargeGraph::getAssignedLoci()
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

void LargeGraph::checkLocusOrderConsistency(vector<string> loci)
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


void LargeGraph::checkConsistency(bool terminalCheck)
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

		string emission = lexical_cast<string>(e->largeEmission);
		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->largeEmission << "(end)\n";
			errEx("Unsigned Char Conversion problem!");
		}
		
		assert(e->count >= 0);
		if(!(CODE.knowCode(e->locus_id, e->largeEmission)))
		{
			cerr << "Coding problem: at " << e->locus_id << ", do not know what the code " << e->largeEmission << " stands for\n";
		}
		assert(CODE.knowCode(e->locus_id, e->largeEmission));

		if(!(CODE.knowCode(e->locus_id, e->largeEmission2)))
		{
			cerr << "Coding problem: at " << e->locus_id << ", do not know what the code " << e->largeEmission2 << " stands for\n";
		}
		assert(CODE.knowCode(e->locus_id, e->largeEmission2));

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

void LargeGraph::removeStarPaths()
{

	// Find edges that are the start of paths that we potentially may want to remove
	// Such paths are defined as: there is an edge with a star symbol, and there are
	// other edges with non-star symbols emanating from the same node

	map<Edge*, bool> edgesPotentialPathStarts;
	vector<int> starSymbols;
	vector<int> gapSymbols;

	for(int level = 0; level < (int)NodesPerLevel.size(); level++)
	{
		int star_symbol;
		int gap_symbol;

		string locusID = "";
		Node* firstNode = *(NodesPerLevel.at(level).begin());

		if(firstNode->Outgoing_Edges.size() > 0)
		{
			Edge* firstEdge = *(firstNode->Outgoing_Edges.begin());
			locusID = firstEdge->locus_id;
			star_symbol = CODE.doCode(locusID, "*");
			gap_symbol = CODE.doCode(locusID, "_");
			starSymbols.push_back(star_symbol);
			gapSymbols.push_back(gap_symbol);
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
				if(e->largeEmission == star_symbol)
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
					edgesPotentialPathStarts[e] = true;
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
		// the path criterion is: on the way only starred edges (or gaps) + ending in a
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
						int gap_symbol = gapSymbols.at(lastEdgeTo->level);

						if((e->largeEmission != star_symbol) && (e->largeEmission != gap_symbol))
						{
							assert(CODE.doCode(e->locus_id, "*") == star_symbol);
							assert(CODE.doCode(e->locus_id, "_") == gap_symbol);
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
		assert((CODE.deCode(e->locus_id, e->largeEmission) == "*") || (CODE.deCode(e->locus_id, e->largeEmission) == "_"));
		deletedKmers_cache[e->To].push_back(CODE.deCode(e->locus_id, e->largeEmission));
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

LargeGraph::~LargeGraph()
{
	freeMemory();
}

void LargeGraph::freeMemory()
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

void LargeGraph::registerEdge(Edge* e)
{
	assert(Edges.count(e) == 0);
	Edges.insert(e);
}

void LargeGraph::assignEdgeNodeIndices()
{
	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;

		if(NodeCounter.count(n) == 0)
		{
			int s = NodeCounter.size();
			NodeCounter[n] = s + 1;
			assert(idx2Node.count(s+1) == 0);
			idx2Node[s+1] = n;
		}
	}

	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;

		if(EdgeCounter.count(e) == 0)
		{
			int s = EdgeCounter.size();
			EdgeCounter[e] = s + 1;
			assert(idx2Edge.count(s+1) == 0);
			idx2Edge[s+1] = e;
		}
	}
}

void LargeGraph::writeToFile(string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	checkConsistency(false);

	vector<string> linesFromCode = CODE.serializeIntoVector();
	vector<string> linesForNodes;
	vector<string> linesForEdges;

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;

		if(NodeCounter.count(n) == 0)
		{
			int s = NodeCounter.size();
			NodeCounter[n] = s + 1;
			assert(idx2Node.count(s+1) == 0);
			idx2Node[s+1] = n;
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
			assert(idx2Edge.count(s+1) == 0);
			idx2Edge[s+1] = e;
		}
		int e_index = EdgeCounter[e];

		vector<string> line_fields;
		line_fields.push_back(lexical_cast<string>(e_index));
		line_fields.push_back(lexical_cast<string>(e->locus_id));
		line_fields.push_back(lexical_cast<string>(e->count));
		
		string emission = lexical_cast<string>(e->largeEmission);
		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->largeEmission << "(end)\n";
			errEx("Conversion problem!");
		}
		line_fields.push_back(emission);

		string emission2 = lexical_cast<string>(e->largeEmission2);
		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->largeEmission2 << "(end)\n";
			errEx("Conversion problem!");
		}
		line_fields.push_back(emission2);

		Node* n1 = e->From;
		Node* n2 = e->To;
		assert(NodeCounter.count(n1) > 0);
		assert(NodeCounter.count(n2) > 0);

		int n1_i = NodeCounter[n1];
		int n2_i = NodeCounter[n2];

		line_fields.push_back(lexical_cast<string>(n1_i));
		line_fields.push_back(lexical_cast<string>(n2_i));

		line_fields.push_back(e->label);
		line_fields.push_back(lexical_cast<string>(e->pgf_protect));

		vector<string> levelsNucleotideGraph_strings;
		for(int lI = 0; lI < e->levelsNucleotideGraph.size(); lI++)
		{
			int l = e->levelsNucleotideGraph.at(lI);
			levelsNucleotideGraph_strings.push_back(Utilities::ItoStr(l));
		}
		line_fields.push_back(Utilities::join(levelsNucleotideGraph_strings, ","));

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

vector<levelInfo> LargeGraph::getLevelInfo()
{
	vector<levelInfo> forReturn;

	for(int l = 0; l < (int)NodesPerLevel.size(); l++)
	{
		int nodes = NodesPerLevel.at(l).size();
		int edges = 0;
		int symbols_CODE = 0;

		std::set<std::string> symbols_set;

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
					std::string emission = CODE.deCode((*outgoingIt)->locus_id, (*outgoingIt)->emission);
					symbols_set.insert(emission);
				}

				if(nodeIt == NodesPerLevel.at(l).begin())
				{
					symbols_CODE = CODE.getAlleles(locus).size();
				}
			}
		}

		int symbols = symbols_set.size();

		levelInfo lI;
		lI.nodes = nodes;
		lI.edges = edges;
		lI.symbols_CODE = symbols_CODE;
		lI.symbols = symbols;

		forReturn.push_back(lI);
	}

	return forReturn;
}

void LargeGraph::stats()
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
void LargeGraph::readFromFile(string filename)
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
			if((lineCounter % 200000) == 0)
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
	int expect_length = -1;
	vector<string> loci = CODE.getLoci();
	for(vector<string>::iterator locusIt = loci.begin(); locusIt != loci.end(); locusIt++)
	{
		string locus = *locusIt;
		vector<int> alleles = CODE.getAlleles(locus);
		for(vector<int>::iterator allelesCodeIt = alleles.begin(); allelesCodeIt != alleles.end(); allelesCodeIt++)
		{
			string kMer = CODE.deCode(locus, *allelesCodeIt);

			if((kMer != "_") && (kMer != "*"))
			{
				if(expect_length == -1)
					expect_length = kMer.length();

				if(expect_length !=  (int)kMer.length())
				{
					errEx("Assume that I am reading in a kMer graph. Edge output does not conform to guessed kMer size of "+Utilities::ItoStr(expect_length)+". Output string is "+kMer+", length "+Utilities::ItoStr(kMer.length()));
				}
			}
		}
	}
	kMerSize = expect_length;

	EdgeCounter.clear();
	NodeCounter.clear();
	idx2Node.clear();
	idx2Edge.clear();

	for(unsigned int i = 0; i < linesForNodes.size(); i++)
	{
		if((i % 200000) == 0)
			cout << "\r" << i << "/" << linesForNodes.size() << "[N]" << flush;
		
		string line = linesForNodes.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if(fields.size() != 3)
		{
			errEx("Cannot node-parse this line (expext 3 fields): "+line);
		}
		int n_idx = lexical_cast<int>(fields.at(0));
		int level = lexical_cast<int>(fields.at(1));
		bool terminal = lexical_cast<bool>(fields.at(2));

		Node* n = new Node();
		n->level = level;
		n->terminal = terminal;
		registerNode(n, level);

		idx2Node[n_idx] = n;
		NodeCounter[n] = n_idx;
	}
	
	bool encounteredErr = false;
	for(unsigned int i = 0; i < linesForEdges.size(); i++)
	{
		if((i % 200000) == 0)
			cout << "\r" << i << "/" << linesForEdges.size() << "[E]" << flush;
		
		string line = linesForEdges.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if(fields.size() != 10)
		{
			errEx("Cannot edge-parse this line 8: "+line+"\n\nIs this a not-kmerified graph or an old version graph file?\n\n"+Utilities::ItoStr(fields.size()));
		}
		int e_idx = lexical_cast<int>(fields.at(0));
		string locusID = fields.at(1);
		double count = lexical_cast<double>(fields.at(2));
		if(fields.at(3) == substitue_indicator)
		{
			fields.at(3) = "|";
		}
		
		int emission;
        try
        {
        	emission = lexical_cast<int>(fields.at(3));
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
        
		int emission2;
        try
        {
        	emission2 = lexical_cast<int>(fields.at(4));
        }
        catch(bad_lexical_cast &)
        {
            cerr << "Cannot cast to unsigned char: " << fields.at(4) << "--" << line << "\n";
            for(int i = 0; i < (int)line.length(); i++)
            {
            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
            }
            encounteredErr = true;

            exit(1);
        }



		int from_idx = lexical_cast<int>(fields.at(5));
		int to_idx = lexical_cast<int>(fields.at(6));

		string label = "";
		label = fields.at(7);

		bool protected_pgf;
        try
        {
        	protected_pgf = lexical_cast<bool>(fields.at(8));
        }
        catch(bad_lexical_cast &)
        {
            cerr << "Cannot cast to bool: " << fields.at(8) << "--" << line << "\n";
            for(int i = 0; i < (int)line.length(); i++)
            {
            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
            }
            encounteredErr = true;

            exit(1);
        }

        string lNucGraph_joined = fields.at(9);
        vector<int> levelsNucleotideGraph;

		vector<string> levelsNucleotideGraph_strings = Utilities::split(lNucGraph_joined, ",");
		for(int lI = 0; lI < levelsNucleotideGraph_strings.size(); lI++)
		{
			string l = levelsNucleotideGraph_strings.at(lI);
			int pushI = lexical_cast<int>(l);
			if(!(pushI >= 0))
			{
				cout << "pushI: " << pushI << "\n";
				cout << "lNucGraph_joined: " << lNucGraph_joined << "\n" << flush;
				assert(1 == 0);
			}
			assert(pushI >= 0);
			levelsNucleotideGraph.push_back(pushI);
		}

		Edge* e = new Edge();
		e->count = count;
		e->largeEmission = emission;
		e->largeEmission2 = emission2;
		e->locus_id = locusID;
		e->label = label;
		e->pgf_protect = protected_pgf;
		e->levelsNucleotideGraph = levelsNucleotideGraph;

		idx2Edge[e_idx] = e;
		EdgeCounter[e] = e_idx;

		registerEdge(e);

		assert(idx2Node.count(from_idx) > 0);
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
	
	cout << "\rCONSISTENCY...             " << flush;
	checkConsistency(false);
	cout << "\n\n" << flush;
}


vector<int> LargeGraph::edgePointerPathToEdgeIndexPath(vector<Edge*>& edgePointerPath)
{
	vector<int> forReturn;
	int path_length = edgePointerPath.size();
	for(int path_pos = 0; path_pos < path_length; path_pos++)
	{
		Edge* e1 = edgePointerPath.at(path_pos);
		assert(EdgeCounter.count(e1) > 0);
		int runningEdgeIndex = EdgeCounter[e1];
		forReturn.push_back(runningEdgeIndex);
	}
	return forReturn;
}

diploidEdgePath LargeGraph::edgePointerPathToEdgeIndexPath(diploidEdgePointerPath& edgePointerPath)
{
	diploidEdgePath forReturn;
	assert(edgePointerPath.h1.size() == edgePointerPath.h2.size());
	forReturn.h1 = edgePointerPathToEdgeIndexPath(edgePointerPath.h1);
	forReturn.h2 = edgePointerPathToEdgeIndexPath(edgePointerPath.h2);
	assert(forReturn.h1.size() == forReturn.h2.size());
	return forReturn;

}

vector<Edge*> LargeGraph::edgeIndexPathToEdgePointerPath(vector<int>& edgeIndexPath)
{
	vector<Edge*> forReturn;
	int path_length = edgeIndexPath.size();
	for(int path_pos = 0; path_pos < path_length; path_pos++)
	{
		int eI = edgeIndexPath.at(path_pos);
		assert(idx2Edge.count(eI) > 0);
		forReturn.push_back(idx2Edge[eI]);
	}
	return forReturn;
}

diploidEdgePointerPath LargeGraph::edgeIndexPathToEdgePointerPath(diploidEdgePath& edgeIndexPath)
{
	diploidEdgePointerPath forReturn;
	assert(edgeIndexPath.h1.size() == edgeIndexPath.h2.size());
	forReturn.h1 = edgeIndexPathToEdgePointerPath(edgeIndexPath.h1);
	forReturn.h2 = edgeIndexPathToEdgePointerPath(edgeIndexPath.h2);
	assert(forReturn.h1.size() == forReturn.h2.size());
	return forReturn;
}

diploidNucleotidePath LargeGraph::diploidPathToAlignedNucleotides(diploidEdgePointerPath& path)
{
	diploidNucleotidePath forReturn;

	assert(path.h1.size() == path.h2.size());
	forReturn.h1.resize(path.h1.size()+kMerSize, "_");
	forReturn.h2.resize(path.h1.size()+kMerSize, "_");

	for(int pos = 0; pos < (int)path.h1.size(); pos++)
	{
		Edge* e1 = path.h1.at(pos);
		string locusID1 = e1->locus_id;
		string kMer1 = CODE.deCode(locusID1, e1->largeEmission2);

		Edge* e2 = path.h2.at(pos);
		string locusID2 = e2->locus_id;
		string kMer2 = CODE.deCode(locusID2, e2->largeEmission2);

		if(pos == 0)
		{
			assert(e1->levelsNucleotideGraph.size() == kMer1.size());
			for(int cI = 0; cI < kMer1.size(); cI++)
			{
				int correspondingPosition = e1->levelsNucleotideGraph.at(cI);
				forReturn.h1.at(correspondingPosition) = kMer1.substr(cI, 1);
			}

			assert(e2->levelsNucleotideGraph.size() == kMer2.size());
			for(int cI = 0; cI < kMer2.size(); cI++)
			{
				int correspondingPosition = e2->levelsNucleotideGraph.at(cI);
				forReturn.h2.at(correspondingPosition) = kMer2.substr(cI, 1);
			}
 		}
		else
		{
			assert((e1->levelsNucleotideGraph.size() == 0) || (e1->levelsNucleotideGraph.size() == 1));
			if(e1->levelsNucleotideGraph.size() > 0)
			{
				string kMer1LastCharacter = kMer1.substr(kMer1.size() - 1, 1);
				int correspondingPosition1 = e1->levelsNucleotideGraph.at(0);
				assert(correspondingPosition1 >= 0);
				assert(correspondingPosition1 < forReturn.h1.size());
				forReturn.h1.at(correspondingPosition1) = kMer1LastCharacter;
			}

			assert((e2->levelsNucleotideGraph.size() == 0) || (e2->levelsNucleotideGraph.size() == 1));
			if(e2->levelsNucleotideGraph.size() > 0)
			{
				string kMer2LastCharacter = kMer2.substr(kMer2.size() - 1, 1);
				int correspondingPosition2 = e2->levelsNucleotideGraph.at(0);
				assert(correspondingPosition2 >= 0);
				assert(correspondingPosition2 < forReturn.h2.size());
				forReturn.h2.at(correspondingPosition2) = kMer2LastCharacter;
			}
		}
	}

	return forReturn;
}

diploidGenomeString LargeGraph::diploidPathToGenomeString(diploidEdgePath& path)
{
	diploidEdgePointerPath pathEdgePointers = edgeIndexPathToEdgePointerPath(path);
	return diploidPathToGenomeString(path);
}

diploidGenomeString LargeGraph::diploidPathToGenomeString(diploidEdgePointerPath& path)
{
	diploidNucleotidePath diploidNucleotidePath;

	assert(path.h1.size() == path.h2.size());
	std::vector<bool> phaseIsLost;
	phaseIsLost.resize(path.h1.size()+kMerSize, false);

	diploidNucleotidePath.h1.resize(path.h1.size()+kMerSize, "_");
	diploidNucleotidePath.h2.resize(path.h1.size()+kMerSize, "_");

	for(int pos = 0; pos < (int)path.h1.size(); pos++)
	{
		Edge* e1 = path.h1.at(pos);
		string locusID1 = e1->locus_id;
		string kMer1 = CODE.deCode(locusID1, e1->largeEmission2);

		Edge* e2 = path.h2.at(pos);
		string locusID2 = e2->locus_id;
		string kMer2 = CODE.deCode(locusID2, e2->largeEmission2);

		if(pos == 0)
		{
			assert(e1->levelsNucleotideGraph.size() == kMer1.size());
			for(int cI = 0; cI < kMer1.size(); cI++)
			{
				int correspondingPosition = e1->levelsNucleotideGraph.at(cI);
				diploidNucleotidePath.h1.at(correspondingPosition) = kMer1.substr(cI, 1);
				if(e1 == e2)
				{
					phaseIsLost.at(correspondingPosition) = true;
				}
			}

			assert(e2->levelsNucleotideGraph.size() == kMer2.size());
			for(int cI = 0; cI < kMer2.size(); cI++)
			{
				int correspondingPosition = e2->levelsNucleotideGraph.at(cI);
				diploidNucleotidePath.h2.at(correspondingPosition) = kMer2.substr(cI, 1);
				if(e1 == e2)
				{
					phaseIsLost.at(correspondingPosition) = true;
				}
			}
 		}
		else
		{
			assert((e1->levelsNucleotideGraph.size() == 0) || (e1->levelsNucleotideGraph.size() == 1));
			if(e1->levelsNucleotideGraph.size() > 0)
			{
				string kMer1LastCharacter = kMer1.substr(kMer1.size() - 1, 1);
				int correspondingPosition1 = e1->levelsNucleotideGraph.at(0);
				assert(correspondingPosition1 >= 0);
				assert(correspondingPosition1 < diploidNucleotidePath.h1.size());
				diploidNucleotidePath.h1.at(correspondingPosition1) = kMer1LastCharacter;
				if(e1 == e2)
				{
					phaseIsLost.at(correspondingPosition1) = true;
				}
			}

			assert((e2->levelsNucleotideGraph.size() == 0) || (e2->levelsNucleotideGraph.size() == 1));
			if(e2->levelsNucleotideGraph.size() > 0)
			{
				string kMer2LastCharacter = kMer2.substr(kMer2.size() - 1, 1);
				int correspondingPosition2 = e2->levelsNucleotideGraph.at(0);
				assert(correspondingPosition2 >= 0);
				assert(correspondingPosition2 < diploidNucleotidePath.h2.size());
				diploidNucleotidePath.h2.at(correspondingPosition2) = kMer2LastCharacter;
				if(e1 == e2)
				{
					phaseIsLost.at(correspondingPosition2) = true;
				}
			}
		}
	}

	diploidGenomeString forReturn;

	std::pair<std::vector<std::string>, std::vector<std::string> > runningExtractedCompartment;

	auto clearRunningCompartment = [&]() {
		assert(runningExtractedCompartment.first.size() == runningExtractedCompartment.second.size());

		if(runningExtractedCompartment.first.size() > 0)
		{
			std::string h1 = Utilities::join(runningExtractedCompartment.first, "");
			std::string h2 = Utilities::join(runningExtractedCompartment.second, "");

			std::vector<std::string> thisCompartmentParts;
			thisCompartmentParts.push_back(h1);

			if(h1 != h2)
				thisCompartmentParts.push_back(h2);

			forReturn.push_back(thisCompartmentParts);
		}

		runningExtractedCompartment.first.clear();
		runningExtractedCompartment.second.clear();
	};

	for(unsigned int pI = 0; pI < diploidNucleotidePath.h1.size(); pI++)
	{
		if(phaseIsLost.at(pI))
		{
			clearRunningCompartment();
		}

		runningExtractedCompartment.first.push_back(diploidNucleotidePath.h1.at(pI));
		runningExtractedCompartment.second.push_back(diploidNucleotidePath.h2.at(pI));
	}
	clearRunningCompartment();

	size_t accumulatedCompartmentsSize = 0;
	for(unsigned int i = 0; i < forReturn.size(); i++)
	{
		assert(forReturn.at(i).size() > 0);
		if(forReturn.at(i).size() > 1)
		{
			assert(forReturn.at(i).size() == 2);
			assert(forReturn.at(i).at(0).size() == forReturn.at(i).at(1).size());
		}
		accumulatedCompartmentsSize += forReturn.at(i).at(0).size();
	}

	if(!(diploidNucleotidePath.h1.size() == accumulatedCompartmentsSize))
	{
		std::cerr << "!(diploidNucleotidePath.h1.size() == accumulatedCompartmentsSize);\n";
		std::cerr << "diploidNucleotidePath.h1.size(): " << diploidNucleotidePath.h1.size() << "\n";
		std::cerr << "accumulatedCompartmentsSize: " << accumulatedCompartmentsSize << "\n";
		std::cerr << "forReturn.size(): " << forReturn.size() << "\n";
		std::cerr << std::flush;
	}
	assert(diploidNucleotidePath.h1.size() == accumulatedCompartmentsSize);

	return forReturn;
}



diploidNucleotidePath LargeGraph::diploidPathToAlignedNucleotides(diploidEdgePath& path)
{
	diploidEdgePointerPath pathEdgePointers = edgeIndexPathToEdgePointerPath(path);
	return diploidPathToAlignedNucleotides(path);
}


diploidNucleotidePath LargeGraph::diploidPathToNucleotides(diploidEdgePath& path)
{
	diploidNucleotidePath forReturn;

	assert(path.h1.size() == path.h2.size());

	for(int pos = 0; pos < (int)path.h1.size(); pos++)
	{
		int edgeIndex = path.h1.at(pos);
		assert(idx2Edge.count(edgeIndex) > 0);
		Edge* e = idx2Edge[edgeIndex];
		string locusID = e->locus_id;
		string kMer = CODE.deCode(locusID, e->largeEmission2);
		string kMerFirst = kMer.substr(0, 1);
		forReturn.h1.push_back(kMerFirst);

		edgeIndex = path.h2.at(pos);
		assert(idx2Edge.count(edgeIndex) > 0);
		e = idx2Edge[edgeIndex];
		locusID = e->locus_id;
		kMer = CODE.deCode(locusID, e->largeEmission2);
		kMerFirst = kMer.substr(0, 1);
		forReturn.h2.push_back(kMerFirst);
	}

	assert(forReturn.h1.size() == path.h1.size());
	assert(forReturn.h2.size() == forReturn.h1.size());

	return forReturn;
}

vector<string> LargeGraph::haploidPathToNucleotides(vector<int> path)
{
	vector<string> forReturn;

	for(int pos = 0; pos < (int)path.size(); pos++)
	{
		int edgeIndex = path.at(pos);
		assert(idx2Edge.count(edgeIndex) > 0);
		Edge* e = idx2Edge[edgeIndex];
		string locusID = e->locus_id;
		string kMer = CODE.deCode(locusID, e->largeEmission2);
		string kMerFirst = kMer.substr(0, 1);
		forReturn.push_back(kMerFirst);
	}

	assert(forReturn.size() == path.size());

	return forReturn;
}


readSimulationResults LargeGraph::simulateReadsForPath(diploidEdgePath& path)
{
	readSimulationResults results;

	map<string, long long> readCount;
	map<string, long long> kMerOccurence;
	map<string, set<int> > kMerOccurenceLocalization;
	vector< vector<int> > haplotypeStartPositions;

	vector<string> reqKMers = requiredKMers();
	for(int i = 0; i < (int)reqKMers.size(); i++)
	{
		string kMer = reqKMers.at(i);
		readCount[kMer] = 0;
		kMerOccurence[kMer] = 0;
	}

	int haploidCoverage = 15;
	int readLength = 85;

	readCount["TotalSeq"] = 3200000000.0 * 2 * 15;
	readCount["MeanReadLen"] = (double)readLength;
	readCount["TotalKMerCoverage"] = 0;

	double onePositionStartReadProbability = (double) haploidCoverage / (double)readLength;

    boost::mt19937 rnd_gen;   //Mersenne Twister generator
    typedef boost::variate_generator< boost::mt19937, boost::poisson_distribution<> > rnd_poisson_t;
	rnd_poisson_t rnd_poisson( rnd_gen, boost::poisson_distribution<>( onePositionStartReadProbability ));
    rnd_poisson = rnd_poisson_t( rnd_gen, boost::poisson_distribution<>( onePositionStartReadProbability ));

    assert(kMerSize > 0);

	for(int haplotype = 0; haplotype < 2; haplotype++)
	{
		vector<int> kMerHaplotype;
		vector<int> haplotypeStart;

		if(haplotype == 0)
		{
			kMerHaplotype = path.h1;
		}
		else if(haplotype == 1)
		{
			kMerHaplotype = path.h2;
		}
		else
		{
			assert( 1 == 0 );
		}

		string haplotype_string = Utilities::join(haploidPathToNucleotides(kMerHaplotype), "");
		results.haplotypes_untrimmed.push_back(haplotype_string);

		for(int pI = 0; pI < (int)haplotype_string.size(); pI++)
		{
			if(haplotype_string.substr(pI, 1) != "_")
			{
				haplotypeStart.push_back(pI);
			}
		}
		erase_all(haplotype_string, "_");

		results.haplotypes.push_back(haplotype_string);
		results.haplotypeStartPositions.push_back(haplotypeStart);

		for(int pos = 0; pos <= (int)(haplotype_string.length() - kMerSize); pos++)
		{
			string kMer = haplotype_string.substr(pos, kMerSize);
			if(kMer.find("*") == string::npos)
			{
				assert((int)kMer.size() == kMerSize);
				assert(kMerOccurence.count(kMer) > 0);
				kMerOccurence[kMer]++;
				kMerOccurenceLocalization[kMer].insert(pos);
			}
		}

		for(int pos = 0; pos < (int)haplotype_string.length(); pos++)
		{
			int starting_reads = rnd_poisson();
			if(starting_reads > 0)
			{
				int nCharacters_Read = readLength;
				if((pos + nCharacters_Read) > haplotype_string.size())
				{
					nCharacters_Read = haplotype_string.size() - pos;
					assert(nCharacters_Read >= 1);
				}
				string proposedRead = haplotype_string.substr(pos, nCharacters_Read);
				assert((int)proposedRead.length() == nCharacters_Read);

				if(((int)proposedRead.length() - kMerSize) >= 0)
				{
					for(int kMerI = 0; kMerI <= ((int)proposedRead.length() - kMerSize); kMerI++)
					{
						if(!((int)proposedRead.length() >= kMerSize))
						{
							cout << "proposedRead.length(): " << proposedRead.length() << "\n";
							cout << "proposedRead.length() - kMerSize: " << (int)proposedRead.length() - kMerSize << "\n";
							cout << "kMerSize: " << kMerSize << "\n";
							cout << "nCharacters_Read: " << nCharacters_Read << "\n";


						}
						assert((int)proposedRead.length() >= kMerSize);

						string proposedkMer = proposedRead.substr(kMerI, kMerSize);
						assert((int)proposedkMer.length() == kMerSize);

						if(proposedkMer.find("*") == string::npos)
						{
							assert(readCount.count(proposedkMer) > 0);
							readCount[proposedkMer] += starting_reads;
						}
					}
				}
			}
		}
	}

	results.kMerOccurenceLocalization = kMerOccurenceLocalization;
	results.kMerOccurence = kMerOccurence;
	assert(results.haplotypes.size() == 2);
	results.simulatedReads = readCount;
	return results;
}

vector< vector<string> > LargeGraph::categorizeEdgeLevels()
{
	vector< vector<string> > forReturn;
	int levels = NodesPerLevel.size();
	for(int level = 0; level < (levels - 1); level++)
	{
		set<string> seenKMers;
		bool isINDEL = false;
		bool haveStar = false;
		set<string> seenNucleotides;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* node = *nodeIt;


			for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				assert(Edges.count(e) > 0);
				string kMer = CODE.deCode(e->locus_id, e->largeEmission);
				if(kMer == "_")
				{
					isINDEL = true;
				}
				else if(kMer == "*")
				{
					haveStar = true;
				}
				else
				{
					seenKMers.insert(kMer);
					seenNucleotides.insert(kMer.substr(0, 1));

				}
			}
		}

		vector<string> categories;
		if(isINDEL == true)
		{
			categories.push_back("INDEL");
		}
		if((! haveStar) && ((seenNucleotides.size() > 1) || ((seenNucleotides.size() > 0) && (isINDEL == true))))
		{
			categories.push_back("PolyNoStar");
		}
		if(seenNucleotides.size() > 1)
		{
			categories.push_back("SNP");
			if((! isINDEL) && (! haveStar))
			{
				categories.push_back("CleanSNP");
			}
			if(! isINDEL)
			{
				categories.push_back("SNPnoINDEL");
			}
		}

		forReturn.push_back(categories);
	}

	return forReturn;
}

diploidEdgePath LargeGraph::simulateRandomPath()
{

	vector< vector<int> > edgePaths;
	for(int pI = 0; pI < 2; pI++)
	{
		vector<int> currentEdgePath;

		Node* currentNode;
		Node* n0 = *(NodesPerLevel.at(0).begin());

		currentNode = n0;
		while(currentNode->Outgoing_Edges.size() != 0)
		{
			int n_edges = currentNode->Outgoing_Edges.size();
			vector<Edge*> currentEdges (currentNode->Outgoing_Edges.begin(), currentNode->Outgoing_Edges.end());

			double f = (double)rand() / RAND_MAX;

			f = f * n_edges;
			int selected_edge = (int) f;
			if(selected_edge == n_edges)
			{
				selected_edge = n_edges - 1;
			}

			assert(selected_edge >= 0);
			assert(selected_edge < n_edges);

			Edge* selectedEdge = currentEdges.at(selected_edge);

			assert(EdgeCounter.count(selectedEdge) > 0);
			int selectedEdgeIndex = EdgeCounter[selectedEdge];

			currentEdgePath.push_back(selectedEdgeIndex);
			currentNode = selectedEdge->To;
		}

		assert(currentEdgePath.size() == (NodesPerLevel.size() - 1));

		edgePaths.push_back(currentEdgePath);

	}

	diploidEdgePath forReturn;
	forReturn.h1 = edgePaths.at(0);
	forReturn.h2 = edgePaths.at(1);

	return forReturn;
}


vector<string> LargeGraph::requiredKMers()
{
	int levels = NodesPerLevel.size();
	set<string> requiredKMers;
	for(int level = 0; level < (levels-1); level++)
	{
		vector<Node*> nodes(NodesPerLevel.at(level).begin(), NodesPerLevel.at(level).end());
		for(unsigned int nI = 0; nI < NodesPerLevel.at(level).size(); nI++)
		{
			Node* n = nodes.at(nI);
			for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
			{
				Edge* e = (*outgoingIt);
				string kMer = CODE.deCode(e->locus_id, e->largeEmission);
				if((kMer != "_") && (kMer != "*"))
				{
					requiredKMers.insert(kMer);
				}
			}
		}
	}

	vector<string> forReturn(requiredKMers.begin(), requiredKMers.end());

	return forReturn;
}


