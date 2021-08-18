#include "query_graph.h"
#include <fstream>
#include <iostream>


Query :: Query(const string iNodeFile, const string iEdgeFile)
{
	graph.clear();

	string id, label, edge;  // Elements to be read from the files
	id = label = edge = "";

	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
	{
		cout<<"Unable to open input [ID -> label] file for query"<<endl;
		exit(1);
	}
	

	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		Vertex *node = new Vertex(id, label);
		graph[id] = node;

		iNF >> id >> label;
	}

	iNF.close();


	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	if(!iEF.is_open())
	{
		cout<<"Unable to open input [ID -> neighbour] file for query"<<endl;
		exit(1);
	}

	iEF >> id >> edge;  // Reads a new node ID and its label
	while(!iEF.eof())  // Read the neighbours for vertices
	{
		Vertex *node = graph.at(id);
		node->addNeighbours(edge);

		// Add the reverse edges for undirected graph
		node = graph.at(edge);
		node->addNeighbours(id);

		iEF >> id >> edge;  // Reads a new node ID and its label
	}

	iEF.close();
}

void Query :: createIndex(const unordered_map<string, unsigned int>& uniqLabel)
{
    // Construct the neighbouring labels for the vertices in the input graph
    unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
        (it->second)->createNeighbourLabels(graph);

	// Construct the neighbour label vector for the vertices
	it = graph.begin();
    for(; it!=graph.end(); it++)
    	(it->second)->createNeighLabelVec(uniqLabel);
}


Query :: ~Query()
{
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();
}

const unsigned Query :: getGraphSize(void) const
{
	return graph.size();
}


const set<string> Query :: getVertexIDs(void) const
{
        set<string> vertex;
        unordered_map<string, Vertex*>::const_iterator it = graph.begin();

        for(; it!=graph.end(); it++)
                vertex.insert(it->first);

        return vertex;
}

const unordered_map<string, Vertex*>& Query :: getGraph(void) const
{
	return graph;
}


void Query :: printVertex(void) const
{
	cout<<"Query graph contains the following vertices with the associated labels:"<<endl;
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();

	for(; it!=graph.end(); it++)
		cout<<it->first<<" ("<<(it->second)->getLabel()<<"), ";

	cout<<endl;
}


void Query :: printGraph(void) const
{
	cout<<"Query graph contains "<<graph.size()<<" vertices."<<endl;

	
	printVertex();
	cout<<"The query graph contains the following vertices:"<<endl;

	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();
	
}

