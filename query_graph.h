/*** Defines the query graph as a collection of vertices. Also keeps an 
 *** index mapping the vertex labels to neighbouring labels.
***/


#ifndef QUERY_GRAPH_H_
 #define QUERY_GRAPH_H_


#include "vertex.h"
#include <vector>


class Query
{
   private:
	unordered_map<string, Vertex*> graph;  // Stores the query graph as a map of vertex ID to the corresponding object

	

   public:
	// Creates the query graph from two files - (1) mapping of node ID to label, (2) list of neighbour labels
	Query(const string, const string);

	~Query();  // Deallocates the graph


	const unsigned getGraphSize(void) const;  // Returns the number of nodes in the query graph
	const set<string> getVertexIDs(void) const;  // Returns the vertex IDs in the query graph
	const unordered_map<string, Vertex*>& getGraph(void) const;  // Returns the query graph
	void createIndex(const unordered_map<string, unsigned int>&);  // Creates the indices of the vertices of the graph

	void printVertex(void) const;  // Prints the vertex IDs and corresponding labels of the query graph
	void printGraph(void) const;
};


#endif
