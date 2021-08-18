/*** Defines the basic vertex structure of a graph comprising an unique ID, corresponding label,
 *** and a set of neighbouring vertex labels. These vertices define the Query graph, and will be
 *** extended for representing the input graph vertices.
***/


#ifndef VERTEX_H_
 #define VERTEX_H_

#include "const.h"
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

using namespace std;


class Vertex
{
   protected:
	const string ID;  // Stores the unique ID of the vertex
	const string label;  // Stores the associated label with the vertex
	unsigned degree;  // Stores the degree of the vertex

	set<string> neighbourIDs;  // Stores the ID of the adjacent vertices
	vector<string> neighbourLabels;  // Stores the Lables of the adjacent vertices
	vector<unsigned int> neigLabelVec;  // Stores a vector map of the number of label occurrences in the neighbourhood


   public:
	// Constructs a vertex with the corresponding characteristics
	Vertex(const string, const string);  // Sets the ID and the label of a vertex

	~Vertex();  // Deallocates space

	const string getID(void) const;  // Returns vertex ID
	const string getLabel(void) const;  // Returns vertex label
	const unsigned getDegree(void) const;  // Returns vertex degree

	const set<string>& getNeighbourIDs(void) const;  // Returns the ID list of the neighbouring vertices
	const vector<string>& getNeighbourLabels(void) const;  // Returns the label list of the neighbouring vertices

	void addNeighbours(string);  // Add the ID of an adjacent vertex, one at a time
	void createNeighbourLabels(const unordered_map<string, Vertex*>&);  // Populates the neighbourLabels index to find the triplets
	void createNeighLabelVec(const unordered_map<string, unsigned int>&);  // Sorts the labels in the neighbourhood

	friend double findVertexSimi(const Vertex*, const Vertex*);  // Computes the neighbourhood label similarity between two vertices
	friend double findVertexSimi(const Vertex*, const Vertex*, double&);  // Computes bi-directional neighbourhood label similarity between two vertices
	void print(void) const;  // Prints the vertex characteristics
};

#endif
