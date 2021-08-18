/*** Defines the input graph as a collection of vertices (of vertex.h). Also computes 
 *** the output approximate sub-graph matching to the query provided as input.
***/


#ifndef INPUT_GRAPH_H_
 #define INPUT_GRAPH_H_


#include "vertex.h"
#include "query_graph.h"
#include <queue>


using namespace std;


class Triple
{
	public:
	   Vertex* v;
	   Vertex* q;
	   double simi;

		Triple()
		{
			v = q = NULL;
			simi = 0.0;
		}

		Triple(Vertex* ver, Vertex* qry, double sim)
		{
			v = ver;
			q = qry;
			simi = sim;
		}

		~Triple()
		{
			v = q = NULL;
			simi = 0.0;
		}

		void setVal(Vertex* ver, Vertex* qry, double sim)
		{
			v = ver;
			q = qry;
			simi = sim;
		}
};


/*

// Compare operator for min-heap operation
class Compare_min
{
   public:
	bool operator() (Triple a, Triple b)
	{
		double a_val = a.simi;
		double b_val = b.simi;

		return (a_val > b_val);
	}
};
*/

class Compare_max
{
   public:
	bool operator() (Triple a, Triple b)
	{
		double a_val = a.simi;
		double b_val = b.simi;

		return (a_val < b_val);
	}
};

class Input
{
   private:
	unsigned long numVert;  // Stores number of vertices in the graph
	unsigned long numLabel;  // Stores the number of unique labels in the graph

	unordered_map<string, Vertex*> graph;  // Stores the input graph as a map of vertex ID to the corresponding object
	unordered_map<string, vector<Vertex*> > vertexLabel;  // Stores the IDs of vertices having the same labels
	unordered_map<string, unsigned int> uniqLabel;  // Stores the index map for the unique labels

	double avgSim, stdSim;  // Computed average similarity and standard deviation between neighbour label set of vertex pairs
	unsigned long numSym;  // Computes the number of match symbols that can be generated
	double *symProb;
	

	//priority_queue<Vertex*, vector<Vertex*>, Compare_min> heap;  // Min-heap structure to keep the vertices with top-k chi-sq values (for greedy approach)

	void computeSymbolProb(double);  // Computes the number of symbols and their occurrence probabilites

   public:
	// Creates the input graph from two input files - (1) mapping of node ID to label, (2) list of neighbour labels
	Input(const string, const string);

	~Input();  // Deallocates the constructed graph

	void computeGraphCharacteristics(string = "");  // Computes the similarity characteristics and symbols of the graph

	// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
	vector<vector<Triple> > getMatchingSubgraphs(const Query&) const;

	// Returns the number of vertices in the graph
	const unsigned long getGraphSize(void) const;

	// Returns the index map of unique labels in the graph
	const unordered_map<string, unsigned int>& getUniqLabel(void) const;

	// Returns the vertex label index
	const unordered_map<string, vector<Vertex*> >& getVertexLabel(void) const;

	// Returns the symbol probabilities of the graph
	const double* getSymProb(void) const;

	// Computes mapping between neighbouring vertices after primary vertex match between input and query graph
	vector<Triple> findNeighbourMapping(const Vertex*, const Vertex*, const unordered_map<string, Vertex*>&) const;
	vector<Triple> findNeighbourMapping(const Vertex*, const Vertex*, const unordered_map<string, Vertex*>&, const set<string>&) const;

	// Computes the chi-square of a vertex matched to that of a query
	double computeChiSquare(const Vertex*, const Vertex*, double, const vector<Triple>&, const Query&) const;

	void printGraph(void) const;
};


#endif

