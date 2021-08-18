/*** Implements the input graph and its functionalities
***/


#include "input_graph.h"
#include "query_graph.h"
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>


Input :: Input(const string iNodeFile, const string iEdgeFile)
{
	string id, label, edge;  // Elements to be read from the files
	id = label = edge = "";

	numVert = 0UL;
	numLabel = 0UL;
	numSym = 0U;
	symProb = NULL;
	avgSim = stdSim = 0.0;

	graph.clear();
	vertexLabel.clear();
	uniqLabel.clear();

	clock_t begin, end;

	begin = clock();
	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
	{
		cout<<"Unable to open input [ID -> label] file for input"<<endl;
		exit(1);
	}


	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		Vertex *ver = new Vertex(id, label);
		graph[id] = ver;

		if(vertexLabel.find(label) == vertexLabel.end())
		{
			vector<Vertex*> v;
			v.push_back(ver);
			vertexLabel[label] = v;
		}
		else
			(vertexLabel.at(label)).push_back(ver);


		if(uniqLabel.find(label) == uniqLabel.end())
		{
			uniqLabel[label] = numLabel;
			numLabel ++;
		}

		iNF >> id >> label;
	}

	iNF.close();

	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	if(!iEF.is_open())
	{
		cout<<"Unable to open input [ID -> neighbour] file for input"<<endl;
		exit(1);
	}

	iEF >> id >> edge;  // Reads a new node ID and its neighbour
	while(!iEF.eof())  // Read the neighbours for vertices
	{
		Vertex *ver = graph.at(id);
		ver->addNeighbours(edge);

		// Add the reverse edges for undirected graph
		ver = graph.at(edge);
		ver->addNeighbours(id);

		iEF >> id >> edge;  // Reads a new node ID and its label
	}

	iEF.close();

	end = clock();
	cout<<"Reading time of Input Graph = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;  // Input Graph Read Time

	// Construct the neighbouring labels for the vertices in the input graph
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();

	for(; it!=graph.end(); it++)
		(it->second)->createNeighbourLabels(graph);

	numVert = graph.size();
}


Input :: ~Input()
{
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();
	vertexLabel.clear();
	uniqLabel.clear();
	
	numVert = 0;
	numSym = 0;
	numLabel = 0;
	avgSim = stdSim = 0.0;

	delete symProb;
	symProb = NULL;
}


const unsigned long Input :: getGraphSize(void) const
{
	return numVert;
}

const unordered_map<string, unsigned int>& Input :: getUniqLabel(void) const
{
	return uniqLabel;
}

const double* Input :: getSymProb(void) const
{
	return symProb;
}

const unordered_map<string, vector<Vertex*> >& Input :: getVertexLabel(void) const
{
	return vertexLabel;
}

void Input :: printGraph(void) const
{
	cout<<"Input graph contains "<<numVert<<" vertices with "<<numLabel<<" unique labels."<<endl;
	cout<<"The number of symbols computed are: "<<numSym<<endl;


/*
	cout<<"Vertex label index constructed: ";
	unordered_map<string, unsigned int>::const_iterator i = uniqLabel.begin();
	for(; i!=uniqLabel.end(); i++)
		cout<<i->first<<":"<<i->second<<" , ";

	cout<<endl;

	cout<<"The vertices present are as:"<<endl;

	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();

	cout<<endl<<"The vertices associated with the labels are:"<<endl;
	unordered_map<string, vector<Vertex*> >::const_iterator it1 = vertexLabel.begin();

	for(; it1!=vertexLabel.end(); it1++)
	{
		cout<<it1->first<<": ";
		const vector<Vertex*> vert = it1->second;

		for(unsigned j=0; j<vert.size(); j++)
			cout<<(vert.at(j))->getID()<<" ";

		cout<<endl;
	}

	cout<<endl<<"The number of symbols computed are: "<<numSym<<" and the associated probabilities are:"<<endl;
	for(unsigned long i=0; i<numSym; i++)
		cout<<"P[s_"<<i<<"] = "<<symProb[i]<<endl;
*/
}

void Input :: computeGraphCharacteristics(string iCharFile)		// Default value of iCharFile=""
{
	// Create neighbour label vector for each of the input graoh vertices
	unordered_map<string, Vertex*>::const_iterator it1 = graph.begin();
	for(; it1!=graph.end(); it1++)
		(it1->second)->createNeighLabelVec(uniqLabel);
	
	if(iCharFile=="")
	{
		// Information not pre-computed, need to compute
		
		bool one_pass = true;

		if(one_pass)
		{
			it1 = graph.begin();
			unsigned long long numSimi = 0ULL;
			long double avg = 0.0L;
			long double M2 = 0.0L;
			double d, d2;

			unsigned long count = 0UL;
			cout<<"Computing Input Graph Characteristics...";
			
			for(; it1!=graph.end(); it1++)
			{
				unordered_map<string, Vertex*>::const_iterator it2 = it1;
				numSimi ++;
				d = 1.0 - avg;
				avg += (d / numSimi);
				d2 = 1.0 - avg;
				M2 += (d * d2);
				count ++;
				it2++;
				for(; it2!=graph.end(); it2++)
				{
					double rev_simi;
					double simi = findVertexSimi(it1->second, it2->second, rev_simi);
					numSimi ++; 
					d = simi - avg;
					avg += (d / numSimi);
					d2 = simi - avg;
					M2 += (d * d2);

					numSimi ++; 
					d = rev_simi - avg;
					avg += (d / numSimi);
					d2 = rev_simi - avg;
					M2 += (d * d2);
				}
				if(count%400 == 0)
				{
					cout<<count<<"...";
					fflush(stdout);
				}
			}
			cout<<endl;

			M2 /= (numSimi-1.0);  // Bessel Correction performed
			M2 = sqrt(M2);

			avgSim = avg;
			stdSim = M2;
		}
		else
		{
			long double avg = 0.0L;
			unsigned long long numSimi = 0ULL;
			unsigned long count = 0UL;

			it1 = graph.begin();
			cout<<"Computing Input Graph Average Characteristics...";
			for(; it1!=graph.end(); it1++)
			{
				unordered_map<string, Vertex*>::const_iterator it2 = it1;
				avg += 1.0;
				numSimi ++;
				it2++;
				count ++;
				for(; it2!=graph.end(); it2++)
				{
					double rev_simi;
					double simi = findVertexSimi(it1->second, it2->second, rev_simi);
					avg += simi;
					avg += rev_simi;
					numSimi += 2;
				}
				if(count%400 == 0)
				{
					cout<<count<<"...";
					fflush(stdout);
				}
			}
			cout<<endl;

			avg /= numSimi;


			long double std = 0.0L;
			count = 0UL;
			it1 = graph.begin();
			cout<<"Computing Input Graph Standard Deviation Characteristics...";
			for(; it1!=graph.end(); it1++)
			{
				unordered_map<string, Vertex*>::const_iterator it2 = it1;
				double dev = 1.0 - avg;
				std += pow(dev, 2.0);
				it2++;
				count ++;
				for(; it2!=graph.end(); it2++)
				{
					double rev_simi;
					double simi = findVertexSimi(it1->second, it2->second, rev_simi);
					dev = simi - avg;
					std += pow(dev, 2.0);
					dev = rev_simi - avg;
					std += pow(dev, 2.0);
				}
				if(count%400 == 0)
				{
					cout<<count<<"...";
					fflush(stdout);
				}
			}
			cout<<endl;


			std /= (numSimi-1.0);  // Bessel Correction performed
			std = sqrt(std);

			avgSim = avg;
			stdSim = std;
		}
	}
	else
	{
		// Information already available
		ifstream iCF;
		iCF.open(iCharFile.c_str());  // Contains previously average and std deviation in that order

		if(!iCF.is_open())
		{
			cout<<"Unable to open "<<iCharFile<<" file for input"<<endl;
			exit(1);
		}

		iCF >> avgSim >> stdSim;
		iCF.close();

		cout<<"Read Input Graph Characteristics from "<<iCharFile<<" as, avgSim = "<<avgSim<<" and stdSim = "<<stdSim<<endl;
	}

	if(stdSim == 0.0)
	{
		cout<<"Standard Deviation found to be 0."<<endl;
		stdSim = 0.00001;
	}

	double max_dev = max(abs(1.0 - avgSim), abs(avgSim - 0)) / stdSim;	

	cerr<<avgSim<<"\t"<<stdSim<<endl;
	
	computeSymbolProb(max_dev);
}

void Input :: computeSymbolProb(double max_dev)
{
	if(max_dev <= 1.0)
	{
		cout<<"Input graph is either a clique or too dense. Incompatible with current algorithm. Aborting!"<<endl;
		exit(1);
	}

	double max_dev_pos = abs(1.0 - avgSim) / stdSim;
	double max_dev_neg = abs(avgSim - 0.0) / stdSim;
	// cout<<max_dev_pos<<" "<<max_dev_neg<<endl;

	long numSym_pos = static_cast<int>(ceil((max_dev_pos-1) / STEP_SIZE));
	long numSym_neg = static_cast<int>(ceil((max_dev_neg-1) / STEP_SIZE));
	// cout<<numSym_pos<<" "<<numSym_neg<<endl;

	numSym = std::max(numSym_pos, numSym_neg); // include <algorithm>

	symProb = new double[numSym];
	double total_prob = 0.0;
	for(unsigned int i=1; i<numSym; i++)
	{
		double a = 1 + i*STEP_SIZE;
		double b = 1 + (i+1)*STEP_SIZE;

		double factor = 0;
		if ((avgSim + a*stdSim) <= 1.0)
			factor = factor + 0.5;
		if ((avgSim - a*stdSim) >= 0.0)
			factor = factor + 0.5;

		symProb[i] = (1.0 / pow(a,2.0)) - (1.0 / pow(b,2.0));
		symProb[i] *= factor;
		total_prob += symProb[i];
	}
	symProb[0] = 1.0 - total_prob;
}

vector<vector<Triple> > Input :: getMatchingSubgraphs(const Query& qry) const
{
	// Generate vertex pairs
	const unordered_map<string, Vertex*> qryGraph = qry.getGraph();
	unordered_map<string, Vertex*>::const_iterator it = qryGraph.begin();
	// Primary Heap - <input vertex, query vertex, chiSq value>
	priority_queue<Triple, vector<Triple>, Compare_max> prim_heap;
	// Create primary heap between query and input graph
	for(; it!=qryGraph.end(); it++)
	{
		Vertex* q = it->second;
		string qryLabel = q->getLabel();
		const vector<Vertex*> candVert = vertexLabel.at(qryLabel);
		
		for(unsigned long i=0; i<candVert.size(); i++)
		{
			Vertex* v = candVert.at(i);
			double cand_simi = findVertexSimi(v, q);
		
			// Compute the mapping between the neighbours of the candidate vertex and the query vertex
			// Computes <input neighbour vertex, query neighbour vertex, similarity>
			const vector<Triple> neig_pair = findNeighbourMapping(v, q, qryGraph);
			// Compute chi-sq of this pair and push to heap
			double chiVal = computeChiSquare(v, q, cand_simi, neig_pair, qry);
			
			// Push to primary heap
			prim_heap.push(Triple(v, q, chiVal));		// TODO: CHANGE chiVal*cand_simi
		}
	}

	vector< vector<Triple> > answer;

	for(unsigned short i=0; i<TOPK; i++)
	{
		if(prim_heap.empty())
			break;

		vector<Triple> subgraph;
		Triple tr = prim_heap.top();
		prim_heap.pop();
		subgraph.push_back(tr);  // Stores the vertex mapping and the computed chi square of the answer

		set<string> qry_duplicate, ver_duplicate;  // Stores vertices already mapped to remove duplicates
		unsigned long num_qry_ver = qryGraph.size();

		priority_queue<Triple, vector<Triple>, Compare_max> scnd_heap;

		do{
			Vertex* v = tr.v;
			Vertex* q = tr.q;
			qry_duplicate.insert(q->getID());
			ver_duplicate.insert(v->getID());

			const vector<Triple> neig_pair = findNeighbourMapping(v, q, qryGraph, ver_duplicate);
			for(unsigned long j=0; j<neig_pair.size(); j++)
			{
				tr = neig_pair.at(j);
				v = tr.v;
				q = tr.q;
				if( (qry_duplicate.find(q->getID()) != qry_duplicate.end()) || (ver_duplicate.find(v->getID()) != ver_duplicate.end()) )
					continue;

				double cand_simi = findVertexSimi(v, q);
				const vector<Triple> n_p = findNeighbourMapping(v, q, qryGraph); 
				double chiVal = computeChiSquare(v, q, cand_simi, n_p, qry);

				scnd_heap.push(Triple(v, q, chiVal));		// TODO: CHANGE chiVal*cand_simi
				qry_duplicate.insert(q->getID());
				ver_duplicate.insert(v->getID());
			}

			if(!scnd_heap.empty())
			{
				tr = scnd_heap.top();
				scnd_heap.pop();
				subgraph.push_back(tr);
			}
		}while( (!scnd_heap.empty()) && (subgraph.size() < num_qry_ver) );

		answer.push_back(subgraph);
	}

	return answer;
}


double Input :: computeChiSquare(const Vertex* v, const Vertex* q, double cand_simi, const vector<Triple>& neig_pair, const Query& qry) const
{
	const set<string> q_neig = q->getNeighbourIDs();
	const unsigned long seq_len = q_neig.size() + 1;

	double* exp_symb_count = new double[numSym];
	for(unsigned long i=0; i<numSym; i++)
		exp_symb_count[i] = symProb[i] * seq_len;

	double* occ_symb_count = new double[numSym];
	for(unsigned long i=0; i<numSym; i++)
		occ_symb_count[i] = 0.0;

	unsigned long sym;

	for(unsigned long i=0; i<neig_pair.size(); i++)
	{
		Triple tr = neig_pair.at(i);
		const double simi = tr.simi;

		// Compute neighbour symbol
		const double k = (simi - avgSim) / stdSim;			// TODO: TAKE ABSOLUTE OF NUMERATOR FOR NO SHORTING
		if( (k-1) <= 0)
			sym = 0;
		else
			sym = ceil((k - 1.0) / STEP_SIZE) - 1;
		occ_symb_count[sym] ++;
	}
	
	// Compute symbol
	const double k = (cand_simi - avgSim) / stdSim;			// TODO: CHANGE ABSOLUTE OF NUMERATOR FOR NO SHORTING
	if( (k-1) <= 0)
		sym = 0;
	else
		sym = ceil((k - 1.0) / STEP_SIZE) - 1;
	occ_symb_count[sym] ++;

	// Compute Chi Square value
	double chiSq = 0.0;
	for(unsigned long i=0; i<numSym; i++)
	{
		double numer = pow((occ_symb_count[i] - exp_symb_count[i]), 2.0);
		double denom = exp_symb_count[i];
		chiSq += (numer / denom);
	}

	delete exp_symb_count;
	delete occ_symb_count;

	return chiSq;
}


/* Finds the best possible mapping between the neighbours of the primarily matched input-query vertex pair.
 * Implements a greedy algorithm to find best match, which is dependent on order of search. For non-matching 
 * query neighbours, all non-mapped input neighbour is considered for possible approximate matches.
*/
vector<Triple> Input :: findNeighbourMapping(const Vertex* v, const Vertex* q, const unordered_map<string, Vertex*>& qryGraph) const
{
	const set<string> v_neig = v->getNeighbourIDs();
	const set<string> q_neig = q->getNeighbourIDs();

	unordered_map<string, vector<Vertex*> > qry_map, ver_map;  // Stores the vertices grouped according to labels

	// Create map for collecting vertices with the same label in the query neighbourhood
	set<string>::const_iterator it = q_neig.begin();
	for(; it!=q_neig.end(); it++)
	{
		Vertex* q = qryGraph.find(*it)->second;
		string lab = q->getLabel();

		if(qry_map.find(lab) == qry_map.end())
			qry_map[lab] = vector<Vertex*>{q};
		else
			qry_map.at(lab).push_back(q);
	}

	// Create map for collecting vertices with the same label in the input neighbourhood
	unordered_map<string, Vertex*> add_ver_map;  // Stores the mapping between vertex ID and the object
	it = v_neig.begin();
	for(; it!=v_neig.end(); it++)
	{
		Vertex* v = graph.at(*it);
		string lab = v->getLabel();
		string id = v->getID();

		add_ver_map[id] = v;

		if(ver_map.find(lab) == ver_map.end())
			ver_map[lab] = vector<Vertex*>{v};
		else
			ver_map.at(lab).push_back(v);
	}

	vector<Vertex*> add_ver;  // Stores unmapped vertices in the input vertex neighbourhood
	vector<Vertex*> add_qry;  // Stores unmapped query vertices
	vector<Triple> neig_pairs;
	neig_pairs.clear();

	unordered_map<string, vector<Vertex*> >::const_iterator it1 = qry_map.begin();
	for(; it1!=qry_map.end(); it1++)
	{
		string lab = it1->first;
		const vector<Vertex*> qry = it1->second;

		if(ver_map.find(lab) != ver_map.end())
		{
			vector<Vertex*> ver = ver_map.at(lab);  // Candidate vectors with the same label
			for(unsigned long i=0; i<qry.size(); i++)
			{
				Vertex* q = qry.at(i);
				double max_simi = -1.0;
				Triple best_pair;  // Stores input and query vertex best mapping with similarity score
				unsigned long pos;

				for(unsigned long j=0; j<ver.size(); j++)
				{
					Vertex* v = ver.at(j);
					double simi = findVertexSimi(v, q);

					if(simi > max_simi)
					{
						max_simi = simi;
						best_pair.setVal(v, q, simi);
						pos = j;
					}
				}

				if(max_simi != -1.0)
				{
					add_ver_map.erase(ver.at(pos)->getID()); // Remove the mapped candidate from list of unmapped vertices
					ver.erase(ver.begin()+pos);  // Remove the candidate vertex from the input graph
					neig_pairs.push_back(best_pair);
				}
				else
					add_qry.push_back(q);  // Add query vertex as an unmapped neighbour
			}
		}
		else
			add_qry.insert(add_qry.end(), qry.begin(), qry.end());  // Store unmapped query vertices
	}


	if(add_qry.size() > 0)
	{
		// For query vertices which do not have a vertex in the input graph with the same label, we 
		// consider all unmapped vertices of the neighbourhood as possible candidates. This enables
		// approximate matching of vertices between input and query graph in case of mis-matches.
		unordered_map<string, Vertex*>::const_iterator it = add_ver_map.begin();
		for(; it!=add_ver_map.end(); it++)
			add_ver.push_back(it->second);

		for(unsigned long i=0; i<add_qry.size(); i++)
		{
			Vertex* q = add_qry.at(i);
			double max_simi = -1.0;
			Triple best_pair;
			unsigned long pos;

			for(unsigned long j=0; j<add_ver.size(); j++)
			{
				Vertex* v = add_ver.at(j);
				double simi = findVertexSimi(v, q);

				if(simi > max_simi)
				{
					max_simi = simi;
					best_pair.setVal(v, q, simi);
					pos = j;
				}
			}

			if(max_simi != -1.0)
			{
				add_ver.erase(add_ver.begin()+pos);
				neig_pairs.push_back(best_pair);
			}
		}
	}

	return neig_pairs;
}


vector<Triple> Input :: findNeighbourMapping(const Vertex* v, const Vertex* q, const unordered_map<string, Vertex*>& qryGraph, const set<string>& ver_dup) const
{
	const set<string> v_neig = v->getNeighbourIDs();
	const set<string> q_neig = q->getNeighbourIDs();

	unordered_map<string, vector<Vertex*> > qry_map, ver_map;  // Stores the vertices grouped according to labels

	// Create map for collecting vertices with the same label in the query neighbourhood
	set<string>::const_iterator it = q_neig.begin();
	for(; it!=q_neig.end(); it++)
	{
		Vertex* q = qryGraph.find(*it)->second;
		string lab = q->getLabel();

		if(qry_map.find(lab) == qry_map.end())
			qry_map[lab] = vector<Vertex*>{q};
		else
			qry_map.at(lab).push_back(q);
	}

	// Create map for collecting vertices with the same label in the input neighbourhood
	unordered_map<string, Vertex*> add_ver_map;  // Stores the mapping between vertex ID and the object
	it = v_neig.begin();
	for(; it!=v_neig.end(); it++)
	{
		Vertex* v = graph.at(*it);
		string lab = v->getLabel();
		string id = v->getID();

		add_ver_map[id] = v;

		if(ver_map.find(lab) == ver_map.end())
			ver_map[lab] = vector<Vertex*>{v};
		else
			ver_map.at(lab).push_back(v);
	}

	vector<Vertex*> add_ver;  // Stores unmapped vertices in the input vertex neighbourhood
	vector<Vertex*> add_qry;  // Stores unmapped query vertices
	vector<Triple> neig_pairs;
	neig_pairs.clear();

	unordered_map<string, vector<Vertex*> >::const_iterator it1 = qry_map.begin();
	for(; it1!=qry_map.end(); it1++)
	{
		string lab = it1->first;
		const vector<Vertex*> qry = it1->second;

		if(ver_map.find(lab) != ver_map.end())
		{
			vector<Vertex*> ver = ver_map.at(lab);  // Candidate vectors with the same label
			for(unsigned long i=0; i<qry.size(); i++)
			{
				Vertex* q = qry.at(i);
				double max_simi = -1.0;
				Triple best_pair;  // Stores input and query vertex best mapping with similarity score
				unsigned long pos;

				for(unsigned long j=0; j<ver.size(); j++)
				{
					Vertex* v = ver.at(j);
					double simi = findVertexSimi(v, q);

					if( (simi > max_simi) && (ver_dup.find(v->getID()) == ver_dup.end()) )
					{
						max_simi = simi;
						best_pair.setVal(v, q, simi);
						pos = j;
					}
				}

				if(max_simi != -1.0)
				{
					add_ver_map.erase(ver.at(pos)->getID()); // Remove the mapped candidate from list of unmapped vertices
					ver.erase(ver.begin()+pos);  // Remove the candidate vertex from the input graph
					neig_pairs.push_back(best_pair);
				}
				else
					add_qry.push_back(q);  // Add query vertex as an unmapped neighbour
			}
		}
		else
			add_qry.insert(add_qry.end(), qry.begin(), qry.end());  // Store unmapped query vertices
	}


	if(add_qry.size() > 0)
	{
		// For query vertices which do not have a vertex in the input graph with the same label, we 
		// consider all unmapped vertices of the neighbourhood as possible candidates. This enables
		// approximate matching of vertices between input and query graph in case of mis-matches.
		unordered_map<string, Vertex*>::const_iterator it = add_ver_map.begin();
		for(; it!=add_ver_map.end(); it++)
			add_ver.push_back(it->second);

		for(unsigned long i=0; i<add_qry.size(); i++)
		{
			Vertex* q = add_qry.at(i);
			double max_simi = -1.0;
			Triple best_pair;
			unsigned long pos;

			for(unsigned long j=0; j<add_ver.size(); j++)
			{
				Vertex* v = add_ver.at(j);
				double simi = findVertexSimi(v, q);

				if( (simi > max_simi) && (ver_dup.find(v->getID()) == ver_dup.end()) )
				{
					max_simi = simi;
					best_pair.setVal(v, q, simi);
					pos = j;
				}
			}

			if(max_simi != -1.0)
			{
				add_ver.erase(add_ver.begin()+pos);
				neig_pairs.push_back(best_pair);
			}
		}
	}

	return neig_pairs;
}
