/*** The main file to create the input and query graphs and subsequently obtain the 
 *** approximate matching sub-graphs in the input graph to the query graph provided
***/


#include "const.h"
#include "query_graph.h"
#include "input_graph.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <fstream>


bool desc_sort(std::pair<vector<Vertex*>, double> a, std::pair<vector<Vertex*>, double> b)
{
	return (a.second > b.second);
}


bool vert_sort(Vertex *a, Vertex *b)
{
	return (a->getID() < b->getID());
}


int main(int argc, char *argv[])
{
	clock_t begin, end;

	if(argc!=4 && argc!=5)
	{
		//cout<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4];
		cout<<"USAGE: ./subgraph <input_graph_vertex_label_file> <input_graph_edge_file> <arg file> [input_characteristic_file]";
        cout<<"\n"<<"Please go through the readMe file\n";
        return 0;
	}

	// Create the input graph
	string i_label_file = argv[1];		// Input Vertex Label File
	string i_edge_file = argv[2];		// Input Edge File
	string i_char_file = "";			// Input Characteristics File
	if(argc == 5)
		i_char_file = argv[4];

	cout<<"Reading Input Graph files ("<<i_label_file<<", "<<i_edge_file<<") and Indexing Input Structure with STEP_SIZE "<<STEP_SIZE<<" ...\n";
	fflush(stdout);
	begin = clock();
	Input inp(i_label_file, i_edge_file);
	end = clock();
	cout<<"Reading + Indexing time of Input Graph = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;  // Input Graph Indexing Time

	begin = clock();
	inp.computeGraphCharacteristics(i_char_file);
	end = clock();
	cout<<"Time to compute characteristics of Input Graph = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;  // Input Graph Characterizing Time

	cout<<"Size of:\t"<<sizeof(inp)<<endl;
	inp.printGraph();


	// Create the query graph
	ifstream qfile;
	
	// Combined arg file time
	clock_t arg_file_time = clock();	
	
	qfile.open(argv[3]);

	string f1, f2;
	qfile>>f1>>f2;
	
	while(!qfile.eof())
	{
		string q_label_file = f1;	// Query Vertex Label File
		string q_edge_file = f2;	// Query Edge File
		
		// Check for query file existence
		ifstream iNQ(q_label_file.c_str()), iEQ(q_edge_file.c_str());
		
		if(!(iNQ.good() && iEQ.good()))
		{
			// Some file missing
			iNQ.close();
			iEQ.close();
			// Read next set of query files
			qfile>>f1>>f2;
			continue;
		}

		iNQ.close();
		iEQ.close();
		cout<<"\nQuery files:\t"<<f1<<"\t"<<f2<<endl;

		begin = clock();
		Query qry(q_label_file, q_edge_file);
		qry.createIndex(inp.getUniqLabel());
		end = clock();
		cout<<endl<<endl;
		cout<<"Reading + Pre-processing time of Query Graph = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;  // Input Graph Indexing Time

		qry.printGraph();

		begin = clock();
		const vector<vector<Triple> > subGraphs = inp.getMatchingSubgraphs(qry);
		end = clock();

		cout<<endl<<"Time to find top-"<<TOPK<<" approx. subgraph matches to Query Graph = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;  // Total run-time for finding matching sub-graph
		cout<<"The vertex mapping and chi-squared values of the obtained subgraph are:"<<endl;

		for(unsigned long i=0; i<subGraphs.size(); i++)
		{
			const vector<Triple> subGraph = subGraphs.at(i);
			for(unsigned j=0; j<subGraph.size(); j++)
			{
				Triple tr = subGraph.at(j);
				cout<<tr.v->getID()<<"("<<tr.v->getLabel()<<"):"<<tr.q->getID()<<"("<<tr.q->getLabel()<<") ["<<tr.simi<<"];\t";
			}
			cout<<endl;
		}

		cout<<"\n**************************************************************************"<<endl;
		qfile>>f1>>f2;
	}
	qfile.close();
	
	// Combined arg file time
	arg_file_time = ((clock() - arg_file_time))*1.0/CLOCKS_PER_SEC;
	
	cout<<"Combined arg file time processing time:\t"<<arg_file_time<<" sec."<<endl;

	return 0;
}
