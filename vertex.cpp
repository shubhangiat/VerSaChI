/*** Implements the vertex structure for the input graph (in vertex.h) and 
 *** computes the chi-square value of the vertices based on the occurrence
 *** and expected counts of the symbols
***/


#include "vertex.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>


Vertex :: Vertex(const string id, const string lab) : ID(id), label(lab)
{
	degree = 0;

	neigLabelVec.clear();
	neighbourIDs.clear();
	neighbourLabels.clear();
}


Vertex :: ~Vertex()
{
	degree = 0;

	neighbourIDs.clear();
	neighbourLabels.clear();
	neigLabelVec.clear();
}

const string Vertex :: getLabel(void) const
{
	return label;
}

const string Vertex :: getID(void) const
{
	return ID;
}

const unsigned Vertex :: getDegree(void) const
{
	return degree;
}

const set<string>& Vertex :: getNeighbourIDs(void) const
{
	return neighbourIDs;
}

const vector<string>& Vertex :: getNeighbourLabels(void) const
{
	return neighbourLabels;
}

void Vertex :: createNeighbourLabels(const unordered_map<string, Vertex*>& graph)
{
	set<string>::const_iterator it = neighbourIDs.begin();
	degree = neighbourIDs.size();

	for(; it!=neighbourIDs.end(); it++)
		neighbourLabels.push_back( (graph.find(*it))->second->getLabel() );
}

void Vertex :: addNeighbours(string neigh)
{
	neighbourIDs.insert(neigh);
}

void Vertex :: print(void) const
{
	cout<<"ID: "<<ID<<"\tLabel: "<<label<<"\tDegree: "<<degree<<"\tNeighbourIDs: ";

	set<string>::const_iterator it = neighbourIDs.begin();
	for(; it!=neighbourIDs.end(); it++)
			cout<<*it<<" ";

	cout<<"\tNeighbourLabels: ";
	for(unsigned i=0; i<neighbourLabels.size(); i++)
		cout<<neighbourLabels.at(i)<<" ";

	cout<<"\tNeighbourLabelVector: ";
	for(unsigned long i=0; i<neigLabelVec.size(); i++)
		cout<<neigLabelVec.at(i)<<" ";
	cout<<endl;	
}


void Vertex :: createNeighLabelVec(const unordered_map<string, unsigned int>& uniqLabel)
{
	vector<unsigned int> tmp_v(uniqLabel.size(), 0U);
	neigLabelVec = tmp_v;

	unsigned long indx = uniqLabel.find(this->label)->second;
	neigLabelVec[indx] = 1U;

	for(unsigned long i=0; i<neighbourLabels.size(); i++)
	{
		unsigned long indx = uniqLabel.find(neighbourLabels.at(i))->second;
		neigLabelVec[indx] ++;
	}
}


double findVertexSimi(const Vertex* v1, const Vertex* v2)  // Computes similarity between vertices: v1 from target graph and v2 from query graph
{
	const vector<unsigned int> src_neigLabVec = v1->neigLabelVec;
	const vector<unsigned int> trg_neigLabVec = v2->neigLabelVec;
	const unsigned long l_size = src_neigLabVec.size();

	if(SIMILARITY == 0)  // Overlap coefficient or Szymkiewicz–Simpson coefficient
	{
		unsigned long overlap, src_size, trg_size;
		overlap = src_size = trg_size = 0UL;

		for(unsigned long i=0; i<l_size; i++)
		{
			const unsigned int src_count = src_neigLabVec.at(i);
			const unsigned int trg_count = trg_neigLabVec.at(i);

			src_size += src_count;
			trg_size += trg_count;

			overlap += src_count < trg_count? src_count : trg_count;
		}
				
		const unsigned long denom = src_size > trg_size? trg_size : src_size;

		return ((1.0*overlap) / denom);
	}
	
	if(SIMILARITY == 1)  // Tversky index (alpha = beta = 1 -> Tanimoto coefficient, alpha = beta = 0.5 -> Sorensen-Dice coefficient / F1 measure)
	{
		unsigned long overlap, xydiff, yxdiff;
		overlap = xydiff = yxdiff = 0UL;

		for(unsigned long i=0; i<l_size; i++)
		{
			const unsigned int src_count = src_neigLabVec.at(i);
			const unsigned int trg_count = trg_neigLabVec.at(i);

			overlap += src_count < trg_count? src_count : trg_count;
			if(src_count > trg_count)
				xydiff += (src_count - trg_count);
			else if(trg_count > src_count)
				yxdiff += (trg_count - src_count);
		}

		//double alpha = 0.5, beta = 0.5;
		//double denom = overlap + alpha*xydiff + beta*yxdiff;
		const double denom = overlap + pow(yxdiff, 3.0);

		return ((1.0*overlap) / denom);
	}
}


double findVertexSimi(const Vertex* v1, const Vertex* v2, double& rev_simi)  // Bi-Directional similarity between the vertices
{
	const vector<unsigned int> src_neigLabVec = v1->neigLabelVec;
	const vector<unsigned int> trg_neigLabVec = v2->neigLabelVec;
	const unsigned long l_size = src_neigLabVec.size();

	if(SIMILARITY == 0)  // Overlap coefficient or Szymkiewicz–Simpson coefficient
	{
		unsigned long overlap, src_size, trg_size;
		overlap = src_size = trg_size = 0UL;

		for(unsigned long i=0; i<l_size; i++)
		{
			const unsigned int src_count = src_neigLabVec.at(i);
			const unsigned int trg_count = trg_neigLabVec.at(i);

			src_size += src_count;
			trg_size += trg_count;

			overlap += src_count < trg_count? src_count : trg_count;
		}
				
		const unsigned long denom = src_size > trg_size? trg_size : src_size;
		rev_simi = ((1.0*overlap) / denom);

		return rev_simi;
	}
	
	if(SIMILARITY == 1)  // Tversky index (alpha = beta = 1 -> Tanimoto coefficient, alpha = beta = 0.5 -> Sorensen-Dice coefficient / F1 measure)
	{
		unsigned long overlap, xydiff, yxdiff;
		overlap = xydiff = yxdiff = 0UL;

		for(unsigned long i=0; i<l_size; i++)
		{
			const unsigned int src_count = src_neigLabVec.at(i);
			const unsigned int trg_count = trg_neigLabVec.at(i);

			overlap += src_count < trg_count? src_count : trg_count;
			if(src_count > trg_count)
				xydiff += (src_count - trg_count);
			else if(trg_count > src_count)
				yxdiff += (trg_count - src_count);
		}

		//double alpha = 0.5, beta = 0.5;
		//double denom = overlap + alpha*xydiff + beta*yxdiff;
		const double denom = overlap + pow(yxdiff, 3.0);
		rev_simi = (1.0*overlap) / (overlap + pow(xydiff, 3.0));

		return ((1.0*overlap) / denom);
	}
}
