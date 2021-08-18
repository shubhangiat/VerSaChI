#ifndef CONST_H_
 #define CONST_H_


// Define the operational constants of the system
#define ORDER_CONSTANT 10  // The factor of the maximum size of the internal heap that is considered

#define STEP_SIZE 4  // Defines the deviation buckets for symbol probability

#define TOPK 1  // Defines the number of best matched subgraphs to be returned

#define SIMILARITY 1 // Vertex label similarity measure used; 0 -> Overlap Coefficient, 1 -> Tversky Index

/*
// Define the input files for the query graph
#define QUERY_VERTEX_FILE "./Data/qry_node.txt"
#define QUERY_EDGE_FILE "./Data/qry_edge.txt"

// Define the input files for the input graph
#define INPUT_VERTEX_FILE "./Data/inp_node.txt"
#define INPUT_EDGE_FILE "./Data/inp_edge.txt"
*/

/*
// Define the input files for the query graph
#define QUERY_VERTEX_FILE "../../dataset/q_v5_label.txt"
#define QUERY_EDGE_FILE "../../dataset/q_v5_edge.txt"

// Define the input files for the input graph
#define INPUT_VERTEX_FILE "../../dataset/ip_v100_label.txt"
#define INPUT_EDGE_FILE "../../dataset/ip_v100_edge.txt"
*/

#endif
