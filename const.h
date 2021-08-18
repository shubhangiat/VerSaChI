#ifndef CONST_H_
 #define CONST_H_


// Define the operational constants of the system
#define ORDER_CONSTANT 10  // The factor of the maximum size of the internal heap that is considered

#define STEP_SIZE 0.001  // Defines the deviation buckets for symbol probability

#define TOPK 1  // Defines the number of best matched subgraphs to be returned

#define SIMILARITY 1 // Vertex label similarity measure used; 0 -> Overlap Coefficient, 1 -> Tversky Index

#endif
