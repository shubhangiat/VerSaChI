# *VerSaChI*

## Introduction
*VerSaChI* is an approximate subgraph matching framework based on 2-hop label
and structural overlap similarity with the query. The similarity is
characterized using Chebyshev’s inequality to compute the chi-square statistical
significance for measuring the degree of matching of the subgraphs.

Please cite our paper, if you use our source code.
* "VerSaChI: Finding Statistically Significant Subgraph Matches using Chebyshev's Inequality. CIKM'21"

## How to run the binary file?

Compile all the header files and dependent cpp files through make command and then run the binary file with appropriate arguments.

```
make  
./subgraph <input graph vertex-label file> <input graph edge file> <list of query graph files>
```

For instance,

```
./subgraph vlabel vedge arg_qfile
```

## Parameters

The _k_ for top-_k_ matching subgraphs and the step size parameter _κ_ (used for
discretizing standard deviations) can be set in the _const.h_ header file.

## Argument and Graph files format:

### Vertex file 
File format: vid label

e.g.
```
v1 l1  
v2 l2  
v3 l1  
```

### Edge file

File format: vid1 vid2

e.g.
```
v1 v2
v2 v3
v3 v1
```

This format is followed for both input target graph.

### Format of third argument (list of query graph files)

```
<query graph vertex-label file 1> <query graph edge file 1>
<query graph vertex-label file 2> <query graph edge file 2>
```
