CFLAGS = -Wall -O3 -std=c++11
CC = -c

.PHONY = all

all: subgraph clean

subgraph: main.o query_graph.o input_graph.o vertex.o vertex.o
	g++ -o subgraph $(CFLAGS) main.o query_graph.o input_graph.o vertex.o

vertex.o: vertex.cpp vertex.h const.h query_graph.h
	g++ $(CC) $(CFLAGS) vertex.cpp

query_graph.o: query_graph.cpp query_graph.h vertex.h
	g++ $(CC) $(CFLAGS) query_graph.cpp

input_graph.o: input_graph.cpp input_graph.h vertex.h const.h
	g++ $(CC) $(CFLAGS) input_graph.cpp

main.o: main.cpp query_graph.h input_graph.h vertex.h const.h
	g++ $(CC) $(CFLAGS) main.cpp


clean:
	rm *.o
