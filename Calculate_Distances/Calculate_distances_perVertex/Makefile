
# *****************************************************
# Variables to control Makefile operation

CC = g++
CCFLAGS = -Wall -O2 -std=c++11 

# ****************************************************
# 

make_calculate_distances_i: makeObj/calculate_distances_i.o makeObj/common.o
	$(CC) $(CCFLAGS) -o make_calculate_distances_i makeObj/calculate_distances_i.o makeObj/common.o -s


makeObj/calculate_distances_i.o: common.h geodesic_algorithm_dijkstra.h geodesic_algorithm_subdivision.h geodesic_algorithm_exact.h
	$(CC) $(CCFLAGS) -c calculate_distances_i.cpp -o makeObj/calculate_distances_i.o

makeObj/common.o: common.h
	$(CC) $(CCFLAGS) -c common.cpp -o makeObj/common.o


