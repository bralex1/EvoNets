# EvoNets

EvoNets.cpp is used to simulate an evolving regulatory network starting with N nodes and no edges, 
adding one edge at a time following a specified selection rule (chosen by "mode"). The size of the 
giant components of the main network, the size of the giant components of the damaged subnetwork, 
assortativity values, and degree-sensitivity correlations are tracked and output to data files. These
files are turned into graphs with the MATLAB script perc_processing.m.

For a specified value of edge-density, a network structure file and can be used to produce a network
image with a circular layout, ordered by node sensitivity.

Algorithm.pdf contains information about the algorithm used in EvoNets.cpp.

The following versions were used for the relevant files:
* C++14
* Python 3.7
* Matlab 2017a
