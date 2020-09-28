#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <deque>
#include <stack>
#include <time.h>
#include <random>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <string>
#include <string.h>
//#include <cstring>
#include <math.h>
#include <chrono>
#include <ctime>
#include <sstream>
#include <cmath>

using namespace std;


#define N_NODE 100 // number of nodes
#define D 6 // maximum edge density
#define N_EDGE (D*N_NODE) // maximum number of edges
#define N_POINTS (D*N_NODE) // number of data points to record
#define ASSORT_STOP -1 	// number of steps after which process is stopped and
						// assortativity values are printed (-1 = no stopping)
#define DIRECTED 1	// boolean marker for directed (1) or undirected (0)
#define ALPHA 0 // constant used in certain selection methods
#define TAG 5 // number to add to end of save files
#define N_RUNS 1 // number of full realizations to run through
#define CP 5.75 // edge density at which network structure files are saved
/*
 * MODES:
 * 	- 1: edge competition
 * 	- 2: local
 * 	- 3: quasi-local -- **disabled**
 * 	- 4: locally tree-like -- **disabled**
 * 	- 5: competitive direction selection
 * 	- 6: edge competition with damage reset (annealed)
 * 	- 7: AB hybrid node selection
 * 	- 8: AB size node selection
 * 	- 9: AB sensitivity node selection
 * 	- 10: TRRUST (needs N = 603, D = 5)
 * 	- 11: ENCODE (needs N = 116, D = 5)
 * 	- 12: Tissue (needs N = 644, D = 71)
 * 	- other: random edge selection
 */
#define MODE 7 // edge selection method to use (see above)
#define M 3 // number of edges to consider in competitive selection processes
#define L 0 // depth of components in local rules
#define DAMAGE_FLAG 0 // used in output files
#define FILL_IN 0 // prioritize connecting nodes in same component


typedef unordered_set<int> NODES;
typedef chrono::high_resolution_clock clock_type;

struct EDGE {
	int i; //source
	int j; //destination
};

vector<EDGE> gene_edges;
vector<EDGE> gene_edges_copy;
unordered_map<int,int> gene_map;

NODES empty_set;

int inmax[N_NODE];
int outmax[N_NODE];


// graph information
NODES graph[N_NODE]; // graph
NODES r_graph[N_NODE]; // reversed graph

char in_gin[N_NODE]; // tracks which nodes belong to GIN
char in_gout[N_NODE]; // tracks which nodes belong to GOUT

NODES gin; // set representation of GIN
NODES gout; // set representation of GOUT

int gscc_size;
int gout_size;
int gin_size;

int din[N_NODE]; // in-degree of each node
int dout[N_NODE]; // out-degree of each node

// variables for computing correlations/assortativity
double Ts[2];
double Td[2];
double Tss[2];
double Tdd[2];
double Tsd[4];
double assort[4];

double Tq;
double Tqq;
double Td2[2];
double Tdd2[2];
double Tqd[2];
double corrqd[2];

double iso_Tq;
double gout_Tq;
double gout_star_Tq;
//-----------------------------

NODES sources;
NODES dests;


// damaged subgraph information
int NSTAR;
float b[N_NODE];
float q[N_NODE];
bool is_damaged[N_NODE];
NODES damaged;

unsigned long outstar[N_NODE];
unsigned long totoutstar;
unsigned long maxoutstar;


NODES graph_star[N_NODE];
NODES r_graph_star[N_NODE];

char in_gin_star[N_NODE];
char in_gout_star[N_NODE];

NODES gin_star;
NODES gout_star;

int gscc_star_size;
int gout_star_size;
int gin_star_size;

int din_star[N_NODE];
int dout_star[N_NODE];

bool active[N_NODE];
int ntilde;

double Ts_star[2];
double Td_star[2];
double Tss_star[2];
double Tdd_star[2];
double Tsd_star[4];
double assort_star[4];

NODES sources_star;
NODES dests_star;


// BFS and Tarjan alg variables
stack<int> st;
int low[N_NODE];
int depth[N_NODE];
bool on_stack[N_NODE];
int tarj_index;

int node_list[M][2];

int ptr[N_NODE];


// random number generator
mt19937 gen((int) time(0));
//mt19937 gen(8);
uniform_int_distribution<int> dist(0,N_NODE-1);
mt19937 gen2((int) time(0));
//mt19937 gen2(3);
uniform_real_distribution<float> b_dist(0.0,1.0);

//void test();

int main();
void main2();
void clear_vars();

void percolate();

void add_edge(int, int);
bool is_adjacent(int, int);

void new_edge(int &, int &);
bool select_edge(int &, int &);
bool select_edge_er(int &, int &);
bool select_edge_ap(int &, int &, int);
bool select_edge_AB_hybrid(int &, int &, int);
bool select_edge_AB_size(int &, int &, int);
bool select_edge_AB_sensitivity(int &, int &, int);
bool select_edge_cds(int &, int &, int);
bool select_edge_local(int &, int &, int, int);
bool select_edge_quasilocal(int &, int &, int, int);
bool select_edge_tl(int &, int &, int);

bool bfs(int, int, NODES *, NODES *, NODES &, int &, int);
bool bfs(int, int, NODES *, NODES *, NODES &, int &, int, bool, char *);

bool check_path(int, int);
bool check_path_star(int, int);
double get_prod(int, int);
double get_in(int);
double get_out(int);

void reset();

void scc_util(int, NODES *, int &);
//void get_scc(int, NODES &);

int find_root(int);

void us_union(NODES &, NODES &, NODES &);
void us_intersect(NODES &, NODES &, NODES &);

void mark(int, char *, NODES &, int &, NODES *, bool, NODES &, bool);
void mark(int, char *, NODES &, int &, NODES *, bool, NODES &, bool, double &);
void new_gscc(int, int);
void new_gscc_star(int, int);

void update_assort(int, int, int);
void update_assort_star(int, int, int);
void update_assort(int, int, int, double *, double *, double *, double *, double *, double *, int *, int *);
void update_corr(int, int, int);

void load_gene_files();
void randomize_edges();

int main() {

	// load data files if needed
	if (MODE == 10 || MODE == 11 || MODE == 12) {
		load_gene_files();
	}

//	for (int i = 0; i < N_NODE; i++) {
//		printf("%d: %f\n", i, q[i]);
//	}

	for (int n = 0; n < N_RUNS; n ++) {
		cout << "Run " << n+1 << ": ";

		if (MODE == 10 || MODE == 11 || MODE == 12) {
			randomize_edges();
			for (EDGE e : gene_edges) {
				gene_edges_copy.push_back(e);
			}
			reverse(gene_edges_copy.begin(),gene_edges_copy.end());
		}

		main2(); // run a realization
		clear_vars(); // reset variables
	}

//	Beep(523,200);

	return 0;
}



void main2() {

//	LARGE_INTEGER start_time;
	clock_type::time_point start_time;
//	LARGE_INTEGER stop_time;
	clock_type::time_point stop_time;
//	LARGE_INTEGER frequency;

//	QueryPerformanceFrequency(&frequency);
	start_time = clock_type::now();
//	QueryPerformanceCounter(&start_time);

	ofstream runtime_file;
	runtime_file.open("runtimes.csv", ofstream::out | ofstream::app);

	percolate(); // main algorithm


//	QueryPerformanceCounter(&stop_time);
	stop_time = clock_type::now();
//	timer = (stop_time.QuadPart - start_time.QuadPart) * 1.0 / frequency.QuadPart;
	chrono::duration<double> timer = stop_time - start_time;
	cout << "Done! Time elapsed: " << timer.count() << "s.\n";
	cout << flush;


	runtime_file << N_NODE << ", " << D << ", " << DIRECTED << ", " << DAMAGE_FLAG << ", ";
	runtime_file << MODE << ", " << M << ", " << L << ", " << timer.count() << "\n";
	runtime_file.close();

}

// reset variables
void clear_vars() {

	for (int i = 0; i < N_NODE; i++) {
		inmax[i] = 0;
		outmax[i] = 0;

		graph[i].clear();
		r_graph[i].clear();

		in_gin[i] = 0;
		in_gout[i] = 0;

		din[i] = 0;
		dout[i] = 0;

		//b[i] = 0;
		//q[i] = 0;
		is_damaged[i] = false;
		outstar[i] = 0;

		graph_star[i].clear();
		r_graph_star[i].clear();

		in_gin_star[i] = 0;
		in_gout_star[i] = 0;

		din_star[i] = 0;
		dout_star[i] = 0;

		active[i] = false;

		low[i] = 0;
		depth[i] = 0;
		on_stack[i] = false;
	}


	gin.clear();
	gout.clear();

	gscc_size = 0;
	gout_size = 0;
	gin_size = 0;

	iso_Tq = 0;
	gout_Tq = 0;
	gout_star_Tq = 0;

	Ts[0] = 0;
	Ts[1] = 0;
	Td[0] = 0;
	Td[1] = 0;

	Tss[0] = 0;
	Tss[1] = 0;
	Tdd[0] = 0;
	Tdd[1] = 0;
	Tsd[0] = 0;
	Tsd[1] = 0;
	Tsd[2] = 0;
	Tsd[3] = 0;

	assort[0] = 0;
	assort[1] = 0;
	assort[2] = 0;
	assort[3] = 0;

	sources.clear();
	dests.clear();


	// damaged subgraph information
	NSTAR = 0;
	damaged.clear();

	totoutstar = 0;
	maxoutstar = 0;

	gin_star.clear();
	gout_star.clear();

	gscc_star_size = 0;
	gout_star_size = 0;
	gin_star_size = 0;


	ntilde = 0;

	Ts_star[0] = 0;
	Ts_star[1] = 0;
	Td_star[0] = 0;
	Td_star[1] = 0;

	Tss_star[0] = 0;
	Tss_star[1] = 0;
	Tdd_star[0] = 0;
	Tdd_star[1] = 0;
	Tsd_star[0] = 0;
	Tsd_star[1] = 0;
	Tsd_star[2] = 0;
	Tsd_star[3] = 0;

	assort_star[0] = 0;
	assort_star[1] = 0;
	assort_star[2] = 0;
	assort_star[3] = 0;

	sources_star.clear();
	dests_star.clear();

	tarj_index = 0;

	Tq = 0;
	Tqq = 0;
	Td2[0] = 0;
	Td2[1] = 0;
	Tdd2[0] = 0;
	Tdd2[1] = 0;
	Tqd[0] = 0;
	Tqd[1] = 0;
	corrqd[0] = 0;
	corrqd[1] = 0;

}

// generate a network using the selection rule set by MODE and generates files
void percolate() {
	int i,j,k;

	stringstream prestream, prestream2;
	string prefix;
	prestream.str(prefix);
	prestream2.str(prefix);

	prestream << "data_dir_" << MODE << "_m" << M << "_l" << L;
	prestream2 << "data_dir_" << MODE << "_m" << M << "_l" << L;

	if (MODE == 2) {
		prestream << "_dmg" << DAMAGE_FLAG;
		prestream2 << "_dmg" << DAMAGE_FLAG;
	} else if (MODE == 1 || MODE == 5) {
		prestream << "_dmg" << DAMAGE_FLAG;
		prestream2 << "_dmg" << DAMAGE_FLAG;
	} else {
		prestream << "_dmg0";
		prestream2 << "_dmg0";
	}

	if (MODE == 5) {
		prestream << "_alp" << ALPHA;
		prestream2 << "_alp" << ALPHA;
	}

	prestream << "_tag" << TAG;

	prefix = prestream2.str();


	ofstream main_file;
	main_file.open(prefix + ".csv", ofstream::out | ofstream::app);

	prefix = prestream.str();

	ofstream outfile;
	outfile.open(prefix + ".csv");

	ofstream outfile_star;
	outfile_star.open(prefix + "_star.csv");

	ofstream assortfile;
	assortfile.open(prefix + "_assort.csv");

	ofstream assortfile_star;
	assortfile_star.open(prefix + "_assort_star.csv");

	ofstream structure_file;
	structure_file.open(prefix + "_struct.csv");

//	structure_file << "SOURCE,TARGET,INTERACTION\n";

	ofstream node_file;
	node_file.open(prefix + "_node.csv");

	node_file << "NODE,B,Q,DAMAGED,DIN,DOUT\n";

	ofstream comp_file;
	comp_file.open(prefix + "_comp.csv");


	ofstream corr_file;
	corr_file.open(prefix + "_corrqd.csv");

	int skip = N_EDGE / N_POINTS;

	k = dist(gen); // random starting node

	in_gin[k] = 1;
	in_gout[k] = 1;

	gout.insert(k);
	gin.insert(k);

	gscc_size = 1;
	gout_size = 1;
	gin_size = 1;

	bool star_set = false;
	int nstar = 0;

	int new_size = 1;

	bool is_path;

	NODES visited;


//	LARGE_INTEGER frequency;
//
//
//	QueryPerformanceFrequency(&frequency);


	NODES temp_nodes;
	NODES temp_nodes2;
	EDGE temp_edge;

	NODES temp_out;
	NODES temp_in;
	int temp_out_size;
	int temp_in_size;


	for (int n = 0; n < N_NODE; n++) {

		inmax[n]++;
		outmax[n]++;

		// assign sensitivities
		if (!(MODE == 10 || MODE == 11 || MODE == 12)) {
			b[n] = b_dist(gen2);
			q[n] = 2*b[n]*(1-b[n]);
		}

		iso_Tq += q[n];

		// assign damaged status
		if (b_dist(gen2) < q[n]) {
			is_damaged[n] = true;
			NSTAR++;

			if (!star_set) {
				in_gin_star[n] = 1;
				in_gout_star[n] = 1;

				gout_star.insert(n);
				gin_star.insert(n);

				gscc_star_size = 1;
				gout_star_size = 1;
				gin_star_size = 1;

				star_set = true;
			}


		}


//		structure_file << n << "," << n << ",s\n";

	}


	for (int n = 0; n < N_EDGE; n++) {

		if (MODE == 10 || MODE == 11 || MODE == 12) {
			if (gene_edges_copy.empty()) {
				break;
			} else {
				temp_edge = gene_edges_copy.back();
				gene_edges_copy.pop_back();
				i = temp_edge.i;
				j = temp_edge.j;
				is_path = check_path(i, j);
			}

		// select and add new edge
		} else {
			do {
				is_path = select_edge(i, j);
			} while (i == j || graph[i].count(j));
		}

		add_edge(i, j);
//		printf("Adding edge %d -> %d\n", i+1, j+1);

		if (n < (int) N_NODE*CP) {
			structure_file << i << "," << j << "\n";
		}

		// update corr/assort values
		if (!active[i]) {
			ntilde++;
			active[i] = true;
			Tq += q[i];
			Tqq += q[i]*q[i];
			iso_Tq -= q[i];
		}

		if (!active[j]) {
			ntilde++;
			active[j] = true;
			Tq += q[j];
			Tqq += q[j]*q[j];
			iso_Tq -= q[j];
		}

		din[j]++;
		dout[i]++;

		sources.insert(i);
		dests.insert(j);

		inmax[j] += inmax[i];
		outmax[i] += outmax[j];


		update_assort(i, j, n+1);
		update_corr(i, j, ntilde);

		// determine method to use based on membership of i,j to GIN and GOUT
		k = 8*in_gin[i] + 4*in_gout[i] + 2*in_gin[j] + in_gout[j];

		switch (k) {

			case 0:
			case 5:
				if (!is_path) {
					scc_util(j, graph, new_size);
					reset();
				}
				break;


			case 2:
				scc_util(i, r_graph, new_size);
				reset();
				break;


			case 4:
				scc_util(j, graph, new_size);
				reset();
				break;


			case 10:
				if (!is_path) {
					scc_util(i, r_graph, new_size);
					reset();
				}
				break;
		}

		tarj_index = 0;

		// update if new GSCC forms
		if (new_size > gscc_size) {
			new_gscc(i, new_size);

		// update GIN, GOUT if no new GSCC; use method based on membership of i,j to GIN and GOUT
		} else {
			switch (k) {
				case 2:
				case 3:
					mark(i, in_gin, gin, gin_size, r_graph, false, empty_set, false);
					break;
				case 4:
				case 12:
					mark(j, in_gout, gout, gout_size, graph, false, empty_set, false, gout_Tq);
					break;
				case 6:
					mark(i, in_gin, gin, gin_size, r_graph, true, gout, false);
					mark(j, in_gout, gout, gout_size, graph, true, gin, false, gout_Tq);
					break;
				case 7:
					mark(i, in_gin, gin, gin_size, r_graph, true, gout, false);
					break;
				case 14:
					mark(j, in_gout, gout, gout_size, graph, true, gin, false, gout_Tq);
					break;
			}
		}

		new_size = 1;

		// update damaged subgraph using similar approach as above
		if (is_damaged[i] && is_damaged[j]) {
			nstar++;

			is_path = check_path_star(i, j);

			graph_star[i].insert(j);
			r_graph_star[j].insert(i);

			din_star[j]++;
			dout_star[i]++;

			sources_star.insert(i);
			dests_star.insert(j);

			update_assort_star(i, j, nstar);


			k = 8*in_gin_star[i] + 4*in_gout_star[i] + 2*in_gin_star[j] + in_gout_star[j];

			switch (k) {

				case 0:
				case 5:
					if (!is_path) {
						scc_util(j, graph_star, new_size);
						reset();
					}
					break;


				case 2:
					scc_util(i, r_graph_star, new_size);
					reset();
					break;


				case 4:
					scc_util(j, graph_star, new_size);
					reset();
					break;


				case 10:
					if (!is_path) {
						scc_util(i, r_graph_star, new_size);
						reset();
					}
					break;
			}

			tarj_index = 0;


			if (new_size > gscc_star_size) {
				new_gscc_star(i, new_size);

			} else {
				switch (k) {
					case 2:
					case 3:
						mark(i, in_gin_star, gin_star, gin_star_size, r_graph_star, false, empty_set, true);
						break;
					case 4:
					case 12:
						mark(j, in_gout_star, gout_star, gout_star_size, graph_star, false, empty_set, true, gout_star_Tq);
						break;
					case 6:
						mark(i, in_gin_star, gin_star, gin_star_size, r_graph_star, true, gout_star, true);
						mark(j, in_gout_star, gout_star, gout_star_size, graph_star, true, gin_star, true, gout_star_Tq);
						break;
					case 7:
						mark(i, in_gin_star, gin_star, gin_star_size, r_graph_star, true, gout_star, true);
						break;
					case 14:
						mark(j, in_gout_star, gout_star, gout_star_size, graph_star, true, gin_star, true, gout_star_Tq);
						break;
				}
			}

			new_size = 1;

		}

		// output files
		if (n % skip == 0) {
			assortfile << (float) n / N_NODE << ", ";
			assortfile << assort[0] << ", " << assort[1] << ", ";
			assortfile << assort[2] << ", " << assort[3] << "\n";

			// in-in, in-out, out-in, out-out
			assortfile_star << (float) n / N_NODE << ", " << (float) nstar / NSTAR << ", ";
			assortfile_star << assort_star[0] << ", " << assort_star[1] << ", ";
			assortfile_star << assort_star[2] << ", " << assort_star[3] << "\n";

			// q*din, q*dout
			corr_file << (float) n/N_NODE << ", " << corrqd[0] << ", " << corrqd[1] << "\n";


	//		printf("Adding edge %d -> %d: New Size = %d\n", i+1, j+1, gout_alt_size);

	//		output_nodes(outfile);
			outfile << (double) n /N_NODE << ", " << ((double) gscc_size) / N_NODE << ", ";
			outfile << ((double) gin_size) / N_NODE << ", " << ((double) gout_size) / N_NODE << ", ";
			outfile << ((double) (gin_size + gout_size - gscc_size)) / N_NODE << "\n";
	//		print_graph();

			outfile_star << (double) n /N_NODE << ", " << (float) nstar / NSTAR << ", " << ((double) gscc_star_size) / NSTAR << ", ";
			outfile_star << ((double) gin_star_size) / NSTAR << ", " << ((double) gout_star_size) / NSTAR << ", ";
			outfile_star << ((double) (gin_star_size + gout_star_size - gscc_star_size)) / NSTAR << ",";
			//outfile_star << ((double) totoutstar) / (NSTAR) << ", " << ((double) maxoutstar) << "\n";
			outfile_star << (double) gout_star_size / N_NODE << ", " << (double) ntilde / N_NODE << "\n";

			if (n != 0) {
				main_file << ", ";
			}

			main_file << (double) n /N_NODE << ", " << (double) ntilde / N_NODE << ", " << (double) gout_star_size / N_NODE;
			main_file << ", " << (double) gout_size / N_NODE;
//			main_file << ", ";
//			if (ntilde < N_NODE) {
//				main_file << iso_Tq / (N_NODE - ntilde);
//			} else {
//				main_file << (int) 0;
//			}
			//main_file << ", " << gout_star_Tq / gout_star_size << ", " << gout_Tq / gout_size;
			main_file << ", " << assort[1] << ", " << assort_star[1];
			main_file << ", " << (double) gscc_size / N_NODE << ", " << (double) gscc_star_size / N_NODE;
			main_file << ", " << (double) (gout_size + gin_size - gscc_size) / N_NODE;
			main_file << ", " << (double) (gout_star_size + gin_star_size - gscc_star_size) / N_NODE;
		}


		if (n == (int) CP*N_NODE) {
			for (int m = 0; m < N_NODE; m++) {
				node_file << m << "," << b[m] << "," << q[m] << ",";
				if (is_damaged[m]) {
					node_file << "true,";
				} else {
					node_file << "false,";
				}
				node_file << din[m] << "," << dout[m] << "\n";
			}
		}


//		if (n == N_NODE || n == 3*N_NODE || n == 5*N_NODE) {
//			for (int m = 0; m < N_NODE; m++) {
//				temp_out.clear();
//				temp_out.insert(m);
//				temp_out_size = 1;
//				bfs(m, -1, graph, r_graph, temp_out, temp_out_size, -1);
//
//				temp_in.clear();
//				temp_in.insert(m);
//				temp_in_size = 1;
//				bfs(m, -1, r_graph, graph, temp_in, temp_in_size, -1);
//
//				temp_nodes.clear();
//				us_intersect(temp_in, temp_out, temp_nodes);
//
//				comp_file << (double) temp_nodes.size() / N_NODE << "," << (double) temp_out_size / N_NODE << "," << (double) temp_in_size / N_NODE;
//
//				if (m == N_NODE - 1) {
//					comp_file << "\n";
//				} else {
//					comp_file << ",";
//				}
//			}
//		}

	}

//	main_file << "\n";
	main_file << ",";


	outfile.close();
	assortfile.close();
	outfile_star.close();
	assortfile_star.close();
	main_file.close();
	structure_file.close();
	node_file.close();
	comp_file.close();
	corr_file.close();

}

// add edge i->j to the graph
void add_edge(int i, int j) {
	if (i != j) {
		graph[i].insert(j);
		r_graph[j].insert(i);
	}
}

// determines if edge s->d exists
bool is_adjacent(int s, int d) {
	return graph[s].count(d);
}

// randomly select a new edge i->j that is not yet in the graph
void new_edge(int &i, int &j) {
	do {
		i = dist(gen);
		j = dist(gen);
	} while (i == j || graph[i].count(j));
}

// select an edge using the method specified by MODE
bool select_edge(int &i, int &j) {

	switch (MODE) {
	case 1:
	case 6:
		// Achlioptas-like selection
		return select_edge_ap(i,j,M);
		break;

	case 2:
		// Local
		return select_edge_local(i,j,M,L);
		break;

//	case 3:
//		// Quasi-local
//		return select_edge_quasilocal(i,j,M,L);
//		break;

//	case 4:
//		// Local Tree-like
//		return select_edge_tl(i,j,M);
//		break;

	case 5:
		// Competitive Direction Selection
		return select_edge_cds(i,j,M);
		break;

	case 7:
		// AB hybrid
		return select_edge_AB_hybrid(i,j,M);
		break;

	case 8:
		// AB size
		return select_edge_AB_size(i,j,M);
		break;

	case 9:
		// AB sensitivity
		return select_edge_AB_sensitivity(i,j,M);
		break;

	default:
		// Erdos-Renyi
		return select_edge_er(i,j);
		break;
	}


}


// select new edge uniformly randomly
bool select_edge_er(int &i, int &j) {
	new_edge(i, j);
	return check_path(i, j);
}

// DISABLED
bool select_edge_tl(int &i, int &j, int m) {
	int e[m][2];

	int itemp, jtemp;

	for (int k = 0; k < m; k++) {
		new_edge(itemp, jtemp);

		if (check_path(itemp, jtemp)) {
			i = itemp;
			j = jtemp;
			return true;
		}

		e[k][0] = itemp;
		e[k][1] = jtemp;
	}


	double minprod, tempprod;

	i = e[0][0];
	j = e[0][1];
	minprod = ((double) inmax[i] / N_NODE) * ((double) outmax[j] / N_NODE);

	for (int k = 1; k < m; k++) {
		itemp = e[k][0];
		jtemp = e[k][1];

		tempprod = ((double) inmax[itemp] / N_NODE) * ((double) outmax[jtemp] / N_NODE);

		if (tempprod < minprod) {
			i = itemp;
			j = jtemp;
			minprod = tempprod;
		}

	}

	return false;
}

// Achloptas-like process (competitive directed selection)
bool select_edge_ap(int &i, int &j, int m) {
	int e[m][2];

	int itemp, jtemp;

	for (int k = 0; k < m; k++) {
		new_edge(itemp, jtemp);

		if (check_path(itemp, jtemp)) {
			i = itemp;
			j = jtemp;
			return true;
		}

		e[k][0] = itemp;
		e[k][1] = jtemp;
	}


	double minprod, tempprod;

	i = e[0][0];
	j = e[0][1];
	minprod = get_prod(i, j);

	for (int k = 1; k < m; k++) {
		itemp = e[k][0];
		jtemp = e[k][1];

		tempprod = get_prod(itemp, jtemp);

		if (tempprod < minprod) {
			i = itemp;
			j = jtemp;
			minprod = tempprod;
		}

	}

	return false;
}

// select edges, direction always from low sensitivity to high (with probability of flipping)
// and choose the edge with the smallest product IN(i)*OUT(j)
bool select_edge_cds(int &i, int &j, int m) {
	int e[m][2];

	int itemp, jtemp, temp;

	for (int k = 0; k < m; k++) {

		do {
			itemp = dist(gen);
			jtemp = dist(gen);

			// make sure edges go from low sensitivity to high
			if (q[itemp] > q[jtemp]) {
				temp = itemp;
				itemp = jtemp;
				jtemp = temp;
			}

			// possibility for edge to go in reverse (allow for cycles)
			if (b_dist(gen2) > 0.5*pow(1 - 2*abs(q[itemp] - q[jtemp]), ALPHA)) {
				temp = itemp;
				itemp = jtemp;
				jtemp = temp;
			}
		} while (itemp == jtemp || graph[itemp].count(jtemp));

		if (check_path(itemp, jtemp)) {
			i = itemp;
			j = jtemp;
			return true;
		}

		e[k][0] = itemp;
		e[k][1] = jtemp;
	}


	double minprod, tempprod;

	i = e[0][0];
	j = e[0][1];
	minprod = get_prod(i, j);

	for (int k = 1; k < m; k++) {
		itemp = e[k][0];
		jtemp = e[k][1];

		tempprod = get_prod(itemp, jtemp);

		if (tempprod < minprod) {
			i = itemp;
			j = jtemp;
			minprod = tempprod;
		}

	}

	return false;
}

// competitively select i with lowest sensitivity*IN and j with lowest sensitivity*OUT
// from m options each
bool select_edge_AB_hybrid(int &i, int &j, int m) {


	int itemp, jtemp;
	double mina, minb, tempa, tempb;

	i = dist(gen);
	mina = q[i]*get_in(i);

	j = dist(gen);
	minb = q[j]*get_out(j);

	node_list[0][0] = i;
	node_list[0][1] = j;

	for (int k = 1; k < m; k++) {
		itemp = dist(gen);
		tempa = q[itemp]*get_in(itemp);

		if (tempa < mina) {
			i = itemp;
			mina = tempa;
		}

		jtemp = dist(gen);
		tempb = q[jtemp]*get_out(jtemp);

		if (tempb < minb) {
			j = jtemp;
			minb = tempb;
		}

		node_list[k][0] = itemp;
		node_list[k][1] = jtemp;

	}

	if (FILL_IN) {
		for (int k1 = 0; k1 < m; k1++) {
			for (int k2 = 0; k2 < m; k2++) {
				itemp = node_list[k1][0];
				jtemp = node_list[k2][1];

				// a path exists
				if (check_path(itemp,jtemp)) {
					i = itemp;
					j = jtemp;
					return true;
				}
			}
		}
	}


	return check_path(i,j);
}

// competitively select i with smallest IN and j with smallest OUT from m options each
bool select_edge_AB_size(int &i, int &j, int m) {

	int itemp, jtemp;
	double mina, minb, tempa, tempb;

	i = dist(gen);
	mina = get_in(i);

	j = dist(gen);
	minb = get_out(j);

	node_list[0][0] = i;
	node_list[0][1] = j;

	for (int k = 1; k < m; k++) {
		itemp = dist(gen);
		tempa = get_in(itemp);

		if (tempa < mina) {
			i = itemp;
			mina = tempa;
		}

		jtemp = dist(gen);
		tempb = get_out(jtemp);

		if (tempb < minb) {
			j = jtemp;
			minb = tempb;
		}

		node_list[k][0] = itemp;
		node_list[k][1] = jtemp;

	}

	if (FILL_IN) {
		for (int k1 = 0; k1 < m; k1++) {
			for (int k2 = 0; k2 < m; k2++) {
				itemp = node_list[k1][0];
				jtemp = node_list[k2][1];

				// a path exists
				if (check_path(itemp,jtemp)) {
					i = itemp;
					j = jtemp;
					return true;
				}
			}
		}
	}

	return check_path(i,j);
}

// competitively select nodes i and j independently with lowest sensitivity from m options
bool select_edge_AB_sensitivity(int &i, int &j, int m) {

	int itemp, jtemp;
	double mina, minb, tempa, tempb;

	i = dist(gen);
	mina = q[i];

	j = dist(gen);
	minb = q[j];

	node_list[0][0] = i;
	node_list[0][1] = j;

	for (int k = 1; k < m; k++) {
		itemp = dist(gen);
		tempa = q[itemp];

		if (tempa < mina) {
			i = itemp;
			mina = tempa;
		}

		jtemp = dist(gen);
		tempb = q[jtemp];

		if (tempb < minb) {
			j = jtemp;
			minb = tempb;
		}

		node_list[k][0] = itemp;
		node_list[k][1] = jtemp;

	}

	if (FILL_IN) {
		for (int k1 = 0; k1 < m; k1++) {
			for (int k2 = 0; k2 < m; k2++) {
				itemp = node_list[k1][0];
				jtemp = node_list[k2][1];

				// a path exists
				if (check_path(itemp,jtemp)) {
					i = itemp;
					j = jtemp;
					return true;
				}
			}
		}
	}

	return check_path(i,j);
}

// Achlioptas-like selection with maximum depth when counting size of IN(i) and OUT(j)
bool select_edge_local(int &i, int &j, int m, int max_depth) {

	int e[m][2];

	int itemp, jtemp;
	NODES visited;
	int size;

	for (int k = 0; k < m; k++) {
		new_edge(itemp, jtemp);

		visited.clear();
		visited.insert(itemp);
		size = 1;

		if (bfs(itemp, jtemp, graph, r_graph, visited, size, max_depth)) {
			i = itemp;
			j = jtemp;
			return true;
		}

		e[k][0] = itemp;
		e[k][1] = jtemp;
	}


	double minprod, tempprod;
	int tempin, tempout;

	i = e[0][0];
	j = e[0][1];

	visited.clear();
	visited.insert(i);
	tempin = 1;
	bfs(i, -1, r_graph, graph, visited, tempin, max_depth);

	visited.clear();
	visited.insert(j);
	tempout = 1;
	bfs(j, -1, graph, r_graph, visited, tempout, max_depth);

	minprod = (((double) tempin) / N_NODE) * (((double) tempout) / N_NODE);

	for (int k = 1; k < m; k++) {
		itemp = e[k][0];
		jtemp = e[k][1];

		if (L == 1) {
			tempin = din[itemp];
			tempout = dout[jtemp];

		} else {

			visited.clear();
			visited.insert(itemp);
			tempin = 1;
			bfs(itemp, -1, r_graph, graph, visited, tempin, max_depth);

			visited.clear();
			visited.insert(jtemp);
			tempout = 1;
			bfs(jtemp, -1, graph, r_graph, visited, tempout, max_depth);
		}


		tempprod = (((double) tempin) / N_NODE) * (((double) tempout) / N_NODE);

		if (tempprod < minprod) {
			i = itemp;
			j = jtemp;
			minprod = tempprod;
		}

	}

	return check_path(i,j);
}

// DISABLED
bool select_edge_quasilocal(int &i, int &j, int m, int max_depth) {
	int e[m][2];

	int itemp, jtemp;
	NODES visited;

	for (int k = 0; k < m; k++) {
		new_edge(itemp, jtemp);


		if (check_path(itemp, jtemp)) {
			i = itemp;
			j = jtemp;
			return true;
		}

		e[k][0] = itemp;
		e[k][1] = jtemp;
	}


	double minprod, tempprod;
	int tempin, tempout;

	i = e[0][0];
	j = e[0][1];

	visited.clear();
	visited.insert(i);
	tempin = 1;
	bfs(i, -1, r_graph, graph, visited, tempin, max_depth);

	visited.clear();
	visited.insert(j);
	tempout = 1;
	bfs(j, -1, graph, r_graph, visited, tempout, max_depth);

	minprod = (((double) tempin) / N_NODE) * (((double) tempout) / N_NODE);

	for (int k = 1; k < m; k++) {
		itemp = e[k][0];
		jtemp = e[k][1];

		visited.clear();
		visited.insert(itemp);
		tempin = 1;
		bfs(itemp, -1, r_graph, graph, visited, tempin, max_depth);

		visited.clear();
		visited.insert(jtemp);
		tempout = 1;
		bfs(jtemp, -1, graph, r_graph, visited, tempout, max_depth);

		tempprod = (((double) tempin) / N_NODE) * (((double) tempout) / N_NODE);

		if (tempprod < minprod) {
			i = itemp;
			j = jtemp;
			minprod = tempprod;
		}

	}

	return false;
}

// utility function for Tarjan's algorithm
void scc_util(int i, NODES * to, int &scc_size) {

	tarj_index++;
	depth[i] = tarj_index;
	low[i] = tarj_index;

	st.push(i);
	on_stack[i] = true;


	for (int k : to[i]) {

		if (depth[k] == 0) {
			scc_util(k, to, scc_size);
			low[i] = min(low[i], low[k]);
		} else if (on_stack[k]) {
			low[i] = min(low[i], depth[k]);
		}

	}


	int j;
	if (low[i] == depth[i]) {
		if (depth[i] == 1) {
			scc_size = st.size();
		}

		do {
			j = st.top();
			st.pop();
			on_stack[j] = false;

		} while (j != i);

	}


}


// determines if a directed path exists from i to j
bool check_path(int i, int j) {

	if (DIRECTED == 0) {
		return find_root(i) == find_root(j);
	}

	int k = 8*in_gin[i] + 4*in_gout[i] + 2*in_gin[j] + in_gout[j];
	bool is_path;
	NODES visited;
	int size = 1;

	switch (k) {
		// i,j not in GBT or
		// j in gout, i not in gin
		// search for path from i to j
		case 0:
		case 1:
		case 5:
			visited.insert(i);
			is_path = bfs(i, j, graph, r_graph, visited, size, -1);
			break;

		// i in gin, j not in gout
		// search for reverse path from j to i
		case 8:
		case 10:
			visited.insert(j);
			is_path = bfs(j, i, r_graph, graph, visited, size, -1);
			break;

		// i in gin, j in gout
		// path must exist
		case 9:
		case 11:
		case 13:
		case 15:
			is_path =  true;
			break;

		// all other cases cannot have a path
		default:
			is_path = false;
	}

	return is_path;
}

// determines if a directed path exists in damaged subgraph from i to j
bool check_path_star(int i, int j) {


	int k = 8*in_gin_star[i] + 4*in_gout_star[i] + 2*in_gin_star[j] + in_gout_star[j];
	bool is_path;
	NODES visited;
	int size = 1;

	switch (k) {
		// i,j not in GBT or
		// j in gout, i not in gin
		// search for path from i to j
		case 0:
		case 1:
		case 5:
			visited.insert(i);
			is_path = bfs(i, j, graph_star, r_graph_star, visited, size, -1);
			break;

		// i in gin, j not in gout
		// search for reverse path from j to i
		case 8:
		case 10:
			visited.insert(j);
			is_path = bfs(j, i, r_graph_star, graph_star, visited, size, -1);
			break;

		// i in gin, j in gout
		// path must exist
		case 9:
		case 11:
		case 13:
		case 15:
			is_path =  true;
			break;

		// all other cases cannot have a path
		default:
			is_path = false;
	}

	return is_path;
}

// returns the product of size(IN(i)) and size(OUT(j)) using info about giant components
double get_prod(int i, int j) {

	if (DIRECTED == 0) {
		int r1, r2;
		r1 = find_root(i);
		r2 = find_root(j);
		return (((double) -ptr[r1]) / N_NODE) * (((double) -ptr[r2]) / N_NODE);
	}

	int in_size;
	int out_size;
	NODES visited;
	visited.clear();

	// find |IN(i)|
	if (in_gout[i]) {
		in_size = gin_size;

		// if i is in gin as well, then IN(i) = GIN
		// otherwise, add all missed nodes
		if (!in_gin[i]) {
			visited.insert(i);
			in_size++;
			bfs(i, -1, r_graph, graph, visited, in_size, -1, true, in_gin);
		}

	//perform full BFS if i not in GOUT
	} else {
		in_size = 1;
		visited.insert(i);
		bfs(i, -1, r_graph, graph, visited, in_size, -1);
	}

	// reset list of visited nodes
	visited.clear();

	// find |OUT(j)|
	// case if j is in GIN
	if (in_gin[j]) {
		out_size = gout_size;

		// if j is also in GOUT, then OUT(j) = GOUT
		// otherwise add missing nodes
		if (!in_gout[j]) {
			visited.insert(j);
			out_size++;
			bfs(j, -1, graph, r_graph, visited, out_size, -1, true, in_gout);
		}

	// otherwise, full BFS
	} else {
		out_size = 1;
		visited.insert(j);
		bfs(j, -1, graph, r_graph, visited, out_size, -1);
	}

	return (((double) in_size) / N_NODE) * (((double) out_size) / N_NODE);
}

// get the size of IN(i)
double get_in(int i) {

	if (DIRECTED == 0) {
		int r1;
		r1 = find_root(i);
		return (((double) -ptr[r1]) / N_NODE);
	}

	int in_size;
	NODES visited;
	visited.clear();

	// find |IN(i)|
	if (in_gout[i]) {
		in_size = gin_size;

		// if i is in gin as well, then IN(i) = GIN
		// otherwise, add all missed nodes
		if (!in_gin[i]) {
			visited.insert(i);
			in_size++;
			bfs(i, -1, r_graph, graph, visited, in_size, -1, true, in_gin);
		}

	//perform full BFS if i not in GOUT
	} else {
		in_size = 1;
		visited.insert(i);
		bfs(i, -1, r_graph, graph, visited, in_size, -1);
	}

	return (((double) in_size) / N_NODE);
}

// get the size of OUT(j)
double get_out(int j) {

	if (DIRECTED == 0) {
		int r2;
		r2 = find_root(j);
		return (((double) -ptr[r2]) / N_NODE);
	}

	int out_size;
	NODES visited;
	visited.clear();

	// find |OUT(j)|
	// case if j is in GIN
	if (in_gin[j]) {
		out_size = gout_size;

		// if j is also in GOUT, then OUT(j) = GOUT
		// otherwise add missing nodes
		if (!in_gout[j]) {
			visited.insert(j);
			out_size++;
			bfs(j, -1, graph, r_graph, visited, out_size, -1, true, in_gout);
		}

	// otherwise, full BFS
	} else {
		out_size = 1;
		visited.insert(j);
		bfs(j, -1, graph, r_graph, visited, out_size, -1);
	}

	return (((double) out_size) / N_NODE);
}

// BFS (see below)
bool bfs(int s, int d, NODES * to, NODES * from, NODES &visited, int &out_size, int max_depth) {
	return bfs(s, d, to, from, visited, out_size, max_depth, false, NULL);
}


/**
 * Performs a BFS starting at node s, ending at node d using the graph "to"
 * and its reverse graph "from", storing all visited nodes in "visited" and number of
 * nodes in the out component in "out_size".
 *
 * BFS does not expand any nodes in "visited". Make sure "visited" and "out_size" are set
 * to their initial values.
 *
 * Returns whether there is a path from s to d, and stops the BFS early.
 * Set d to -1 to perform a full BFS.
 *
 * BFS does not expand nodes at "max_depth".
 * Set "max_depth" to -1 to perform a full search.
 */
bool bfs(int s, int d, NODES * to, NODES * from, NODES &visited, int &out_size, int max_depth, bool flag, char * ignore) {


	if (to[s].empty()) {
		return false;
	}

	// BREADTH FIRST SEARCH

	reset();



	// create queue for BFS
	deque<int> queue;

	// mark starting node as visited and enqueue it
	queue.push_back(s);
	depth[s] = 0;

	int curr_depth;

	while (!queue.empty()) {
		//dequeue a vertex from queue
		s = queue.front();
		queue.pop_front();
		curr_depth = depth[s];

		if (max_depth == -1 || curr_depth < max_depth) {

			// get all adjacent vertices of s
			for (int v : to[s]) {

				// if adjacent node is destination, the return true
				if (d != -1 && v == d) {
					return true;
				}


				// skip if in ignore list
				if (!flag || !ignore[v]) {
					// if root has not been visited, mark visited and enqueue
					if (!visited.count(v)) {
						visited.insert(v);
						queue.push_back(v);
						depth[v] = curr_depth + 1;
						out_size++;
					}
				}
			}
		}

	}


	// completed search and could not find path
	return false;
}



// calculate the intersection of two sets "a" and "b"
void us_intersect(NODES &a, NODES &b, NODES &aib) {

	// b is smaller, so use as iterator
	if (a.size() > b.size()) {
		us_intersect(b, a, aib);
		return;

	// a is smaller
	} else {
		NODES::const_iterator iter;
		int v;

		aib.clear();
		for (iter = a.begin(); iter != a.end(); iter++) {
			v = *iter;
			if (b.count(v))
				aib.insert(v);
		}
	}
}

// calculate the union of two sets "a" and "b"
void us_union(NODES &a, NODES &b, NODES &aub) {
	aub = a;
	aub.insert(b.begin(), b.end());
}

// resets depth variable
void reset() {
	memset(depth, 0, sizeof(depth));
}


void mark(int i, char * in_g, NODES &g_set, int &g_size, NODES * to, bool comp, NODES &g_comp, bool star) {
	double temp = 0;
	mark(i, in_g, g_set, g_size, to, comp, g_comp, star, temp);
}

// follows edges in graph "to" starting from node i, marking nodes as beloning to g_set
// Ex: if "to" is the standard graph, then identifies GOUT; if "to" is reversed graph, then identifies GIN
void mark(int i, char * in_g, NODES &g_set, int &g_size, NODES * to, bool comp, NODES &g_comp, bool star, double &tot_q) {
	in_g[i] = 1;
	g_set.insert(i);
	g_size++;
	tot_q += q[i];

	if (comp && g_comp.count(i)) {
		if (star) {
			gscc_star_size++;
		} else {
			gscc_size++;
		}
	}

	for (int k : to[i]) {
		if (!in_g[k]) {
			mark(k, in_g, g_set, g_size, to, comp, g_comp, star, tot_q);
		}
	}
}

// sets GIN and GOUT membership for new GSCC containing node i
void new_gscc(int i, int new_size) {
	memset(in_gin, 0, sizeof(in_gin));
	memset(in_gout, 0, sizeof(in_gout));
	gin.clear();
	gout.clear();

	gin_size = 0;
	gout_size = 0;
	gscc_size = new_size;
	gout_Tq = 0;

	mark(i, in_gin, gin, gin_size, r_graph, false, empty_set, false);
	mark(i, in_gout, gout, gout_size, graph, false, empty_set, false, gout_Tq);
}

// sets GIN* and HOUT* for new GSCC* containing node i
void new_gscc_star(int i, int new_size) {
	memset(in_gin_star, 0, sizeof(in_gin_star));
	memset(in_gout_star, 0, sizeof(in_gout_star));
	gin_star.clear();
	gout_star.clear();

	gin_star_size = 0;
	gout_star_size = 0;
	gscc_star_size = new_size;
	gout_star_Tq = 0;

	mark(i, in_gin_star, gin_star, gin_star_size, r_graph_star, false, empty_set, true);
	mark(i, in_gout_star, gout_star, gout_star_size, graph_star, false, empty_set, true, gout_star_Tq);
}


// calculates the new assortativity based on new edge i->j
void update_assort(int i, int j, int n) {
	update_assort(i, j, n, Ts, Td, Tss, Tdd, Tsd, assort, din, dout);
}

// calculates the new assortativity in damaged subnetwork with new edge i->j
void update_assort_star(int i, int j, int n) {
	update_assort(i, j, n, Ts_star, Td_star, Tss_star, Tdd_star, Tsd_star, assort_star, din_star, dout_star);
}

// updates assortativity for specified graph with the specified new edge
void update_assort(int i, int j, int n, double *Ts, double *Td, double *Tss, double *Tdd, double *Tsd, double *assort, int *din, int *dout) {

	double vars[2];
	double vard[2];
	double covsd[4];


	Ts[0] += dout[j] + din[i];
	Ts[1] += 2*dout[i] - 1;
	Td[0] += 2*din[j] - 1;
	Td[1] += din[i] + dout[j];

//	if (dout[i] == 1 || din[j] == 1) {
//		Ts[0] += din[i];
//		Td[1] += dout[j];
//	}

	Tss[0] += dout[j]*(2*din[j] - 1) + din[i]*din[i];
	Tss[1] += (dout[i] - 1)*(dout[i] - 1) + dout[i]*(2*dout[i] - 1);
	Tdd[0] += (din[j] - 1)*(din[j] - 1) + din[j]*(2*din[j] - 1);
	Tdd[1] += din[i]*(2*dout[i] - 1) + dout[j]*dout[j];

	Tsd[0] += din[i]*din[j] - din[i];
	Tsd[1] += din[i]*dout[j];
	Tsd[2] += dout[i]*din[j] - din[j] - dout[i];
	Tsd[3] += dout[i]*dout[j] - dout[j];


	for (int k : graph[i]) {
		Tsd[2] += din[k];
		Tsd[3] += dout[k];
	}

	for (int k : graph[j]) {
		Tsd[0] += din[k];
		Tsd[1] += dout[k];
	}

	for (int k : r_graph[i]) {
		Tsd[1] += din[k];
		Tsd[3] += dout[k];
	}

	for (int k : r_graph[j]) {
		Tsd[0] += din[k];
		Tsd[2] += dout[k];
	}


	vars[0] = Tss[0] - Ts[0]*Ts[0]/n;
	vars[1] = Tss[1] - Ts[1]*Ts[1]/n;

	vard[0] = Tdd[0] - Td[0]*Td[0]/n;
	vard[1] = Tdd[1] - Td[1]*Td[1]/n;

	covsd[0] = Tsd[0] - Ts[0]*Td[0]/n;
	covsd[1] = Tsd[1] - Ts[0]*Td[1]/n;
	covsd[2] = Tsd[2] - Ts[1]*Td[0]/n;
	covsd[3] = Tsd[3] - Ts[1]*Td[1]/n;

	assort[0] = covsd[0]/(sqrt(vars[0])*sqrt(vard[0]));
	assort[1] = covsd[1]/(sqrt(vars[0])*sqrt(vard[1]));
	assort[2] = covsd[2]/(sqrt(vars[1])*sqrt(vard[0]));
	assort[3] = covsd[3]/(sqrt(vars[1])*sqrt(vard[1]));

	for (int k = 0; k < 4; k++) {
		if (!isfinite(assort[k])) {
			assort[k] = 0;
		}
	}

	if (n == ASSORT_STOP) {
		printf("i -> j : %d -> %d, %lu\n", i, j, graph[j].count(i));
		printf("din[i] = %d, dout[i] = %d\n", din[i], dout[i]);
		printf("din[j] = %d, dout[j] = %d\n\n", din[j], dout[j]);
		printf("Fast:\n");
		printf("E[sin] = %f\n", Ts[0]);
		printf("E[sout] = %f\n", Ts[1]);
		printf("E[din] = %f\n", Td[0]);
		printf("E[dout] = %f\n", Td[1]);
		printf("Var[sin] = %f\n", vars[0]/n);
		printf("Var[sout] = %f\n", vars[1]/n);
		printf("Var[din] = %f\n", vard[0]/n);
		printf("Var[dout] = %f\n", vard[1]/n);
		printf("Cov[sin,din] = %f\n", covsd[0]/n);
		printf("Cov[sin,dout] = %f\n", covsd[1]/n);
		printf("Cov[sout,din] = %f\n", covsd[2]/n);
		printf("Cov[sout,dout] = %f\n", covsd[3]/n);
		printf("r[in,in] = %f\n", assort[0]);
		printf("r[in,out] = %f\n", assort[1]);
		printf("r[out,in] = %f\n", assort[2]);
		printf("r[out,out] = %f\n", assort[3]);
	}

}

// updates degree-sensitivity correlation based on new edge i->j
void update_corr(int i, int j, int n) {

	double varq = Tqq - Tq*Tq/n;
	double vard[2];
	double covqd[2];

	Td2[0]++;
	Td2[1]++;

	Tdd2[0] += 2*din[j] - 1;
	Tdd2[1] += 2*dout[i] - 1;

	Tqd[0] += q[j];
	Tqd[1] += q[i];



	vard[0] = Tdd2[0] - Td2[0]*Td2[0]/n;
	vard[1] = Tdd2[1] - Td2[1]*Td2[1]/n;

	covqd[0] = Tqd[0] - Tq*Td2[0]/n;
	covqd[1] = Tqd[1] - Tq*Td2[1]/n;

	corrqd[0] = covqd[0]/(sqrt(varq)*sqrt(vard[0]));
	corrqd[1] = covqd[1]/(sqrt(varq)*sqrt(vard[1]));

	for (int k = 0; k < 2; k++) {
		if (!isfinite(corrqd[k])) {
			corrqd[k] = 0;
		}
	}

}

// used in undirected percolation
int find_root(int i) {
	if (ptr[i] < 0) {
		return i;
	}

	return ptr[i] = find_root(ptr[i]);
}

// load files containing human gene network information
void load_gene_files() {
	string edge_path = "";
	if (MODE == 10) {
		edge_path = "Human-Tissue-Data/edges_trrust.csv";
	} else if (MODE == 11) {
		edge_path = "Human-Tissue-Data/edges_encode.csv";
	} else if (MODE == 12) {
		edge_path = "Human-Tissue-Data/edges_tissue.csv";
	}

	string bias_path = "Human-Tissue-Data/bias_all.csv";

	ifstream edge_file;
	edge_file.open(edge_path, ifstream::in);
	printf("Edge file opened.\n");

	ifstream bias_file;
	bias_file.open(bias_path, ifstream::in);
	printf("Bias file opened.\n");

	string line;
	int num;
	int ind = 0;
	int node_num;
	float q_val;

	while (getline(edge_file, line)) {
		stringstream ss(line);
		string field;

		EDGE new_edge;
		getline(ss, field, ',');
		stringstream s_num1(field);
		s_num1 >> num;

		if (gene_map.count(num)) {
			node_num = gene_map[num];
		} else {
			std::pair<int,int> newpair (num,ind);
			gene_map.insert(newpair);
			node_num = ind;
			//printf("ind = %d\n", ind);
			ind++;

		}

		new_edge.i = node_num;

		getline(ss, field, ',');
		stringstream s_num2(field);
		s_num2 >> num;

		if (gene_map.count(num)) {
			node_num = gene_map[num];
		} else {
			std::pair<int,int> newpair (num,ind);
			gene_map.insert(newpair);
			node_num = ind;
			ind++;
		}

		new_edge.j = node_num;
		gene_edges.push_back(new_edge);

	}

	printf("Edges loaded.\n");

	while (getline(bias_file, line)) {
		stringstream ss(line);
		string field;

		getline(ss, field, ',');
		stringstream s_num1(field);
		s_num1 >> num;

		if (gene_map.count(num)) {
			node_num = gene_map[num];
		} else {
			continue;
		}


		getline(ss, field, ',');
		stringstream s_num2(field);
		s_num2 >> q_val;

		q[node_num] = q_val;

	}

	printf("Sensitivities loaded.\n");

	edge_file.close();
	bias_file.close();

}

// randomizes edges for human gene network
void randomize_edges() {
	random_shuffle(gene_edges.begin(), gene_edges.end());
}
