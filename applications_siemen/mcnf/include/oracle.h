/*
 *  PDCGM  (Primal-Dual Column Generation Method)
 *  Copyright:     Jacek Gondzio, 1996, 2013.
 *
 *  COLUMN GENERATION FOR MULTICOMMODITY NETWORK FLOW PROBLEMS
 *  USING TWO DIFFERRENT FORMULATIONS: AGGREGATED AND DISAGGREGATED
 * 
 *  Authors: Pedro Munari [munari@dep.ufscar.br], 
 *           Pablo González-Brevis, 
 *           Jacek Gondzio
 *  
 *  PARAMETERS AND FUNCTION HEADERS
 *  Last Modified: July, 2013.
 * 
 *  On top of HOPDM  (Higher Order Primal-Dual Method)
 *  Copyright:     Jacek Gondzio, 1990, 2010.
 *
 *  Using the Dijkstra algorithm implemented by:
 *      Jacek Gondzio and Robert Sarkissian
 *
 */
 
 #ifndef ORACLE_H
#define ORACLE_H

#include "../pdcgm_env.h"
#include "../pdcgm_SMatrix.h"
#include <climits>
#include <cfloat>

 #ifdef __cplusplus
extern "C" {
#endif

typedef struct Graph_struct 
{
    int    n;
    int    m;
    int    oriented;
    int   *org;
    int   *ex;
    int   *la;
    int   *lp;
    int   *ls;
    double *K;
    double *cost;
} Graph;

typedef struct 
{
    int s; /* source de la demande */
    int t; /* destination de la demande */
    double val; /* quantite de flot */
} 
Demande;

typedef int Vertex;
typedef int Arc;
typedef double Flow;

#define EdgeEq(N,i,j)     ((N.org[i] == N.org[j] && N.ex[i] == N.ex[j])  ||   \
                 (N.org[i] == N.ex[j] && N.ex[i] == N.org[j]) )

#define EdgeEqp(N,i,j)    EdgeEq((*N),i,j)

#define max_edge_graph(n) ((n)*((n)-1)/2)

#define VECT(type, dim)  (type*) calloc((unsigned) dim, (unsigned) sizeof(type))

#define min(a,b)       ((a<b) ? a : b)
#define max(a,b)       ((a>b) ? a : b)

#define show(expr)     printf(#expr " = %g\n", expr)


static Graph *Network;
static Demande *demande;
static int nb_demande;
static double *d_tableau, *net_cost;
static int *prec_tableau;
static short int *dijkstra_from_node;

static int m_graphe, 
           n_graphe, 
           *or_graphe, 
           *ex_graphe,
           nbr_demande_graphe;

static double *capac_graphe, *cout_graphe;

static Demande *dem_graphe;

static void adjacency_list(Graph*);
static void succesor_list(Graph*);

static int *active_set;
static int active_set_n;

static int *source_nodes;
static int *source_nodes_len;
static int *source_nodes_commodities;
static int source_nodes_n;

double max_cap;
double max_cost;

/* Parse the input given in the text line */
static void parse_input (int argc, char *argv[], int *aggregated, int *oriented, 
    char DefltFlnme[], int *begin_zero, int *with_cost);

/* Read the instance data from an input file */
void read_instance(char file_name[], int zero, int aggregated);

/* Header functions of the Dijkstra algorithm */
void set_network_Dijkstra(Graph *G, Demande **dem, int *nb_dem, int oriented);
void init_Dijkstra (Graph * G, Demande * dem, int nb_dem);
void Reset_Dijkstra (Graph * H);
void Dijkstra_solve (Graph * H, Vertex org, Vertex ex, double *cost);

/* Generate a *sparse* column using the flow of a shortest path found by the Dijkstra algorithm */
/* It is used in the disaggregated oracle */
int get_solution_from_Dijkstra(Graph * G, PDCGM_SMatrix_CW * M, int orig, int extr, double val);

/* Add to a dense column the flow of a shortest path obtained by the Dijkstra algorithm */
/* It is used in the aggregated oracle */
int add_solution_from_Dijkstra(Graph * G, double *vector, int orig, int extr, double val);

#ifdef __cplusplus
}
#endif
#endif
