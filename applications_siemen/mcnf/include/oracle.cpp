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
 *  AUXILIARY FUNCTIONS FOR THE ORACLE PROCEDURE
 *  Last Modified: July, 2013.
 * 
 *  On top of HOPDM  (Higher Order Primal-Dual Method)
 *  Copyright:     Jacek Gondzio, 1990, 2010.
 *
 *  Using the Dijkstra algorithm implemented by:
 *      Jacek Gondzio and Robert Sarkissian
 *
 */

/* Bring in auxiliary headers */
#include "oracle.h"

/* Parse the input given in the text line */
void parse_input (int argc, char *argv[], int *aggregated, int *oriented, 
         char DefltFlnme[], int *begin_zero, int *with_cost)
{
    /* Default values of control parameters */
    *aggregated = 0;
    *oriented   = 1;
    *begin_zero = 0;
    *with_cost  = 1;

    if (argc > 1) 
    {
        /* Get the name of the instance file */
        strcpy (DefltFlnme, argv[1]);
    }
    else
    {
        printf("\nError: Instance file name not specified.");
        exit(1);
    }
    
    if(argc > 2)
    {
        /* Check if aggregation is required */
        if(argv[2][0] == '1') 
        {
            *aggregated = 1;
        }
        else if(argv[2][0] == '2') 
        {
            *aggregated = 2;
        }
    } 

    if(argc > 3)
    {
        /* Check if aggregation is required */
        if(argv[3][0] == '0') 
        {
            USE_ACTIVE_SET = 0;
        }
    } 
}

/* Read the instance data from an input file */
void read_instance(char file_name[], int zero, int aggregated)
{   
    FILE *input_file;
    int    i, p;
    char   *ret, firstline[1000];// = "#Instance: \n";
    double flot, scale_factor, aux_double;

    /* Open the instance file */
    if (!(input_file = fopen(file_name, "r+t"))) 
    {
        printf("\nError: Could not open input file %s!\n", file_name);
        exit(1);
    }

    /* Ignore the first line of the file */
    ret = fgets(firstline, 1000, input_file);
    
    /* Read the dimensions of the instance */
    ret = fgets(firstline, 1000, input_file);
    p = sscanf(firstline, "%d %d %d %lf", &n_graphe, &m_graphe, &nbr_demande_graphe, &scale_factor);
    if(p < 3) 
    {
        printf("\nError: the input file is not in the required format (1).\n");
        exit(1);
    }
    else if(p < 4) scale_factor = 1.0; 
    
    /* Allocate memory */
    or_graphe = VECT(int, m_graphe);
    ex_graphe = VECT(int, m_graphe);
    capac_graphe = VECT(double, m_graphe);
    cout_graphe = VECT(double, m_graphe);
    dem_graphe = VECT(Demande, nbr_demande_graphe);
    
    if(aggregated == 2)
    {
        source_nodes = VECT(int, n_graphe);
        source_nodes_len = VECT(int, n_graphe);
        source_nodes_commodities = VECT(int, nbr_demande_graphe);
        source_nodes_n = 0;
    }
    
    /* Initialize the scale factor */
    max_cap = -1.0;
    max_cost = -1.0;

    /* Read information about arcs */
    for (i = 0; i < m_graphe; i++) 
    {
        ret = fgets(firstline, 1000, input_file);
        p = sscanf(firstline,  "%d %d %lf %lf", or_graphe + i, ex_graphe + i, cout_graphe + i, capac_graphe + i);
        if(p < 3) 
        {
            printf("\nError: the input file is not in the required format (1).\n");
            exit(1);
        }
        else if(p < 4) 
        {
            capac_graphe[i] = cout_graphe[i];
            cout_graphe[i] = 100.0;
        }
        
        /* Check if it is the largest cost */
        if(max_cost < cout_graphe[i]) max_cost = cout_graphe[i]; 
        
        /* Check if it is the largest capacity */
        if(capac_graphe[i] > max_cap) max_cap = capac_graphe[i];

        if (zero == 0) 
        {
            (or_graphe[i])--;
            (ex_graphe[i])--;
        }
    }
    
    /* "Round up" max_cap */
    aux_double = max_cap;
    max_cap = 1.0;
    if(aggregated != 0) while(max_cap < aux_double) max_cap *= 10.0;

    /* SCALE THE CAPACITY OF ARCS */
    if(max_cap > 1.0) for (i = 0; i < m_graphe; i++) capac_graphe[i] /= max_cap; 
    
    /* "Round up" max_cost */
    aux_double = max_cost;
    max_cost = 1.0;
    while(max_cost < aux_double) max_cost *= 10.0;

    /* SCALE THE LINEAR COST */
    if(max_cost > 1.0) for (i = 0; i < m_graphe; i++) cout_graphe[i] /= max_cost;
    
    /* Read information about commodities */
    for (i = 0; i < nbr_demande_graphe; i++) 
    {
        ret = fgets(firstline, 1000, input_file);
        p = sscanf(firstline,  "%d %d %lf", &((dem_graphe + i)->s), &((dem_graphe + i)->t), &flot);
        if(p < 3) 
        {
            printf("\nError: the input file is not in the required format (2).\n");
            exit(1);
        }

        /* Apply the scaling factor read from the input file */
        (dem_graphe + i)->val = scale_factor * flot;
        
        /* SCALE THE DEMAND OF COMMODITIES */
        if(max_cap > 1.0) (dem_graphe + i)->val /= max_cap; 
        
        if (zero == 0) 
        {
            ((dem_graphe + i)->s)--;
            ((dem_graphe + i)->t)--;
        }

        if(aggregated == 2)
        {
            /* Update the number of commodities associated to the corresponding source node */
            if(source_nodes_len[(dem_graphe + i)->s] == 0) source_nodes_n++;
            source_nodes_len[(dem_graphe + i)->s]++;
        }
    }
    
    if(aggregated == 2)
    {
        /* Assign comodities to their corresponding source nodes */
        if(source_nodes_n > 0)
        {
            /* Set the header of each group */
            source_nodes[0] = 0;
            for (i = 1; i < n_graphe; i++) 
            {
                source_nodes[i] = source_nodes[i-1] + source_nodes_len[i-1];
                source_nodes_len[i-1] = 0;
            }
            source_nodes_len[n_graphe-1] = 0;
            
            /* Set commodities in their groups */
            for (i = 0; i < nbr_demande_graphe; i++) 
            {
                source_nodes_commodities[ source_nodes[(dem_graphe + i)->s] + source_nodes_len[(dem_graphe + i)->s] ] = i;
                source_nodes_len[(dem_graphe + i)->s]++;
            }
        }
    }
    
    /* Close the input file */    
    fclose(input_file);
}

/*---------------------------- adjacency_list ----------------------------*/
void adjacency_list(Graph * N)
{
  /* calcule la liste adjacence la (arcs adjacents a un sommet)  du graphe */
  int    i, j;
  int    nb = 0;

  N->lp[0] = 0;

  if (N->oriented == 1) {
    for (i = 0; i < N->n; i++) {
      for (j = 0; j < N->m; j++) {
        if (N->org[j] == i) {
          N->la[nb] = j;
          nb++;
        }
      }
      N->lp[i + 1] = nb;
    }
  }
  else {
    for (i = 0; i < N->n; i++) {
      for (j = 0; j < N->m; j++) {
        if (N->org[j] == i || N->ex[j] == i) {
          N->la[nb] = j;
          nb++;
        }
      }
      N->lp[i + 1] = nb;
    }
  }
}

/*---------------------------- succesor_list ----------------------------*/
void succesor_list(Graph * N)
{
  /* calcule la liste adjacence ls (sommets adjacents a un sommet)  du graphe */
  int    i, k, u;

  for (k = 0; k < N->n; k++) {
    for (i = N->lp[k]; i < N->lp[k + 1]; i++) {
      u = N->la[i];
      if ((N->ls[i] = N->org[u]) == k)
        N->ls[i] = N->ex[u];
    }
  }
}

void set_network_Dijkstra(Graph *G, Demande **dem, int *nb_dem, int oriented)
{
    G->oriented = oriented;
    G->org = or_graphe;
    G->ex = ex_graphe;
    G->n = n_graphe;
    G->m = m_graphe;
    G->cost = cout_graphe;
    G->K= capac_graphe;
    G->la = VECT(int, 2 * m_graphe);
    G->ls = VECT(int, 2 * m_graphe);
    G->lp = VECT(int, n_graphe + 1);
    adjacency_list(G);
    succesor_list(G); 

    *dem = dem_graphe;
    *nb_dem = nbr_demande_graphe;
}

/* Initialize the data structure used by the Dijkstra algorithm */
void init_Dijkstra (Graph * G, Demande * dem, int nb_dem)
{
    Network = G;
    demande = dem;
    nb_demande = nb_dem;
    net_cost = VECT (double, G->m);
    prec_tableau = VECT (int, G->n * G->n);
    d_tableau = VECT (double, G->n * G->n);
    dijkstra_from_node = VECT (short int, G->n);
    active_set = NULL;
    if(USE_ACTIVE_SET) active_set = VECT (int, G->m);
    active_set_n = 0; 
}

/*---------------------------- Dijkstra ----------------------------*/
/* Find the shortest path in a graph using Dijkstra algorithm.      */
/* The adjacency lists H->la and H-​>ls are used.                    */
/*------------------------------------------------------------------*/
void ErrChem (int is, int it)
{
    printf ("There is no path between %i and %i\n", is, it);
    exit (1);
}

void Reset_Dijkstra (Graph * H)
{
  int i;
  for (i = 0; i < H->n; i++)
    dijkstra_from_node[i] = 0;
}

void Dijkstra_solve (Graph * H, Vertex org, Vertex ex, double *cost)
{
  int i, j, k;
  Vertex is, it;
  Arc u;
  double z;

  int *visited;
    
  is = org;
  it = ex;
  int *prec = prec_tableau + is * H->n;
  double *d = d_tableau + is * H->n;
  if (!dijkstra_from_node[is])
  {
     dijkstra_from_node[is] = 1;

    visited = VECT (int, H->n);

      for (i = 0; i < H->n; i++)
    {
      d[i] = DBL_MAX;
      prec[i] = INT_MIN;
      visited[i] = 0;
    }

      d[is] = 0;
      for (k = H->lp[is]; k < H->lp[is + 1]; k++)
    {
      u = H->la[k];
      j = H->ls[k];
      prec[j] = u;
      d[j] = cost[u];
    }
      visited[is] = 1;
      i = is;           /* d(is)=0, S={is}, i=is */

      do
    {
      z = DBL_MAX;
      i = INT_MIN;
      for (j = 0; j < H->n; j++)
        {
          if (visited[j] == 0 && d[j] <= z)
        {
          i = j;
          z = d[j];
        }
        }
      if (i < 0) break;

      for (k = H->lp[i]; k < H->lp[i + 1]; k++)
        {
          u = H->la[k];
          j = H->ls[k];
          if (visited[j] == 0)
        {
          z = d[i] + cost[u];
          if (z < d[j])
            {
              d[j] = z;
              prec[j] = u;
            }
        }
        }
      visited[i] = 1;
      }
      while (1);
      
      if (!visited[it]) ErrChem (is, it);
    
        free(visited);
    }
    //std::cout << "finished solving Dijkstra\n";
}

/* Generate a *sparse* column using the flow of a shortest path found by the Dijkstra algorithm */
/* It is used in the disaggregated oracle */
int get_solution_from_Dijkstra(Graph * G, PDCGM_SMatrix_CW * M, int orig, int extr, double val)
{
    int u, som, *prec;

    som = extr;
    prec = prec_tableau + orig * G->n;
    if (prec[extr] == INT_MIN) 
    {
        /* No path was found by the Dijkstra algorithm */
        printf("\n### Error: No path was found by the Dijkstra algorithm.\n");
        fflush(stdout);
        return 0;
    }
    
    while (som != orig)
    {
        u = prec[som];
        M->m_rwnmbs.push_back( u );
        //printf("-val=%lf\n",-val);
        M->m_coeff.push_back( -val);
        //printf("-val=%lf\n",-val);
        (M->m_nz)++;
        if (som == G->ex[u]) som = G->org[u];
        else som = G->ex[u];
    }
    //printf("In Dijkstra: \n");
    //(*M).print_SMatrix_CW();
    
    return 1;
}

/* Add to a dense column the flow of a shortest path obtained by the Dijkstra algorithm */
/* It is used in the aggregated oracle */
int add_solution_from_Dijkstra(Graph * G, double *vector, int orig, int extr, double val)
{
    int u, som, *prec;

    som = extr;
    prec = prec_tableau + orig * G->n;
    if (prec[extr] == INT_MIN)
    {
        /* No path was found by the Dijkstra algorithm */
        printf("\n### Error: No path was found by the Dijkstra algorithm.\n");
        fflush(stdout);
        return 0;
    }
    
    while (som != orig)
    {
        u = prec[som];
        vector[u] -= val;
        
        if (som == G->ex[u]) som = G->org[u];
        else som = G->ex[u];
    }
    
    return 1;
}
