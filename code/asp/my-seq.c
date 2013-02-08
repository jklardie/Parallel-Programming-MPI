#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>

//#define VERBOSE

#define MAX_DISTANCE 256    // max distance between places randomly generated

/**
 * Initialize a random adjacency matrix of size n x n.
 */
void init_random_adjacency_matrix(int n, int *num_edges, int ***matrix, int oriented){
    int **am;   // adjacency matrix
    int i, j, m = n * n;

    am = (int **) malloc(n * sizeof(int *));
    if(am == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<n; i++){
        am[i] = (int *) malloc(n * sizeof(int));
        if(am[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix\n");
            exit(EXIT_FAILURE);
        }

        am[i][i] = 0;   // distance from one place to itself is 0

        // assign random distances to all edges
        for(j=0; j<i; j++){
            am[i][j] = 1 + (int)((double)MAX_DISTANCE * rand() / (RAND_MAX + 1.0));
            if(oriented){
                // oriented graph, so i to j can have a different distance than j to i
                am[j][i] = 1 + (int)((double)MAX_DISTANCE * rand() / (RAND_MAX + 1.0));
            } else {
                // not oriented so i to j has the same distance as j to i
                am[j][i] = am[i][j];
            }

            if(am[i][j] == MAX_DISTANCE) m--;
            if(am[j][i] == MAX_DISTANCE) m--;
        }
    }

    *matrix = am;
    *num_edges = m;
}

/**
 * Read the adjacency matrix from a file.
 *
 * Returns: the number of edges that are "incorrect" in the file. That is, in case
 * the graph is not oriented, but there are different entries for symmetrical pairs
 * of edges, the second such edge is ignored, yet counted for statistics reasons.
 *
 * File structure:
 *  - first line: [num vertices] [num edges] [oriented (0/1)]
 *  - following [num edges] lines: [source] [destination] [distance]
 */
int read_adjacency_matrix(char *filename, int *num_vertices, int *num_edges, int ***matrix, int *oriented){
    int **am;
    int i, j;
    int n;
    int source, dest, dist;
    int num_bad_edges = 0;
    FILE *fp;

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("Error opening file '%s': %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d %d %d \n", num_vertices, num_edges, oriented);

    #ifdef VERBOSE
        printf("%d %d %d\n", *num_vertices, *num_edges, *oriented);
    #endif

    n = *num_vertices;
    am = (int **) malloc(n * sizeof(int *));
    if(am == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<n; i++){
        am[i] = (int *) malloc(n * sizeof(int));
        if(am[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix\n");
            exit(EXIT_FAILURE);
        }

        // init each distance to 0 or infinity (max distance in this case)
        for(j=0; j<n; j++){
            am[i][j] = (i == j) ? 0 : MAX_DISTANCE;
        }
    }

    while(!feof(fp)){
        fscanf(fp, "%d %d %d \n", &source, &dest, &dist);
        if(!(*oriented)){
            // distance i to j is the same as j to i
            if(am[source-1][dest-1] < MAX_DISTANCE){
                // verify that we did not set this distance yet.
                // if we did, we increment the number of bad edges, and ignore the distance
                num_bad_edges++;
            } else {
                am[source-1][dest-1] = dist;
                am[dest-1][source-1] = dist;
            }
        } else {
            am[source-1][dest-1] = dist;
        }
    }

    fclose(fp);

    #ifdef VERBOSE
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++){
                printf("%5d", am[i][j]);
            }
            printf("\n");
        }
    #endif

    *matrix = am;

    return num_bad_edges;
}

void usage() {
    printf ("Run the asp program with the following parameters. \n");
    printf (" -read filename :: reads the graph from a file.\n");
    printf (" -random N 0/1  :: generates a NxN graph, randomly. \n");
    printf ("                :: if 1, the graph is oriented, otherwise it is not oriented\n");
    return ;
}

int main(int argc, char **argv){
    int num_vertices = 0;
    int num_edges = 0;
    int i;                  // counter to loop through parameters
    int print = 0;          // print resulting adjacency matrix?
    int oriented = 0;       // whether or not the random matrix is oriented
    int **matrix;           // adjacency matrix
    int num_bad_edges = 0;  // number of bad edges in file
    char filename[100];     // filename if given as param



    // print usage
    usage();

    // read params
    for(i=1; i<argc; i++){
        if(strcmp(argv[i], "-print") == 0){
            print = 1;
        } else if(strcmp(argv[i], "-read") == 0){
            // assume filename is not longer than 100 chars
            strcpy(filename, argv[++i]);
        } else if(strcmp(argv[i], "-random") == 0){
            num_vertices = atoi(argv[++i]);
            oriented = atoi(argv[++i]);
        } else {
            num_vertices = 4000;
            oriented = 1;
        }
    }

    if(num_vertices > 0){
        // init random adjacency matrix for testing
        init_random_adjacency_matrix(num_vertices, &num_edges, &matrix, oriented);
    } else {
        // generate adjacency matrix from file
        num_bad_edges = read_adjacency_matrix(filename, &num_vertices, &num_edges, &matrix, &oriented);
    }

    fprintf(stderr, "Running ASP with %d rows and %d edges (%d are bad)\n", num_vertices, num_edges, num_bad_edges);

    return 0;
}
