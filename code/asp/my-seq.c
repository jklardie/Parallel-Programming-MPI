#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <limits.h>

//#define VERBOSE

// use half of max integer value because we add two values later on (and we don't want overflows)
// note that the max_distance in the example (256) was incorrect, because some distances where > 256
#define INFINITY INT_MAX/2

#define MAX_RANDOM_DIST 256

/**
 * Apply the Floyd-Warshall algorithm to the adjacency matrix,
 * and keep track of paths in the parent matrix.
 */
void floyd_warshall(int **matrix, int **parent_matrix, int n){
    int i, j, k;
    int new_dist;

    for(k=0; k<n; k++){
        for(i=0; i<n; i++){
            if(i != k){
                for(j=0; j<n; j++){
                    new_dist = matrix[i][k] + matrix[k][j];
                    if(new_dist < matrix[i][j]){
                        matrix[i][j] = new_dist;

                        // we just added a new node to the path between i and j (k),
                        // so add k as the parent of j. So to go from i to j, k will
                        // be visited directly before j.
                        parent_matrix[i][j] = parent_matrix[k][j];
                    }
                }
            }
        }
    }


}

int calculate_diameter(int **matrix, int n){
    int diameter = 0, i, j;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(matrix[i][j] > diameter && matrix[i][j] < INFINITY){
                diameter = matrix[i][j];
            }
        }
    }

    return diameter;
}

/**
 * Recursively print path between from i to j
 */
void print_path(int **parent_matrix, int i, int j){
    if(j == -1){
        return;
    }

    print_path(parent_matrix, i, parent_matrix[i][j]);
    printf("%2d, ", j+1);
}

/**
 * Prints paths between each combination of i and j, if there is any
 */
void print_paths(int **matrix, int **parent_matrix, int n){
    int i, j;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            // no path between i and j
            if(matrix[i][j] == INFINITY || i == j) continue;

            printf("%2d to %2d : %2d : ", i+1, j+1, matrix[i][j]);
            print_path(parent_matrix, i, j);

            printf("\n");
        }
    }
}

/**
 * Print a given n x n matrix
 */
void print_matrix(int **matrix, int n){
    int i, j;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%2d ", matrix[i][j]);
        }
        printf("\n");
    }
}

/**
 * Free all memory used by the n x n matrix
 */
void free_matrix(int **matrix, int n){
    int i;

    for(i=0; i<n; i++){
        free(matrix[i]);
    }

    free(matrix);
}

/**
 * Initialize a random adjacency matrix of size n x n.
 *
 * Return: total distance of usable roads
 */
int init_random_adjacency_matrix(int n, int *num_edges, int ***matrix, int ***parent_matrix, int oriented){
    int **am;   // adjacency matrix
    int **pm;   // parent matrix
    int i, j, m = n * n;
    int total_distance = 0;

    am = (int **) malloc(n * sizeof(int *));
    pm = (int **) malloc(n * sizeof(int *));
    if(am == (int **)0 || pm == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<n; i++){
        am[i] = (int *) malloc(n * sizeof(int));
        pm[i] = (int *) malloc(n * sizeof(int));
        if(am[i] == (int *)0 || pm[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
            exit(EXIT_FAILURE);
        }

        am[i][i] = 0;   // distance from one place to itself is 0

        // assign random distances to all edges
        for(j=0; j<i; j++){
            am[i][j] = 1 + (int)((double)MAX_RANDOM_DIST * rand() / (RAND_MAX + 1.0));

            if(am[i][j] < MAX_RANDOM_DIST){
                total_distance += am[i][j];
            }

            if(oriented){
                // oriented graph, so i to j can have a different distance than j to i
                am[j][i] = 1 + (int)((double)MAX_RANDOM_DIST * rand() / (RAND_MAX + 1.0));

                if(am[j][i] < MAX_RANDOM_DIST){
                    total_distance += am[j][i];
                }
            } else {
                // not oriented so i to j has the same distance as j to i
                am[j][i] = am[i][j];
            }

            // if a distance is the max distance, then assume the road does not exist,
            // and set the distance to infinity
            if(am[i][j] == MAX_RANDOM_DIST) {
                am[i][j] = INFINITY;
                m--;
            }
            if(am[j][i] == MAX_RANDOM_DIST) {
                am[j][i] = INFINITY;
                m--;
            }

            // store parents, used to restore path later on
            if(i == j || am[i][j] == INFINITY){
                pm[i][j] = -1;
            } else {
                pm[i][j] = i;
            }
        }
    }

    *matrix = am;
    *parent_matrix = pm;
    *num_edges = m;

    return total_distance;
}

/**
 * Read the adjacency matrix from a file.
 *
 *
 * Return: total distance of usable roads
 *
 * num_bad_edges: the number of edges that are "incorrect" in the file. That is, in case
 * the graph is not oriented, but there are different entries for symmetrical pairs
 * of edges, the second such edge is ignored, yet counted for statistics reasons.
 *
 * File structure:
 *  - first line: [num vertices] [num edges] [oriented (0/1)]
 *  - following [num edges] lines: [source] [destination] [distance]
 */
int read_adjacency_matrix(char *filename, int *num_vertices, int *num_edges, int ***matrix, int ***parent_matrix, int *oriented, int *num_bad_edges){
    int **am;
    int **pm;
    int i, j;
    int n;
    int source, dest, dist;
    int total_distance = 0;

    *num_bad_edges = 0;
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
    pm = (int **) malloc(n * sizeof(int *));
    if(am == (int **)0 || pm == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<n; i++){
        am[i] = (int *) malloc(n * sizeof(int));
        pm[i] = (int *) malloc(n * sizeof(int));
        if(am[i] == (int *)0 || pm[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
            exit(EXIT_FAILURE);
        }

        // init each distance to 0 or infinity (max distance in this case)
        // init each parent to -1.
        for(j=0; j<n; j++){
            am[i][j] = (i == j) ? 0 : INFINITY;
            pm[i][j] = -1;
        }
    }

    while(!feof(fp)){
        fscanf(fp, "%d %d %d \n", &source, &dest, &dist);
        i = source-1;
        j = dest-1;

        if(!(*oriented)){
            // distance i to j is the same as j to i
            if(am[i][j] < INFINITY){
                // verify that we did not set this distance yet.
                // if we did, we increment the number of bad edges, and ignore the distance
                (*num_bad_edges)++;
            } else {
                am[i][j] = dist;
                am[j][i] = dist;

                pm[i][j] = i;
                pm[j][i] = j;

                total_distance += dist;
            }
        } else {
            am[i][j] = dist;
            pm[i][j] = i;

            total_distance += dist;
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
    *parent_matrix = pm;

    return total_distance;
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
    int **parent_matrix;    // matrix used to restore paths
    int num_bad_edges = 0;  // number of bad edges in file
    int total_distance;
    int diameter;
    char filename[100];     // filename if given as param

    struct timeval start, end;
    double time;


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
        total_distance = init_random_adjacency_matrix(num_vertices, &num_edges, &matrix, &parent_matrix, oriented);
    } else {
        // generate adjacency matrix from file
        total_distance = read_adjacency_matrix(filename, &num_vertices, &num_edges, &matrix, &parent_matrix, &oriented, &num_bad_edges);
    }

    fprintf(stderr, "Running ASP with %d rows and %d edges (%d are bad)\n", num_vertices, num_edges, num_bad_edges);

    if(gettimeofday(&start, 0) != 0){
        fprintf(stderr, "Error starting timer\n");
        exit(EXIT_FAILURE);
    }

    floyd_warshall(matrix, parent_matrix, num_vertices);

    diameter = calculate_diameter(matrix, num_vertices);

    if(gettimeofday(&end, 0) != 0){
        fprintf(stderr, "Error stopping timer\n");
        exit(EXIT_FAILURE);
    }

    time = (end.tv_sec + end.tv_usec / 1000000.0) -
            (start.tv_sec + start.tv_usec / 1000000.0);


    fprintf(stderr, "Total distance: %d\n", total_distance);
    fprintf(stderr, "Diameter: %d\n", diameter);
    fprintf(stderr, "ASP took %10.3f seconds\n", time);

    if(print){
        print_paths(matrix, parent_matrix, num_vertices);
    }

    free_matrix(matrix, num_vertices);
    free_matrix(parent_matrix, num_vertices);

    return 0;
}
