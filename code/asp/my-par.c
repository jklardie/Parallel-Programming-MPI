#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <limits.h>
#include "mpi.h"

//#define VERBOSE

// use half of max integer value because we add two values later on (and we don't want overflows)
// note that the max_distance in the example (256) was incorrect, because some distances where > 256
#define INFINITY            INT_MAX/2
#define MAX_RANDOM_DIST     256
#define MASTER_PROC_ID      0
#define ROWS_TAG            10000
#define PARENT_ROWS_TAG     20000
#define NUM_ROWS_TAG        3
#define START_ROW_TAG       4
#define NUM_VERTICES_TAG    5

int my_proc_id;


/**
 * Broadcast a row from the adjacency matrix and the parent matrix to all nodes.
 *
 * Note that we don't send the size; all nodes already know the
 * size of the matrix (n).
 *
 * We also don't send the row index; all nodes know which index
 * they are receiving.
 */
void broadcast(int **row, int **parent_matrix_row, int proc_id, int n){
	int mpi_error; 
	
    if(proc_id != my_proc_id){
    	// I'm receiving, so allocate memory for the row
        int *r = (int *) malloc(n * sizeof(int));
        int *pr = (int *) malloc(n * sizeof(int));
    	
        // receive adjacency matrix row
		mpi_error = MPI_Bcast(r, n, MPI_INT, proc_id, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		// receive parent matrix row
		mpi_error = MPI_Bcast(pr, n, MPI_INT, proc_id, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		*row = r;
		*parent_matrix_row = pr;
    } else {
    	// broadcast adjacency matrix row
    	mpi_error = MPI_Bcast(*row, n, MPI_INT, proc_id, MPI_COMM_WORLD);
    	if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
    	
    	// broadcast parent matrix row
    	mpi_error = MPI_Bcast(*parent_matrix_row, n, MPI_INT, proc_id, MPI_COMM_WORLD);
    	if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
    }

}

/**
 * Apply the Floyd-Warshall algorithm to the adjacency matrix,
 * and keep track of paths in the parent matrix.
 */
void floyd_warshall(int **matrix, int **parent_matrix, int start_row, int end_row, int n, int num_procs){
    int i, j, k;
    int new_dist;
    int responsible_proc = 0;
    int avg_rows_per_proc = n / num_procs;

    for(k=0; k<n; k++){
		if(k > 0 && k % avg_rows_per_proc == 0 && responsible_proc < num_procs-1){
			responsible_proc++;
		}
		
        // broadcast row k in the adjacency matrix and parent matrix to others, or receive it
        broadcast(&(matrix[k]), &(parent_matrix[k]), responsible_proc, n);

        for(i=start_row; i<end_row; i++){
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

/**
 * Receive the calculated adjacency rows and parent matrix rows from slaves. 
 * 
 * Return: the diameter in the received rows (longest shortest path between two places)
 */
int receive_final_rows(int **matrix, int **parent_matrix, int num_procs, int n, int my_start, int my_end){
	int num_rows_to_receive, start_row, end_row, i, j, k;
	int mpi_error;
	int diameter = 0;
	
	// note start loop at i=1, because we don't need to receive it from ourselve
	for(i=1; i<num_procs; i++){
		// receive int telling how much rows we will receive
		mpi_error = MPI_Recv(&num_rows_to_receive, 1, MPI_INT, i, NUM_ROWS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		// receive int telling what our first row is we'll receive
		mpi_error = MPI_Recv(&start_row, 1, MPI_INT, i, START_ROW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		end_row = start_row + num_rows_to_receive;

		// receive the actual rows
		for(j=start_row; j<end_row; j++){
			// receive adjacency matrix row
			mpi_error = MPI_Recv(matrix[j], n, MPI_INT, i, ROWS_TAG+j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(mpi_error != MPI_SUCCESS){
				exit(EXIT_FAILURE);
			}
			
			// receive parent matrix row
			mpi_error = MPI_Recv(parent_matrix[j], n, MPI_INT, i, PARENT_ROWS_TAG+j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(mpi_error != MPI_SUCCESS){
				exit(EXIT_FAILURE);
			}
			
			for(k=0; k<n; k++){
				if(matrix[j][k] > diameter && matrix[j][k] < INFINITY){
					diameter = matrix[j][k];
				}
			}
		}
	}
	
	// calculate actual diameter (need to look at my own rows)
	for(i=my_start; i<my_end; i++){
	    for(j=0; j<n; j++){
            if(matrix[i][j] > diameter && matrix[i][j] < INFINITY){
                diameter = matrix[i][j];
            }
        }
	}

	return diameter;
}

/**
 * Send my adjacency rows and parent matrix rows to the master
 */
void send_rows_to_master(int **matrix, int **parent_matrix, int start_row, int end_row, int n){
    int i, num_rows_to_send, mpi_error;

	num_rows_to_send = end_row - start_row;

	// first send number of rows that will follow
	mpi_error = MPI_Send(&num_rows_to_send, 1, MPI_INT, MASTER_PROC_ID, NUM_ROWS_TAG, MPI_COMM_WORLD);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}

	// send the start_row to the worker
	mpi_error = MPI_Send(&start_row, 1, MPI_INT, MASTER_PROC_ID, START_ROW_TAG, MPI_COMM_WORLD);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
	
	// note: we don't have to send the number of vertices, the master already knows that

	// then send the actual rows
	for(i=start_row; i<end_row; i++){
		// send adjacency matrix row
		mpi_error = MPI_Send(matrix[i], n, MPI_INT, MASTER_PROC_ID, ROWS_TAG+i, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		// send parent matrix row
		mpi_error = MPI_Send(parent_matrix[i], n, MPI_INT, MASTER_PROC_ID, PARENT_ROWS_TAG+i, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
	}
}

/**
 * Receive part of the adjacency matrix from the master
 *
 * Return: the number of rows
 */
void receive_rows(int ***matrix, int ***parent_matrix, int *n, int *start_row, int *end_row){
    int num_rows_to_receive, i;
    int mpi_error;
    int **am;	// received adjacency matrix
    int **pm;	// received parent matrix

    // receive int telling how much rows we will receive
    mpi_error = MPI_Recv(&num_rows_to_receive, 1, MPI_INT, MASTER_PROC_ID, NUM_ROWS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}

    // receive int telling what our first row is we'll receive
    mpi_error = MPI_Recv(start_row, 1, MPI_INT, MASTER_PROC_ID, START_ROW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
	
    *end_row = *start_row + num_rows_to_receive;

    // receive int telling how many vertices there are
    mpi_error = MPI_Recv(n, 1, MPI_INT, MASTER_PROC_ID, NUM_VERTICES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
	
    am = (int **) malloc((*n) * sizeof(int *));
    pm = (int **) malloc((*n) * sizeof(int *));
    if(am == (int **)0 || am == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix (while receiving)\n");
        exit(EXIT_FAILURE);
    }

    // receive the actual rows
    for(i=(*start_row); i<(*end_row); i++){
        am[i] = (int *) malloc((*n) * sizeof(int));
        pm[i] = (int *) malloc((*n) * sizeof(int));

        if(am[i] == (int *)0 || pm[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix (while receiving)\n");
            exit(EXIT_FAILURE);
        }
        
		// receive adjacency matrix row
        mpi_error = MPI_Recv(am[i], *n, MPI_INT, MASTER_PROC_ID, ROWS_TAG+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
		
		// receive parent matrix row
		mpi_error = MPI_Recv(pm[i], *n, MPI_INT, MASTER_PROC_ID, PARENT_ROWS_TAG+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}
    }

    *matrix = am;
    *parent_matrix = pm;
}

/**
 * Distribute contiguous rows to processes. Note that this approach
 * might be improved by using a single array for the matrix
 * representation. That way, only one send can be done to send multiple
 * rows, because the memory is adjacent.
 *
 * Right now, with a multidimensional array, the arrays are not adjacent,
 * and therefore have to be sent independently.
 *
 */
void distribute_rows(int n, int num_procs, int **matrix, int **parent_matrix, int *start_row, int *end_row){
    int i, j, start, end, num_rows_to_send, mpi_error;
    int avg_rows_per_proc = n / num_procs;

    for(i=0; i<num_procs; i++){
        start = i * avg_rows_per_proc;
        end = (i+1) * avg_rows_per_proc;

        if((n - end) < avg_rows_per_proc){
            end = n;
        }

        if(i == 0){
            // for the master we only want to know the start and end rows
            *start_row = start;
            *end_row = end;
            continue;
        }

        num_rows_to_send = end - start;

        // first send number of rows that will follow
        mpi_error = MPI_Send(&num_rows_to_send, 1, MPI_INT, i, NUM_ROWS_TAG, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}

        // send the start_row to the worker
        mpi_error = MPI_Send(&start, 1, MPI_INT, i, START_ROW_TAG, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}

        // then send the number of ints per row (eg number of vertices)
        mpi_error = MPI_Send(&n, 1, MPI_INT, i, NUM_VERTICES_TAG, MPI_COMM_WORLD);
		if(mpi_error != MPI_SUCCESS){
			exit(EXIT_FAILURE);
		}

        // then send the actual rows
        for(j=start; j<end; j++){
        	// send adjacency matrix row
            mpi_error = MPI_Send(matrix[j], n, MPI_INT, i, ROWS_TAG+j, MPI_COMM_WORLD);
    		if(mpi_error != MPI_SUCCESS){
    			exit(EXIT_FAILURE);
    		}
    		
    		// send parent matrix row
    		mpi_error = MPI_Send(parent_matrix[j], n, MPI_INT, i, PARENT_ROWS_TAG+j, MPI_COMM_WORLD);
			if(mpi_error != MPI_SUCCESS){
				exit(EXIT_FAILURE);
			}
        }
    }
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
int read_adjacency_matrix(char *filename, int *n, int *num_edges, int ***matrix, int ***parent_matrix, int *oriented, int *num_bad_edges){
    int **am;
    int **pm;
    int i, j;
    int source, dest, dist;
    int total_distance = 0;

    *num_bad_edges = 0;
    FILE *fp;

    fp = fopen(filename, "r");
    if(fp == NULL){
        fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    if(fscanf(fp, "%d %d %d \n", n, num_edges, oriented) == EOF){
    	fprintf(stderr, "Error reading file '%s': %s\n", filename, strerror(errno));
    	exit(EXIT_FAILURE);
    }

    #ifdef VERBOSE
        printf("%d %d %d\n", *n, *num_edges, *oriented);
    #endif

    am = (int **) malloc((*n) * sizeof(int *));
    pm = (int **) malloc((*n) * sizeof(int *));
    if(am == (int **)0 || pm == (int **)0){
        fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<(*n); i++){
        am[i] = (int *) malloc((*n) * sizeof(int));
        pm[i] = (int *) malloc((*n) * sizeof(int));
        if(am[i] == (int *)0 || pm[i] == (int *)0){
            fprintf(stderr, "Error while allocating memory for adjacency matrix or parent matrix\n");
            exit(EXIT_FAILURE);
        }

        // init each distance to 0 or infinity (max distance in this case)
        // init each parent to -1.
        for(j=0; j<(*n); j++){
            am[i][j] = (i == j) ? 0 : INFINITY;
            pm[i][j] = -1;
        }
    }

    while(!feof(fp)){
        if(fscanf(fp, "%d %d %d \n", &source, &dest, &dist) == EOF){
			fprintf(stderr, "Error rows from file '%s': %s\n", filename, strerror(errno));
			exit(EXIT_FAILURE);
		}
        
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
    int num_vertices = 0;	// note: in other methods referred to as n
    int num_edges = 0;
    int i;                  // counter to loop through parameters and mpi node ids
    int print = 0;          // print resulting adjacency matrix?
    int oriented = 0;       // whether or not the random matrix is oriented
    int **matrix;           // adjacency matrix
    int **parent_matrix;    // matrix used to restore paths
    int num_bad_edges = 0;  // number of bad edges in file
    int total_distance = 0; // total distance of available roads
    int diameter = 0;       // diameter of the resulting adjacency matrix (longest shortest path)
    char filename[100];     // filename if given as param
    int mpi_error, num_procs;   // used for mpi communication
    int start_row = 0, end_row = 0;

    double wtime = 0;


    // initialize mpi, get my process id, get total number of processes
    mpi_error = MPI_Init(&argc, &argv);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
	
    mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &my_proc_id);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
	
    mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}


    // print usage
    if(my_proc_id == MASTER_PROC_ID){
        usage();
    }

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

    // let process 0 create the adjacency matrix, and send it to the slaves
    if(my_proc_id == MASTER_PROC_ID){
        // create matrix
        if(num_vertices > 0){
            total_distance = init_random_adjacency_matrix(num_vertices, &num_edges, &matrix, &parent_matrix, oriented);
        } else {
            total_distance = read_adjacency_matrix(filename, &num_vertices, &num_edges, &matrix, &parent_matrix, &oriented, &num_bad_edges);
        }

        fprintf(stderr, "Running ASP with %d rows and %d edges (%d are bad)\n", num_vertices, num_edges, num_bad_edges);

        // start timer
        wtime = MPI_Wtime();

        // distribute rows in matrix
        distribute_rows(num_vertices, num_procs, matrix, parent_matrix, &start_row, &end_row);
    } else {
        // slave, so wait until the rows from the matrix are received
        receive_rows(&matrix, &parent_matrix, &num_vertices, &start_row, &end_row);
    }

    // at this moment all workers have their working matrix in my_matrix.
    // this matrix contains <num_vertices> columns, and <num_rows> rows.

    floyd_warshall(matrix, parent_matrix, start_row, end_row, num_vertices, num_procs);
    
    if(my_proc_id == MASTER_PROC_ID){
    	// master needs to receive all other rows
    	diameter = receive_final_rows(matrix, parent_matrix, num_procs, num_vertices, start_row, end_row);
    } else {
    	// slave needs to send his rows to the master
    	send_rows_to_master(matrix, parent_matrix, start_row, end_row, num_vertices);
    }

    // let master print runtime, distance, diameter, and paths
    if(my_proc_id == MASTER_PROC_ID){
        wtime = MPI_Wtime() - wtime;

        fprintf(stderr, "Total distance: %d\n", total_distance);
        fprintf(stderr, "Diameter: %d\n", diameter);
        fprintf(stderr, "ASP took %10.3f seconds\n", wtime);

        // please note: printing paths might be very slow. Especially chain_10K.gr.new!
        if(print){
            print_paths(matrix, parent_matrix, num_vertices);
        }
    }

    free_matrix(matrix, num_vertices);

    if(my_proc_id == MASTER_PROC_ID){
        free_matrix(parent_matrix, num_vertices);
    }
    
    mpi_error = MPI_Finalize();
	if(mpi_error != MPI_SUCCESS){
		exit(EXIT_FAILURE);
	}
    
    return 0;
}
