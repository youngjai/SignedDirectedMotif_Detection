/*
	This is a collection of fucntions for coding.

	It was written by Youngjai in Nov. 2, 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// from mt19937ar.c
void init_genrand(unsigned long s);
double genrand_real2(void);

clock_t time_begin(char *path)
{	
	clock_t begin = clock();
	FILE *log = fopen(path, "wt");
	fprintf(log, "Log file is opened by youngjai.\r\n");
	fclose(log);

	return begin;
}

void time_end(clock_t begin, char *path)
{
	FILE *log = fopen(path, "at");
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; /* unit : second */
	int dd, hh, mm;
	double ss;
	dd = time_spent/86400;
	hh = (time_spent-dd*86400)/3600;
	mm = (time_spent-dd*86400-hh*3600)/60;
	ss = time_spent-dd*86400-hh*3600-mm*60;

	fprintf(log, "\r\nThe runing time is %d days %02d hours %02d minutes %02.4lf seconds.\r\n", dd, hh, mm, ss);
	fclose(log);
}

void swap_two_integers(int *num_i, int *num_j)
{
	int temp = *num_i;	/* temporary storage */
	*num_i = *num_j; 
	*num_j = temp;
}

void swap_two_floats(double *num_i, double *num_j)
{
	double temp = *num_i;	/* temporary storage */
	*num_i = *num_j; 
	*num_j = temp;
}

int call = 0;
double normal_dist(double avg, double stdev, unsigned long seed)
{
	init_genrand(seed);
	/* generate a number from the normal distribution. */
	double v1, v2, s, temp;
	static double x1, x2;
	// static int call = 0;

	if(call == 1){
		call = !call;
		return (avg + stdev * (double) x2);
	}

	do{
		v1 = 2. * genrand_real2() - 1.;
		v2 = 2. * genrand_real2() - 1.;
		s = v1 * v1 + v2 * v2;
	} while(s >= 1. || s == 0.);

	temp = sqrt((-2. * log(s)) / s);
	x1 = v1 * temp;
	x2 = v2 * temp;

	call = !call;

	return (avg + stdev * (double) x1);
}

int *ascending_bubble_sorting_int(int *array, int length)
{
	int i, j;

	for(i=0; i<length; i++) for(j=i+1; j<length; j++)
		if(array[j] < array[i])
			swap_two_integers(&array[j], &array[i]);

	return array;
}

int *shuffle_array_int(int length, unsigned long seed)
{
	init_genrand(seed);
	int i, j;
	int *shuffle = malloc(sizeof(int)*(length+1));
	for(i=0; i<length; i++) shuffle[i] = i;

	for(i=length-1; i>0; i--){
		j = (int) ((i+1)*genrand_real2());
		swap_two_integers(&shuffle[i], &shuffle[j]);
	}
	return shuffle;
}

double Euclidean_distance(int dx, int dy, int dz)
{
	return sqrt(dx*dx + dy*dy + dz*dz);
}

int ***initialize_arr3d_int(int *size, int initial_value)
{
	int i, j, k;

	// size[0] is the range of x axis.
	// size[1] is the range of y axis.
	// size[2] is the range of z axis.
	int ***arr3d = malloc(sizeof(int **) * size[0]);
	for(i=0; i<size[0]; i++){
		arr3d[i] = malloc(sizeof(int *) * size[1]);
		for(j=0; j<size[1]; j++){
			arr3d[i][j] = malloc(sizeof(int) * size[2]);
			for(k=0; k<size[2]; k++) arr3d[i][j][k] = initial_value;
		}
	}
	return arr3d;
}

double ***initialize_arr3d_float(int *size, double initial_value)
{
	int i, j, k;

	// size[0] is the range of x axis.
	// size[1] is the range of y axis.
	// size[2] is the range of z axis.
	double ***arr3d = malloc(sizeof(double **) * size[0]);
	for(i=0; i<size[0]; i++){
		arr3d[i] = malloc(sizeof(double *) * size[1]);
		for(j=0; j<size[1]; j++){
			arr3d[i][j] = malloc(sizeof(double) * size[2]);
			for(k=0; k<size[2]; k++) arr3d[i][j][k] = initial_value;
		}
	}
	return arr3d;
}

void free_arr3d_int(int ***arr3d, int *size)
{
	int i, j;

	// size[0] is the range of x axis.
	// size[1] is the range of y axis.
	// size[2] is the range of z axis.
	for(i=0; i<size[0]; i++){
		for(j=0; j<size[1]; j++) free(arr3d[i][j]);
		free(arr3d[i]);
	} free(arr3d);
}

void free_arr3d_float(double ***arr3d, int *size)
{
	int i, j;

	// size[0] is the range of x axis.
	// size[1] is the range of y axis.
	// size[2] is the range of z axis.
	for(i=0; i<size[0]; i++){
		for(j=0; j<size[1]; j++) free(arr3d[i][j]);
		free(arr3d[i]);
	} free(arr3d);
}

void gcc_3D_square_lattice(int ***arr3d, int from, int to, int *size)	// 3-Dimension
{
	// size[0] is the range of x axis.
	// size[1] is the range of y axis.
	// size[2] is the range of z axis.
	int i, j, k;
	int x, y, z;
	int gnum, key, csize;
	int gcs=0;

	int ***compare = initialize_arr3d_int(size, 0);
	// Cluster Grouping
	// X and Y are applied for the periodic boundary conditions.
	gnum=1;
	for(i=0; i<size[0]; i++) for(j=0; j<size[1]; j++) for(k=0; k<size[2]; k++){
		if(arr3d[i][j][k] == from){
			if(i==0 && j==0 && k==0){			// First site
				compare[i][j][k] = gnum++;
			}
			else if(j==0 && k==0){					// First x axis
				if(arr3d[i-1][j][k] == from)
					compare[i][j][k] = compare[i-1][j][k];
				else compare[i][j][k] = gnum++;
			}
			else if(k==0 && i==0){					// First y axis
				if(arr3d[i][j-1][k] == from)
					compare[i][j][k] = compare[i][j-1][k];
				else compare[i][j][k] = gnum++;
			}
			else if(i==0 && j==0){					// First z axis
				if(arr3d[i][j][k-1] == from)
					compare[i][j][k] = compare[i][j][k-1];
				else compare[i][j][k] = gnum++;
			}
			else if(k == 0){						// First xy plane
				if(arr3d[i-1][j][k]==from && arr3d[i][j-1][k]==from){
					if(compare[i-1][j][k] > compare[i][j-1][k]){
						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
					else if(compare[i-1][j][k] == compare[i][j-1][k])
						compare[i][j][k] = compare[i-1][j][k];
					else{
						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
				}
				else if(arr3d[i-1][j][k] == from)
					compare[i][j][k] = compare[i-1][j][k];
				else if(arr3d[i][j-1][k] == from)
					compare[i][j][k] = compare[i][j-1][k];
				else compare[i][j][k] = gnum++;
			}
			else if(i == 0){						// First yz plane
				if(arr3d[i][j-1][k]==from && arr3d[i][j][k-1]==from){
					if(compare[i][j-1][k] > compare[i][j][k-1]){
						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
					else if(compare[i][j-1][k] == compare[i][j][k-1])
						compare[i][j][k] = compare[i][j-1][k];
					else{
						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
				}
				else if(arr3d[i][j-1][k] == from)
					compare[i][j][k] = compare[i][j-1][k];
				else if(arr3d[i][j][k-1] == from)
					compare[i][j][k] = compare[i][j][k-1];
				else compare[i][j][k] = gnum++;
			}
			else if(j == 0){						// First zx plane
				if(arr3d[i][j][k-1]==from && arr3d[i-1][j][k]==from){
					if(compare[i][j][k-1] > compare[i-1][j][k]){
						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
					else if(compare[i][j][k-1] == compare[i-1][j][k])
						compare[i][j][k] = compare[i][j][k-1];
					else{
						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
				}
				else if(arr3d[i][j][k-1] == from)
					compare[i][j][k] = compare[i][j][k-1];
				else if(arr3d[i-1][j][k] == from)
					compare[i][j][k] = compare[i-1][j][k];
				else compare[i][j][k] = gnum++;
			}
			else{									// Generally
				if(arr3d[i-1][j][k]==from && \
					arr3d[i][j-1][k]==from && \
					arr3d[i][j][k-1]==from){

					if(compare[i][j-1][k]>compare[i-1][j][k] && \
						compare[i][j][k-1]>compare[i-1][j][k]){

						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
					else if(compare[i][j][k-1]>compare[i][j-1][k] && \
						compare[i-1][j][k]>compare[i][j-1][k]){

						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
					else if(compare[i-1][j][k]>compare[i][j][k-1] && \
						compare[i][j-1][k]>compare[i][j][k-1]){

						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
					else if(compare[i-1][j][k]==compare[i][j-1][k] && \
						compare[i][j][k-1]>compare[i-1][j][k]){

						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
					else if(compare[i][j-1][k]==compare[i][j][k-1] && \
						compare[i-1][j][k]>compare[i][j-1][k]){

						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
					else if(compare[i][j][k-1]==compare[i-1][j][k] && \
						compare[i][j-1][k]>compare[i][j][k-1]){

						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
					else if(compare[i-1][j][k]==compare[i][j-1][k] && \
						compare[i-1][j][k]==compare[i][j][k-1])

						compare[i][j][k] = compare[i-1][j][k];
				}
				else if(arr3d[i-1][j][k]==from && arr3d[i][j-1][k]==from){
					if(compare[i-1][j][k] > compare[i][j-1][k]){
						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
					else if(compare[i-1][j][k] == compare[i][j-1][k])
						compare[i][j][k] = compare[i-1][j][k];
					else{
						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
				}
				else if(arr3d[i][j-1][k]==from && arr3d[i][j][k-1]==from){
					if(compare[i][j-1][k] > compare[i][j][k-1]){
						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i][j-1][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
					else if(compare[i][j-1][k] == compare[i][j][k-1])
						compare[i][j][k] = compare[i][j-1][k];
					else{
						compare[i][j][k] = compare[i][j-1][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j-1][k];
					}
				}
				else if(arr3d[i][j][k-1]==from && arr3d[i-1][j][k]==from){
					if(compare[i][j][k-1] > compare[i-1][j][k]){
						compare[i][j][k] = compare[i-1][j][k];
						key = compare[i][j][k-1];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i-1][j][k];
					}
					else if(compare[i][j][k-1] == compare[i-1][j][k])
						compare[i][j][k] = compare[i][j][k-1];
					else{
						compare[i][j][k] = compare[i][j][k-1];
						key = compare[i-1][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k-1];
					}
				}
				else if(arr3d[i-1][j][k]==from)
					compare[i][j][k] = compare[i-1][j][k];
				else if(arr3d[i][j-1][k]==from)
					compare[i][j][k] = compare[i][j-1][k];
				else if(arr3d[i][j][k-1]==from)
					compare[i][j][k] = compare[i][j][k-1];
				else compare[i][j][k] = gnum++;
			}

			// X and Y are applied for the periodic boundary conditions.
			if(i==size[0]-1 && j==size[1]-1){
				if(arr3d[0][j][k]==from && arr3d[i][0][k]==from){
					if(compare[i][0][k]>compare[i][j][k] && \
						compare[i][0][k]>compare[0][j][k]){

						key = compare[i][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][0][k];
						key = compare[0][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][0][k];
					}
					else if(compare[0][j][k]>compare[i][j][k] && \
						compare[0][j][k]>compare[i][0][k]){

						key = compare[i][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[0][j][k];
						key = compare[i][0][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[0][j][k];
					}
					else if(compare[i][j][k]>compare[i][0][k] && \
						compare[i][j][k]>compare[0][j][k]){

						key = compare[i][0][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k];
						key = compare[0][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][j][k];
					}
				}
			}
			else if(i == size[0]-1){
				if(arr3d[0][j][k] == from){
					if(compare[i][j][k] > compare[0][j][k]){
						key = compare[i][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[0][j][k];
					}
					else if(compare[i][j][k] == compare[0][j][k]);
					else{
						key = compare[0][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[0][j][k];
					}
				}
			}
			else if(j == size[1]-1){
				if(arr3d[i][0][k] == from){
					if(compare[i][j][k] > compare[i][0][k]){
						key = compare[i][j][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][0][k];
					}
					else if(compare[i][j][k] == compare[i][0][k]);
					else{
						key = compare[i][0][k];
						for(x=0; x<size[0]; x++) for(y=0; y<size[1]; y++) for(z=0; z<size[2]; z++)
							if(compare[x][y][z] == key)
								compare[x][y][z] = compare[i][0][k];
					}
				}
			}
		}
	}
	
	// Cluster Profiling
	for(x=1; x<gnum+1; x++){
		csize=0;
		for(i=0; i<size[0]; i++) for(j=0; j<size[1]; j++) for(k=0; k<size[2]; k++)
			if(compare[i][j][k] == x) csize++;

		if(csize){
			if(csize > gcs){
				gcs = csize;
				key = x;
			}
		}
	}

	for(i=0; i<size[0]; i++) for(j=0; j<size[1]; j++) for(k=0; k<size[2]; k++)
		if(compare[i][j][k] == key) arr3d[i][j][k] = to;
		else arr3d[i][j][k] = from;

	free_arr3d_int(compare, size);
}

void gcc_2D_square_lattice(int **arr2d, int from, int to, int x_range, int y_range)
{
	int i, j, x, y;
	int gnum, key, key1;
	int gcs = 0;
	int csize;

	int comp[x_range][y_range];
	for(i=0; i<x_range; i++) for(j=0; j<y_range; j++) comp[i][j] = 0;

	// Cluster Grouping
	gnum=1;
	for(i=0; i<x_range; i++){
		for(j=0; j<y_range; j++){
			if(arr2d[i][j] == from){
				if(i==0 && j==0) comp[i][j] = gnum++;
				else if(i==0)   //First line
					if(arr2d[i][j-1] == from) comp[i][j] = comp[i][j-1];
					else comp[i][j] = gnum++;
				else if(j==0)   //First low
					if(arr2d[i-1][j] == from) comp[i][j] = comp[i-1][j];
					else comp[i][j] = gnum++;

				else{      //Generally
					//case that 'left and upper site is occufied'
					if(arr2d[i][j-1]  == from && arr2d[i-1][j] == from){
						if(comp[i][j-1] > comp[i-1][j]){
							comp[i][j] = comp[i-1][j];
							key = comp[i][j-1];
							for(x=0; x<x_range; x++) for(y=0; y<y_range; y++)
								if(comp[x][y] == key) comp[x][y] = comp[i-1][j];
						}
						else if(comp[i][j-1] == comp[i-1][j])
							comp[i][j] = comp[i-1][j];
						else{
							comp[i][j] = comp[i][j-1];
							key = comp[i-1][j];
							for(x=0; x<x_range; x++) for(y=0; y<y_range; y++)
								if(comp[x][y] == key) comp[x][y] = comp[i][j-1];
						}   
					}
					else if(arr2d[i][j-1] == from)
						comp[i][j] = comp[i][j-1];   //case that 'left site is occufied'
					else if(arr2d[i-1][j] == from)
						comp[i][j] = comp[i-1][j];   //case that 'upper site is occufied'
					else comp[i][j] = gnum++;
				}

				//periodic boundary only x axis
				if(i==x_range-1 && comp[0][j]){
					if(comp[i][j] > comp[0][j]){
						key = comp[i][j];
						for(x=0; x<x_range; x++) for(y=0; y<y_range; y++)
							if(comp[x][y]==key) comp[x][y] = comp[0][j];
					}
					else if(comp[i][j] == comp[0][j]);
					else{
						key = comp[0][j];
						for(x=0; x<x_range; x++) for(y=0; y<y_range; y++)
							if(comp[x][y]==key) comp[x][y] = comp[i][j];
					}
				}
			}
		}
	}
	
	// Cluster Profiling
	for(x=1; x<gnum+1; x++){
		csize=0;
		for(i=0; i<x_range; i++) for(j=0; j<y_range; j++)
			if(comp[i][j] == x) csize++;
		if(csize && csize>gcs){
			gcs = csize;
			key = x;
		}
	}

	for(i=0; i<x_range; i++) for(j=0; j<y_range; j++)
		if(comp[i][j] != key) arr2d[i][j] = to;
}

double RK4(double x, double t, double h, double *values, double (*f)(double, double, double*))
{
	double h2 = 0.5*h;
	double k1 = f(x, t, values);
	double k2 = f(x+h2*k1, t+h2, values);
	double k3 = f(x+h2*k2, t+h2, values);
	double k4 = f(x+h*k3, t+h, values);

	return x + h/6.*(k1+2.*k2+2.*k3+k4);
}

double Adaptive_RK4_h(double x, double t, double h, double *values, double (*f)(double, double, double*))
{
	double x_next1, x_next2;
	double rho = 0.;
	double accuracy = 1e-12; // target accuracy
	// If the difference between two results is less than 'precision', these are same.
	double precision = 1e-9; // 소수 9번째 자리까지 정도 확인

	while(rho<1.){
		x_next1 = RK4(RK4(x, t, h, values, f), t, h, values, f);
		x_next2 = RK4(x, t, 2.*h, values, f);

		if(fabs(x_next1-x_next2) < precision){
			h *= 2.;
			rho = 10.;
		}
		else{
			rho = 30.*h*accuracy/fabs(x_next1-x_next2);

			if(rho<1.) h *= pow(rho, 0.25);
			else h *= 2.;
		}
	}

	return h;
}