/*
	This code is a collection of fucntions
	for calculating the standard score of motifs.

	It was written by Youngjai in Nov. 19, 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define Motifs 27

typedef struct __list {
	struct __node *head;
	struct __node *tail;
} Linked_list;

typedef struct __node {
	int node_label;
	int degree;
	int *neighbor;
	int *direction;
	double *sign;
	double fitness;
	struct __node *next;
} Node;

typedef struct {
	int nodes;
	int links;
	Linked_list *node_list;
} Net_info;

typedef struct {
	Net_info net_info;
	int node_label;
} Load_profile;

/* Load some packages from the 'mt19937.c'. */
void init_genrand(unsigned long s);

/* Load some packages from the 'youngjai_pacakeges.c'. */
void swap_two_integers(int *num_i, int *num_j);
double normal_dist(double avg, double stdev);
clock_t time_begin(char *path);
void time_end(clock_t begin, char *path);

/* Load some packages from the 'measure_motif.c'. */
double *calculate_motifs(Linked_list *node_list, double *motif);
double *measure_z_score(Net_info net_info, int node_index, double *motif, int ensemble, int rand_level);
void calculate_fitness_all(Linked_list *node_list);
void print_network_profile(char *path, Net_info net_info, int rand_level);

/* Load some packages from the 'make_networks.c'. */
Load_profile load_networks(char *filename);
Net_info test_network(int nodes, int links);
void print_test_temp(Net_info net_info);
void free_Linked_list(Linked_list *node_list);

void print_motif(double *motif, char *path);
void print_z_score(double *z_score, char *path, int rand_level);

int main(int argc, char* argv[])
{
	if(argc != 4) {
		printf("Error! wrong number of arguments.\r\n");
		return 1;
	}
	int ensemble = atoi(argv[1]);
	char *path = malloc(strlen(argv[2])+1);	strcpy(path, argv[2]);
	char *filename = malloc(strlen(argv[3])+1);	strcpy(filename, argv[3]);
	char *str = malloc(strlen(path)+10);	strcpy(str, path);
	mkdir(path, 0755);

	clock_t begin = time_begin(strcat(str, "log.dat"));

	// unsigned long seed = time(NULL);	init_genrand(seed);
	unsigned long seed = 1;	init_genrand(seed);
	strcpy(str, path);
	FILE *log = fopen(strcat(str, "log.dat"), "at");
	fprintf(log, "random seed : %ld\r\n", seed);
	fclose(log);

	int i;
	Load_profile load_profile;

	/* Load the network information. */
	load_profile = load_networks(filename);
	// load_profile.net_info = test_network(10, 20); load_profile.node_label = load_profile.net_info.nodes;
	// print_network_profile(path, load_profile.net_info, -1);
	double *motif = malloc(sizeof(double) * (Motifs+1));
	motif = calculate_motifs(load_profile.net_info.node_list, motif);
	print_motif(motif, path);
	// print_test_temp(load_profile.net_info);

	strcpy(str, path);
	log = fopen(strcat(str, "log.dat"), "at");
	fprintf(log, "node %d\tlinks %d\r\n", load_profile.net_info.nodes, load_profile.net_info.links);
	fclose(log);

	/* Calculate the standard score (z score) comparing with null model network. */
	double *z_score = malloc(sizeof(double) * (Motifs+1));
	for(i=0; i<4; i++){
		z_score = measure_z_score(load_profile.net_info, load_profile.node_label, motif, ensemble, i);
		print_z_score(z_score, path, i);
	}

	free_Linked_list(load_profile.net_info.node_list);
	free(filename); free(motif); free(z_score);
	strcpy(str, path);
	time_end(begin, strcat(str, "log.dat"));
	free(path); free(str);

	return 0;
}

void print_motif(double *motif, char *path)
{
	int i;
	double all_tripets = 0.;
	double closed_tripets = 0.;
	double gcc;	// global clustering coefficient
	char *str = malloc(strlen(path)+100); strcpy(str, path);

	sprintf(str, "%smotif.dat", str);
	FILE *motif_dist = fopen(str, "wt"); free(str);

	for(i=0; i<22; i++){
		if(i > 10) closed_tripets += motif[i];
		all_tripets += motif[i];
	}
	if(all_tripets == 0) gcc = 0.;
	else gcc = closed_tripets/all_tripets;
	
	for(i=1; i<=10; i++) fprintf(motif_dist, "M(3,2)-%d\t%.0lf\r\n", i, motif[i-1]);
	for(i=11; i<=22; i++) fprintf(motif_dist, "M(3,3)-%d\t%.0lf\r\n", i, motif[i-1]);
	for(i=23; i<=Motifs; i++) fprintf(motif_dist, "M(unsigned)-%d\t%.0lf\r\n", i-22, motif[i-1]);
	fprintf(motif_dist, "gcc(global clustering coefficient)\t%lf\r\n", gcc);
	fclose(motif_dist);
}

void print_z_score(double *z_score, char *path, int rand_level)
{
	int i;
	char *str = malloc(strlen(path)+100); strcpy(str, path);

	sprintf(str, "%sz_score_%d.dat", str, rand_level);
	FILE *z_score_dist = fopen(str, "wt"); free(str);

	for(i=1; i<=10; i++) fprintf(z_score_dist, "M(3,2)-%d\t%.4lf\r\n", i, z_score[i-1]);
	for(i=11; i<=22; i++) fprintf(z_score_dist, "M(3,3)-%d\t%.4lf\r\n", i, z_score[i-1]);
	for(i=23; i<=Motifs; i++) fprintf(z_score_dist, "M(unsigned)-%d\t%.4lf\r\n", i-22, z_score[i-1]);
	fclose(z_score_dist);
}
