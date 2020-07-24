/*
	This code is a collection of fucntions
	for calculating the standard score of motifs.

	It was written by Youngjai in Nov. 19, 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

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
	int *ref_array;
	int *ref_links;
	int **ref_indices;
} Ref_info;

/* Load some packages from the 'youngjai_pacakeges.c'. */
void swap_two_integers(int *num_i, int *num_j);
double normal_dist(double avg, double stdev);

/* Load some packages from the 'mt19937.c'. */
double genrand_real2(void);

/* Load some packages from the 'make_network.c'. */
void free_Linked_list(Linked_list *node_list);
void create_node(Linked_list *node_list, int node_label);
Ref_info reference_links(Net_info net_info, int node_index, int rand_level);
Linked_list *make_ER_network(Net_info net_info, int ref_links);
Linked_list *make_rewired_network(Net_info net_info, Ref_info ref_info, int node_index, int rand_level);
void print_test_temp(Net_info net_info);

void calculate_fitness_all(Linked_list *node_list)
{
	int i;
	double fitness;
	Node *cursor = node_list->head;

	while(cursor->next != NULL){
		cursor = cursor->next;
		fitness = 0.;
		for(i=0; i<cursor->degree; i++) if(cursor->direction[i] < 0) fitness += cursor->sign[i];
		cursor->fitness = fitness;
	}
}
void calculate_fitness(Linked_list *node_list, int name_i)
{
	int i;
	double fitness = 0.;
	Node *cursor = node_list->head;

	while(cursor->node_label != name_i) cursor = cursor->next;
	for(i=0; i<cursor->degree; i++) if(cursor->direction[i] < 0) fitness += cursor->sign[i];
	cursor->fitness = fitness;
}

void print_network_profile(char *path, Net_info net_info, int rand_level)
{
	int i;
	int incoming, outgoing;
	double fitness;
	Node *cursor = net_info.node_list->head;
	char *str = malloc(strlen(path)+100);
	if(rand_level < 0) sprintf(str, "%sdegree_profile.dat", path);
	else sprintf(str, "%sdegree_%d_profile.dat", path, rand_level);
	FILE *degree_profile = fopen(str, "at"); free(str);

	while(cursor->next != NULL){
		cursor = cursor->next;
		incoming = 0; outgoing = 0; fitness = 0.;
		for(i=0; i<cursor->degree; i++)
			if(cursor->direction[i] < 0){
				incoming++; fitness += cursor->sign[i];
			}
			else outgoing++;
		cursor->fitness = fitness;
		fprintf(degree_profile, "%d\t%d\t%d\t%d\t%d\t%lf\r\n", \
			net_info.nodes, net_info.links, incoming+outgoing, incoming, outgoing, fitness);
	}
	fclose(degree_profile);
}

Net_info deep_copy_network(Net_info to, Net_info from)
{
	int i;

	to.node_list = malloc(sizeof(Linked_list));
	to.node_list->head = NULL; to.node_list->tail = NULL;
	create_node(to.node_list, 0);
	for(i=0; i<from.nodes; i++) create_node(to.node_list, 0);

	Node *cursor, *cursor1;
	cursor = from.node_list->head;
	cursor1 = to.node_list->head;
	
	while(cursor->next != NULL){
		cursor = cursor->next;
		cursor1 = cursor1->next;

		cursor1->node_label = cursor->node_label;
		cursor1->degree = cursor->degree;

		cursor1->neighbor = malloc(sizeof(int) * (cursor->degree+1));
		cursor1->direction = malloc(sizeof(int) * (cursor->degree+1));
		cursor1->sign = malloc(sizeof(double) * (cursor->degree+1));
		memcpy(cursor1->neighbor, cursor->neighbor, sizeof(int) * (cursor->degree+1));
		memcpy(cursor1->direction, cursor->direction, sizeof(int) * (cursor->degree+1));
		memcpy(cursor1->sign, cursor->sign, sizeof(double) * (cursor->degree+1));

		cursor1->fitness = cursor->fitness;
	}
	to.nodes = from.nodes;
	to.links = from.links;

	return to;
}	

double *calculate_motifs(Linked_list *node_list, double *motif)
{
	int i, j, k;
	int count;

	for(i=0; i<Motifs; i++) motif[i] = 0;

	Node *cursor = node_list->head;
	Node *cursor1, *cursor2;

	while(cursor->next != NULL){
		cursor = cursor->next;
		for(i=0; i<cursor->degree; i++){
			if(cursor->neighbor[i] != cursor->node_label){ // Avoid self-loops.
				cursor1 = node_list->head;
				while(cursor1->node_label != cursor->neighbor[i]) cursor1 = cursor1->next;
				for(j=0; j<cursor1->degree; j++){
					// Avoid mult-edges and self-loops.
					if(cursor1->neighbor[j] != cursor->node_label && cursor1->neighbor[j] != cursor1->node_label){
						cursor2 = node_list->head; count = 0;
						while(cursor2->node_label != cursor1->neighbor[j]) cursor2 = cursor2->next;
						for(k=0; k<cursor2->degree; k++){
							if(cursor2->neighbor[k] == cursor->node_label){
								/* M(3,3) motifs (11~22) */
								count++;

								if(cursor->direction[i] > 0){
									if(cursor->sign[i] > 0){
										if(cursor1->direction[j] > 0){
											if(cursor1->sign[j] > 0){
												if(cursor2->direction[k] > 0){
													if(cursor2->sign[k] > 0) motif[10] += 2.;
													else motif[11] += 6.;
												}
												else{
													if(cursor2->sign[k] > 0) motif[12] += 6.;
													else motif[13] += 6.;
												}
											}
											else{
												if(cursor2->direction[k] > 0){
													if(cursor2->sign[k] < 0) motif[14] += 6.;
												}
												else{
													if(cursor2->sign[k] > 0) motif[15] += 6.;
													else motif[16] += 6.;
												}
											}
										}
										else{
											if(cursor1->sign[j] > 0){
												if(cursor2->direction[k] < 0){
													if(cursor2->sign[k] < 0) motif[17] += 6.;
												}
											}
											else{
												if(cursor2->direction[k] > 0){
													if(cursor2->sign[k] < 0) motif[18] += 6.;
												}
												else{
													if(cursor2->sign[k] < 0) motif[19] += 6.;
												}
											}
										}
									}
									else{
										if(cursor1->direction[j] > 0){
											if(cursor1->sign[j] < 0){
												if(cursor2->direction[k] > 0){
													if(cursor2->sign[k] < 0) motif[20] += 2.;
												}
												else{
													if(cursor2->sign[k] < 0) motif[21] += 6.;
												}
											}
										}
									}
								}
							}
						}
						/* M(3,2) motifs (1~10) */
						if(count == 0){
							if(cursor->direction[i] > 0){
								if(cursor->sign[i] > 0){
									if(cursor1->direction[j] > 0){
										if(cursor1->sign[j] > 0) motif[0] += 6.;
										else motif[1] += 6.;
									}
									else{
										if(cursor1->sign[j] > 0) motif[2] += 3.;
										else motif[3] += 6.;
									}
								}
								else{
									if(cursor1->direction[j] > 0){
										if(cursor1->sign[j] > 0) motif[4] += 6.;
										else motif[5] += 6.;
									}
									else{
										if(cursor1->sign[j] < 0) motif[6] += 3.;
									}
								}
							}
							else{
								if(cursor->sign[i] > 0){
									if(cursor1->direction[j] > 0){
										if(cursor1->sign[j] > 0) motif[7] += 3.;
										else motif[8] += 6.;
									}
								}
								else{
									if(cursor1->direction[j] > 0){
										if(cursor1->sign[j] < 0) motif[9] += 3.;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	motif[22] = motif[0]+motif[1]+motif[4]+motif[5];
	motif[23] = motif[2]+motif[3]+motif[6];
	motif[24] = motif[7]+motif[8]+motif[9];
	motif[25] = motif[10]+motif[11]+motif[14]+motif[20];
	motif[26] = motif[12]+motif[13]+motif[15]+motif[16]+motif[17]+motif[18]+motif[19]+motif[21];
	for(i=0; i<Motifs; i++) motif[i] /= 6.;

	return motif;
}

double *measure_z_score(Net_info net_info, int node_index, double *motif, int ensemble, int rand_level)
{
	Net_info net_info_rewired;
	net_info_rewired = deep_copy_network(net_info_rewired, net_info);

	double *z_score = malloc(sizeof(double) * (Motifs+1));
	double *motif_rewired = malloc(sizeof(double) * (Motifs+1));

	int i, j, ens;
	double avg[Motifs], stdev[Motifs];

	for(i=0; i<Motifs; i++){
		avg[i] = 0.; stdev[i] = 0.;
	}

	Ref_info ref_info = reference_links(net_info_rewired, node_index, rand_level);
	if((rand_level == 3) && (ref_info.ref_links[0] < 2 || ref_info.ref_links[1] < 2 || net_info_rewired.nodes < 4)) rand_level = 0;
	else if((rand_level != 3) && (ref_info.ref_links[0]+ref_info.ref_links[1] < 2 || net_info_rewired.nodes < 4)) rand_level = 0;

	if(rand_level == 0){
		for(ens=0; ens<ensemble; ens++){
			net_info_rewired.node_list = make_ER_network(net_info_rewired, ref_info.ref_links[0]);
			motif_rewired = calculate_motifs(net_info_rewired.node_list, motif_rewired);
			free_Linked_list(net_info_rewired.node_list);

			for(i=0; i<Motifs; i++){
				avg[i] += motif_rewired[i];
				stdev[i] += motif_rewired[i]*motif_rewired[i];
			}
		}
	}
	else{
		for(ens=0; ens<ensemble; ens++){
			net_info_rewired.node_list = make_rewired_network(net_info_rewired, ref_info, node_index, rand_level);
			motif_rewired = calculate_motifs(net_info_rewired.node_list, motif_rewired);

			for(i=0; i<Motifs; i++){
				avg[i] += motif_rewired[i];
				stdev[i] += motif_rewired[i]*motif_rewired[i];
			}
		}
		free_Linked_list(net_info_rewired.node_list);
	}

	free(ref_info.ref_array); free(ref_info.ref_links);
	for(i=0; i<3; i++) free(ref_info.ref_indices[i]);
	free(ref_info.ref_indices);

	for(i=0; i<Motifs; i++){
		avg[i] /= ensemble;	stdev[i] = sqrt(stdev[i]/ensemble - avg[i]*avg[i]) + 1e-3;
		z_score[i] = (motif[i]-avg[i])/stdev[i];
		// printf("motif : %2d\t%6.3lf\t%6.3lf\t%d\r\n", i+1, avg[i], stdev[i], (int)(avg[i]*ensemble)%1);
	}
	return z_score;
}