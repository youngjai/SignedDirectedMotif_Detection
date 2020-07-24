/*
	This code is for loading pajex files.

	It was written by Youngjai in Nov. 19, 2018.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

typedef struct {
	Net_info net_info;
	int node_label;
} Load_profile;

/* Load some packages from the 'mt19937.c'. */
double genrand_real2(void);

/* Load some packages from the 'youngjai_pacakeges.c'. */
double normal_dist(double avg, double stdev);
void swap_two_integers(int *num_i, int *num_j);
void swap_two_floats(double *num_i, double *num_j);
int *ascending_bubble_sorting_int(int *array, int length);

void print_test_temp(Net_info net_info)
{
	int i;
	int positive, negative;	positive = 0;	negative = 0;
	Node *cursor = net_info.node_list->head;
	printf("(Nodes : %d\tLinks : %d)\r\n", net_info.nodes, net_info.links);
	while(cursor->next != NULL){
		cursor = cursor->next;
		printf("%2d\t%2d\t%.2lf:\t", cursor->node_label, cursor->degree, cursor->fitness);
		for(i=0; i<cursor->degree; i++){
			printf("%2d(%2d, %.2lf)\t", cursor->neighbor[i], cursor->direction[i], cursor->sign[i]);
			if(cursor->sign[i] > 0) positive++;
			else negative++;
		}
		printf("\r\n");
	}
	printf("+ : %d\t- : %d\r\n", positive/2, negative/2);
}

void free_Node(Node *cursor)
{
	free(cursor->neighbor);
	free(cursor->direction);
	free(cursor->sign);
	free(cursor);
}

void free_Linked_list(Linked_list *node_list)
{
	Node *cursor = node_list->head;
	Node *cursor1;
	while(cursor != NULL){
		cursor1 = cursor;
		cursor = cursor->next;
		free_Node(cursor1);
	}
}

void create_node(Linked_list *node_list, int node_label)
{
	Node *new_node = malloc(sizeof(Node));
	new_node->node_label = node_label;
	new_node->degree = 0;
	new_node->neighbor = malloc(sizeof(int));
	new_node->direction = malloc(sizeof(int));
	new_node->sign = malloc(sizeof(double));
	new_node->fitness = 0.;
	new_node->next = NULL;

	if(node_list->head == NULL && node_list->tail == NULL)
		node_list->head = node_list->tail = new_node;
	else{
		node_list->tail->next = new_node;
		node_list->tail = new_node;
	}
}

/* Make Erdős–Rényi network. */
void create_ER_connection(Linked_list *node_list, int name_i, int name_j, int sign_ij)
{
	Node *cursor = node_list->head;
	Node *node_i, *node_j;
	int count = 0;
	int direction_ij;
	/* '1' means an outgoing degree, '-1' is an incoming degree. */
	if(genrand_real2() < 0.5) direction_ij = -1;
	else direction_ij = 1;

	while(cursor->next != NULL){
		cursor = cursor->next;
		if(cursor->node_label == name_i){
			node_i = cursor;
			count++;
		}
		if(cursor->node_label == name_j){
			node_j = cursor;
			count++;
		}
		if(count == 2) break;
	}

	node_i->neighbor = realloc(node_i->neighbor, sizeof(int) * (node_i->degree+2));
	node_i->direction = realloc(node_i->direction, sizeof(int) * (node_i->degree+2));
	node_i->sign = realloc(node_i->sign, sizeof(double) * (node_i->degree+2));

	node_j->neighbor = realloc(node_j->neighbor, sizeof(int) * (node_j->degree+2));
	node_j->direction = realloc(node_j->direction, sizeof(int) * (node_j->degree+2));
	node_j->sign = realloc(node_j->sign, sizeof(double) * (node_j->degree+2));

	node_i->neighbor[node_i->degree] = node_j->node_label;
	node_j->neighbor[node_j->degree] = node_i->node_label;

	node_i->direction[node_i->degree] = -direction_ij;
	node_j->direction[node_j->degree] = direction_ij;

	node_i->sign[node_i->degree] = sign_ij;
	node_j->sign[node_j->degree] = sign_ij;

	node_i->degree++;
	node_j->degree++;
}
Linked_list *make_ER_network(Net_info net_info, int ref_links)
{
	int i, j;
	int node_label = 0;
	int count = 0;
	int name_i, name_j;

	int *shuffled = malloc(sizeof(int) * (net_info.nodes*net_info.nodes+1));
	for(i=0; i<net_info.nodes*net_info.nodes; i++) shuffled[i] = i;
	for(i=net_info.nodes*net_info.nodes-1; i>0; i--){
		j = (int) ((i+1) * genrand_real2());
		swap_two_integers(&shuffled[i], &shuffled[j]);
	}

	//Linked_list pointer define start
	net_info.node_list = malloc(sizeof(Linked_list));
	net_info.node_list->head = NULL; net_info.node_list->tail = NULL;
	create_node(net_info.node_list, node_label++);		// define head pointer
	//Linked_list pointer define end

	for(i=0; i<net_info.nodes; i++) create_node(net_info.node_list, node_label++);

	i = 0;
	while(count < net_info.links && i < net_info.nodes*net_info.nodes){
		name_i = shuffled[i]/net_info.nodes%net_info.nodes;
		name_j = shuffled[i++]%net_info.nodes;

		if(name_i < name_j){	// Avoid selfloop and multi-edges.
			/* ++ is because the index start to 0. */
			if(count < ref_links) create_ER_connection(net_info.node_list, ++name_i, ++name_j, 1);
			else create_ER_connection(net_info.node_list, ++name_i, ++name_j, -1);
			count++;
		}
	}
	free(shuffled);
	return net_info.node_list;
}

/* Make random network conserved information of degree, sign, and direction. */
int check_double_links(Ref_info ref_info, int nodes, int node_index, int name_i, int name_j, int name_k, int name_l, int rand_level)
{
	int i, j;
	int count = 0;
	int multi_count = 0;
	int index_i, index_j, index_k, index_l;

	/* Avoid multi-edges. */
	if(rand_level == 1){
		for(i=0; i<ref_info.ref_links[0]+ref_info.ref_links[1]; i++){
			if(ref_info.ref_array[i] == name_i*node_index + name_l || \
				ref_info.ref_array[i] == name_l*node_index + name_i || \
				ref_info.ref_array[i] == name_k*node_index + name_j || \
				ref_info.ref_array[i] == name_j*node_index + name_k){
				multi_count++; break;
			}
		}
	}
	else{
		for(i=0; i<nodes+1; i++){
			if(ref_info.ref_indices[0][i]-1 == name_i){ index_i = i; count++;}
			else if(ref_info.ref_indices[0][i]-1 == name_j){ index_j = i; count++;}
			if(ref_info.ref_indices[0][i]-1 == name_k){ index_k = i; count++;}
			else if(ref_info.ref_indices[0][i]-1 == name_l){ index_l = i; count++;}

			if(count == 4) break;
		}

		for(i=0; i<2; i++){
			for(j=i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_i-1]; \
				j<i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_i]; j++){
				if(ref_info.ref_array[j]%node_index == name_l){ multi_count++; break;}
			}
			if(multi_count != 0) break;
			for(j=i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_j-1]; \
				j<i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_j]; j++){
				if(ref_info.ref_array[j]%node_index == name_k){ multi_count++; break;}
			}
			if(multi_count != 0) break;
			for(j=i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_k-1]; \
				j<i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_k]; j++){
				if(ref_info.ref_array[j]%node_index == name_j){ multi_count++; break;}
			}
			if(multi_count != 0) break;
			for(j=i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_l-1]; \
				j<i*ref_info.ref_links[0]+ref_info.ref_indices[i+1][index_l]; j++){
				if(ref_info.ref_array[j]%node_index == name_i){ multi_count++; break;}
			}
			if(multi_count != 0) break;
		}
	}

	return multi_count;
}
Ref_info reference_links(Net_info net_info, int node_index, int rand_level)
{
	int i, j, k;
	int ref_index;
	Ref_info ref_info;
	ref_info.ref_array = malloc(sizeof(int) * (net_info.links+1));
	ref_info.ref_links = malloc(sizeof(int) * (2+1));
	ref_info.ref_indices = malloc(sizeof(int *) * (3+1));
	int **array_temp = malloc(sizeof(int *) * (2+1));
	int *ref_degree = malloc(sizeof(int) * (2+1));
	for(i=0; i<2; i++){
		ref_info.ref_links[i] = 0;
		array_temp[i] = malloc(sizeof(int) * (net_info.links+1));
	}
	for(i=0; i<3; i++){
		ref_info.ref_indices[i] = malloc(sizeof(int) * (net_info.nodes+2));
		ref_info.ref_indices[i][0] = 0;
	}
	Node *cursor = net_info.node_list->head; i = 0;
	while(cursor->next != NULL){
		cursor = cursor->next;

		ref_degree[0] = 0; ref_degree[1] = 0;
		for(j=0; j<cursor->degree; j++){
			if(cursor->direction[j] > 0){ 
				if(cursor->sign[j] > 0) ref_index = 0;
				else ref_index = 1;
				array_temp[ref_index][ref_info.ref_links[ref_index]++] = \
					(cursor->node_label-1)*node_index + (cursor->neighbor[j]-1);
				ref_degree[ref_index]++;
			}
		}
		ref_info.ref_indices[0][i] = cursor->node_label;
		ref_info.ref_indices[1][i] = ref_degree[0];
		ref_info.ref_indices[2][i++] = ref_degree[1];
	}

	for(i=0; i<net_info.nodes; i++) for(j=0; j<(net_info.nodes-i-1); j++)
		if(ref_info.ref_indices[0][j] > ref_info.ref_indices[0][j+1]){
			for(k=0; k<3; k++) swap_two_integers(&ref_info.ref_indices[k][j], &ref_info.ref_indices[k][j+1]);
	}

	if(rand_level == 1){
		for(i=net_info.nodes-1; i>=0; i--){
			for(j=0; j<i; j++){
				ref_info.ref_indices[1][i] += ref_info.ref_indices[1][j] + ref_info.ref_indices[2][j];
			}
			ref_info.ref_indices[1][i] += ref_info.ref_indices[2][i];
			ref_info.ref_indices[2][i] = 0;
		}
		for(i=0; i<2; i++) for(j=0; j<ref_info.ref_links[i]; j++)
				ref_info.ref_array[i*ref_info.ref_links[0]+j] = array_temp[i][j];
		ref_info.ref_array = ascending_bubble_sorting_int(ref_info.ref_array, net_info.links);
	}
	else if(rand_level == 2){
		for(i=net_info.nodes-1; i>=0; i--){
			for(j=0; j<i; j++){
				ref_info.ref_indices[1][i] += ref_info.ref_indices[1][j] + ref_info.ref_indices[2][j];
			}
			ref_info.ref_indices[1][i] += ref_info.ref_indices[2][i];
			ref_info.ref_indices[2][i] = 0;
		}
		for(i=0; i<2; i++) for(j=0; j<ref_info.ref_links[i]; j++)
				ref_info.ref_array[i*ref_info.ref_links[0]+j] = array_temp[i][j];
		ref_info.ref_array = ascending_bubble_sorting_int(ref_info.ref_array, net_info.links);
	}
	else if(rand_level == 3){
		for(i=net_info.nodes-1; i>=0; i--){
			for(j=0; j<i; j++){
				ref_info.ref_indices[1][i] += ref_info.ref_indices[1][j];
				ref_info.ref_indices[2][i] += ref_info.ref_indices[2][j];
			}
		}
		for(i=0; i<2; i++){
			array_temp[i] = ascending_bubble_sorting_int(array_temp[i], ref_info.ref_links[i]);
			for(j=0; j<ref_info.ref_links[i]; j++)
				ref_info.ref_array[i*ref_info.ref_links[0]+j] = array_temp[i][j];
		}
	}
	for(i=0; i<2; i++) free(array_temp[i]);
	free(array_temp); free(ref_degree);

	return ref_info;
}
void create_rewired_interaction(Net_info net_info, int name_i, int name_j, int name_k, int name_l, int rand_level)
{
	int index_i, index_j, index_k, index_l;
	Node *cursor = net_info.node_list->head;
	Node *node_i, *node_j, *node_k, *node_l;
	int count = 0;

	while(cursor->next != NULL){
		cursor = cursor->next;
		if(cursor->node_label == name_i){
			node_i = cursor;
			count++;
		}
		else if(cursor->node_label == name_j){
			node_j = cursor;
			count++;
		}
		if(cursor->node_label == name_k){
			node_k = cursor;
			count++;
		}
		else if(cursor->node_label == name_l){
			node_l = cursor;
			count++;
		}
		if(count == 4) break;
	}
	for(index_i=0; index_i<node_i->degree; index_i++)
		if(node_i->neighbor[index_i] == name_j) break;
	for(index_j=0; index_j<node_j->degree; index_j++)
		if(node_j->neighbor[index_j] == name_i) break;
	for(index_k=0; index_k<node_k->degree; index_k++)
		if(node_k->neighbor[index_k] == name_l) break;
	for(index_l=0; index_l<node_l->degree; index_l++)
		if(node_l->neighbor[index_l] == name_k) break;

	swap_two_integers(&node_i->neighbor[index_i], &node_k->neighbor[index_k]);
	swap_two_integers(&node_j->neighbor[index_j], &node_l->neighbor[index_l]);
	if(rand_level < 3){
		swap_two_floats(&node_j->sign[index_j], &node_l->sign[index_l]);
		if(rand_level < 2)
			swap_two_integers(&node_j->direction[index_j], &node_l->direction[index_l]);
	}
}
Linked_list *make_rewired_network(Net_info net_info, Ref_info ref_info, int node_index, int rand_level)
{
	int i;
	int limit = net_info.links*net_info.links;
	int count = 0;
	int name_i, name_j, name_k, name_l;
	int index[2];
	int trial = 0;

	while(count < net_info.links && trial++ < limit){
		if(rand_level == 3){
			if(genrand_real2() < 0.5) i = 0;
			else i = 1;

			do{
				index[0] = (int) (i*ref_info.ref_links[0] + ref_info.ref_links[i]*genrand_real2());
				index[1] = (int) (i*ref_info.ref_links[0] + ref_info.ref_links[i]*genrand_real2());
			} while(index[0] == index[1]);
		}
		else{
			do{
				index[0] = (int) (net_info.links*genrand_real2());
				index[1] = (int) (net_info.links*genrand_real2());
			} while(index[0] == index[1]);
		}
		name_i = ref_info.ref_array[index[0]]/node_index%node_index;
		name_j = ref_info.ref_array[index[0]]%node_index;
		name_k = ref_info.ref_array[index[1]]/node_index%node_index;
		name_l = ref_info.ref_array[index[1]]%node_index;

		if(rand_level == 1){
			if(genrand_real2() > 0.5) swap_two_integers(&name_i, &name_j);
			if(genrand_real2() > 0.5) swap_two_integers(&name_k, &name_l);
		}

		if(name_i != name_l && name_j != name_k){ // Avoid selfloop. 
			/* ++ is because the node label start to 1. */
			if(check_double_links(ref_info, net_info.nodes, node_index, name_i, name_j, name_k, name_l, rand_level) == 0){
				ref_info.ref_array[index[0]] = name_i*node_index + name_l;
				ref_info.ref_array[index[1]] = name_k*node_index + name_j;
				create_rewired_interaction(net_info, ++name_i, ++name_j, ++name_k, ++name_l, rand_level);
				count++; trial = 0;
			}
		}
	}
	return net_info.node_list;
}

/* Make a loaded network. */
Net_info create_loaded_connection(Net_info net_info, int **links_profile)
{
	int i, j;
	int links = net_info.links;
	int count;
	int sign_ij;
	int name_i, name_j;
	Node *cursor, *node_i, *node_j;

	for(i=0; i<links; i++){
		cursor = net_info.node_list->head;
		count = 0;
		name_i = links_profile[i][0];
		name_j = links_profile[i][1];
		sign_ij = links_profile[i][2];
		// '1' means an outgoing degree, '-1' is an incoming degree. 

		if(name_i != name_j){
			while(cursor->next != NULL){
				cursor = cursor->next;
				if(cursor->node_label == name_i){
					node_i = cursor;
					count++;
				}
				if(cursor->node_label == name_j){
					node_j = cursor;
					count++;
				}
				if(count == 2) break;
			}

			count = 0;
			for(j=0; j<node_i->degree; j++) if(node_i->neighbor[j] == name_j){ count++; break;}

			if(count == 0){
				node_i->neighbor = realloc(node_i->neighbor, sizeof(int) * (node_i->degree+2));
				node_i->direction = realloc(node_i->direction, sizeof(int) * (node_i->degree+2));
				node_i->sign = realloc(node_i->sign, sizeof(double) * (node_i->degree+2));

				node_j->neighbor = realloc(node_j->neighbor, sizeof(int) * (node_j->degree+2));
				node_j->direction = realloc(node_j->direction, sizeof(int) * (node_j->degree+2));
				node_j->sign = realloc(node_j->sign, sizeof(double) * (node_j->degree+2));

				node_i->neighbor[node_i->degree] = name_j;
				node_j->neighbor[node_j->degree] = name_i;

				node_i->direction[node_i->degree] = 1;
				node_j->direction[node_j->degree] = -1;

				node_i->sign[node_i->degree] = sign_ij;
				node_j->sign[node_j->degree] = sign_ij;

				node_i->degree++;
				node_j->degree++;
			}
			else{
				for(j=0; j<node_i->degree; j++) if(node_i->neighbor[j] == name_j){
					node_i->neighbor[j] = name_j;
					node_i->direction[j] = 1;
					node_i->sign[j] = sign_ij;
				}
				for(j=0; j<node_j->degree; j++) if(node_j->neighbor[j] == name_i){
					node_j->neighbor[j] = name_i;
					node_j->direction[j] = -1;
					node_j->sign[j] = sign_ij;
				}
				net_info.links--;
			}
		}
		else net_info.links--;
	}
	return net_info;
}
Net_info make_loaded_network(Net_info net_info, int *nodes_profile, int **links_profile)
{
	int i, j;

	//Linked_list pointer define start
	net_info.node_list = malloc(sizeof(Linked_list));
	net_info.node_list->head = NULL;	net_info.node_list->tail = NULL;
	create_node(net_info.node_list, 0);		// define head pointer
	//Linked_list pointer define end

	for(i=0; i<net_info.nodes; i++) create_node(net_info.node_list, nodes_profile[i]);
	net_info = create_loaded_connection(net_info, links_profile);

	return net_info;
}
Load_profile load_networks(char *filename)
{
	Load_profile load_profile;
	load_profile.node_label = 0;
	int *nodes_profile;
	int **links_profile;

	int i, index;
	char *tmp = malloc(sizeof(char) * 10);

	FILE *network = fopen(filename, "rt");

	fscanf(network, "%s %d\r\n", tmp, &load_profile.net_info.nodes);
	nodes_profile = malloc(sizeof(int) * (load_profile.net_info.nodes+1));
	for(i=0; i<load_profile.net_info.nodes; i++){
		fscanf(network, "%d %d\r\n", &index, &nodes_profile[i]);
		if(nodes_profile[i] > load_profile.node_label) load_profile.node_label = nodes_profile[i];
	}

	fscanf(network, "%s %d\r\n", tmp, &load_profile.net_info.links);
	links_profile = malloc(sizeof(int *) * (load_profile.net_info.links+1));
	for(i=0; i<load_profile.net_info.links; i++){
		links_profile[i] = malloc(sizeof(int) * (4+1));
		fscanf(network, "%d\t%d\t%d\t%d\r\n", &links_profile[i][0], &links_profile[i][1], \
			&links_profile[i][2], &links_profile[i][3]);
		// links_profile[i] = malloc(sizeof(int) * (3+1));
		// fscanf(network, "%d\t%d\t%d\r\n", &links_profile[i][0], &links_profile[i][1], \
			&links_profile[i][2]);
	}
	fclose(network);

	index = load_profile.net_info.links;
	load_profile.net_info = make_loaded_network(load_profile.net_info, nodes_profile, links_profile);
	for(i=0; i<index+1; i++) free(links_profile[i]);
	free(links_profile); free(nodes_profile); free(tmp);

	return load_profile;
}

/* Make a test ER network. */
Net_info test_network(int nodes, int links)
{
	Net_info net_info;
	net_info.nodes = nodes;
	net_info.links = links;
	net_info.node_list = make_ER_network(net_info, links/2);
	return net_info;
}