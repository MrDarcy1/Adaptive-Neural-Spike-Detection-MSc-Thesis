/*From https://github.com/dxtr/list*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"

/* Creates a list (node) and returns it
 * Arguments: The data the list will contain or NULL to create an empty
 * list/node
 */
list_node* list_create(int data)
{
	list_node *l = malloc(sizeof(list_node));
	if (l != NULL) {
		l->next = NULL;
		l->data = data;
	}

	return l;
}

/* Completely destroys a list
 * Arguments: A pointer to a pointer to a list
 */
void list_destroy(list_node **list)
{
	if (list == NULL) return;
	while (*list != NULL) {
		list_remove(list, *list);
	}
}

/* Creates a list node and inserts it after the specified node
 * Arguments: A node to insert after and the data the new node will contain
 */
list_node* list_insert_after(list_node *node, int data)
{
	list_node *new_node = list_create(data);
	if (new_node) {
		new_node->next = node->next;
		node->next = new_node;
	}
	return new_node;
}

/* Creates a new list node and inserts it in the beginning of the list
 * Arguments: The list the node will be inserted to and the data the node will
 * contain
 */
list_node* list_insert_beginning(list_node *list, int data)
{
	list_node *new_node = list_create(data);
	if (new_node != NULL) { new_node->next = list; }
	return new_node;
}

/* Creates a new list node and inserts it at the end of the list
 * Arguments: The list the node will be inserted to and the data the node will
 * contain
 */
list_node* list_insert_end(list_node *list, int data)
{
	list_node *new_node = list_create(data);
	if (new_node != NULL) {
		for(list_node *it = list; it != NULL; it = it->next) {
			if (it->next == NULL) {
				it->next = new_node;
				break;
			}
		}
	}
	return new_node;
}

/* Removes a node from the list
 * Arguments: The list and the node that will be removed
 */
void list_remove(list_node **list, list_node *node)
{
	list_node *tmp = NULL;
	if (list == NULL || *list == NULL || node == NULL) return;

	if (*list == node) {
		*list = (*list)->next;
		free(node);
		node = NULL;
	} else {
		tmp = *list;
		while (tmp->next && tmp->next != node) tmp = tmp->next;
		if (tmp->next) {
			tmp->next = node->next;
			free(node);
			node = NULL;
		}
	}
}


/* Find an element in a list by the pointer to the element
 * Arguments: A pointer to a list and a pointer to the node/element
 */
list_node* list_find_node(list_node *list, list_node *node)
{
	while (list) {
		if (list == node) break;
		list = list->next;
	}
	return list;
}

/* Finds an elemt in a list by the data pointer
 * Arguments: A pointer to a list and a pointer to the data
 */
list_node* list_find_by_data(list_node *list, int data)
{
	while (list) {
		if (list->data == data) break;
		list = list->next;
	}
	return list;
}


/**************************************************/
/*wrote @ 06/18/2019 
  wrote by Zheng Zhang*/


void list_remove_by_index(list_node **list, int idx){
/*list_remove_by_index - remove the indexed node
parameters:
	list: address of the list head
	idx: index of the node to be removed */

	list_node *tmp1 = NULL;
	list_node *tmp2 = NULL;
	int count = 0;
	if (list == NULL || *list == NULL)
		return;

	if (idx == 0) {
		*list = (*list)->next;
	} else {
		tmp1 = *list;
		while (tmp1->next && count < idx){
			tmp2 = tmp1;
			tmp1 = tmp1->next;
			count++;
		}
		if (tmp1) {
			tmp2->next = tmp1->next;
			free(tmp1);
			tmp1 = NULL;
		}
	}
}

int list_indexing(list_node *list, int idx){
/*list_indexing -return the value of indexed node
parameters:
	list: list to be indexed
	idx: index of the desired node*/
	list_node *p = malloc(sizeof(list_node));
    if(list==NULL){
		printf("\nlist is empty");
	}
	else{
		p = list;
		for (int i = 0; i < idx; i++){
			if (p==NULL){
				printf("not enough data");
				return -1;
			}
			else{
				p=p->next;
			}
        }
    }
	return p->data;
}

void list_print(list_node *list, int iterms){
/*lisy_print - print the values stored in the lists
parameters: 
	list: list to be printed
	iterms: number of iterms to be printed */
	struct list_node *p;
    if(list==NULL){
		printf("\nlist is empty");
	}
	else{
		p = list;
		for (int i = 0; i < iterms;i++){
			if (p==NULL){
				printf("not enough data");
				return;
			}
			else{
				printf("%d ",p->data);
				p=p->next;
			}
        }
    }
}

