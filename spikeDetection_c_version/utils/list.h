/*From https://github.com/dxtr/list*/

#ifndef _LIST_H
#define _LIST_H

typedef struct list_node {
	int data;
	struct list_node *next;
} list_node;

/* linked list */
list_node* list_create(int data);
void list_destroy(list_node **list);
list_node* list_insert_after(list_node *node, int data);
list_node* list_insert_beginning(list_node *list, int data);
list_node* list_insert_end(list_node *list, int data);
void list_remove(list_node **list, list_node *node);
void list_remove_by_index(list_node **list, int idx);

list_node *list_find_node(list_node *list, list_node *node);
list_node* list_find_by_data(list_node *list, int data);
list_node* list_find(list_node *list, int(*func)(list_node*,void*), int data);

int list_indexing(list_node *list, int idx);

void list_print(list_node *list, int iterms);

#endif

