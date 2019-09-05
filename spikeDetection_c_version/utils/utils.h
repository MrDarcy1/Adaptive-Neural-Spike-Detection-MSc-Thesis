
#ifndef UTILS_H
#define UTILS_H

#include "list.h"


int loadData(long start, int length, char *name, double *data);
list_node* loadLocation(int start, int end, list_node *loc, char *name, double *data,int *N);
int loadNoise(char *name, double *data);

double expDistribRV(double lambda,int r);

double max(double *arr, int size);
int maxIndex(double *arr, int size);
double mean(double *arr, int L);

int arrayCopy(double *ori, double *cop, int size);
int arrayCopyConst(double ori, double* cop, int size);
int arrAbs(double *abs_data, double *data, int length);

int *locationCompare(list_node *locations, int N, list_node *interval, list_node *spike_detected, int M, list_node *FP, list_node *FN, list_node *TP);

void writeData(double *data, int length, char *name);

#endif