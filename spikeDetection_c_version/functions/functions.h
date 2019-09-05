

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "../utils/list.h"

int addNoise(double *noise_data, double *data, int N, double *noise_base, double SNR, double lambda, int cells, int fs);
int meanRemovalFilter(double* mean_removed_data, double *data, int N, int L);

int ASO(double *preprocessed_data, double *data, int k, int N);
int NEO(double *preprocessed_data, double *data, int k, int N);

int thresholding(list_node *spike_detected,list_node *interval, double *threshold, double *data, double c, int L, int N, int update_freq);

#endif