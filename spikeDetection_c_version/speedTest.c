/*
Wrote @  13/06/2019
Wrote by Zheng Zhang
program to detect the spikes with aso operator and improved mean thresholding
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utils/utils.h"
#include "utils/list.h"
#include "functions/functions.h"
int main(){
    //load data
    clock_t start, finish;
    double duration;
    char loc_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/spike_location_1.txt";
    char data_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/realDataWithLFP_1.txt";
    char noise_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/noise_base.txt";
    int data_start = 0;
    int num_of_data = 40000;
    int N = 0;
    double data[num_of_data];
	list_node* loc = (list_node*)malloc(sizeof(list_node));
    double* noise_base=(double*)malloc(__SIZEOF_DOUBLE__*1000*64);
    loadData(data_start, num_of_data, data_file_name, data);
    loadLocation(data_start, num_of_data+data_start, loc, loc_file_name,data,&N);
    loadNoise(noise_file_name, noise_base);

    /*add noise*/
    double noise_data[num_of_data];
    double SNR = __DBL_MAX__;
    // double SNR = 20;
    double lambda = 20;
    double cells = 3;
    double fs = 24144;
    addNoise(noise_data, data, num_of_data, noise_base, SNR, lambda, cells, fs);
    free(noise_base);

    int window_length = 16;
    double mean_removed_data[num_of_data];

    int resolution = 2;
    double preprocessed_data[num_of_data];

    list_node* spike_detected = (list_node*)malloc(sizeof(list_node));
    list_node* interval = (list_node*)malloc(sizeof(list_node));

    double c = 25;
    int L = 64;
    int update_freq = 10000;
    int M = 0;
    double threshold[num_of_data];


    /*subtract mean*/
    start = clock();
    for (int i = 0; i < 100000;i++){
        meanRemovalFilter(mean_removed_data, noise_data, num_of_data, window_length);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of filtering: %lf sec\n", duration);

    /* spike emphasis */
    //neo
    start = clock();
    for (int i = 0; i < 800000;i++){
        NEO(preprocessed_data, mean_removed_data, resolution, num_of_data);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of NEO emphasising: %lf sec\n", duration);

    //aso
    start = clock();
    for (int i = 0; i < 800000;i++){
        ASO(preprocessed_data, mean_removed_data, resolution, num_of_data);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of ASO emphasising: %lf sec\n", duration);

    /*thresholding*/
    update_freq = 0;
    start = clock();
    for (int i = 0; i < 1000;i++){
        M = thresholding(spike_detected, interval, threshold, preprocessed_data, c, L, num_of_data, update_freq);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of thresholding with update freqeuncy %d: %lf sec\n", update_freq,duration);


    update_freq = 100;
    start = clock();
    for (int i = 0; i < 1000;i++){
        M = thresholding(spike_detected, interval, threshold, preprocessed_data, c, L, num_of_data, update_freq);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of thresholding with update freqeuncy %d: %lf sec\n", update_freq,duration);


    update_freq = 1000;
    start = clock();
    for (int i = 0; i < 1000;i++){
        M = thresholding(spike_detected, interval, threshold, preprocessed_data, c, L, num_of_data, update_freq);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of thresholding with update freqeuncy %d: %lf sec\n", update_freq,duration);

    update_freq = 1000;
    start = clock();
    for (int i = 0; i < 1000;i++){
        M = thresholding(spike_detected, interval, threshold, preprocessed_data, c, L, num_of_data, update_freq);
    }
    finish= clock();
    duration = (double)(finish- start) / CLOCKS_PER_SEC;
    printf("duration of thresholding with update freqeuncy %d: %lf sec\n", update_freq,duration);
}
