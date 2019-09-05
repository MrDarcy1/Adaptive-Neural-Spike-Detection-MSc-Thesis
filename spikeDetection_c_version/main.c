/*
Wrote @  13/06/2019
Wrote by Zheng Zhang
program to detect the spikes with aso operator and improved mean thresholding
*/
#include <stdio.h>
#include <stdlib.h>


#include "utils/utils.h"
#include "utils/list.h"
#include "functions/functions.h"
int main(){
    //load data

    char loc_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/spike_location_1.txt";
    char data_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/realDataWithLFP_1.txt";
    char noise_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/noise_base.txt";
    int data_start = 1000;
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

    /*subtract mean*/
    int window_length = 16;
    double mean_removed_data[num_of_data];
    meanRemovalFilter(mean_removed_data, noise_data, num_of_data, window_length);

    //spike emphasis
    int resolution = 2;
    double preprocessed_data[num_of_data];
    ASO(preprocessed_data, mean_removed_data, resolution, num_of_data);
    // NEO(preprocessed_data, mean_removed_data, resolution, num_of_data);
   
    /*thresholding*/
    list_node* spike_detected = (list_node*)malloc(sizeof(list_node));
    list_node* interval = (list_node*)malloc(sizeof(list_node));

    double c = 25;
    int L = 64;
    int update_freq = 1000;
    int M = 0;
    double threshold[num_of_data];
    M = thresholding(spike_detected, interval, threshold, preprocessed_data, c, L, num_of_data, update_freq);
    
    /*evaluation*/
    list_node *FP = (list_node *)malloc(sizeof(list_node));
    list_node *TP = (list_node *)malloc(sizeof(list_node));
    list_node *FN = (list_node *)malloc(sizeof(list_node));
    int *count = (int *)malloc(__SIZEOF_INT__ * 3); //TP,FP,FN
    count=locationCompare(loc, N, interval, spike_detected, M, FP, FN, TP);
    printf("TP: %d\n", count[0]);
    printf("FP: %d\n", count[1]);
    printf("FN: %d\n", count[2]);

    float Sens = (float)count[0] / (float)(count[0] + count[2]);
    float FDR = (float)count[2] / (float)(count[0] + count[1]);
    float Acc = (float)count[0] / (float)(count[0] + count[1] + count[2]);
    printf("Sens: %f, FDR: %f, Acc: %f", Sens, FDR, Acc);

    
    char thr_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/threshold_1.txt";
    char processed_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/Processed_Data_1.txt";

    writeData(threshold, num_of_data, thr_file_name);
    writeData(preprocessed_data, num_of_data, processed_file_name);
    
    // printf("\n");

    // for (int i = 0; i < 10; i++){ //print  data
    //     printf("%lf ", data[i]);
    // }

    // printf("\n");

    // list_print(spike_detected, 20);
    // printf("\n");
    // list_print(loc, 20);
    // 

}
