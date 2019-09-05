/*
Wrote @  14/06/2019
Wrote by Zheng Zhang
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include "../utils/utils.h"
#include "../utils/list.h" 

int addNoise(double* noise_data, double* data, int N,double* noise_base, double SNR, double lambda,int cells, int fs){
    /*addNoise- add background activity (spikes) from neighbor cells whose arrival follows poisson distribution
    parameters：    
        noise_data:hold the noise added signal
        data:original signal
        N:length of signals
        noise_base: an 2D array of 1000 different noise spike samples with normalised peak value (peak at 1)
        SNR:desired samples
        lambda: average number of arrivals per sec, typically 20Hz
        cell:number of cells to simulate, typicall 3 cells
        fs:sampling frequency, in this case 24414Hz*/
    int D = 64; //duration of noise base
    int num_of_samples = 1000;
    int index_of_arrival = 0;
    lambda = lambda / fs; //20 spikes/sec -> 20/fs spikes/data points
    int index_interval = 0;
    double *noise_sel;

    //cal noise power
    double signal_power=0;
    double noise_power;
    for (int i = 0; i < N; i++){
        signal_power += pow(data[i] , 2);
    }
    signal_power /= N;
    noise_power = sqrt(signal_power / pow(10.0 ,SNR / 10));

    arrayCopy(data, noise_data,N);
    
    if (SNR != __DBL_MIN__){
        for (int i = 0; i < cells; i++){
            index_of_arrival = 0;
            while (index_of_arrival < N){
                index_interval = expDistribRV(lambda,(i+1)*(1+index_of_arrival)); //the interval of arrival follows exp distrib
                // printf("%d ", index_of_arrival);

                index_of_arrival += index_interval;    // the accumulation of arrival interval is the absolute index
                noise_sel = noise_base + D * (rand() % num_of_samples); //in total 1000 different noise bases
                for (int j = 0; j < D; j++){
                    if(index_of_arrival + j>N) // could out of range in the end
                        break;
                    noise_data[index_of_arrival + j] += noise_power * noise_sel[j]; //add noise to original data
                }
            }
        }
    }
    return 0;
}

int meanRemovalFilter(double* mean_removed_data, double *data, int N, int L){
    /*meanRemovalFilter-filter removes the local field potential by subtract
     the mean calculated by a moving average filer.
     parameters:
        mean_removed_data:resulatant signal
        data:singal to be filtered
        N:length of the signal
        L:length of the filter, 64 is chosen.
    */
    double MEAN;
    for (int i = 0; i < N;i++){
        if(i<N-L){  //when filter is full
            MEAN = mean(&data[i], L);
            mean_removed_data[i] = data[i] - MEAN;
        }
        else{       //when data is exiting filter
            MEAN = mean(&data[i], N-i);
            mean_removed_data[i] = data[i] - MEAN;
        }
    }
    return 0;
}

int ASO(double *preprocessed_data, double *data, int k, int N){
    /*ASO- Amplitude Scope Operator emphasises data with equation
    y(n)=x(n)(x(n)-x(n-1))
    parameters:
        preprocessed_data:resultant data
        data: signal to be processed
        k: resolution
        N: length of the signal
    */
    for (int i = 0; i < N; i++){
        if (i >= k || i <= N - k)
            preprocessed_data[i] = (data[i] - data[i - k]) * data[i];
        else //end and beginning
            preprocessed_data[i] = data[i];
    }
    return 0;
}

int NEO(double *preprocessed_data, double *data, int k, int N){
    /*ASO- Amplitude Scope Operator emphasises data with equation
    y(n)=x(n)(x(n)-x(n-1))
    parameters:
        preprocessed_data:resultant data
        data: signal to be processed
        k: resolution
        N: length of the signal
    */
    for (int i = 0; i < N; i++){
        if (i >= k || i <= N - k)
            preprocessed_data[i] = data[i] * data[i] - data[i + k] * data[i - k];
        else //end and beginning
            preprocessed_data[i] = data[i];
    }
    return 0;
}

double update(double *data, double c, int size, double* pre_mean){
    /*update-update function of thresholding, used by function 'thresholding'
    out=c*mean(data)
    parameters:
        data:input date
        c:user chosen parameters
        size:length of the window
        pre_mean: mean of the data
    */
    *pre_mean = mean(data, size);
    return c * *pre_mean;
}

int thresholding(list_node *spike_detected, list_node *interval, double *threshold, double *data, double c, int L, int N,int update_freq){
    /*Thresholding - adaptive mean thresholding according to thr=c*mean(data) data is the signal inside the adaptive window
    several tricks are applied
    1. Update frequency can be set to reduce computation
    2. When a spike is detected, spikes with not be replaced with the previous local mean 
    to avoid the abroupt change of threshold after a spike and also avoid redetecting the same spike.
    3. A subthreshold is set to be half of the threshold and if the signal outnumber this thr. they will 
    also be replaced by the previous mean but not count as a spike

    parameters:
        spike_detected: detected spike locations
        threshold: adaptive threshold calculated with mean of the signal
        data: signal to be thresholded 
        c: user chosen parameter
        L: adaptive window length
        N: signal length
        update_freq: frequency of the thresholding
    */
    double data_for_thr[N];
    arrAbs(data_for_thr, data, N); //absolute value is used to calculate the threshold
    
    int before = 5;
    int after = 10;
    int count = 0;
    int idx = 0;
    int M = 0; //number of spikes detected
    int indicator = 0;
    double pre_mean;  // contain the mean calculate from last updating, avoid division
    double init_threshold;
    init_threshold=update(data_for_thr, 3, L, &pre_mean); //initial threshold. c is set to be 3

    arrayCopyConst(init_threshold, threshold, L);

    for (int i = L; i < N; i++){
        double abs_data = fabs(data[i]);
        if ( abs_data< threshold[i - 1]){ // not detected
            if (abs_data >= 0.5 * threshold[i - 1])  //trick 3
                arrayCopyConst(pre_mean, &data_for_thr[i-before], before+after); //replace the original data with previous mean

            if (count == update_freq){ // if it's time to update, calculate the new threshold
            //Data used for calculating the threshold is from L+'before' to 'before' before the current location
            //Avoid using the raised data in the begining of a spike. Does NOT make much difference using 'before' or not.
                if (i - L - before >= 1) 
                    threshold[i] = update(&data_for_thr[i - L - before], c, L, &pre_mean); 
                else
                    threshold[i] = update(&data_for_thr[1], c, i - before, &pre_mean);
                count = 0; 
            }
            else { // otherise use the previous threshold
                threshold[i] = threshold[i - 1];
                count++; 
            }
        }
        else{ // spike detected
            M++;
            // this step can be reduced in parctice
            if (i + after < N - 1) // find the maximum from the current location to 'after' data afterwards
                idx = maxIndex(&data[i], after);
            else
                idx = maxIndex(&data[i], N-i-1);

            // this list data structure could be inefficient, it will iterate the hold list to insert in the end.
            if (indicator == 0){ // the fitst time, create  list
                *spike_detected = *list_create(i + idx);
                *interval = *list_create(i);
                indicator = 1;
            }
            else{ // store the detected location
                list_insert_end(spike_detected, i + idx);
                list_insert_end(interval, i);
            }

            arrayCopyConst(pre_mean,&data_for_thr[i-before], before+after);//replace the original data with previous mean
            arrayCopyConst(threshold[i-1],&threshold[i], after);// hold the threshold

 
            count=update_freq; //must update the threshold after spike detected
            
            i+= (after-1); //current location jumps to current+after
        }
    }

    /*naive*/
    // int last_spike=-10;
    // for (int i = L; i < N; i++){
    //     if (i - L - before >= 1)
    //         threshold[i] = update(&data_for_thr[i - L - before], c, L, &pre_mean);
    //     else
    //         threshold[i] = update(&data_for_thr[1], c, i - before, &pre_mean);
    //     if (data[i] > threshold[i - 1]){
    //         if (i + after < N - 1)
    //             idx = maxIndex(&data[i], after);
    //         else
    //             idx = maxIndex(&data[i], N-i-1);
    //         if (indicator == 0)
    //             *spike_detected = *list_create(i + idx);
    //         else if(i+idx>last_spike+after)
    //             list_insert_end(spike_detected, i + idx);
    //         indicator = 1;
    //         last_spike = i + idx;
    //     }      
    // }

    /*standard*/
    // for (int i = L; i < N; i++){
    //     // if (i>hold){
    //         if(data[i] < threshold[i - 1]){
    //             if (i - L - before >= 1)
    //                 threshold[i] = update(&data_for_thr[i - L - before], c, L, &pre_mean);
    //             else
    //                 threshold[i] = update(&data_for_thr[1], c, i - before, &pre_mean);
    //         }
    //         else {
    //             if (i + after < N - 1)
    //                 idx = maxIndex(&data[i], after);
    //             else
    //                 idx = maxIndex(&data[i], N-i-1);

    //             if (indicator == 0)
    //                 *spike_detected = *list_create(i + idx);
    //             else
    //                 list_insert_end(spike_detected, i + idx);
    //             indicator = 1;
    //             arrayCopyConst(pre_mean,&data_for_thr[i-before], before+after);
    //             arrayCopyConst(threshold[i-1],&threshold[i], after);
    //             i += (after - 1);
    //         }    
    //     }  
    // // }
    writeData(data_for_thr, N, "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/data_for_thr_1.txt");

    return M; //return the number of detected spikes 
}


