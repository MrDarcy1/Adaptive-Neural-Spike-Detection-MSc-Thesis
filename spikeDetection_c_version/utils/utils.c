/*
Wrote @  13/06/2019
Wrote by Zheng Zhang
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "list.h"

//#include <float.h>
int loadData(long start, int length, char* name,double* data){
    /*loadData- load the amplitude of signals from target files
    parameters:
        start: start index
        length: number of data points to load
        name: file name
        data: pointer for loaded data*/
    FILE *fp;

    if((fp =fopen(name,"r")) == NULL)
	{
        printf(" %s file does not exist\n",name);
        return 0;    
	}
    double temp;
    for (int i = 0; i < start;i++){
        fscanf(fp, "%lf ", &temp ); //load data from the correct location for length data points

    }
        // fseek(fp, 10 * start, 0); //move to start index; 1 double+' ' = 10bytes
    for(int i=0; i < length; i++)
	{
        fscanf(fp, "%lf ", data+i ); //load data from the correct location for length data points
    }
    fclose(fp);
    return 1;
}

list_node* loadLocation(int start, int end,list_node *loc, char* name,double *data, int* N){
    /*loadLocation- load the index of spikes
    parameters:
        start:index to start
        end: index to end
        name: file name
        data: pointer for loaded spike locations*/
    FILE *fp;
    int temp_loc;
    if((fp =fopen(name,"r")) == NULL)
	{
        printf(" %s file does not exist\n",name);
        return 0;    
	}
    int indicator = 0;

    while(1){
        fscanf(fp,"%d ",&temp_loc); //load the current index
        if(temp_loc>=start){// if the index is within desired range, load to loc.
            if (indicator==0){
                *loc = *list_create(temp_loc-start);
                (*N)++;
                indicator = 1;
            }
            else{
                if(temp_loc > end)
                    break;
                list_insert_end(loc, temp_loc-start);
                (*N)++;
            }
        }
    } 

    fclose(fp);
    return loc;
}
int loadNoise(char* name, double* data){
    /*loadData- load the amplitude of signals from target files
    parameters:
        start: start index
        length: number of data points to load
        name: file name
        data: pointer for loaded data*/
    FILE *fp;

    if((fp =fopen(name,"r")) == NULL)
	{
        printf(" %s file does not exist\n",name);
        return 0;    
	}
    for(int i=0; i < 1000*64; i++)
	{
        fscanf(fp, "%lf ", data+i ); //load data from the correct location for length data points
    }
    fclose(fp);
    return 1;
}

double expDistribRV(double lambda,int r){
    /*expDistribRV- Generate exponential random variable with average time of arrival lambda
    f(x)=lambda*exp^(-lambda*x); E(x)=1/lambda
    -> x=-log(1-(f(x)/lambda))/(lambda)
    parameters:
        lambda:average time of arrival
    */
   srand((int)time(NULL)*(int)time(NULL)*r); 
    double pV; //prob density
    while(1){ // pV cannot be 1
		pV = (double)rand()/(double)RAND_MAX;
        if (pV != 1)
			break;
    }
    pV = -log(1 - pV) / lambda;
    return pV;
}


double max(double* arr,int size){
    /*max-find the max of an array
    parameters:
        arr: target array
        size: length of arr*/
    double MAX = __DBL_MIN__;
    for (int i = 0; i < size;i++){
        if(*(arr + i) > MAX)
            MAX = *(arr + i);
    }
    return MAX;
}

int maxIndex(double* arr,int size){
    /*maxIndex-find the index of the maximum of an array
    parameters:
        arr: target array
        size: length of arr*/
    double MAX = __DBL_MIN__;
    int max_index;
    for (int i = 0; i < size;i++){
        if(*(arr + i) > MAX){
            MAX = *(arr + i);
            max_index = i;
        }
    }
    return max_index;
}

double mean(double* arr, int size){
    /*mean- calculate the mean of an array
    parameters:
        arr:input array
        size:length of the array
    */
    double sum=0;
    for (int i = 0; i < size; i++)
        sum += arr[i];
    return sum / size;
}


int arrayCopy(double* ori, double* cop,int size){
    /*arrCopy-copy array ori to cop
    parameters:
        ori:original array
        cop: resultant array
        size:length of ori*/
    for (int i = 0; i < size;i++)
        *cop++ = *ori++;
    return 0;
}

int arrayCopyConst(double ori, double* cop,int size){
    /*arrCopy-copy array ori to cop
    parameters:
        ori:original array
        cop: resultant array
        size:length of ori*/
    for (int i = 0; i < size;i++)
        *cop++ = ori;
    return 0;
}

int arrAbs(double* abs_data, double *data, int length){
    /*arrCopy-copy array ori to cop
    parameters:
        ori:original array
        cop: resultant array
        size:length of ori*/
    for (int i = 0; i < length;i++)
        *abs_data++ = fabs(*(data++));
    return 0;
}

void writeData(double *data, int length,char* name){
    /*writeData - write data to a file
    patameters:
        data:data to be written
        length:number of data to be written
        name: out file name*/
    FILE *fp;
    fp=fopen(name,"w");
    for (int i = 0; i < length; i++)
        fprintf(fp,"%lf ", data[i]);
    fclose(fp);
}

int* locationCompare(list_node *locations, int N, list_node *interval, list_node *spike_detected, int M, list_node *FP, list_node *FN, list_node *TP){
    /*locationCompare-compare the detected spike location with the ground truth and return the number of TP, FP, FN
    the if the detected spikes are within the location-10~locations+5, they will be count as a TP, 
    otherwise FP. If one location is not detected, then it is a FN
    parameters:
        locations:ground truth spike location
        N:number of real spikes
        interval:detected locations
        spike_detected: detected spike peaks
        M:number of detected locations
        FP: list to hold FP locations
        FN: list to hold FN locations
        TP: list to hold TP locatopns
    Return value:
        count: an address of an array that holds the number of TP, FP, FN*/
    list_node *loc = locations;
    list_node *intv = interval;
    list_node *idx = list_create(0);
    int TP_count = 0;
    int FN_count = M;
    int FP_count = 0;

    for (int i = 1; i < M;i++) // for indexing the list
        list_insert_end(idx, i);

    for (int i = 0; i < N; i++){ //iterate all locations
        int detected = 0;
        intv = interval;
        for (int j = 0; j < FN_count; j++){ //checking every unmatcahed detected spikes
            if (loc->data >= intv->data-5 && loc->data <= intv->data + 10){ // if matched
                if (TP_count == 0) //the first time, create list
                    *TP = *list_create(list_indexing(spike_detected,list_indexing(idx, j)));
                else //add the peak location into TP list
                    list_insert_end(TP, list_indexing(spike_detected,list_indexing(idx, j)));
                list_remove(&interval, intv);//remove the spike from the lists
                list_remove_by_index(&idx, j); //remove the detected index


                detected = 1;
                TP_count++;
                FN_count--; //number of rest detected spikes
                break;
            }
            intv = intv->next;
        }        
        if (!detected){ // is no matching to the current location, then FP
            if (FP_count == 0)
                *FP = *list_create(list_indexing(locations, i));
            else
                list_insert_end(FP, list_indexing(locations, i));
            FP_count++;
        }
        loc = loc->next;
    }
    if(idx->next){// all rest idx are FN
        *FN = *list_create(list_indexing(spike_detected, idx->data));
        while (idx->next)
        {
            list_insert_end(TP, list_indexing(spike_detected, idx->data));
            idx = idx->next;
        }
    }

    int *toReturn=(int*)malloc(__SIZEOF_INT__*3);
    toReturn[0] = TP_count;
    toReturn[1] = FP_count;
    toReturn[2] = FN_count;
    return toReturn;
}

//TEST
// int main(){
//     char loc_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/spike_location_1.txt";
//     char data_file_name[] = "D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/realDataWithLFP_1.txt";

//     int loc_start = 0;
//     int data_start = 0;
//     int num_of_data = 10000;
//     int loc[num_of_data];
//     double data[num_of_data];

//     loadData(data_start, num_of_data, data_file_name, data);
//     loadLocation(loc_start, loc_start+num_of_data, loc_file_name, loc);

//     for (int i = 0; i < 10; i++){ //print  data
//         printf("%lf ", data[i]);
//     }

// //     printf("\n");

// //     for (int i = 0; i < 10; i++){ //print locations
// //         printf("%d ", loc[i]);
// //     }
// //     return 0;
//  }