/*wrote @ 06/29/2019
wrote by Zheng Zhang
This program is capable for visualisation with ociliscope, if one only wants spike locations computation can be even reduced
i.e. keep slience when holding data*/

/*Notes: mini out ~3.3mV */
/*Sampling freq of source 24414 Hz, Vpp: 500mV~2V, offset 500mV -> offset is important*/

#include "mbed.h"

Serial pc(USBTX, USBRX);
AnalogIn raw(A0);
AnalogOut aso_out(DAC0_OUT);
PwmOut thr(A4);
DigitalOut indicator(A2);
//DigitalOut red(LED1);
//DigitalOut green(LED2);
//DigitalOut blue(LED3);



Ticker Sampling;
float update = 1;
float hold_data = 0;
float detected = 0;
float adding = 0;
float Sum = 0;


    int mean_buffer_end = 0;
    float mean_buffer_mean = 0;
    float previous_demean = 0;
    float mean_buffer[16]; // a ring buffer
    
    float aso = 0;
    
    int thr_buffer_end = 0;
    float thr_buffer_mean = 0;
    float threshold = 0;
    
//    float hold_time = 0.0005;
//    float update_period = 0.1;
//    float detected_time = 0.0005;
    
float Median(float a[],int N)
{
    int i,j;
    int max;
    int t;
    for(i=0;i<N-1;i++)
    {
        max=i;
        for(j=i+1;j<N;j++)
        if(a[j]>a[max]) max=j;
        t=a[i];a[i]=a[max];a[max]=t;
    }
    return a[(N-1)/2];
} 
    void iter(){
        float temp1 =  mean_buffer[mean_buffer_end];
        float temp2 = raw.read();
        mean_buffer[mean_buffer_end] = temp2;
        float demean = mean_buffer[mean_buffer_end]-mean_buffer_mean;
        mean_buffer_mean = mean_buffer_mean - (temp1 - temp2)/16;

        

        if(mean_buffer_end + 1 == 16)
            mean_buffer_end = 0;
        else 
            mean_buffer_end++;
        /*ASO*/
        aso = abs(demean - previous_demean)*demean;   
        previous_demean = demean;    
        if(detected == 6){
            indicator = 0;
            
            if(aso > threshold){
                detected = 0;
                indicator = 1;
            }
                
            if(update < 4300)
                update ++;
            else if(adding < 64){
                if(aso<=threshold/2){
                    Sum += aso;
                    adding++;
                }
            }
            else{
                if(aso <= threshold/2){
                    Sum += aso;
                    threshold = Sum*0.357; //24/64 = 0.357
                    thr.write((float)threshold * (1.0f / (float)0x7FFF));
                    Sum = 0;
                    adding = 0;
                    update = 0;
                }
            }
        }
        else{
            detected ++;
            update ++;
        }

        aso_out.write((float)aso * (1.0f / (float)0xFFFF));


        }  
int main(){
//    red = 1;
//    green = 1;
    thr.period(0.0005); //pwm period
    float sum;
    float thr_buffer[64]; // a ring buffer

//fill buffers
    for(int i = 0; i < 16; i++){
        mean_buffer[mean_buffer_end] = raw.read_u16();
        sum += mean_buffer[i];
        wait(0.0005);
    }
    mean_buffer_mean = sum/16;
    sum = 0;
    for(int j = 1; j < 64; j++){
        /*subtract mean*/
        float temp1 =  mean_buffer[mean_buffer_end];
        float temp2 = raw.read();
        mean_buffer[mean_buffer_end] = temp2;
        float demean = mean_buffer[mean_buffer_end]-mean_buffer_mean;
        mean_buffer_mean = mean_buffer_mean - (temp1 - temp2)/16;

        if(mean_buffer_end + 1 == 16)
            mean_buffer_end = 0;
        else 
            mean_buffer_end++;
        /*ASO*/
        aso = abs(demean - previous_demean)*demean;   
        aso_out.write((float)(aso) * (1.0f / (float)0xFFFF));
        thr_buffer[j] = aso;
        wait(0.0005);
    }
    float median = Median(thr_buffer, 64);

    threshold =  median*6; 

    thr.write(threshold);

    Sampling.attach(&iter, 0.00014); // sample in every 0.4ms
    
    while(1){}
}


