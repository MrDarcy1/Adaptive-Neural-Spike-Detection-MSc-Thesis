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
unsigned short Count = 0;
short update = 1;
short hold_data = 0;
short detected = 0;
short adding = 0;
short Sum = 0;


    unsigned short mean_buffer_end = 0;
    unsigned short mean_buffer_mean = 0;
    int previous_demean = 0;
    unsigned short mean_buffer[16]; // a ring buffer
    
    unsigned short aso = 0;
    
    unsigned short thr_buffer_end = 0;
    unsigned short thr_buffer_mean = 0;
    unsigned short threshold = 0;
    
//    short hold_time = 0.0005;
//    short update_period = 0.1;
//    short detected_time = 0.0005;
    
unsigned short Median(unsigned short a[],int N)
{
    int i,j;
    unsigned short max;
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
        Count++;    
        unsigned short temp1 =  mean_buffer[mean_buffer_end];
        unsigned short temp2 = raw.read_u16();
        mean_buffer[mean_buffer_end] = temp2;
        int demean = (int)mean_buffer[mean_buffer_end]-(int)mean_buffer_mean;
        int tp = (int)temp1 - (int)temp2;
        if(tp > 0)
            mean_buffer_mean = (unsigned short)((int)mean_buffer_mean - (tp>>4));
         else
            mean_buffer_mean = (unsigned short)((int)mean_buffer_mean + (abs(tp)>>4));

        if(Count == 3500){
            int sum = 0;
            for(int i = 0; i < 16; i++)
                sum+=(int)mean_buffer[i];
            mean_buffer_mean = sum>>4;
            sum = 0;
        }
        

        if(mean_buffer_end + 1 == 16)
            mean_buffer_end = 0;
        else 
            mean_buffer_end++;
        /*ASO*/
        int temp3 = previous_demean;
        previous_demean = demean;
        unsigned short c = 0;
        unsigned short v= abs(demean);
        __asm__{
            CLZ c, v
        }
        c = 32-c;
        aso = ((abs(demean - temp3)<<c)>>9);      
        if(detected == 6){
            indicator = 0;
            
            if(aso > threshold){
                detected = 0;
                indicator = 1;
            }
                
            if(update < 4300)
                update ++;
            else if(adding < 64){
                if(aso<=threshold >> 1){
                    Sum += aso >> 6;
                    adding++;
                }
            }
            else{
                if(aso <= threshold >> 1){
                    Sum += aso>>6;
                    threshold = Sum<<4+Sum<<3;
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
    int sum;
    unsigned short thr_buffer[64]; // a ring buffer

//fill buffers
    for(int i = 0; i < 16; i++){
        mean_buffer[mean_buffer_end] = raw.read_u16();
        sum += mean_buffer[i];
        wait(0.0005);
    }
    mean_buffer_mean = sum>>4;
    sum = 0;
    for(int j = 1; j < 64; j++){
        /*subtract mean*/
        unsigned short temp1 =  mean_buffer[mean_buffer_end];
        unsigned short temp2 = raw.read_u16();
        mean_buffer[mean_buffer_end] = temp2;
//        int demean = mean_buffer[mean_buffer_end]-mean_buffer_mean;
//        mean_buffer_mean = mean_buffer_mean - ((temp1-temp2)>>4);
        int demean = (int)mean_buffer[mean_buffer_end]-(int)mean_buffer_mean;

        int tp = (int)temp1 - (int)temp2;
        if(tp > 0)
            mean_buffer_mean = (unsigned short)((int)mean_buffer_mean - (tp>>4));
         else
            mean_buffer_mean = (unsigned short)((int)mean_buffer_mean + (abs(tp)>>4));
        if(mean_buffer_end + 1 == 16)
            mean_buffer_end = 0;
        else 
            mean_buffer_end++;
        /*ASO*/
        int temp3 = previous_demean;
        previous_demean = demean;
        unsigned short c = 0;
        unsigned short v= abs(demean);
//        v |= (v >> 1);
//        v |= (v >> 2);
//        v |= (v >> 4);
//        v |= (v >> 8);
//        v |= (v >> 16);
//        v= v - (v >> 1);
        __asm__{
            CLZ c, v
        }
        c = 32-c;
        aso = ((abs(demean - temp3)<<c)>>9);     
//        aso = demean * ( demean - previous_demean );
        aso_out.write((float)(aso) * (1.0f / (float)0xFFFF));
        thr_buffer[j] = aso;
        wait(0.0005);
    }
    unsigned short median = Median(thr_buffer, 64);

    threshold =  (median<<2)+median<<1; 

    thr.write(threshold);

    Sampling.attach(&iter, 0.00014); // sample in every 0.4ms
    
    while(1){}
}


