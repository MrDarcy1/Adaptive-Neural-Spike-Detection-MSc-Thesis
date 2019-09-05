
#include <stdio.h>
#include <time.h>



int main(){
//     clock_t start, finish;
//     double duration;

//     double a = 2.1;
//     double b = 4.112;
//     double c = 3.1123;
//     double d;
//     start = clock();
//     for (long long i = 0; i < 10000000000;i++){
//             d = 1*b-(a*c);  
//     }
//     finish= clock();
//     duration = (double)(finish- start)/CLOCKS_PER_SEC;
//     printf("duration 1: %lf sec, b = %lf \n",duration,d);

//     start = clock();
//     for (long long  i = 0; i < 10000000000;i++){
//             d = b-(a*c); 
//     }
//     finish= clock();
//     duration = (double)(finish- start)/CLOCKS_PER_SEC;
//     printf("duration 2: %lf sec, b = %f \n",duration,d);
unsigned short v = 10;
unsigned short c = 0;
__asm__("CLZ c,v");

printf("%d", c);
}
