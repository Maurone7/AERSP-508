#include <stdio.h>

int B_n_loop(n){
    int n = 471;
    int i;
    int sum = 0;
    /* Integral of e^x*/

    for(i = 0; i < n; i++){
        sum += exp(-n);
    }
    return sum;
}


int main(void){
    printf(B_n_loop(10));

    return 0;
}
