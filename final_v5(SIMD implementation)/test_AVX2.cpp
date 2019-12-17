#include <immintrin.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <array>
#include <stdio.h>
#include <emmintrin.h>
using namespace std;

int main(int argc, char** argv){
  __attribute__((aligned (16))) int16_t a[16];
  __attribute__((aligned (16))) int16_t b[16];
  int decr = 1;
  for(int i = 0; i < 16; i++){
    a[i] = decr;
    b[i] = 0;
    if(i == 3 || i == 7 || i == 11)
      decr--;
  }
  __m256i A = _mm256_load_si256((__m256i*) &a);
  __m256i B = _mm256_load_si256((__m256i*) &b);


  __m256i RES = _mm256_cmpgt_epi16(B, A);

  _mm256_store_si256((__m256i*)&a, RES);

  cout << result << endl;
  for(int i = 0; i < 16; i++){
    cout << a[i] << " ";
  }
  cout << endl;
  return EXIT_SUCCESS;
}
