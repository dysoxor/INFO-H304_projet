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

  int16_t decr = 1;

  cout << "matrices initiales: "<< endl;
  for(int i = 0; i < 16; i++){
    a[i] = decr;
    cout << a[i] << " ";
    b[i] = 1;
    if(i == 3 || i == 7 || i == 11)
      decr--;
  }
  cout << endl;
  for(int j = 0; j < 16; j++){
    cout << b[j] << " ";
  }
  cout << endl;

  __m256i A = _mm256_load_si256((__m256i*) &a[0]);
  __m256i B = _mm256_load_si256((__m256i*) &b[0]);



  __m256i RES = _mm256_mullo_epi16(B, A);
  cout << "cmplgt" << endl;
  _mm256_store_si256((__m256i*)&a[0], RES);

  cout << "matrice finale: " << endl;
  for(int i = 0; i < 16; i++){
    cout << a[i] << " ";
  }
  cout << endl;
  return EXIT_SUCCESS;
}
