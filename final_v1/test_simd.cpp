#include <emmintrin.h>
#include <iostream>
#include <math.h>
#include <ctime>

using namespace std;


/*int main(){
  float a[] __attribute__ ((aligned (16))) = { 41982.,  81.5091, 3.14, 42.666 };
  __m128* ptr = (__m128*)a;
  __m128 t = _mm_sqrt_ps(*ptr);
  float *p = (float*) &t;
  cout << p[0] << endl;
  return EXIT_SUCCESS;
}*/
/*union {
    __m128i m128;
    int8_t i8[16];
} residuePtr;
union {
    __m128i m128;
    int8_t i8[16];
} toAddPtr;
union {
    __m128i m128;
    __m128i* m128s;
    int8_t i8[16];
} sumPtr;*/


/*union {
    __m128i m128;
    __m128i* m128s;
    int8_t i8[16];
} H;
union {
    __m128i m128;
    __m128i* m128s;
    int8_t i8[16];
} P;*/
int selfMadeBlosum[12][20];

void matching_SIMD(int8_t* Score, int8_t* diagGap, int8_t* leftGap,
                  int8_t* bestScore, int8_t* upGap, int8_t* gapExtension,
                  int8_t* gapOpen, int8_t* saveScore){

  /*--------------------------------------------------------------------------
  paddsb P[q], H // H = H + P[q]
  pmaxub F, H // H = max(H, F)
  pmaxub E, H // H = max(H, E)
  pmaxub H, S // S = max(S, H)
  psubsb R, F // F = F – R
  psubsb R, E // E = E – R
  movdqa H, N // N = H
  psubsb Q, H // H = H – Q
  pmaxub H, E // E = max(H, E)
  pmaxub H, F // F = max(H, F)
  ----------------------------------------------------------------------------*/
  __m128i H = _mm_load_si128((__m128i*) Score);
  __m128i P = _mm_load_si128((__m128i*) diagGap);
  __m128i E = _mm_load_si128((__m128i*) leftGap);
  __m128i S = _mm_load_si128((__m128i*) bestScore);
  __m128i F = _mm_load_si128((__m128i*) upGap);
  __m128i R = _mm_load_si128((__m128i*) gapExtension);
  __m128i Q = _mm_load_si128((__m128i*) gapOpen);
  __m128i N = _mm_load_si128((__m128i*) saveScore);


  H = _mm_add_epi8(H,P);
  H = _mm_max_epu8(H,F);
  H = _mm_max_epu8(H,E);
  S = _mm_max_epu8(S,H);
  F = _mm_sub_epi8(F,R);
  E = _mm_sub_epi8(E,R);
  N = H;
  H = _mm_sub_epi8(H,Q);
  E = _mm_max_epu8(H,E);
  F = _mm_max_epu8(H,F);
  /*
  residuePtr.m128 = _mm_load_si128((__m128i*) residues);
  toAddPtr.m128 = _mm_load_si128((__m128i*) toAdd);
  sumPtr.m128 = _mm_load_si128((__m128i*) sum);

  // show content of Xi and Ai
  for(int j = 0; j < 16; j++) {
      printf("residuePtr[%d] = %d\t residue[%d] = %d\n", j, residuePtr.i8[j], j, residues[j]);
  }

  //_mm_store_ps(sum, _mm_add_ps(*residuePtr,*addPtr));
  _mm_store_si128((__m128i*)sum, _mm_add_epi8(residuePtr.m128,toAddPtr.m128));
  sumPtr.m128 = _mm_load_si128((__m128i*) sum);
  for(int j = 0; j < 16; j++) {
      printf("sumPtr[%d] = %d\t residue[%d] = %d\n", j, sumPtr.i8[j], j, static_cast<int16_t>(sum[j]));
  }*/

}





int main(int argc, char** argv)
{
  int incrementor = 0;
  for(int i = 0 ; i < 12; i ++){
    for(int j = 0; j < 20; j++){
      selfMadeBlosum[i][j] = incrementor;
      incrementor++;
    }
  }



  /*int sequence[16][6] = {12, 1, 8, 7, 5, 18,
                5, 7, 9, 1, 13, 8,
                14, 13, 11, 7, 9, 8,
                12, 4, 13, 7, 12, 1,
                19, 18, 15, 1, 13, 10,
                12, 1, 8, 7, 5, 18,
                5, 7, 9, 1, 13, 8,
                14, 13, 11, 7, 9, 8,
                12, 4, 13, 7, 12, 1,
                19, 18, 15, 1, 13, 10,
                12, 1, 8, 7, 5, 18,
                5, 7, 9, 1, 13, 8,
                14, 13, 11, 7, 9, 8,
                12, 4, 13, 7, 12, 1,
                19, 18, 15, 1, 13, 10,
                16, 7, 17, 18, 8, 9
              };

  __attribute__((aligned (8))) int8_t residue [16];
  //int32_t residue[16];

  for(int i = 0; i < 16; i++) {
      residue[i] = (int8_t)sequence[i][0];
  }


  int8_t sum[16]__attribute__ ((aligned (8)));
  int8_t toAdd[16]__attribute__ ((aligned (8)));
  for(int i = 0 ; i < 16; i++)
    toAdd[i] = 1;

  toAdd[2] = 5;

  matching_SIMD(residue, sum, toAdd);*/

}
