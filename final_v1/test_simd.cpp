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
union {
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
} sumPtr;


void matching_SIMD(int8_t* residues, int8_t* sum, int8_t* toAdd){

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
  //__m128i Xi;
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
  }

}
void sse(float* a, int N)
{
  // We assume N % 4 == 0.
  int nb_iters = N / 4;
  __m128* ptr = (__m128*)a;

  for (int i = 0; i < nb_iters; ++i, ++ptr, a += 4)
    _mm_store_ps(a, _mm_sqrt_ps(*ptr));
}



int main(int argc, char** argv)
{
  int sequence[16][6] = {12, 1, 8, 7, 5, 18,
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
  /*int residue[16]__attribute__ ((aligned (8)));

  for(int i = 0; i<16; i++){
    residue[i] = sequence[i][0];
    cout << residue[i] << " ";
  }
  cout << endl;*/
  __attribute__((aligned (8))) int8_t residue [16];
  //int32_t residue[16];

  for(int i = 0; i < 16; i++) {
      residue[i] = (int8_t)sequence[i][0];
  }


  int8_t sum[16]__attribute__ ((aligned (8)));
  int8_t toAdd[16]__attribute__ ((aligned (8)));


  matching_SIMD(residue, sum, toAdd);

  cout << endl;

  /*if (argc != 2)
    return 1;
  int N = atoi(argv[1]);

  float* a;
  posix_memalign((void**)&a, 16,  N * sizeof(float));

  for (int i = 0; i < N; ++i){
    float t = 3141592.65358;
    a[i] = t;
  }

  clock_t begin = clock();
  normal(a, N);
  clock_t end = clock();
  double interTime = double(end - begin)/CLOCKS_PER_SEC;
  cout << "The normal exection time : " << interTime << endl;

  for (int i = 0; i < N; ++i){
    float t = 3141592.65358;
    a[i] = t;
  }

  clock_t begin2 = clock();
  sse(a, N);
  clock_t end2 = clock();
  double interTime2 = double(end2 - begin2)/CLOCKS_PER_SEC;
  cout << "The sse exection time : " << interTime2 << endl;
  */
}
