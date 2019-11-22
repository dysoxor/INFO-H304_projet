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
void matching_SIMD(int* residues, int* sum, int* toAdd){

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


  __m128* residuePtr = (__m128*)residues;
  __m128* sumPtr = (__m128*)sum;
  __m128* addPtr = (__m128*)toAdd;

  //_mm_store_ps(sum, _mm_add_ps(*residuePtr,*addPtr));
  _mm_store_ps1(sum, _mm_add_ps(*residuePtr,*addPtr));


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
  int sequence[][] = {12, 1, 8, 7, 5, 18,
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
                }
  int* residue;
  posix_memalign((void**)&residue, 8,  16 * sizeof(int));
  for(int i = 0; i<16; i++){
    residue[i] = sequence[i][0];
  }
  int* sum;
  posix_memalign((void**)&sum, 8,  16 * sizeof(int));
  int* toAdd;
  posix_memalign((void**)&toAdd, 8,  16 * sizeof(int));
  for(int i = 0; i<16; i++){
    toAdd[i] = 2;
  }
  matching_SIMD(residue, sum, toAdd);
  for(int i = 0; i<16; i++){
    cout << sum[i] << endl;
  }

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
