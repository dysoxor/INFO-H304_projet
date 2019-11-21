#include <stdlib.h>
#include <array>
#include <immintrin.h>
#include <stdio.h>
#include <iostream>

using namespace std;

/*int32_t sum_array(const int32_t a[], const int n)
{
    __m128i vsum = _mm_set1_epi32(0);       // initialise vector of four partial 32 bit sums
    int32_t sum;
    int i;

    for (i = 0; i < n; i += 4)
    {
        __m128i v = _mm_load_si128(&a[i]);  // load vector of 4 x 32 bit values
        vsum = _mm_add_epi32(vsum, v);      // accumulate to 32 bit partial sum vector
    }
    // horizontal add of four 32 bit partial sums and return result
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
    sum = _mm_cvtsi128_si32(vsum);
    return sum;
}

int main(){
  int32_t *array[16];
  array[0] = 3;
  array[5] = 5;
  cout << sum_array(array, 16) << endl;
}*/




/*union Line{
    short int* row;
    __m128i sRow;
};

short int* transpose_simd(short int * matrice1, short int* matrice2, int N){
    short int i=0, n=0;

    short int* __attribute__ ((aligned (16))) line1; //ligne de matrice
    short int* __attribute__ ((aligned (16))) line2;

    __m128i sLine1, sLine2; //ligne en version sse
    __m128i sLine3, sLine4;


    //There will be a loop surrounding the following code, but first I kept it simple


            line1 = matrice1 + N*i;
            line2 = matrice1 + N*(N/2 + i);

            sLine1 = _mm_loadu_si128((__m128i*) line1); //charge 1 ligne (8 nombres de 16 bits = 128, "coup de bol")
            sLine2 = _mm_loadu_si128((__m128i*) line2);



            sLine3 = _mm_unpacklo_epi16 ( sLine1, sLine2 ); //shuffle les 4 premiers chiffres de line1, voir p74
            sLine4 = _mm_unpackhi_epi16 ( sLine1, sLine2 ); //shuffle les 4 premiers chiffres de line1, voir p74

            _mm_storeu_si128((__m128i*) (matrice2 + N*2*i), sLine3);
            _mm_storeu_si128((__m128i*) (matrice2 + N*(2*i+1)), sLine4);

            return matrice2;
}*/

short int* waterman_simd(short int * matrice1, short int* matrice2, int N){
  //paddsb P[q], H // H = H + P[q]
  __m128i sLine1 =  add_epi16(mi a,mi b)
  pmaxub F, H // H = max(H, F)
  pmaxub E, H // H = max(H, E)
  pmaxub H, S // S = max(S, H)
  psubsb R, F // F = F – R
  psubsb R, E // E = E – R
  movdqa H, N // N = H
  psubsb Q, H // H = H – Q
  pmaxub H, E // E = max(H, E)
  pmaxub H, F // F = max(H, F)
}



int main(void){
    const int N = 8;
    short int  matrice1[] = {
        10, 11, 12, 13, 14, 15, 16, 17,
        20, 21, 22, 23, 24, 25, 26, 27,
        30, 31, 32, 33, 34, 35, 36, 37,
        40, 41, 42, 43, 44, 45, 46, 47,
        50, 51, 52, 53, 54, 55, 56, 57,
        60, 61, 62, 63, 64, 65, 66, 67,
        70, 71, 72, 73, 74, 75, 76, 77,
        80, 81, 82, 83, 84, 85, 86, 87
    };
    short int matrice2[N];

    short int *matrice3 = transpose_simd(matrice1, matrice2, N);
    for(int i = 0; i<N; i++)
      cout << matrice3[i] << endl;


 return 0;
}
