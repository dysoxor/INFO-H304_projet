#include <emmintrin.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <array>
#include <stdio.h>

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
} H;*/
/*union {
    __m128i m128;
    __m128i* m128s;
    int8_t i8[16];
} S;*/
int selfMadeBlosum[12][20];
bool firstResidu = true;

void matching_SIMD(int16_t* score, int* residue, int16_t* diagGap, int16_t* upGap,
  int16_t* HScore, int16_t* leftScore, int16_t* leftGap, int16_t* diagonalScore, int16_t* maxScore){

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

  for(int i = 0; i < 2; i++){

    __m128i S = _mm_load_si128((__m128i*) &score[i*8]);

    __m128i R = _mm_load_si128((__m128i*) &upGap[i*8]);
    __m128i F = _mm_load_si128((__m128i*) &HScore[i*8]);

    __m128i E = _mm_load_si128((__m128i*) &leftScore[i*8]);
    __m128i GL = _mm_load_si128((__m128i*) &leftGap[i*8]);

    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    __m128i P = _mm_load_si128((__m128i*) &diagGap[i*8]);

    __m128i M = _mm_load_si128((__m128i*) &maxScore[i*8]);

    F = _mm_add_epi16(F,R);
    E = _mm_add_epi16(E,GL);
    S = _mm_add_epi16(D,P);

    S = _mm_max_epi16(S,F);
    S = _mm_max_epi16(S,E);

    _mm_store_si128((__m128i*)&score[i*8], S);
    _mm_store_si128((__m128i*)&leftScore[i*8], S);
    _mm_store_si128((__m128i*)&HScore[i*8], S);
    _mm_store_si128((__m128i*)&maxScore[i*8], M = _mm_max_epi16(S,M));

  }



}





int main(int argc, char** argv)
{
  int incrementor = 0;
  for(int i = 0 ; i < 12; i ++){
    for(int j = 0; j < 20; j++){
      selfMadeBlosum[i][j] = incrementor;
      incrementor++;
      if((i+j)%3 == 0)
        incrementor*=-1;
      if(incrementor == 11)
        incrementor = 0;
    }
  }

  cout << "The blosum : \n";
  for(int i = 0 ; i < 12; i ++){
    for(int j = 0; j < 20; j++){
      if(selfMadeBlosum[i][j] < 10 && selfMadeBlosum[i][j]>= 0)
        cout << "  ";
      else if((selfMadeBlosum[i][j] > -10 && selfMadeBlosum[i][j]< 0) || selfMadeBlosum[i][j] >= 10)
        cout << " ";
      cout << selfMadeBlosum[i][j] <<" | ";

    }
    cout << endl;
  }
  int query[20] = {11, 1, 8, 7, 5, 3,
                5, 7, 9, 1, 10, 8,
                8, 1, 11, 7, 9, 8,
                7, 4};

  vector<vector<int>> seqContainer;
  vector<int> seq;
  int lenSeq[17] = {
    30, 29, 17, 23, 11, 18, 19, 54, 10, 3, 22, 9, 32, 22, 43, 66, 19
  };
  incrementor = 0;
  for(int i = 0; i < 17; i++){
    for(int j = 0; j < lenSeq[i]; j++){
      seq.push_back(incrementor);
      if(incrementor<19)
        incrementor++;
      else
        incrementor = 0;
    }
    seqContainer.push_back(seq);
    seq.clear();
  }

  __attribute__((aligned (16))) int16_t score [16];
  //int8_t saveScore[16];
  //int maxScore[16];
  int maxScoreX[16];
  int maxScoreY[16];

  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t upGap[16];
  __attribute__((aligned (16))) int16_t HScore [16];
  __attribute__((aligned (16))) int16_t leftScore[21][16];
  __attribute__((aligned (16))) int16_t leftGap[16];
  __attribute__((aligned (16))) int16_t maxScore[16];


  vector<vector<int>> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  int residue[16];

  bool done = false;
  bool firstQuery = false;

  for(int i = 0; i < seqContainer.size(); i++){
    if(freePosition != -1){

      analysedSequences[freePosition] = seqContainer[i];
      len[freePosition] = seqContainer[i].size();
      analysedSeqIndex[freePosition] = 0;

      freePosition = -1;

      for(int t = 0; t < 16; t++){

        if(len[t] > 0){
          residue[t] = analysedSequences[t][analysedSeqIndex[t]];
        }
        else{
          residue[t] = -1;
          if(freePosition == -1)
            freePosition = t;
        }
      }

    }
    while(!done && (freePosition == -1 || i == seqContainer.size()-1)){
      done = true;

      for(int l = 0; l < 20; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < 21; i++){
                leftScore[i][j] = 0;
              }
              leftGap[j] = -1;
              maxScore[j] = 0;
            }
            else{
              leftGap[j] = -11;
            }
            HScore[j] = 0;
            upGap[j] = -1;
          }
          else if(l == 1){
            upGap[j] = -11;
          }
          if(l == 0){
            firstQuery = true;
          }

          diagGap[j] = selfMadeBlosum[query[l]][residue[j]];
          if(firstResidu && firstQuery){
            cout << "diagGap " << j << " = " << static_cast<int16_t>(diagGap[j]) << " blosum[" << query[l] << "," << residue[j] << "] = "<< selfMadeBlosum[query[l]][residue[j]];
            printf(" leftScore[%d,%d] = %d leftScore[%d,%d] = %d\n", l+1, j, static_cast<int16_t>(leftScore[l+1][j]), l, j, static_cast<int16_t>(leftScore[l][j]));
          }

        }
        matching_SIMD(score,residue, diagGap, upGap, HScore, leftScore[l+1], leftGap, leftScore[l], maxScore);
        if(firstResidu && firstQuery){
          //printf("[i:%d] maxScore ", i);
          for(int j = 0; j < 16; j++) {
              //printf("[j:%d]%d ", j, static_cast<int16_t>(saveScore[j]));
              printf("maxScore[%d] = %d score[%d] = %d HScore[%d] = %d upGap[%d] = %d leftScore[%d,%d] = %d leftScore[%d,%d] = %d\n", j, static_cast<int16_t>(maxScore[j]), j, static_cast<int16_t>(score[j]), j, static_cast<int16_t>(HScore[j]), j, static_cast<int16_t>(upGap[j]), l+1, j, static_cast<int16_t>(leftScore[l+1][j]), l, j, static_cast<int16_t>(leftScore[l][j]));
          }
          printf("\n");
        }
        firstQuery = false;
      }


      if(firstResidu){
        cout << "residue : ";
        for(int m = 0; m< 16; m++){
          cout << residue[m] << " ";
        }
        cout << endl;
      }
      firstResidu = false;
      for(int k = 0 ; k < 16; k++){
        analysedSeqIndex[k]++;
        if(len[k] > analysedSeqIndex[k]){
          residue[k] = analysedSequences[k][analysedSeqIndex[k]];
          done = false;
        }
        else{
          residue[k] = -1;
          if(freePosition == -1)
            freePosition = k;
        }
      }

    }
	}



  cout << "done all" << endl;

}
