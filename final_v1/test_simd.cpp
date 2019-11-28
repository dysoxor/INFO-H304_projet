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

void matching_SIMD(int8_t* maxScore, int8_t* saveScore, int* residue, int8_t* diagGap, int8_t* Score, int8_t* upGap){

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

  /*__attribute__((aligned (8))) int8_t leftGap [16];
  __attribute__((aligned (8))) int8_t upGap [16];
  __attribute__((aligned (8))) int8_t gapExtension [16];
  __attribute__((aligned (8))) int8_t gapOpen [16];*/


  //__m128i P = _mm_load_si128((__m128i*) diagGap);
  __m128i S = _mm_load_si128((__m128i*) maxScore);
  __m128i H = _mm_load_si128((__m128i*) Score);

  __m128i P = _mm_load_si128((__m128i*) diagGap);
  __m128i F = _mm_load_si128((__m128i*) upGap);
    S = _mm_add_epi8(H,P);
    _mm_store_si128((__m128i*)upGap, _mm_sub_epi8())
    _mm_store_si128((__m128i*)maxScore, _mm_add_epi8(H,P));





  /*
  __m128i E = _mm_load_si128((__m128i*) leftGap);
  __m128i F = _mm_load_si128((__m128i*) upGap);
  __m128i R = _mm_load_si128((__m128i*) gapExtension);
  __m128i Q = _mm_load_si128((__m128i*) gapOpen);
  __m128i N = _mm_load_si128((__m128i*) saveScore);*/


  /*H = _mm_add_epi8(H,P);
  H = _mm_max_epu8(H,F);
  H = _mm_max_epu8(H,E);
  S = _mm_max_epu8(S,H);
  F = _mm_sub_epi8(F,R);
  E = _mm_sub_epi8(E,R);
  N = H;
  H = _mm_sub_epi8(H,Q);
  E = _mm_max_epu8(H,E);
  F = _mm_max_epu8(H,F);*/
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
      if(incrementor == 11)
        incrementor = 0;
    }
  }

  cout << "The blosum : ";
  for(int i = 0 ; i < 12; i ++){
    for(int j = 0; j < 20; j++){
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

  vector<vector<int>> Left(16,vector<int>(21,0));
  __attribute__((aligned (8))) int8_t maxScore [16];
  __attribute__((aligned (8))) int8_t saveScore [16];
  //int8_t saveScore[16];
  //int maxScore[16];
  int maxScoreX[16];
  int maxScoreY[16];

  __attribute__((aligned (8))) int8_t diagGap[16];
  __attribute__((aligned (8))) int8_t upGap[16];
  __attribute__((aligned (8))) int8_t Score[16];


  vector<vector<int>> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  int freePosition = 0;
  int residue[16];

  bool done = false;
  bool firstQuery = true;

  for(int i = 0; i < seqContainer.size(); i++){
    if(freePosition != -1){
      analysedSequences[freePosition] = seqContainer[i];
      len[freePosition] = seqContainer[i].size();
      analysedSeqIndex[freePosition] = 0;
      freePosition = -1;

      for(int k = 0; k < 16; k++){
        if(len[k] != 0){
          residue[k] = analysedSequences[k][analysedSeqIndex[k]];
        }
        else{
          residue[k] = -1;
          if(freePosition == -1)
            freePosition = k;
        }
      }
    }
    while(!done && (freePosition == -1 || i == seqContainer.size()-1)){
      done = true;

      for(int l = 0; l < 20; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            upGap[j] = 1;
          }
          else{
            upGap[j] = 11;
          }
          diagGap[j] = selfMadeBlosum[query[l]][residue[j]];
          if(firstResidu && firstQuery){
            cout << "diagGap " << j << " = " << static_cast<int16_t>(diagGap[j]) << " blosum[" << query[l] << "," << residue[j] << "] = "<< selfMadeBlosum[query[l]][residue[j]] <<endl;
          }

        }
        matching_SIMD(maxScore,saveScore,residue, diagGap, Score, upGap);
        if(firstResidu && firstQuery){
          //printf("[i:%d] maxScore ", i);
          for(int j = 0; j < 16; j++) {
              //printf("[j:%d]%d ", j, static_cast<int16_t>(saveScore[j]));
              printf("maxScore[%d] = %d\n", j, static_cast<int16_t>(maxScore[j]));
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
