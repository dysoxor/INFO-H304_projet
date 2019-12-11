#include <emmintrin.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <array>
#include <stdio.h>
#include <immintrin.h>

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

void offset(int16_t* leftScore, int16_t* diagonalScore, int16_t* zero, bool firstLine){
  for(int i = 0; i < 2; i++){

    __m128i L = _mm_load_si128((__m128i*) &leftScore[i*8]);
    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    _mm_store_si128((__m128i*)&leftScore[i*8], D);
    if(firstLine){
      __m128i Z = _mm_load_si128((__m128i*) &zero[i*8]);
      _mm_store_si128((__m128i*)&diagonalScore[i*8], Z);
    }
  }
}

void matching_SIMD(int16_t* score, int16_t* diagGap, int16_t* opGap,
  int16_t* maxColVal, int16_t* maxPosCol, int16_t* maxLineVal,
  int16_t* maxPosLine, int16_t* exGap, int16_t* diagonalScore,
  int16_t* maxScore, int16_t* zero, int16_t* PosCol, int16_t* PosLine,
  int16_t* colScore, int16_t* lineScore){

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

    __m128i Z = _mm_load_si128((__m128i*) &zero[i*8]);
    __m128i RESC = _mm_load_si128((__m128i*) &colScore[i*8]);
    __m128i RESL = _mm_load_si128((__m128i*) &lineScore[i*8]);
    __m128i S = _mm_load_si128((__m128i*) &score[i*8]);

    __m128i OP = _mm_load_si128((__m128i*) &opGap[i*8]);
    __m128i COL = _mm_load_si128((__m128i*) &maxColVal[i*8]);
    __m128i PMCOL = _mm_load_si128((__m128i*) &maxPosCol[i*8]);
    __m128i PCOL = _mm_load_si128((__m128i*) &PosCol[i*8]);

    __m128i EX = _mm_load_si128((__m128i*) &exGap[i*8]);
    __m128i LIN = _mm_load_si128((__m128i*) &maxLineVal[i*8]);
    __m128i PMLIN = _mm_load_si128((__m128i*) &maxPosLine[i*8]);
    __m128i PLIN = _mm_load_si128((__m128i*) &PosLine[i*8]);

    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    __m128i P = _mm_load_si128((__m128i*) &diagGap[i*8]);

    __m128i M = _mm_load_si128((__m128i*) &maxScore[i*8]);


    RESC = _mm_sub_epi16(_mm_sub_epi16(COL, OP),_mm_mulhi_epi16(EX,_mm_sub_epi16(PCOL,PMCOL)));
    RESL = _mm_sub_epi16(_mm_sub_epi16(LIN, OP),_mm_mulhi_epi16(EX,_mm_sub_epi16(PLIN,PMLIN)));
    S = _mm_add_epi16(D,P);

    S = _mm_max_epi16(S,RESC);
    S = _mm_max_epi16(S,RESL);
    S = _mm_max_epi16(S,Z);

    _mm_store_si128((__m128i*)&score[i*8], S);
    _mm_store_si128((__m128i*)&diagonalScore[i*8], S);
    _mm_store_si128((__m128i*)&maxPosCol[i*8], _mm_mulhi_epi16(PCOL,_mm_cmplt_epi16(COL, S)));
    _mm_store_si128((__m128i*)&maxColVal[i*8], _mm_max_epi16(COL, S));
    _mm_store_si128((__m128i*)&maxPosLine[i*8], _mm_mulhi_epi16(PLIN,_mm_cmplt_epi16(LIN, S)));
    _mm_store_si128((__m128i*)&maxLineVal[i*8], _mm_max_epi16(LIN, S));
    _mm_store_si128((__m128i*)&maxScore[i*8], M = _mm_max_epi16(S,M));
    _mm_store_si128((__m128i*)&colScore[i*8], RESC);
    _mm_store_si128((__m128i*)&lineScore[i*8], RESL);

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
      cout << "[" << i << "-" << j << "]";
      cout << selfMadeBlosum[i][j];

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
  cout << "Matrice to build : " << endl;
  for(int i = 0; i < 17; i++){
    for(int j = 0; j < lenSeq[i]; j++){
      seq.push_back(incrementor);
      if(i == 1){
        if(j == 0)
          cout << " ";
        cout<< " " << incrementor;
      }
      if(incrementor<19)
        incrementor++;
      else
        incrementor = 0;
    }
    seqContainer.push_back(seq);
    seq.clear();
  }
  cout << endl;
  for(int i = 0; i < 20; i++){
    cout << query[i] << endl;
  }
  cout << endl;

  __attribute__((aligned (16))) int16_t score [16];
  //int8_t saveScore[16];
  //int maxScore[16];
  int maxScoreX[16];
  int maxScoreY[16];

  int indSearchSeq = 1;
  int lenSearchSeq = lenSeq[indSearchSeq];

  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t opGap[16];
  __attribute__((aligned (16))) int16_t HScore [16];
  __attribute__((aligned (16))) int16_t leftScore[21][16];
  __attribute__((aligned (16))) int16_t exGap[16];
  __attribute__((aligned (16))) int16_t maxScore[16];
  __attribute__((aligned (16))) int16_t zero[16];

  __attribute__((aligned (16))) int16_t maxPosLine[20][16];
  __attribute__((aligned (16))) int16_t maxPosCol[16];
  __attribute__((aligned (16))) int16_t maxColVal[16];
  __attribute__((aligned (16))) int16_t maxLineVal[20][16];
  __attribute__((aligned (16))) int16_t ColPos[16];
  __attribute__((aligned (16))) int16_t LinePos[16];
  __attribute__((aligned (16))) int16_t colScore[16];
  __attribute__((aligned (16))) int16_t lineScore[16];



  for(int i = 0; i < 16; i++){
    zero[i] = 0;
  }


  vector<vector<int>> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  int residue[16];

  bool done = false;
  bool firstQuery = false;

  int matriceScore[20][lenSearchSeq];
  bool WantSeq = true;
  bool firstLine = false;

  int maxS;

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
      if(analysedSeqIndex[indSearchSeq] > lenSearchSeq)
        WantSeq = false;
      for(int l = 0; l < 20; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < 21; i++){
                leftScore[i][j] = 0;
                maxLineVal[i][j] = 0;
                maxLineVal[i][j] = 0;
              }
              opGap[j] = 0;
              exGap[j] = 1;
              maxScore[j] = 0;
            }
            HScore[j] = 0;
            maxColVal[j] = 0;
            maxPosCol[j] = 0;
            colScore[j] = 0;
            lineScore[j] = 0;
          }
          ColPos[j] = l;
          LinePos[j] = analysedSeqIndex[j];
          diagGap[j] = selfMadeBlosum[query[l]][residue[j]];
          /*if(firstResidu && firstQuery){
            cout << "diagGap " << j << " = " << static_cast<int16_t>(diagGap[j]) << " blosum[" << query[l] << "," << residue[j] << "] = "<< selfMadeBlosum[query[l]][residue[j]];
            printf(" leftScore[%d,%d] = %d leftScore[%d,%d] = %d\n", l+1, j, static_cast<int16_t>(leftScore[l+1][j]), l, j, static_cast<int16_t>(leftScore[l][j]));
          }*/

        }
        if(WantSeq && analysedSeqIndex[indSearchSeq] == 0){
          cout << "blosum[" << query[l] << "," << residue[indSearchSeq]
          << "]" << diagGap[indSearchSeq] << " maxY : "  << maxColVal[indSearchSeq]
          << " maxPosY : " << maxPosCol[indSearchSeq] << " PosY : " <<
          ColPos[indSearchSeq] << " colScore : " << colScore[indSearchSeq] << endl;
        }
        matching_SIMD(score, diagGap, opGap, maxColVal, maxPosCol, maxLineVal[l],
          maxPosLine[l], exGap, leftScore[l], maxScore, zero, ColPos, LinePos, colScore, lineScore);
        if(analysedSeqIndex[indSearchSeq] < lenSearchSeq && WantSeq)
          matriceScore[l][analysedSeqIndex[indSearchSeq]] = score[indSearchSeq];
        if(WantSeq && analysedSeqIndex[indSearchSeq] == 0){
          maxS = maxScore[indSearchSeq];
          //cout << "Blosum[" << query[l] << "-" << residue[0] << "] " << selfMadeBlosum[query[l]][residue[0]]
          //<< " Diagonal[" << l << "] " << leftScore[l][0] << " Hscore : " << HScore[0] << endl;
          cout << "Query[" << l << "] "<< query[l]<< " seq[" << analysedSeqIndex[indSearchSeq] << "] "
          << residue[indSearchSeq] << " score : " << score[indSearchSeq] << " maxScore : " << maxScore[indSearchSeq]
          << " left : " << leftScore[l+1][indSearchSeq] << " diagonal : " << leftScore[l][indSearchSeq] <<  endl << endl;
        }
        if(l == 19){
          if(WantSeq && analysedSeqIndex[indSearchSeq] == lenSearchSeq -1)
            cout << "Décalage de LEFT" << endl;
          for(int f = 19; f >= 0; f--){
            if(f == 0)
              firstLine = true;
            offset(leftScore[f+1], leftScore[f], zero, firstLine);
            firstLine = false;
          }
          if(WantSeq && analysedSeqIndex[indSearchSeq] == 0){
            for(int f = 0; f < 20; f++){
              cout << "L : " << leftScore[f+1][indSearchSeq] << " D : " << leftScore[f][indSearchSeq] << endl;
            }
          }
        }
        /*if(firstResidu && firstQuery){
          //printf("[i:%d] maxScore ", i);
          for(int j = 0; j < 16; j++) {
              //printf("[j:%d]%d ", j, static_cast<int16_t>(saveScore[j]));
              printf("maxScore[%d] = %d score[%d] = %d HScore[%d] = %d upGap[%d] = %d leftScore[%d,%d] = %d leftScore[%d,%d] = %d\n", j, static_cast<int16_t>(maxScore[j]), j, static_cast<int16_t>(score[j]), j, static_cast<int16_t>(HScore[j]), j, static_cast<int16_t>(upGap[j]), l+1, j, static_cast<int16_t>(leftScore[l+1][j]), l, j, static_cast<int16_t>(leftScore[l][j]));
          }
          printf("\n");
        }*/
        //firstQuery = false;
      }


      /*if(firstResidu){
        cout << "residue : ";
        for(int m = 0; m< 16; m++){
          cout << residue[m] << " ";
        }
        cout << endl;
      }*/
      //firstResidu = false;
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


  cout << "matrice score : \n";
  for(int i = 0; i < 20; i++){
    for(int j=0; j < lenSearchSeq; j++){
      if(matriceScore[i][j] == maxS)
       cout << "#";
      cout << matriceScore[i][j] << " ";
    }
    cout << endl;

  }

  cout << "done all" << endl;

}
