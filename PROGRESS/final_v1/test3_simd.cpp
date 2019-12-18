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
  int16_t* colScore, int16_t* lineScore, int16_t* unit, int16_t* match, int16_t* alignement){

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


    __m128i MATCH = _mm_load_si128((__m128i*) &match[i*8]);
    __m128i ALIGN = _mm_load_si128((__m128i*) &alignement[i*8]);

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
    __m128i UNIT = _mm_load_si128((__m128i*) &unit[i*8]);


    __m128i SUBC1 = _mm_sub_epi16(COL, OP);
    __m128i SUBC2 = _mm_sub_epi16(PCOL,PMCOL);
    __m128i MULTC = _mm_mullo_epi16(EX,SUBC2);
    RESC = _mm_sub_epi16(SUBC1,MULTC);
    __m128i SUBL1 = _mm_sub_epi16(LIN, OP);
    __m128i SUBL2 = _mm_sub_epi16(PLIN,PMLIN);
    __m128i MULTL = _mm_mullo_epi16(EX,SUBL2);
    RESL = _mm_sub_epi16(SUBL1,MULTL);

    D = _mm_add_epi16(D,P);

    S = _mm_max_epi16(S,D);
    S = _mm_max_epi16(S,RESC);
    S = _mm_max_epi16(S,RESL);
    S = _mm_max_epi16(S,Z);


    __m128i CMPLTC = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, COL));
    __m128i MULTC2 = _mm_mullo_epi16(PCOL,CMPLTC);
    __m128i NEWMPC = _mm_max_epi16(PMCOL, MULTC2);

    __m128i CMPLTL = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, LIN));
    __m128i MULTL2 = _mm_mullo_epi16(PLIN,CMPLTL);
    __m128i NEWMPL = _mm_max_epi16(PMLIN, MULTL2);

    ALIGN = Z;
    ALIGN = _mm_mullo_epi16(MATCH,_mm_cmpeq_epi16(S, D));
    ALIGN = _mm_mullo_epi16(UNIT,_mm_cmpeq_epi16(S, RESC));
    ALIGN = _mm_mullo_epi16(_mm_add_epi16(UNIT, UNIT),_mm_cmpeq_epi16(S, RESL));



    _mm_store_si128((__m128i*)&score[i*8], S);
    _mm_store_si128((__m128i*)&diagonalScore[i*8], S);
    _mm_store_si128((__m128i*)&maxPosCol[i*8], NEWMPC);
    _mm_store_si128((__m128i*)&maxColVal[i*8], _mm_max_epi16(COL, S));
    _mm_store_si128((__m128i*)&maxPosLine[i*8], NEWMPL);
    _mm_store_si128((__m128i*)&maxLineVal[i*8], _mm_max_epi16(LIN, S));
    _mm_store_si128((__m128i*)&maxScore[i*8], M = _mm_max_epi16(S,M));
    _mm_store_si128((__m128i*)&colScore[i*8], RESC);
    _mm_store_si128((__m128i*)&lineScore[i*8], RESL);
    _mm_store_si128((__m128i*)&alignement[i*8], ALIGN);
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
      if(i == 0){
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
  __attribute__((aligned (16))) int16_t unit[16];
  __attribute__((aligned (16))) int16_t match[16];
  __attribute__((aligned (16))) int16_t alignement[16];


  for(int i = 0; i < 16; i++){
    zero[i] = 0;
    unit[i] = 1;
    match[i] = 0;
    alignement[i] = 0;
  }


  vector<vector<int>> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  int residue[16];

  bool done = false;
  bool firstQuery = false;

  int matriceScore[20][30];
  vector<int> allScores;
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
      if(analysedSeqIndex[0] > 29)
        WantSeq = false;
      for(int l = 0; l < 20; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < 21; i++){
                leftScore[i][j] = 0;
                if(i != 20){
                  maxPosLine[i][j] = 0;
                  maxLineVal[i][j] = 0;
                }
              }
              opGap[j] = 0;
              exGap[j] = 1;
              maxScore[j] = 0;
            }
            maxPosCol[j] = 0;
            maxColVal[j] = 0;
          }
          ColPos[j] = l;
          LinePos[j] = analysedSeqIndex[j];

          colScore[j] = 0;
          lineScore[j] = 0;
          diagGap[j] = selfMadeBlosum[query[l]][residue[j]];
          if(query[l] == residue[j])
            match[j] = 4;
          else if(diagGap[j] > 0)
            match[j] = 5;
          else
            match[j] = 3;


        }
        if(WantSeq && analysedSeqIndex[0] == 1){
          cout << "-------------------- AVANT[" << l << "]-------"<< query[l] <<"---" << residue[0] <<"-------------" << endl;
          cout << "score: " << score[0] << " diagap: " << diagGap[0] << " opGap: " << opGap[0]
          << " maxColVal: " << maxColVal[0] << " maxPosCol: " << maxPosCol[0] << " maxLineVal: "
          << maxLineVal[l][0] << "\nmaxPosLine: " << maxPosLine[l][0] << " exGap: " << exGap[0]
          << " leftScore: " << leftScore[l][0] << " maxScore: " << maxScore[0] << " zero: " << zero[0]
          << " ColPos: " << ColPos[0] << "\nLinePos: " << LinePos[0] << " colScore: " << colScore[0]
          << " lineScore: " << lineScore[0] << " unit: " << unit[0] << endl;
        }
        matching_SIMD(score, diagGap, opGap, maxColVal, maxPosCol, maxLineVal[l],
          maxPosLine[l], exGap, leftScore[l], maxScore, zero, ColPos, LinePos, colScore, lineScore, unit, match, alignement);
        /*vector<int> root;
        fo(int d = 0; d < 16; d++){
          if(score[d] != diagGap[d]){
            if(score[d] == colScore[d])
              match[d] = 1;
            else
              match[d] = 2;
          }
          root
        }*/
        if(analysedSeqIndex[0] < 30 && WantSeq)
          matriceScore[l][analysedSeqIndex[0]] = score[0];
        if(WantSeq && analysedSeqIndex[0] == 1){
          cout << "-------------------- APRES[" << l << "]-------"<< query[l] <<"---" << residue[0] <<"-------------" << endl;
          cout << "score: " << score[0] << " diagap: " << diagGap[0] << " opGap: " << opGap[0]
          << " maxColVal: " << maxColVal[0] << " maxPosCol: " << maxPosCol[0] << " maxLineVal: "
          << maxLineVal[l][0] << "\nmaxPosLine: " << maxPosLine[l][0] << " exGap: " << exGap[0]
          << " leftScore: " << leftScore[l][0] << " maxScore: " << maxScore[0] << " zero: " << zero[0]
          << " ColPos: " << ColPos[0] << "\nLinePos: " << LinePos[0] << " colScore: " << colScore[0]
          << " lineScore: " << lineScore[0] << " unit: " << unit[0] << " alignement: " << alignement[0] << endl << endl;
        }
        if(WantSeq)
          maxS = maxScore[0];
        if(l == 19){
          if(WantSeq && analysedSeqIndex[0] == 0)
            cout << "Décalage de LEFT" << endl;
          for(int f = 19; f >= 0; f--){
            if(f == 0)
              firstLine = true;
            offset(leftScore[f+1], leftScore[f], zero, firstLine);
            firstLine = false;
          }
          if(WantSeq && analysedSeqIndex[0] == 0){
            for(int f = 0; f < 20; f++){
              cout << "L : " << leftScore[f+1][0] << " D : " << leftScore[f][0] << endl;
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
          if(freePosition == -1 && i != seqContainer.size() - 1)
            freePosition = k;
        }
      }

    }
	}


  cout << "matrice score : \n";
  for(int i = 0; i < 20; i++){
    for(int j=0; j < 30; j++){
      if(matriceScore[i][j] == maxS)
       cout << "#";
      cout << matriceScore[i][j] << " ";
    }
    cout << endl;

  }

  cout << "done all" << endl;

}