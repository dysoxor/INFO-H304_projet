#include <emmintrin.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <array>
#include <stdio.h>
#include <immintrin.h>

using namespace std;


int selfMadeBlosum[12][20];

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

  int query[10] = {11, 1, 8, 7, 5, 3,
                5, 7, 9, 1};

  int len1 = 10;

  vector<vector<int>> seqContainer;
  vector<int> seq;
  int lenSeq[17] = {
    7, 5, 10, 8, 11, 4, 7, 9, 10, 3, 6, 4, 12, 3, 4, 6, 9
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
  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t opGap[16];
  __attribute__((aligned (16))) int16_t leftScore[len1+1][16];
  __attribute__((aligned (16))) int16_t exGap[16];
  __attribute__((aligned (16))) int16_t maxScore[16];
  __attribute__((aligned (16))) int16_t zero[16];

  __attribute__((aligned (16))) int16_t maxPosLine[len1][16];
  __attribute__((aligned (16))) int16_t maxPosCol[16];
  __attribute__((aligned (16))) int16_t maxColVal[16];
  __attribute__((aligned (16))) int16_t maxLineVal[len1][16];
  __attribute__((aligned (16))) int16_t ColPos[16];
  __attribute__((aligned (16))) int16_t LinePos[16];
  __attribute__((aligned (16))) int16_t colScore[16];
  __attribute__((aligned (16))) int16_t lineScore[16];
  __attribute__((aligned (16))) int16_t unit[16];
  __attribute__((aligned (16))) int16_t match[16];


  for(int i = 0; i < 16; i++){
    zero[i] = 0;
    unit[i] = 1;
    match[i] = 0;
  }


  vector<int> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  int residue[16];

  bool done = false;

  bool firstLine = false;
  int Root[16][len1+1][lenSeq[12]+1];
  for(int i = 0; i < 16; i++){
    for(int j = 0; j < len1+1; j++)
      Root[i][j][0] = 0;
    for(int k = 0; k < lenSeq[12] + 1; k++)
      Root[i][0][k] = 0;
  }
  int maxX[16];
  int maxY[16];

  int maxS;
  int p = 9;
  for(int i = 0; i < seqContainer.size(); i++){
    if(freePosition != -1){
      analysedSequences[freePosition] = i;
      len[freePosition] = seqContainer[i].size();
      analysedSeqIndex[freePosition] = 0;

      freePosition = -1;

      for(int t = 0; t < 16; t++){

        if(len[t] > analysedSeqIndex[t]){
          residue[t] = seqContainer[analysedSequences[t]][analysedSeqIndex[t]];
          done = false;
        }
        else{
          residue[t] = -1;
          if(freePosition == -1 && i != seqContainer.size() - 1)
            freePosition = t;
        }

      }

    }
    while(!done && (freePosition == -1 || i == seqContainer.size()-1)){
      done = true;
      for(int l = 0; l < len1; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < len1+1; i++){
                leftScore[i][j] = 0;
                if(i != len1){
                  maxPosLine[i][j] = 0;
                  maxLineVal[i][j] = 0;
                }
              }
              opGap[j] = 0;
              exGap[j] = 1;
              maxScore[j] = 0;
              maxX[j] = 0;
              maxY[j] = 0;
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

        for(int w = 0; w < 2; w++){

          __m128i Z = _mm_load_si128((__m128i*) &zero[w*8]);
          __m128i RESC = _mm_load_si128((__m128i*) &colScore[w*8]);
          __m128i RESL = _mm_load_si128((__m128i*) &lineScore[w*8]);
          __m128i S = _mm_load_si128((__m128i*) &score[w*8]);

          __m128i OP = _mm_load_si128((__m128i*) &opGap[w*8]);
          __m128i COL = _mm_load_si128((__m128i*) &maxColVal[w*8]);
          __m128i PMCOL = _mm_load_si128((__m128i*) &maxPosCol[w*8]);
          __m128i PCOL = _mm_load_si128((__m128i*) &ColPos[w*8]);

          __m128i EX = _mm_load_si128((__m128i*) &exGap[w*8]);
          __m128i LIN = _mm_load_si128((__m128i*) &maxLineVal[l][w*8]);
          __m128i PMLIN = _mm_load_si128((__m128i*) &maxPosLine[l][w*8]);
          __m128i PLIN = _mm_load_si128((__m128i*) &LinePos[w*8]);

          __m128i D = _mm_load_si128((__m128i*) &leftScore[l][w*8]);
          __m128i P = _mm_load_si128((__m128i*) &diagGap[w*8]);

          __m128i M = _mm_load_si128((__m128i*) &maxScore[w*8]);
          __m128i UNIT = _mm_load_si128((__m128i*) &unit[w*8]);



          __m128i SUBC2 = _mm_sub_epi16(PCOL,PMCOL);
          __m128i MULTC = _mm_mullo_epi16(EX,SUBC2);
          RESC = _mm_sub_epi16(COL,MULTC);

          __m128i SUBL2 = _mm_sub_epi16(PLIN,PMLIN);
          __m128i MULTL = _mm_mullo_epi16(EX,SUBL2);
          RESL = _mm_sub_epi16(LIN,MULTL);

          S = _mm_add_epi16(D,P);

          S = _mm_max_epi16(S,RESC);
          S = _mm_max_epi16(S,RESL);
          S = _mm_max_epi16(S,Z);

          __m128i SUB1 = _mm_sub_epi16(S, OP);

          __m128i CMPLTC = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, RESC));
          __m128i MULTC2 = _mm_mullo_epi16(PCOL,CMPLTC);
          __m128i NEWMPC = _mm_max_epi16(PMCOL, MULTC2);
          __m128i NEWMC = _mm_max_epi16(_mm_mullo_epi16(_mm_sub_epi16(UNIT,CMPLTC), COL),_mm_mullo_epi16(CMPLTC, SUB1));

          __m128i CMPLTL = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, RESL));
          __m128i MULTL2 = _mm_mullo_epi16(PLIN,CMPLTL);
          __m128i NEWMPL = _mm_max_epi16(PMLIN, MULTL2);
          __m128i NEWML = _mm_max_epi16(_mm_mullo_epi16(_mm_sub_epi16(UNIT,CMPLTL), LIN),_mm_mullo_epi16(CMPLTL, SUB1));

          _mm_store_si128((__m128i*)&score[w*8], S);
          _mm_store_si128((__m128i*)&leftScore[l][w*8], S);
          _mm_store_si128((__m128i*)&maxPosCol[w*8], NEWMPC);
          _mm_store_si128((__m128i*)&maxColVal[w*8], NEWMC);
          _mm_store_si128((__m128i*)&maxPosLine[l][w*8], NEWMPL);
          _mm_store_si128((__m128i*)&maxLineVal[l][w*8], NEWML);
          _mm_store_si128((__m128i*)&maxScore[w*8], M = _mm_max_epi16(S,M));
          _mm_store_si128((__m128i*)&colScore[w*8], RESC);
          _mm_store_si128((__m128i*)&lineScore[w*8], RESL);
        }
        for(int x = 0; x < 16; x++){
          if(score[x] == colScore[x])
            match[x] = 1;
          else if (score[x] == lineScore[x])
            match[x] = 2;
          if(score[x] == 0)
            match[x] = 0;
          if(score[x] == maxScore[x]){
            maxX[x] = analysedSeqIndex[x];
            maxY[x] = l;
          }
          Root[x][l+1][analysedSeqIndex[x]+1] = match[x];
        }



        if(l == len1-1){
          for(int f = len1-1; f >= 0; f--){
            if(f == 0)
              firstLine = true;
            offset(leftScore[f+1], leftScore[f], zero, firstLine);
            firstLine = false;
          }
        }
      }

      for(int k = 0 ; k < 16; k++){
        analysedSeqIndex[k]++;
        if(len[k] > analysedSeqIndex[k]){
          residue[k] = seqContainer[analysedSequences[k]][analysedSeqIndex[k]];
          done = false;
        }
        else{
          residue[k] = -1;
          if(freePosition == -1 && i != seqContainer.size() - 1){
            freePosition = k;
          }
        }
      }

    }
	}


  cout << "matrice score: "<< maxScore[p] << " pos: " << maxX[p] << " " << maxY[p] << endl;


  cout << "Root : \n";
  for(int i = 0; i < len1+1; i++){
    for(int j=0; j < len[p]+1; j++){
      if(i == maxY[p]+1 && j == maxX[p]+1)
       cout << "#";
      cout << Root[p][i][j] << " ";
    }
    cout << endl;

  }
  cout << "done all" << endl;

}
