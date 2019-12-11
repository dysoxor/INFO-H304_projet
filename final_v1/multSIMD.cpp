#include <emmintrin.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <array>
#include <stdio.h>
#include <immintrin.h>

using namespace std;

int main (int argc, char** argv){
  __attribute__((aligned (16))) int16_t score[16];
  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t opGap[16];
  __attribute__((aligned (16))) int16_t HScore [16];
  __attribute__((aligned (16))) int16_t diagonalScore[16];
  __attribute__((aligned (16))) int16_t exGap[16];
  __attribute__((aligned (16))) int16_t maxScore[16];
  __attribute__((aligned (16))) int16_t zero[16];

  __attribute__((aligned (16))) int16_t maxPosLine[16];
  __attribute__((aligned (16))) int16_t maxPosCol[16];
  __attribute__((aligned (16))) int16_t maxColVal[16];
  __attribute__((aligned (16))) int16_t maxLineVal[16];
  __attribute__((aligned (16))) int16_t PosCol[16];
  __attribute__((aligned (16))) int16_t PosLine[16];
  __attribute__((aligned (16))) int16_t colScore[16];
  __attribute__((aligned (16))) int16_t lineScore[16];

  __attribute__((aligned (16))) int16_t test[16];
  __attribute__((aligned (16))) int16_t unit[16];

  for(int i = 0; i < 16 ; i++){
    score[i] = 0;
    diagGap[i] = 3;
    opGap[i] = 0;
    HScore[i] = 0;
    diagonalScore[i] = 0;
    exGap[i] = 1;
    maxScore[i] = 0;
    zero[i] = 0;
    maxPosLine[i] = 0;
    maxPosCol[i] = 0;
    maxColVal[i] = 4;
    maxLineVal[i] = 0;
    PosCol[i] = 10;
    PosLine[i] = 0;
    colScore[i] = 0;
    lineScore[i] = 0;
    test[i] = 0;
    unit[i] = 1;
  }

  __m128i Z = _mm_load_si128((__m128i*) &zero);
  __m128i RESC = _mm_load_si128((__m128i*) &colScore);
  __m128i RESL = _mm_load_si128((__m128i*) &lineScore);
  __m128i S = _mm_load_si128((__m128i*) &score);

  __m128i OP = _mm_load_si128((__m128i*) &opGap);
  __m128i COL = _mm_load_si128((__m128i*) &maxColVal);
  __m128i PMCOL = _mm_load_si128((__m128i*) &maxPosCol);
  __m128i PCOL = _mm_load_si128((__m128i*) &PosCol);

  __m128i EX = _mm_load_si128((__m128i*) &exGap);
  __m128i LIN = _mm_load_si128((__m128i*) &maxLineVal);
  __m128i PMLIN = _mm_load_si128((__m128i*) &maxPosLine);
  __m128i PLIN = _mm_load_si128((__m128i*) &PosLine);

  __m128i D = _mm_load_si128((__m128i*) &diagonalScore);
  __m128i P = _mm_load_si128((__m128i*) &diagGap);

  __m128i M = _mm_load_si128((__m128i*) &maxScore);

  __m128i TEST = _mm_load_si128((__m128i*) &test);
  __m128i UNIT = _mm_load_si128((__m128i*) &unit);


  __m128i SUB1 = _mm_sub_epi16(COL, OP);
  __m128i SUB2 = _mm_sub_epi16(PCOL,PMCOL);
  __m128i MULT = _mm_mullo_epi16(EX,SUB2);
  RESC = _mm_sub_epi16(SUB1,MULT);
  RESL = _mm_sub_epi16(_mm_sub_epi16(LIN, OP),_mm_mulhi_epi16(EX,_mm_sub_epi16(PLIN,PMLIN)));
  S = _mm_add_epi16(D,P);

  S = _mm_max_epi16(S,RESC);
  S = _mm_max_epi16(S,RESL);
  S = _mm_max_epi16(S,Z);

  TEST = _mm_add_epi16(UNIT,_mm_cmplt_epi16(RESC, S));
  __m128i MULT2 = _mm_mullo_epi16(PCOL,TEST);

  _mm_store_si128((__m128i*)&test, MULT);
  _mm_store_si128((__m128i*)&score, S);
  _mm_store_si128((__m128i*)&diagonalScore, S);
  _mm_store_si128((__m128i*)&maxPosCol, MULT2);
  _mm_store_si128((__m128i*)&maxColVal, _mm_max_epi16(COL, S));
  _mm_store_si128((__m128i*)&maxPosLine, _mm_mulhi_epi16(PLIN,_mm_cmplt_epi16(LIN, S)));
  _mm_store_si128((__m128i*)&maxLineVal, _mm_max_epi16(LIN, S));
  _mm_store_si128((__m128i*)&maxScore, M = _mm_max_epi16(S,M));
  _mm_store_si128((__m128i*)&colScore, RESC);
  _mm_store_si128((__m128i*)&lineScore, RESL);

  cout << "score : " << score[0] << " maxY : " << maxColVal[0] << " maxPosY : "
  << maxPosCol[0] << " colScore : " << colScore[0] << " test : " << test[0] << endl;
}
