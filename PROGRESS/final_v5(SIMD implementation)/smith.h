#include "binary.h"
#include <immintrin.h>
#include <emmintrin.h>
//#include <avxintrin.h>

using namespace std;

void offset(int16_t* leftScore, int16_t* diagonalScore, int16_t* zero, bool firstLine);
int matching1(int seq1[], int seq2[], int len1, int len2);
int matching2(int seq1[], int seq2[], int len1, int len2);
int matching_SIMD(vector<int> seq1, vector<int> residue, int len1);
int findMax(int tableau[], int size);
int findMax(vector<int> tableau, int size);
void dbAlignmentTest(string s1, string s2);
void setupBlosumMatrix(string pathToBlosumMatrix);
vector<vector<int>> dbAlignment(string db, string query, PSQ* filePSQ, string smMatrix, int gpo, int gpe, int nbResults);
void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right);
void insertion_sortmerge(vector<int> & scorev, vector<int> &indexv,int left, int right);
void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right);
