#include "binary.h"
#include <thread>
#include <mutex>
#include <immintrin.h>
#include <emmintrin.h>

using namespace std;

void offset(int16_t* leftScore, int16_t* diagonalScore, int16_t* zero, bool firstLine);
int matching(int seq1[], int seq2[], int len1, int len2);
int findMax(int tableau[], int size);
int setupBlosumMatrix(string pathToBlosumMatrix);
vector<vector<int>> dbAlignment(string db, string query, PSQ* filePSQ, string smMatrix, int gpo, int gpe, int nbResults);
void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right);
void insertion_sortmerge(vector<int> & scorev, vector<int> &indexv,int left, int right);
void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right);
