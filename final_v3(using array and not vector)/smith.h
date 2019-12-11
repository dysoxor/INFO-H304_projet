#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <byteswap.h>
#include <iomanip>
#include <vector>
#include <bitset>
#include <sstream>
#include <climits>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include <math.h>
#include <vector>

//#include "ScoringMatrix.h"
//#include "BlosumMatrix.h"
#include "binary.h"

using namespace std;

int matching(int seq1[], int seq2[], int len1, int len2);
int matching_SIMD(vector<int> seq1, vector<int> residue, int len1);
int findMax(int tableau[], int size);
int findMax(vector<int> tableau, int size);
void dbAlignmentTest(string s1, string s2);
void setupBlosumMatrix(string pathToBlosumMatrix);
vector<vector<int>> dbAlignment(string db, string query, PSQ* filePSQ, string smMatrix, int gpo, int gpe, int nbResults);
void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right);
void insertion_sortmerge(vector<int> & scorev, vector<int> &indexv,int left, int right);
void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right);
