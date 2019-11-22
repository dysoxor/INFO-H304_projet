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

//#include "ScoringMatrix.h"
//#include "BlosumMatrix.h"
#include "binary.h"

using namespace std;

int matching(vector<int> seq1, string seq2, int len1);
int findMax(int tableau[], int size);
int findMax(vector<int> tableau, int size);

void setupBlosumMatrix(string pathToBlosumMatrix);
void dbAlignment(string db, string query, PIN* filePIN, PSQ* filePSQ);
void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right);
void insertion_sortmerge(vector<int> & scorev, vector<int> &indexv,int left, int right);
void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right);
