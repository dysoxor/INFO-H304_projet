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

#include "ScoringMatrix.h"
#include "BlosumMatrix.h"
#include "binary.h"

using namespace std;

int matching(string prot1, string prot2);
int findMax(int tableau[], int size);
void setupBlosumMatrix(string pathToBlosumMatrix);
void dbAlignment(string db, string query, PIN* filePIN, PSQ* filePSQ);
