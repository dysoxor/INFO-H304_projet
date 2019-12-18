#include "smith.h"

pair<string, string> readFasta(string file);
string scoreString(int index, int score, string db, PIN* filePIN, PHR* filePHR, int maxLine);
string alignementString(vector<int> result ,string query,string db, PIN* filePIN, PHR* filePHR, char dataBase[] , int maxLine);
void writeOutput(vector<vector<int>> results, string outputFile, string queryFileName, string dataBaseFileName, string queryName, string querySequence, clock_t begin, PIN* filePIN, PSQ* filePSQ);
