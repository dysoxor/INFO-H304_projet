#include "smith.h"
#include <chrono>
#include <ctime>



pair<string, string> readFasta(string file);
string scoreString(int index, int score, string db, PIN* filePIN,
  PHR* filePHR, int maxLine);
string alignementString(vector<int> result ,string query,string db,
  PIN* filePIN, PHR* filePHR, char dataBase[] , int maxLine);
void writeOutput(vector<vector<int>> results, string outputFile,
  string queryFileName, string dataBaseFileName,
  string queryName, string querySequence,
  chrono::time_point<chrono::system_clock> begin,
  PIN* filePIN, PSQ* filePSQ,PHR* filePHR, string smFile,
  int gap_ex, int gap_op);
