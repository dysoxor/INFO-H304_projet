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

using namespace std;

class PIN{
private:
  uint32_t version;
  uint32_t database_type;
  uint32_t title_length;
  char* title;
  uint32_t timestamp_length;
  char* timestamp;
  uint32_t num_of_seq;
  uint64_t residue_count;
  uint32_t max_seq;
  uint32_t* hr_offset_table;
  uint32_t* seq_offset_table;
  uint32_t* ambiguity_offset_table;
public:
  int read(string dataFileName);
  int getNumSeq()const;
  int getHrOffset(int i)const;
  int getSqOffset(int i)const;
};

class PSQ{
private:
  int sequence;
public:
  int read(PIN* filePIN, string query, string dataFileName);
  vector<int> queryToInt(string query);
};

class PHR{
private:
  int binary;
  string hexadecimal;
  string string_length_bits;
  bool significantBitOn;
  unsigned long int string_length;
  string seqTitle;
public:
  int read(PIN* filePIN, int index, string dataFileName);
  string int_to_hex( int i );
  string hex_to_string(const string& in);
  string byteToBits(unsigned int u);
  unsigned long toInt(std::string const &s);
};