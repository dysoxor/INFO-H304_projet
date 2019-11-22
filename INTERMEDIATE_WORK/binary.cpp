#include "binary.h"

//-------------------------- Reading .pin file section -------------------------

int PIN::read(string dataFileName){
  ifstream filePIN;
  filePIN.open(dataFileName+".pin", ios::binary | ios::in);
  //check the size of the file to know if the file is empty or does not exist
  filePIN.seekg(0, ios::end);
  int fileSize=filePIN.tellg();
  filePIN.seekg(0, ios::beg);
  if(fileSize == -1){
    return EXIT_FAILURE;
  }
  //it reads .pin file on uint32_t bytes which are stocked at the adress of
  // 'version'. Then,it has to convert big endian to the little endian.
  // Otherwise, computer is not able to understand the sequence of bytes.
  filePIN.read( (char *)&version, sizeof(uint32_t));
  version = __bswap_32(version);

  filePIN.read( (char*)&database_type, sizeof(uint32_t));
  database_type = __bswap_32(database_type);

  filePIN.read( (char*)&title_length, sizeof(uint32_t));
  title_length = __bswap_32(title_length);

  title = new char[title_length];
  filePIN.read((char*)title, sizeof(char)*title_length);

  filePIN.read( (char*)&timestamp_length, sizeof(uint32_t));
  timestamp_length = __bswap_32(timestamp_length);

  timestamp = new char[timestamp_length];
  filePIN.read((char*)timestamp, sizeof(char)*timestamp_length);

  filePIN.read( (char*)&num_of_seq, sizeof(uint32_t));
  num_of_seq = __bswap_32(num_of_seq);

  // this one is already in little endian
  filePIN.read( (char*)&residue_count, sizeof(uint64_t));

  filePIN.read( (char*)&max_seq, sizeof(uint32_t));
  max_seq = __bswap_32(max_seq);

  hr_offset_table = new uint32_t[num_of_seq+1];
  seq_offset_table = new uint32_t[num_of_seq+1];
  int sizeTable = sizeof(uint32_t)*(num_of_seq+1);
  filePIN.read((char*)hr_offset_table, sizeTable);
  filePIN.read((char*)seq_offset_table, sizeTable);
  for (int i=0; i < num_of_seq+1; i++){
    hr_offset_table[i] = __bswap_32(hr_offset_table[i]);
    seq_offset_table[i] = __bswap_32(seq_offset_table[i]);
  }

  filePIN.close();
  return EXIT_SUCCESS;
}

int PIN::getNumSeq()const{
  return num_of_seq;
}

int PIN::getHrOffset(int i)const{
  return hr_offset_table[i];
}

int PIN::getSqOffset(int i)const{
  return seq_offset_table[i];
}

//-------------------------- Reading .psq files section ------------------------

//convert query's sequence into a list of integer
vector<int> PSQ::queryToInt(string query){
  vector<int> res;
  map<char,int> intToChar {
    {'-',0},{'A',1},{'B',2},{'C',3},{'D',4},
    {'E',5},{'F',6},{'G',7},{'H',8},{'I',9},
    {'K',10},{'L',11},{'M',12},{'N',13},{'P',14},
    {'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},
    {'W',20},{'X',21},{'Y',22},{'Z',23},{'U',24},
    {'*',25},{'O',26},{'J',27}
  };

  for (int i = 0; i< query.size(); i++){
    int indice = intToChar[query.at(i)];
    res.push_back(indice);
  }

  return res;
}

int PSQ::read(PIN* filePIN, string query, string dataFileName){
  //convert query's sequence into a vector of integers to compare it faster with
  // sequences from data file
  vector<int> table_query = queryToInt(query);

  ifstream filePSQ;
  filePSQ.open(dataFileName+".psq", ios::binary | ios::in);
  //checking file if it is empty or unopennable
  filePSQ.seekg(0, ios::end);
  int fileSize=filePSQ.tellg();
  filePSQ.seekg(0, ios::beg);
  if(fileSize == -1){
    return EXIT_FAILURE;
  }

  //boolean which knows if the sequence was found in the .psq
  bool finded = false;
  int sizeOfSq;
  //need the index to find the offset of sequence in the header which is the same
  // as this
  int index = -1;
  for (int i = 1; i < filePIN->getNumSeq(); i++){
    if (finded){
      break;
    }
    //difference between the current offset and the following one, minus 1,
    //to get the size of the sequence (last bit is separator)
    sizeOfSq = filePIN->getSqOffset(i+1) - filePIN->getSqOffset(i) - 1;
    //the sequence from the data file has the same size as the query
    if (sizeOfSq == table_query.size() ){
      filePSQ.seekg(filePIN->getSqOffset(i));
      for (int j = 0; j <= sizeOfSq; j++){
        filePSQ.read((char*)&sequence, sizeof(uint8_t));
        //while the seperator '0' is not read it compares each letter of the
        //sequence
        if (sequence != 0){
          if (table_query[j] != sequence){
              break; //this is not a perfect match
          }
        }
        else {
            finded = true;//if the separator is reach without break, it means that
            //this is a perfect match
            index = i;
            break;
        }
      }
    }
  }
  filePSQ.close();
  return index;
}

//-------------------------- Reading .phr files section ------------------------

//conversion of an integer to hexadecimal
string PHR::int_to_hex( int i ){
  stringstream stream;
  stream << std::setw(2) << hex << i;//it allocates a space of 2 character for the conversion of the integer to hex
  return stream.str();
}

//convert byte to bits
string PHR::byteToBits(unsigned int u){
  string bits;
  for(int t=128; t>0; t = t/2) {//in binary it means that t begin from 1000 0000 and end at 0000 0001 it is a walking '1'
    if(u & t) bits+='1';//if 'u' bitwise AND 't' is true it means that a 1 is corresponding between t and u and it is added to the string bits
    else bits+='0';
  }
  return bits;
}

//convert hexadecimal to string
string PHR::hex_to_string(const string& in) {
    string output;
    if ((in.length() % 2) != 0) {//an hexadecimal should be pair
        throw runtime_error("String is not valid length ...");
    }
    size_t cnt = in.length() / 2;
    for (size_t i = 0; cnt > i; ++i) {
        uint32_t s = 0;
        stringstream ss;

        ss << hex << in.substr(i * 2, 2);//convert input string to hexadecimal format stocked into ss
        ss >> s;//convert hexadecimal to uint32_t
        output.push_back(static_cast<unsigned char>(s));//add it to output string with conversion of uint32_t to char
    }
    return output;
}

//convert bits to integer
unsigned long PHR::toInt(std::string const &s) {
    static const std::size_t MaxSize = CHAR_BIT*sizeof(unsigned long);
    if (s.size() > MaxSize) return 0; // handle error or just truncate?

    std::bitset<MaxSize> bits;
    std::istringstream is(s);
    is >> bits;
    return bits.to_ulong();
}

int PHR::read(PIN* filePIN, int index, string dataFileName){
  ifstream filePHR;
  filePHR.open(dataFileName+".phr", ios::binary | ios::in);
  if(!filePHR){
    return EXIT_FAILURE;
  }
  //check if the file is correct
  filePHR.seekg(0, ios::end);
  int fileSize = filePHR.tellg();
  filePHR.seekg(0, ios::beg);

  int seqOffset = filePIN->getHrOffset(index);//position in .psq file of the found sequence
  int size = filePIN->getHrOffset(index+1)-seqOffset;//size of the sequence's header
  bool visible_string = false;
  int byteForSize=0;

  filePHR.seekg(seqOffset);

  for(int i=0; i<size ; i++){
    filePHR.read( (char*)&binary, 1);//read byte by byte the file
    hexadecimal = int_to_hex(binary);//convert integer into the hexadecimal
    if(!visible_string && hexadecimal == "1a"){//'1a' says that the following is visible_string
      visible_string = true;
    }
    else if(visible_string){//if it is reading visible string of the file
      if(string_length_bits == ""){//the byte which is next to '1a' is the length of the visible_string
        string_length_bits = byteToBits(binary);//convert the byte responsible to give the size of the string into bits
        significantBitOn = false;//see if the first bit is off or on
        if(string_length_bits[0] == '1'){
          significantBitOn = true;
        }
        //get the string size by getting a look over the sequence of bits
        if(!significantBitOn){
          string_length = binary;
        }
        else{
          byteForSize = toInt(string_length_bits.substr(1,7));
          filePHR.read((char *)&string_length, byteForSize);
          string_length = __bswap_32(string_length);
        }
      }
      else if(string_length!=-1){//add character which is read while the string size is not reduced to -1
        seqTitle+=hex_to_string(hexadecimal);
        string_length--;
      }
    }
  }
  cout << "The title is : " << seqTitle << endl;
  filePHR.close();
  return EXIT_SUCCESS;
}
