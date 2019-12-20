#include "binary.h"


//-------------------------- Reading .pin file section -------------------------

int PIN::read(string dataFileName){
  ifstream filePIN;
  filePIN.open(dataFileName+".pin", ios::binary | ios::in);
  //check the size of the file to know if the file is empty or does not exists
  filePIN.seekg(0, ios::end);
  int fileSize=filePIN.tellg();
  filePIN.seekg(0, ios::beg);
  if(fileSize == -1){
    return EXIT_FAILURE;
  }
  //it reads .pin file on uint32_t bytes which are stocked at the adress of
  // 'version' since there it has to convert big endian to the little endian
  // otherwize the computer is not able to understand the sequence of bytes.
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

void PIN::end(){
  delete hr_offset_table;
  delete seq_offset_table;
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

int PIN::getmaxSeq()const{
  return max_seq;
}
int PIN::getNumRes()const{
  return residue_count;
}

//-------------------------- Reading .psq files section ------------------------


int PSQ::charge(PIN* fPIN, string dataFileName){
  filePIN = fPIN;
  ifstream filePSQ;
  filePSQ.open(dataFileName+".psq", ios::binary | ios::in);
  //checking file if it is empty or unopennable
  filePSQ.seekg(0, ios::end);
  int fileSize=filePSQ.tellg();
  filePSQ.seekg(0, ios::beg);
  if(fileSize == -1){
    return EXIT_FAILURE;
  }
  database = new char[fileSize];
  filePSQ.read(database, fileSize);
  filePSQ.close();
  return EXIT_SUCCESS;
}

char* PSQ::getDatabase(){
  return database;
}

void PSQ::end(){
  delete database;
}

PIN* PSQ::getPIN(){
  return filePIN;
}

void PSQ::copy(char* out){
  for(int i = 0; i < filePIN->getNumRes(); i++){
    out[i] = database[i];
  }
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
  uint32_t s = 0;
  stringstream ss;

  ss << hex << in;//convert input string to hexadecimal format stocked into ss
  ss >> s;//convert hexadecimal to uint32_t
  output.push_back(static_cast<unsigned char>(s));//add it to output string with conversion of uint32_t to char

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

int PHR::charge(PIN* fPIN, string dataFileName){
  //Charge the file in memory
  filePIN = fPIN;
  ifstream filePHR;
  filePHR.open(dataFileName+".phr", ios::binary | ios::in);
  //checking file if it is empty or unopennable
  filePHR.seekg(0, ios::end);
  fileSize=filePHR.tellg();
  filePHR.seekg(0, ios::beg);
  if(fileSize == -1){
    return EXIT_FAILURE;
  }
  file = new char[fileSize];
  filePHR.read(file, fileSize);
  filePHR.close();
  return EXIT_SUCCESS;
}

string PHR::getTitle(int index){
  string seqTitle = "";
  int binary = 0;
  string hexadecimal = "";
  string string_length_bits = "";
  bool significantBitOn;
  unsigned long int string_length;

  int seqOffset = filePIN->getHrOffset(index);//position in .psq file of the found sequence
  int size = filePIN->getHrOffset(index+1)-seqOffset;//size of the sequence's header
  bool visible_string = false;
  int byteForSize=0;

  for(int i=0; i<size ; i++){
    binary = file[seqOffset+i];//filePHR.read( (char*)&binary, 1);//read byte by byte the file
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
          string_length = file[seqOffset+i+byteForSize];//filePHR.read((char *)&string_length, byteForSize);
          i+=byteForSize;
          string_length = __bswap_32(string_length);
        }
      }
      else if(string_length!=-1){//add character which is read while the string size is not reduced to -1
        seqTitle+=hex_to_string(hexadecimal);
        string_length--;
      }
    }
  }
  return seqTitle;
}

void PHR::end(){
  delete file;
}
