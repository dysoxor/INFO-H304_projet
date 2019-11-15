#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <tuple>

using namespace std;

class Protein{
private:
  string name;
  string sequence;
  Protein* next;
public:
  Protein(string n, string s);
  void setNext(Protein* s);
  Protein* getNext() const;
  string getName() const;
  string getSeqence() const;
};

class List{
private:
  int numOfProtein;
  Protein* head;
public:
  List();
  ~List();

  void insert(string n, string s);
  void del();
  void printList();
  int getNumOfProtein() const;
  Protein* getHead();
};

List* readFasta(string file);
