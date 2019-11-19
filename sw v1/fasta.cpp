#include "fasta.h"

// The purpose of this is to read fasta file (database or query file) and it stocks
// the name and content of the sequences into a chained list
List* readFasta(string file){
  ifstream input(file);
  string line, name, content;

  //creating the chained list
  List* listProtein = new List();

  //If it is not able to open the file it returns error and in the main it exits
  // on failure
  if(!input.good()){
    return listProtein;
  }else{
    //while it not in end of file it read the line
    while( getline( input, line ).good() ){
        //in fasta the symbol '>' precedes the name and nexts lines are the sequences
        // of it
        if( line[0] == '>' ){
          if(!name.empty()){
            // it insert the last protein that it read into the list before starting
            // to read the next one
            listProtein->insert(name, content);
            content = "";
          }
          name = line.substr(1);
        } else if( !name.empty() ){
          content+=line;
        }
    }
    // it inserts the last protein of fasta file that it read because it is not
    // covered by the loop above
    listProtein->insert(name, content);
    // the following commented line prints all the fasta file from list that it
    // just built
    //listProtein->printList();
    return listProtein;
  }
}
// Constructor of an item contained in the List
Protein::Protein(string n, string s){
  name = n;
  sequence = s;
  next = NULL;
}
// the list has an order so each protein keep in memory the pointer of the following
// protein
void Protein::setNext(Protein* s){
  next = s;
}

Protein* Protein::getNext() const{
  return next;
}
// name of the protein
string Protein::getName() const{
  return name;
}
// sequence of the protein
string Protein::getSequence() const{
  return sequence;
}
// initialization of the chained list
List::List(){
  numOfProtein = 0;
  head = NULL;
}
// destructor of the chained list
List::~List(){
  while (head != NULL)
    del();
}
// insert a new protein on the first place of the list
void List::insert(string n, string s){
  Protein* p = new Protein(n, s);
  p->setNext(head);
  head=p;
  numOfProtein++;
}
// free the memory used by the first item of the list
void List::del(){
  Protein* t = head;
  // set the following protein as the first item of the list
  head = head->getNext();
  delete t;
  numOfProtein--;
}
// print the list
void List::printList(){
  Protein* t = head;
  for(int i = 0; i < numOfProtein; i++){
    cout << i << "name:  " << t->getName() << endl;
    cout << "sequence:  " << t->getSequence() << endl;
    t = t->getNext();
  }
}
// quantity of protein
int List::getNumOfProtein() const{
  return numOfProtein;
}
// first protein in the list
Protein* List::getHead(){
  return head;
}
