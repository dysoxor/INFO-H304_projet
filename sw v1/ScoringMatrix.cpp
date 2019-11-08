#include "ScoringMatrix.h"
using namespace std;

ScoringMatrix::ScoringMatrix(){
  sizeX=0;
  sizeY=0;
};

Position* ScoringMatrix::getPosition(int x, int y){
  return tableau[pair<int,int>(x,y)];
};

void ScoringMatrix::setValue(int value, int x, int y){
  tableau[pair<int,int>(x,y)]->setValue(value);
};

const int ScoringMatrix::getValue(int x, int y){
  return tableau[pair<int, int>(x,y)]->getValue();
};


void ScoringMatrix::addPosition(Position* pos){
  tableau.insert(pair<pair<int, int>, Position*> (pair<int,int>(pos->getX(),pos->getY()),pos));
  if (pos->getX()+1 > sizeX){
    sizeX = pos->getX() +1;
  }
  if (pos->getY()+1 > sizeY){
    sizeY = pos->getY() +1;
  }
};

void ScoringMatrix::setRootTarget(Position* rootPos, Position* targetPos){
  rootPos->addTarget(targetPos);
  targetPos->addRoot(rootPos);
};

const void ScoringMatrix::print(){
  for (int i = 0; i < sizeY; i++){
		for (int j = 0; j < sizeX; j++){
			cout << getValue(j, i) << " ";
		}
		cout << endl;
	}
}
