#include "ScoringMatrix.h"

ScoringMatrix::ScoringMatrix(){

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
};

void ScoringMatrix::setRootTarget(Position* rootPos, Position* targetPos){
  rootPos->addTarget(targetPos);
  targetPos->addRoot(rootPos);
};
