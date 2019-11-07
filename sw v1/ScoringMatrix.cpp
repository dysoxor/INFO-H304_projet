#include "ScoringMatrix"

ScoringMatrix::ScoringMatrix(){

}

const Position* ScoringMatrix::getPosition(int x, int y){
  return tableau[pair<int,int>(x,y)]
}

void ScoringMatrix::setValue(int value, int x, int y){
  tableau[pair<int,int>(x,y)]->setValue(value);
}

void ScoringMatrix::getValue(int x, int y){
  return tableau[pair<int, int>(x,y)]->getValue();
}


void addPosition(Position* pos){
  tableau.insert(pair<int,int>(pos->getX(),pos->getY()),pos);
}

void setRootTarget(Position* rootPos, Position* targetPos){
  rootPos.addTarget(targetPos);
  targetPos.addRoot(rootPos);
}
