#include "ScoringMatrix.h"
using namespace std;

ScoringMatrix::ScoringMatrix(){
  sizeX=0;
  sizeY=0;
};

Position* ScoringMatrix::getPosition(int x, int y){
  //Return the position pointer of (x,y)
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
    //If we add a new position with a x greater than our actual maximum x,
    //redefine the sizeX
    sizeX = pos->getX() +1;
  }
  if (pos->getY()+1 > sizeY){
    sizeY = pos->getY() +1;
  }
};

/*void ScoringMatrix::setRootTarget(Position* rootPos, Position* targetPos){
  rootPos->addTarget(targetPos);
  targetPos->addRoot(rootPos);
};*/

const void ScoringMatrix::print(){
  for (int i = 0; i < sizeY; i++){
		for (int j = 0; j < sizeX; j++){
			cout << getValue(j, i) << " ";
		}
		cout << endl;
	}
}

void ScoringMatrix::setMaxLine(int line, Position* pos){
  if (pos->getValue() > maxLines[line]->getValue()){
    maxLines[line] = pos;
  }
}

Position* ScoringMatrix::getMaxLine(int line){
  return maxLines[line];
}

void ScoringMatrix::setMaxColumn(int column, Position* pos){
  if (pos->getValue() > maxColumns[column]->getValue()){
    maxColumns[column] = pos;
  }
}

Position* ScoringMatrix::getMaxColumn(int column){
  return maxColumns[column];
}

void ScoringMatrix::setupMax(int len1, int len2){
  maxLines.assign(len2,new Position());
  maxColumns.assign(len1,new Position());
}

int ScoringMatrix::getDistYWithMax(Position* pos){
  return pos->getY()-maxColumns[pos->getX()]->getY();
}

int ScoringMatrix::getDistXWithMax(Position* pos){
  return pos->getX() - maxLines[pos->getY()]->getX();
}