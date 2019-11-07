all: smith;\

smith: smith.o;\
g++ -o smith smith.o
smith.o : smith.cpp BlosumMatrix.o ScoringMatrix.o;\
g++ -c smith.cpp
BlosumMatrix.o : BlosumMatrix.cpp Position.o;\
g++ -c BlosumMatrix.cpp
ScoringMatrix.o : ScoringMatrix.cpp Position.o;\
g++ -c ScoringMatrix.cpp
Position.o : Position.cpp;\
g++ -c Position.cpp
#test: testBlosum.o;\
g++ -o testBlosum testBlosum.o
#testBlosum.o: testBlosum.cpp BlosumMatrix.o;\
g++ -c testBlosum.cpp
clean: ;\
rm *.o
