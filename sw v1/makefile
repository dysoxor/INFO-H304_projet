#all: smith;\

all :;\
g++ -o smith smith.cpp BlosumMatrix.cpp fasta.cpp
#smith: smith.o;\
g++ -o smith smith.cpp BlosumMatrix.cpp ScoringMatrix.cpp Position.cpp
#smith.o : smith.cpp BlosumMatrix.o ScoringMatrix.o;\
g++ -c smith.cpp BlosumMatrix.o ScoringMatrix.o
#BlosumMatrix.o : BlosumMatrix.cpp Position.o;\
g++ -c BlosumMatrix.cpp Position.o
#ScoringMatrix.o : ScoringMatrix.cpp Position.o;\
g++ -c ScoringMatrix.cpp Position.o
#Position.o : Position.cpp;\
g++ -c Position.cpp
#test: testBlosum.o;\
g++ -o testBlosum testBlosum.o
#testBlosum.o: testBlosum.cpp BlosumMatrix.o;\
g++ -c testBlosum.cpp
clean: ;\
rm *.o
