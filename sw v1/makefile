all: smith;\

smith: smith.o BlosumMatrix.o;\
g++ -o smith smith.o BlosumMatrix.o
smith.o : smith.cpp BlosumMatrix.o;\
g++ -c smith.cpp
BlosumMatrix.o : BlosumMatrix.cpp;\
g++ -c BlosumMatrix.cpp
#test: testBlosum.o;\
g++ -o testBlosum testBlosum.o
#testBlosum.o: testBlosum.cpp BlosumMatrix.o;\
g++ -c testBlosum.cpp
clean: ;\
rm *.o
