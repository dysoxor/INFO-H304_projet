CPPFLAG=g++
OPTIFLAG=-O3 -Ofast
C11FLAG=-std=c++11
PTHREADFLAG=-pthread
EXEC=find

all: $(EXEC)

find: find.o binary.o smith.o iomanagement.o
	$(CPPFLAG) $(OPTIFLAG) $(C11FLAG) $(PTHREADFLAG) -o $@ $^
%.o: %.cpp
	$(CPPFLAG) $(OPTIFLAG) $(C11FLAG) -o $@ -c $^
clean:
	rm -rf *.o
mrproper: clean
	rm -rf find
