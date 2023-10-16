LIBS := ../Libs/

Ex9.out: Ex9.o csvwriter.o random.o Functions.o Chromosome.o Population.o
	g++ -g Ex9.o csvwriter.o random.o Chromosome.o Population.o Functions.o -o Ex9.out

Ex9.o: $(LIBS)/CsvWriter/csvwriter.h Ex9.cpp
	g++ -g -c Ex9.cpp -o Ex9.o

Population.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Functions.h Population.h Population.cpp
	g++ -g -c Population.cpp -o Population.o

Chromosome.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Chromosome.h Chromosome.cpp
	g++ -g -c Chromosome.cpp -o Chromosome.o

Functions.o: Functions.cpp Functions.h
	g++ -g -c Functions.cpp -o Functions.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex9.out

clean:
	rm -rf *.o