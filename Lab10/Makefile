LIBS := ../Libs/
LAB9 := ../Lab9/

GPP := mpicxx

Ex10_1.out: Ex10_1.o csvwriter.o random.o Functions.o Chromosome.o Population.o
	$(GPP) -g Ex10_1.o csvwriter.o random.o Chromosome.o Population.o Functions.o -o Ex10_1.out

Ex10_1.o: $(LIBS)/CsvWriter/csvwriter.h Ex10_1.cpp
	$(GPP) -g -c Ex10_1.cpp -o Ex10_1.o

Population.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h $(LAB9)/Functions.h $(LAB9)/Population.h $(LAB9)/Population.cpp
	$(GPP) -g -c $(LAB9)/Population.cpp -o Population.o

Chromosome.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h $(LAB9)/Chromosome.h $(LAB9)/Chromosome.cpp
	$(GPP) -g -c $(LAB9)/Chromosome.cpp -o Chromosome.o

Functions.o: $(LAB9)/Functions.cpp $(LAB9)/Functions.h
	$(GPP) -g -c $(LAB9)/Functions.cpp -o Functions.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	$(GPP) -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	$(GPP) -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex10_1.out

clean:
	rm -rf *.o