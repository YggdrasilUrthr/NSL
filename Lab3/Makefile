LIBS := ../Libs/

Ex3_1_1.out: Ex3_1_1.o csvwriter.o random.o 
	g++ -g Ex3_1_1.o csvwriter.o random.o -o Ex3_1_1.out

Ex3_1_1.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex3_1_1.cpp
	g++ -g -c Ex3_1_1.cpp -o Ex3_1_1.o

Ex3_1_2.out: Ex3_1_2.o csvwriter.o random.o 
	g++ -g Ex3_1_2.o csvwriter.o random.o -o Ex3_1_2.out

Ex3_1_2.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex3_1_2.cpp
	g++ -g -c Ex3_1_2.cpp -o Ex3_1_2.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex3_1_1.out Ex3_1_2.out

clean:
	rm -rf *.o