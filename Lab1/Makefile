LIBS := ../Libs/

Ex1_1.out: Ex1_1.o csvwriter.o random.o 
	g++ -g Ex1_1.o csvwriter.o random.o -o Ex1_1.out

Ex1_1.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex1_1.cpp
	g++ -g -c Ex1_1.cpp -o Ex1_1.o

Ex1_2.out: Ex1_2.o csvwriter.o random.o
	g++ -g Ex1_2.o csvwriter.o random.o -o Ex1_2.out

Ex1_2.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex1_2.cpp
	g++ -g -c Ex1_2.cpp -o Ex1_2.o

Ex1_3.out: Ex1_3.o csvwriter.o random.o
	g++ -g Ex1_3.o csvwriter.o random.o -o Ex1_3.out

Ex1_3.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex1_3.cpp
	g++ -g -c Ex1_3.cpp -o Ex1_3.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex1_1.out Ex1_2.out Ex1_3.out

clean:
	rm -rf *.o