LIBS := ../Libs/

Ex2_1.out: Ex2_1.o csvwriter.o random.o 
	g++ -g Ex2_1.o csvwriter.o random.o -o Ex2_1.out

Ex2_1.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex2_1.cpp
	g++ -g -c Ex2_1.cpp -o Ex2_1.o

Ex2_2.out: Ex2_2.o csvwriter.o random.o
	g++ -g Ex2_2.o csvwriter.o random.o -o Ex2_2.out

Ex2_2.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex2_2.cpp
	g++ -g -c Ex2_2.cpp -o Ex2_2.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex2_1.out Ex2_2.out

clean:
	rm -rf *.o