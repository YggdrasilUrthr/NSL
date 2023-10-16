LIBS := ../Libs/

Ex8.out: Ex8.o csvwriter.o random.o 
	g++ -g Ex8.o csvwriter.o random.o -o Ex8.out

Ex8.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex8.cpp
	g++ -g -c Ex8.cpp -o Ex8.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex8.out

clean:
	rm -rf *.o