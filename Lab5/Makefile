LIBS := ../Libs/

Ex5_1.out: Ex5_1.o csvwriter.o random.o 
	g++ -g Ex5_1.o csvwriter.o random.o -o Ex5_1.out

Ex5_1.o: $(LIBS)/CsvWriter/csvwriter.h $(LIBS)/RanGen/random.h Ex5_1.cpp
	g++ -g -c Ex5_1.cpp -o Ex5_1.o

csvwriter.o: $(LIBS)/CsvWriter/csvwriter.h
	g++ -g -c $(LIBS)/CsvWriter/csvwriter.cpp -o csvwriter.o

random.o: $(LIBS)/RanGen/random.h
	g++ -g -c $(LIBS)/RanGen/random.cpp -o random.o

all: Ex5_1.out Ex2_2.out

clean:
	rm -rf *.o