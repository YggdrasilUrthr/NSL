main.out: csvwriter.o main.o
	g++ main.o csvwriter.o -o main.out

main.o: csvwriter.h
	g++ -c main.cpp -o main.o

csvwriter.o: csvwriter.h
	g++ -c csvwriter.cpp -o csvwriter.o

clean:
	rm -rf *.o