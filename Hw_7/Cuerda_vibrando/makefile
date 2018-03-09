all: string_10.00.pdf string_10.00.dat

string_10.00.pdf:string_10.00.dat
	 python plot_string.py string_10.00.dat 

string_10.00.dat: string.x
	./string.x 10

string.x: string.c
	cc string.c -o string.x -lm

clean:
	rm -f string_*
	rm -f string.x
