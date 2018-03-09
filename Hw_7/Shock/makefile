all: estado_0.01.pdf

estado_0.01.pdf: graficas.py estado_0.01.dat
	python graficas.py estado_0.01.dat

estado_0.01.dat: sod_test.x
	./sod_test.x 0.01

sod_test.x: sod_test.c
	cc sod_test.c -o sod_test.x -lm

clean: 
	rm -f *.dat *.pdf *.x
