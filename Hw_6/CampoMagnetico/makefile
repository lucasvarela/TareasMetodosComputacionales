all: trayectoria_1000000_45.pdf  

trayectoria_1000000_45.pdf: plot.py trayectoria_1000000_45.dat
	python plot.py trayectoria_1000000_45.dat

trayectoria_1000000_45.dat: particle_in_field.x
	./particle_in_field.x 1000000 45

particle_in_field.x: particle_in_field.c
	cc particle_in_field.c -lm -o ./particle_in_field.x