all: asp-seq asp-par my-asp-seq my-asp-par

asp-seq: asp-seq.c
	gcc -Wall -O3 -o asp-seq asp-seq.c

my-asp-seq: my-seq.c
	gcc -Wall -O3 -o my-asp-seq my-seq.c

asp-par: asp-par.c
	mpicc -Wall -O3 -o asp-par asp-par.c

my-asp-par: my-par.c
	mpicc -Wall -O3 -o my-asp-par my-par.c	

clean:
	rm -f *.o *core *~ asp-seq
