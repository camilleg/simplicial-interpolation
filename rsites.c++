#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

// Improved by Camille Goudeseune to use /dev/urandom as seed
// instead of the number of seconds since the epoch
// (which gives the same values if called within a 1-second interval).

double erand48 (unsigned short xsubi[3]);

int main(int argc, char *argv[])

{	int nsites,d,i;
	unsigned short X[3];
	if (argc<3) exit(3);
	sscanf(argv[1], "%d", &nsites);
	sscanf(argv[2], "%d", &d);
	int fd = open("/dev/urandom", O_RDONLY);
	if (sizeof(X) != read(fd, X, sizeof(X)))
	  fprintf(stderr, "rsites warning: failed to read from /dev/urandom.\n");
	close(fd);
	while(nsites>0){
		for (i=0;i<d;i++) printf("%6.0f ", floor(1e6*erand48(X)));
		printf("\n");
		nsites--;
	}
}
