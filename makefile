mam:repeatfinder.o seqtools.o slider.o
	gcc repeatfinder.o seqtools.o slider.o -o mam -g
repeatfinder.o: repeatfinder.c main.h
	gcc -Wall -c repeatfinder.c -g 
seqtools.o: seqtools.c main.h
	gcc -Wall -c seqtools.c -g
slider.o: slider.c
	gcc -Wall -c slider.c -g
clean:
	rm -f *~ *.o \#* .\#* mam
install:
	cp mam /usr/bin
	
