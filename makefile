
OPTION = -I../boole/src -L../boole/src


ifeq ($(DEBUG),1)
	CFLAGS = -Wall -g  $(OPTION) 
else

	CFLAGS = -O3  -march=native   -Wno-unused-result $(OPTION)

endif


all : prange.exe onthefly.exe zipcheck.exe check.exe read.exe write.exe basis.exe repres.exe nlqfast.exe 

debug:
	$(MAKE)  -B DEBUG=1

common.o  : common.c
	gcc $(CFLAGS)  $^  -c  


prange.exe : common.o prange.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -fopenmp -lboole -lgmp

onthefly.exe : common.o onthefly.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -fopenmp -lboole -lgmp


write.exe : common.o  write.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

read.exe : common.o  read.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

check.exe : common.o check.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

zipcheck.exe : common.o zipcheck.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

basis.exe : basis.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

nlqfast.exe : common.o nlqfast.c
	gcc $(CFLAGS) -fopenmp -Wall  $^  -o $@  -lboole  -lgmp 

repres.exe : repres.c
	gcc $(CFLAGS) -fopenmp -Wall  $^  -o $@  -lboole  -lgmp 

infos      :
	./repres.exe  > info-3-5-7.txt
web        :
	cp B357.html  ~/web-docs/project/B357/B357.html 
	cp NQ-38-40-3-5-7.txt  ~/web-docs/project/B357/. 
	cp info-3-5-7.txt  ~/web-docs/project/B357/. 
	mkdir -p ~/web-docs/project/B357/data 
	cp  data/*  ~/web-docs/project/B357/data/. 

clean :
	rm -f *.exe
