
OPTION = -I../boole/src -L../boole/src


ifeq ($(DEBUG),1)
	CFLAGS = -Wall -g  $(OPTION) 
else

	CFLAGS = -O3  -march=native   -Wno-unused-result $(OPTION)

endif


all : basis.exe repres.exe nlqfast.exe 

debug:
	$(MAKE)  -B DEBUG=1


basis.exe : basis.c
	gcc $(CFLAGS) -Wall  $^  -o $@  -lgmp  -lboole -lgmp

nlqfast.exe : nlqfast.c
	gcc $(CFLAGS) -fopenmp -Wall  $^  -o $@  -lboole  -lgmp 

repres.exe : repres.c
	gcc $(CFLAGS) -fopenmp -Wall  $^  -o $@  -lboole  -lgmp 

infos      :
	./repres.exe  > info-3-5-7.txt
web        :
	cp B357.html  ~/web-docs/project/B357/B357.html 
	cp NQ-38-40-3-5-7.txt  ~/web-docs/project/B357/. 
	cp infos-3-5-7.txt  ~/web-docs/project/B357/. 

clean :
	rm -f *.exe
