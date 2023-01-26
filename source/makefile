TARGET = MLS

CC = g++
OBJE = main.o errorfunc.o log.o 
#CFLAGS  = -O3 -Wall
CFLAGS  = -O3 -Wall -Wno-uninitialized
#CFLAGS  = -g -Wall
LDFLAGS =  -lm


$(TARGET): $(OBJE)
	$(CC) $(CFLAGS) -o $@ $(OBJE) $(LDFLAGS) 

.cpp.o:
	$(CC) $(CFLAGS) -c $< 

clean:
	/bin/rm -f *.o *~ $(TARGET)  $(TARGET).exe 


main.o			:log.h
errorfunc.o		:errorfunc.h
log.o			:errorfunc.h log.h
