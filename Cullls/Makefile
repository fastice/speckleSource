CC =		gcc
CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= cullls.o loadLSCullData.o cullLSData.o cullLSStats.o cullLSSmooth.o cullLSIslands.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)



	mv *.o $(MACHTYPE)-$(OSTYPE)
