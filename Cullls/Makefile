CC =		gcc
CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

OBJS= cullls.o loadLSCullData.o cullLSData.o cullLSStats.o cullLSSmooth.o cullLSIslands.o

.KEEP_STATE:
all:	$(OBJS)



	mv *.o $(MACHTYPE)-$(OSTYPE)
