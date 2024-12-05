CC =		gcc
CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= cullst.o loadCullData.o cullSTData.o writeCullData.o cullStats.o \
cullSmooth.o cullIslands.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)
