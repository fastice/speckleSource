CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS =	$(FLAGS) -fopenmp -c -I$(INCLUDEPATH) -D_LARGEFILE64_SOURCE -D_FILE64 -D_LARGEFILE_SOURCE -DHPLARGEFILES -MMD -MP

SRCS =	strackw.c corrTrackFast.c mallocPerThreadArraysW.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
