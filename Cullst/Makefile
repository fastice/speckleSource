CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH) -MMD -MP -fopenmp

SRCS =	cullst.c loadCullData.c cullSTData.c writeCullData.c cullStats.c \
	cullSmooth.c cullIslands.c cullStatsSmooth.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
