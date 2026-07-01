CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH) -MMD -MP

SRCS =	cullls.c loadLSCullData.c cullLSData.c cullLSStats.c cullLSSmooth.c cullLSIslands.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
