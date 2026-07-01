CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS =	$(FLAGS) -c -I$(INCLUDEPATH) -D_LARGEFILE64_SOURCE -D_FILE64 -D_LARGEFILE_SOURCE -DHPLARGEFILES -MMD -MP -fopenmp

SRCS =	strack.c parseTrack.c parseBase.c parseInitialOffsets.c speckleTrack.c \
	sTrackOut.c getInt.c getMask.c parsePar.c writeVrt.c readBothOffsetsStrack.c dopplerEstimation.c \
	mallocPerThreadArrays.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
