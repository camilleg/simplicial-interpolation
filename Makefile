CC	= g++
CFLAGS	= -O3 -W -Wall
#CFLAGS	= -O0 -g -W -Wall
OBJS	= hull.o ch.o io.o rand.o pointops.o fg.o hullmain.o
OBJS_RSITES = rsites.o
OBJS_SI	= si.o sammon.o ga.o gacli.o det.o bary.o edahiro.o

HDRS	= hull.h points.h pointsites.h stormacs.h
HDRS_SI	= si.h sammon.h ga.h gacli.h bary.h det.h util.h edahiro.h
# gacli.h

EXES    = hull rsites si

all	: doit

doit	: $(EXES)
	./si
#	./si | more

$(OBJS) : $(HDRS)
$(OBJS_SI) : $(HDRS_SI)

%.o: %.c++
	$(CC) -c $(CFLAGS) $<

hullmain.o	: $(HDRS)

hull	: $(OBJS)
	$(CC) -o $@ $^ -lm

rsites	: $(OBJS_RSITES)
	$(CC) -o $@ $^ -lm

$(OBJS_SI): si.h

si: $(OBJS_SI)
	$(CC) -o $@ $^ -lglut -lGLU -lGL -lm

clean	:
	-rm -f $(OBJS) $(OBJS_RSITES) $(OBJS_SI) core a.out $(EXES)
