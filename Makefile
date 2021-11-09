UNAME_S := $(shell sh -c 'uname -s 2>/dev/null')
ifeq ($(UNAME_S),Darwin)
  CC = clang++
  GL_LDFLAGS = -framework GLUT -framework OpenGL
  # For OS X 10.3.9:
  # CC = g++
  # CFLAGS += -I/usr/X11R6/include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers
else
  CC = g++
  GL_LDFLAGS = -lglut -lGLU -lGL
endif

CFLAGS := -std=c++17 -O3 -W -Wall -Werror -Weffc++
ifeq ($(UNAME_S),Darwin)
  CFLAGS += -DGL_SILENCE_DEPRECATION # MacOS 11.2
endif

OBJS_HULL = hull.o ch.o io.o rand.o pointops.o fg.o hullmain.o
OBJS_RSITES = rsites.o
OBJS_SI = si.o sammon.o ga.o gacli.o det.o bary.o edahiro.o

HDRS_HULL = hull.h points.h pointsites.h stormacs.h
HDRS_SI = si.h sammon.h ga.h gacli.h bary.h det.h util.h edahiro.h

EXES = hull rsites si

all: $(EXES)
	./si

$(OBJS_HULL): $(HDRS_HULL)
$(OBJS_SI): $(HDRS_SI)

%.o: %.c++
	$(CC) -c $(CFLAGS) $<

hull: $(OBJS_HULL)
	$(CC) -o $@ $^ -lm

rsites: $(OBJS_RSITES)
	$(CC) -o $@ $^ -lm

si: $(OBJS_SI)
	$(CC) -o $@ $^ $(GL_LDFLAGS) -lm

clean:
	-rm -f $(OBJS_HULL) $(OBJS_RSITES) $(OBJS_SI) $(EXES)

.PHONY: all clean
