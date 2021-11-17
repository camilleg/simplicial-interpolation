MAKEFLAGS += --jobs=99
UNAME_S := $(shell sh -c 'uname -s 2>/dev/null')
ifeq ($(UNAME_S),Darwin)
  CC = clang++
  GL_LDFLAGS = -framework GLUT -framework OpenGL
  # For OS X 10.3.9:
  # CC = g++
  # CFLAGS += -I/usr/X11R6/include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers
LIBS :=
else
  CC = g++
  LIBS := -lglut -lGLU -lGL
endif

CFLAGS := -std=c++17 -O3 -W -Wall -Werror -Weffc++
ifeq ($(UNAME_S),Darwin)
  CFLAGS += -DGL_SILENCE_DEPRECATION # MacOS 11.2
endif

# Optional file containing debugging options for CFLAGS.
-include Rules.debug

OBJS_HULL = hull.o ch.o io.o rand.o pointops.o fg.o hullmain.o
OBJS_SI = si.o sammon.o ga.o gacli.o bary.o edahiro.o

HDRS_HULL = hull.h points.h pointsites.h stormacs.h
HDRS_SI = si.h sammon.h ga.h gacli.h bary.h util.h edahiro.h

EXES = hull si

all: $(EXES)
	./si

$(OBJS_HULL): $(HDRS_HULL)
$(OBJS_SI): $(HDRS_SI)

%.o: %.c++
	$(CC) -c $(CFLAGS) $<

hull: $(OBJS_HULL)
	$(CC) -o $@ $^ $(LIBS)

si: $(OBJS_SI)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f $(OBJS_HULL) $(OBJS_SI) $(EXES)

.PHONY: all clean
