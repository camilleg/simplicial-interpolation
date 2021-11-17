MAKEFLAGS += --jobs=99
UNAME_S := $(shell sh -c 'uname -s 2>/dev/null')
ifeq ($(UNAME_S),Darwin)
  CC = clang++
  # For OS X 10.3.9:
  # CC = g++
  # CFLAGS += -I/usr/X11R6/include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers
  LIBS := -framework GLUT -framework OpenGL
else
  CC = g++
  LIBS := -lglut -lGLU -lGL
endif

CFLAGS := -std=c++17 -O3 -W -Wall -Werror -Weffc++
ifeq ($(UNAME_S),Darwin)
  CFLAGS += -DGL_SILENCE_DEPRECATION # macOS 11.2
endif

# Optional file containing debugging options for CFLAGS.
-include Rules.debug

OBJS_HULL = hull.o ch.o io.o rand.o pointops.o fg.o hullmain.o
OBJS_SI = si.o sammon.o ga.o gacli.o bary.o edahiro.o
OBJS = $(OBJS_SI) $(OBJS_HULL)

EXES = hull si

all: $(EXES)
	./si

hull: $(OBJS_HULL)
	$(CC) -o $@ $^ $(LIBS)

si: $(OBJS_SI)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -rf $(EXES) $(OBJS_HULL) $(OBJS_SI) .depend

DEPENDFLAGS = -MMD -MT $@ -MF $(patsubst %.o,.depend/%.d,$@)
%.o: %.c++
	@mkdir -p .depend
	$(CC) -c $(CFLAGS) $(DEPENDFLAGS) $<

# .depend/*.d is built by %.o:%.c++.
-include $(patsubst %.o,.depend/%.d,$(OBJS))

.PHONY: all clean
