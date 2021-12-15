MAKEFLAGS += --jobs=99
UNAME_S := $(shell sh -c 'uname -s 2>/dev/null')

CFLAGS := -std=c++17 -O3 -W -Wall -Werror -Weffc++
LIBS :=
ifeq ($(UNAME_S),Darwin)
  CFLAGS += -DGL_SILENCE_DEPRECATION # macOS 11.2
  CC = clang++
  LIBS_GLUT := -framework GLUT -framework OpenGL
  # For OS X 10.3.9:
  # CC = g++
  # CFLAGS += -I/usr/X11R6/include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers
else
  CC = g++
  LIBS_GLUT := -lglut -lGLU -lGL
endif

# Optional debugging options for CFLAGS.
-include Rules.debug

OBJS_SI = bary.o callhull.o ga.o gacli.o edahiro.o sammon.o si.o # Camille's code.
OBJS_HULL = ch.o fg.o hull.o io.o pointops.o # Ken Clarkson's code.
OBJS_CORE = $(OBJS_SI) $(OBJS_HULL)

OBJS_TEST = $(OBJS_CORE) selftest.o
OBJS_GLUT = $(OBJS_CORE) glut.o
OBJS = $(sort $(OBJS_GLUT) $(OBJS_TEST))

EXES = selftest glut
all: test demo

test: selftest
	./selftest
demo: glut
	./glut 7 25

selftest: $(OBJS_TEST)
	$(CC) -o $@ $^ $(LIBS)
glut: $(OBJS_GLUT)
	$(CC) -o $@ $^ $(LIBS) $(LIBS_GLUT)

clean:
	-rm -rf $(EXES) $(OBJS) .depend

DEPENDFLAGS = -MMD -MT $@ -MF $(patsubst %.o,.depend/%.d,$@)
%.o: %.c++
	@mkdir -p .depend
	$(CC) -c $(CFLAGS) $(DEPENDFLAGS) $<

# .depend/*.d is built by %.o:%.c++.
-include $(patsubst %.o,.depend/%.d,$(OBJS))

.PHONY: all clean test demo
