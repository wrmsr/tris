.DEFAULT_GOAL=bin/tris

CFLAGS=-std=gnu11 -Wall -Werror -D_XOPEN_SOURCE=500 -Os

bin:
	mkdir bin

bin/tris: obj/tris.o obj/vector.o | bin
	$(CC) $^ -o $@ -lSDL2 -lm -lGL
	stat -c %s $@
	strip --strip-unneeded $@
	stat -c %s $@

obj:
	mkdir obj

obj/%.o: %.c | obj
	$(CC) -c $(CFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -r bin/ obj/
