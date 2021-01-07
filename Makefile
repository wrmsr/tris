.DEFAULT_GOAL=bin/tris

CFLAGS=-std=c99 -pedantic -Wall -Werror -D_XOPEN_SOURCE=500 -Os

bin:
	mkdir bin

bin/tris: obj/tris.o | bin
	$(CC) $^ -o $@ -lSDL2 -lm
	stat -c %s $@
	strip --strip-unneeded $@
	stat -c %s $@
	upx --best $@
	stat -c %s $@

obj:
	mkdir obj

obj/%.o: %.c | obj
	$(CC) -c $(CFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -r bin/ obj/
