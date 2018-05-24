# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libcs173h.a libcs173h.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcs173h.so ${prefix}/lib
	cp libcs173h.a ${prefix}/lib
	cp cs173h.h ${prefix}/include

clean:
	rm -rf $(TARGETS)
	rm -rf *.o

libcs173h.a: cs173h_static.o cs173h_gtl_static.o
	$(AR) rcs $@ $^

libcs173h.so: cs173h.o cs173h_gtl.o
	$(CC) -shared $(AM_FCFLAGS) -o libcs173h.so $^ $(AM_LDFLAGS)

cs173h.o: cs173h.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

cs173h_gtl.o: cs173h_gtl.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
cs173h_static.o: cs173h.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

cs173h_gtl_static.o: cs173h_gtl.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
