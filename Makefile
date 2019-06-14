FC=gfortran
LD=gfortran
FFLAGS=-ffree-line-length-0 -O3 -fno-backtrace

.PHONY: clean install uninstall

OS := $(shell uname)
ifeq ($(OS),Darwin)
  PREFIX=/usr/local
else
  PREFIX=/usr
endif

equib06: source/equib06.o
	$(LD) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

clean:
	rm -f equib06 source/*.o source/*.mod

install:
	test -e ${DESTDIR}${PREFIX}/share/equib06 || mkdir -p ${DESTDIR}${PREFIX}/share/equib06
	test -e ${DESTDIR}${PREFIX}/bin || mkdir -p ${DESTDIR}${PREFIX}/bin
	install -m 644 atomic_data06/*.* ${DESTDIR}${PREFIX}/share/equib06
	install equib06 ${DESTDIR}${PREFIX}/bin

uninstall:
	rm -rf ${DESTDIR}${PREFIX}/share/equib06
	rm -f ${DESTDIR}${PREFIX}/bin/equib06
