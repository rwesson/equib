FC=gfortran
LD=gfortran
FFLAGS=-ffree-line-length-0 -O3 -fno-backtrace

.PHONY: clean install uninstall

equib06: source/equib06.o
	$(LD) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

clean:
	rm -f equib06 source/*.o source/*.mod

install:
	test -e ${DESTDIR}/usr/share/equib06 || mkdir -p ${DESTDIR}/usr/share/equib06
	test -e ${DESTDIR}/usr/bin || mkdir -p ${DESTDIR}/usr/bin
	install -m 644 atomic_data06/*.* ${DESTDIR}/usr/share/equib06
	install equib06 ${DESTDIR}/usr/bin

uninstall:
	rm -rf ${DESTDIR}/usr/share/equib06
	rm -f ${DESTDIR}/usr/bin/equib06 
