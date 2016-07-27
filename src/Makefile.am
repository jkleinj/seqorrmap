#===============================================================================
# Makefile.am : autotools
#===============================================================================

bin_PROGRAMS = seqorrmap

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

seqorrmap_SOURCES = \
arg.c arg.h argpdb.h config.h seqorrmap.c seqorrmap.h error.c error.h \
getcorr.c getcorr.h getpdb.c getpdb.h getseq.c getseq.h \
mapcorr.c mapcorr.h matrix.c matrix.h \
modstring.c modstring.h pattern.c pattern.h \
pdb_structure.h putpdb.c putpdb.h \
putmat.c putmat.h puttcl.c puttcl.h \
safe.c safe.h \
select.c select.h \
seq.h \
vector.c vector.h

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

