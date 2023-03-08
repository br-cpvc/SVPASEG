# Makefile for SVPASEG /GAMIXTURE bundle

CC      = g++
INCLUDE = -I/usr/include  
LDFLAGS = -L/usr/lib  -lstdc++ -lm
CFLAGS  = -O3 $(INCLUDE)


OBJSHOWHEADER =  analyze.o znzlib.o nifti1_io.o showheader.o
PROGSHOWHEADER = ../bin/showheader
OBJSVPASEG = analyze.o znzlib.o nifti1_io.o atlasspec.o nhmrf.o svpaseg.o
PROGSVPASEG = ../bin/svpaseg
OBJRELABEL = relabel_analyze.o analyze.o znzlib.o nifti1_io.o
PROGRELABEL = ../bin/relabel_analyze
OBJCMATRIX = cmatrix.o analyze.o znzlib.o nifti1_io.o
PROGCMATRIX = ../bin/cmatrix
OBJGAMIXTURE = RFRandom.o analyze.o znzlib.o nifti1_io.o nhmrf.o parzen.o atlasspec.o gamixture.o gamixturemain.o
PROGGAMIXTURE = ../bin/gamixture
OBJERODEMASK = analyze.o znzlib.o nifti1_io.o erodemask.o
PROGERODEMASK  = ../bin/erodemask
OBJMEDFILT3D = analyze.o znzlib.o nifti1_io.o medfilt3D.o
PROGMEDFILT3D = ../bin/medfilt3D
OBJKNNFILT = analyze.o znzlib.o nifti1_io.o knnfilt.o
PROGKNNFILT = ../bin/knnfilt
OBJTRANSLATE_ANALYZE = analyze.o znzlib.o nifti1_io.o translate_analyze.o
PROGTRANSLATE_ANALYZE = ../bin/translate_analyze
OBJCONFATLAS = atlasspec.o analyze.o znzlib.o nifti1_io.o config_atlas_dirs.o
PROGCONFATLAS = ../bin/conf_atlas_dirs

all  :  svpaseg gamixture relabel_analyze medfilt3D erodemask knnfilt translate_analyze config_atlas_dirs

svpaseg : $(OBJSVPASEG)
	$(CC) $(CFLAGS)  $(OBJSVPASEG) -o $(PROGSVPASEG)  $(LDFLAGS)

relabel_analyze : $(OBJRELABEL)
	$(CC) $(CFLAGS)  $(OBJRELABEL) -o $(PROGRELABEL)  $(LDFLAGS)

cmatrix : $(OBJCMATRIX)
	$(CC) $(CFLAGS)  $(OBJCMATRIX) -o $(PROGCMATRIX)  $(LDFLAGS)

gamixture : $(OBJGAMIXTURE)
	$(CC) $(CFLAGS)  $(OBJGAMIXTURE) -o $(PROGGAMIXTURE)  $(LDFLAGS)

erodemask: $(OBJERODEMASK)
	$(CC) $(CFLAGS)  $(OBJERODEMASK) -o $(PROGERODEMASK)  $(LDFLAGS)

medfilt3D: $(OBJMEDFILT3D)
	$(CC) $(CFLAGS)  $(OBJMEDFILT3D) -o $(PROGMEDFILT3D)  $(LDFLAGS)

knnfilt: $(OBJKNNFILT)
	$(CC) $(CFLAGS)  $(OBJKNNFILT) -o $(PROGKNNFILT)  $(LDFLAGS)

translate_analyze: $(OBJTRANSLATE_ANALYZE)
	$(CC) $(CFLAGS)  $(OBJTRANSLATE_ANALYZE) -o $(PROGTRANSLATE_ANALYZE)  $(LDFLAGS)

config_atlas_dirs: $(OBJCONFATLAS)
	$(CC) $(CFLAGS)  $(OBJCONFATLAS) -o $(PROGCONFATLAS)  $(LDFLAGS)

znzlib.o               : znzlib.c znzlib.h 
RFRandom.o             : RFRandom.c RFRandom.h
nifti1_io.o            : nifti1_io.c nifti1_io.h nifti1.h znzlib.h
parzen.o               : parzen.c parzen.h analyze.h
analyze.o              : analyze.c analyze.h nifti1_io.h
atlasspec.o            : atlasspec.c atlasspec.h analyze.h 
showheader.o           : showheader.c analyze.h 
nhmrf.o                : nhmrf.c nhmrf.h atlasspec.h analyze.h
svpaseg.o              : svpaseg.c nhmrf.h atlasspec.h analyze.h
relabel_analyze.o      : relabel_analyze.c analyze.h
cmatrix.o              : cmatrix.c analyze.h
gamixture.o            : gamixture.c gamixture.h parzen.h analyze.h nhmrf.h atlasspec.h
gamixturemain.o        : gamixturemain.c gamixture.h parzen.h analyze.h nhmrf.h atlasspec.h RFRandom.h
erodemask.o            : analyze.h erodemask.c
medfilt3D.o            : analyze.h medfilt3D.c
knnfilt.o              : analyze.h knnfilt.c
translate_analyze.o    : analyze.h translate_analyze.c
config_atlas_dirs.o    : atlasspec.h analyze.h config_atlas_dirs.c

clean:
	\rm -f  *.o 
