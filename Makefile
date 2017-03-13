COMPILE := x86_64-w64-mingw32-g++
LINK := x86_64-w64-mingw32-g++
COMPILE_FLAGS := -xc++ -fpermissive -g -O3 -ffast-math -Wno-write-strings -pipe -MMD
LINK_FLAGS := -static-libgcc -static-libstdc++ -Wl,--stack,67108864
LIBS := -lgdi32 -lcomdlg32 -lole32 -loleaut32 -lcomctl32 -luuid
WINDRES := x86_64-w64-mingw32-windres

FRAKTAL_SOURCES_CPP = \
fraktal_sft/CDecNumber.cpp \
fraktal_sft/CFixedFloat.cpp \
fraktal_sft/dbl_functions.cpp \
fraktal_sft/exp_functions2.cpp \
fraktal_sft/exp_functions.cpp \
fraktal_sft/fraktal_sft.cpp \
fraktal_sft/listbox.cpp \
fraktal_sft/main.cpp \
fraktal_sft/newton.cpp

FRAKTAL_SOURCES_H = \
fraktal_sft/CDecNumber.h \
fraktal_sft/CFixedFloat.h \
fraktal_sft/complex.h \
fraktal_sft/floatexp.h \
fraktal_sft/fraktal_sft.h \
fraktal_sft/listbox.h \
fraktal_sft/resource.h

LDBL_SOURCES_CPP = \
ldbl64/ldbl.cpp

COMMON_SOURCES_CPP = \
common/FolderBrowser.cpp \
common/getimage.cpp \
common/parallell.cpp \
common/StringVector.cpp \
common/tooltip.cpp

COMMON_SOURCES_H = \
common/FolderBrowser.h \
common/getimage.h \
common/parallell.h \
common/StringVector.h \
common/tooltip.h 

JPEG_SOURCES_CPP = \
jpeg/jpeg_static.cpp

JPEG_SOURCES_C = \
jpeg/cdjpeg.c \
jpeg/jcapimin.c \
jpeg/jcapistd.c \
jpeg/jccoefct.c \
jpeg/jccolor.c \
jpeg/jcdctmgr.c \
jpeg/jchuff.c \
jpeg/jcinit.c \
jpeg/jcmainct.c \
jpeg/jcmarker.c \
jpeg/jcmaster.c \
jpeg/jcomapi.c \
jpeg/jcparam.c \
jpeg/jcphuff.c \
jpeg/jcprepct.c \
jpeg/jcsample.c \
jpeg/jctrans.c \
jpeg/jdapimin.c \
jpeg/jdapistd.c \
jpeg/jdatadst.c \
jpeg/jdatasrc.c \
jpeg/jdcoefct.c \
jpeg/jdcolor.c \
jpeg/jddctmgr.c \
jpeg/jdhuff.c \
jpeg/jdinput.c \
jpeg/jdmainct.c \
jpeg/jdmarker.c \
jpeg/jdmaster.c \
jpeg/jdmerge.c \
jpeg/jdphuff.c \
jpeg/jdpostct.c \
jpeg/jdsample.c \
jpeg/jdtrans.c \
jpeg/jerror.c \
jpeg/jfdctflt.c \
jpeg/jfdctfst.c \
jpeg/jfdctint.c \
jpeg/jidctflt.c \
jpeg/jidctfst.c \
jpeg/jidctint.c \
jpeg/jidctred.c \
jpeg/jmemansi.c \
jpeg/jmemmgr.c \
jpeg/jquant1.c \
jpeg/jquant2.c \
jpeg/jutils.c \
jpeg/rdbmp.c \
jpeg/rdcolmap.c \
jpeg/rdgif.c \
jpeg/rdppm.c \
jpeg/rdrle.c \
jpeg/rdswitch.c \
jpeg/rdtarga.c \
jpeg/transupp.c \
jpeg/wrbmp.c \
jpeg/wrgif.c \
jpeg/wrppm.c \
jpeg/wrrle.c \
jpeg/wrtarga.c

JPEG_SOURCES_H = \
jpeg/cderror.h \
jpeg/cdjpeg.h \
jpeg/jchuff.h \
jpeg/jconfig.h \
jpeg/jdct.h \
jpeg/jdhuff.h \
jpeg/jerror.h \
jpeg/jinclude.h \
jpeg/jmemsys.h \
jpeg/jmorecfg.h \
jpeg/jpegint.h \
jpeg/jpeglib.h \
jpeg/jversion.h \
jpeg/transupp.h

SOURCES_CPP = $(FRAKTAL_SOURCES_CPP) $(COMMON_SOURCES_CPP) $(LDBL_SOURCES_CPP) $(JPEG_SOURCES_CPP)
SOURCES_C = $(JPEG_SOURCES_C)
SOURCES_H = $(FRAKTAL_SOURCES_H) $(COMMON_SOURCES_H) $(JPEG_SOURCES_H)

SOURCES = $(SOURCES_CPP) $(SOURCES_C) $(SOURCES_H)

OBJECTS_CPP := $(patsubst %.cpp,%.o,$(SOURCES_CPP))
OBJECTS_C := $(patsubst %.c,%.o,$(SOURCES_C))
OBJECTS := $(OBJECTS_CPP) $(OBJECTS_C) res.o

DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))

all: fraktal_sft64.exe

clean:
	rm -f $(OBJECTS) $(DEPENDS)

fraktal_sft64.exe: $(OBJECTS)
	$(LINK) -o fraktal_sft64.exe $(OBJECTS) $(LINK_FLAGS) $(LIBS)

res.o: fraktal_sft/fraktal_sft.rc
	$(WINDRES) -i fraktal_sft/fraktal_sft.rc -o res.o

ldbl64.dll: $(LDBL_SOURCES_CPP)
	$(COMPILE) -o ldbl64.dll $(COMPILE_FLAGS) -shared $(LDBL_SOURCES_CPP) $(LINK_FLAGS)

%.o: %.cpp
	$(COMPILE) $(COMPILE_FLAGS) -o $@ -c $<

%.o: %.c
	$(COMPILE) $(COMPILE_FLAGS) -o $@ -c $<

-include $(DEPENDS)
