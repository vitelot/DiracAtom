#################################################################
#								#
# Makefile for the ATOM program.				#
#  to be processes by gmake only				#
# created  8 March 2002						#
#       VDPS							#
#################################################################

SHELL=/bin/tcsh

SYSTEM  = $(shell uname)

BINDIR  =  ../bin
OBJDIR  =  ../objs/$(SYSTEM)
MANDIR  =  ../manual

OUT	=  $(BINDIR)/ATOM.$(SYSTEM)

ifeq "$(SYSTEM)" "Darwin"
CC 	= gcc -DLINUX
FC	= gfortran
LFLAGS	= -lm
endif

ifeq "$(SYSTEM)" "Linux"
CC 	= gcc -DLINUX # -Wall
FC	= gfortran
LFLAGS	= -lm
endif

# Remember: -g with -O not allowed 
CFLAGS  = -O
FFLAGS  = -O

DATE = \"`date "+%d-%b-%Y"`\"

FSOURCES = $(wildcard *.f)
CSOURCES = $(wildcard *.c)
SOURCES  = $(CSOURCES) $(FSOURCES)
ALLH     = $(wildcard *.h)

COBJS = $(patsubst %.c,$(OBJDIR)/%.o,$(CSOURCES))
FOROBJS	=  $(OBJDIR)/rs.o $(OBJDIR)/rsg.o 

INPUTOBJ = $(OBJDIR)/input.o
COBJS := $(filter-out $(INPUTOBJ),$(COBJS))

OBJS = $(COBJS) $(FOROBJS) $(INPUTOBJ)

# the followings depend on relat.h
RELAT = corevnew.o basis_xtrct.o core_util.o sic_nsph.o diag_matrix.o\
	cleanb3.o cmdline.o  coredir.o integ2d.o \
        deltaw_nc.o dentopot_nc.o diagrs.o diagrsg.o exc.o iexcvxc_nc.o \
	input.o magn2_nc.o main.o niocc.o nosphint.o plotfield.o \
	potat.o potcul.o potcul_nsph.o prepat.o printout.o \
	prepare_qn_tables.o
RELAT := $(patsubst %.o,$(OBJDIR)/%.o,$(RELAT))

# the followings depend on globalvar.h
GLOBAL = $(OBJDIR)/main.o

#pippo:
#	@echo $(RELAT)

#all: manual $(OUT)
all: $(OUT)

manual: $(MANDIR)/manual.tex
	cd $(MANDIR) ; $(MAKE) -f Makefile.manual

$(OUT): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LFLAGS) -o $@
	Store

# header file dependencies 
$(RELAT): relat.h
$(GLOBAL): globalvar.h

$(OBJDIR)/cleanb3.o:     alloc_cleanb.h  free_cleanb.h
$(OBJDIR)/cmdline.o:     manual.h
$(OBJDIR)/corevnew.o:    alloc_core.h  free_core.h core.h
$(OBJDIR)/basis_xtrct.o: alloc_core.h  free_core.h core.h
$(OBJDIR)/deltaw_nc.o:   alloc_deltaw.h free_deltaw.h
$(OBJDIR)/diagrs.o:      my_f2c.h
$(OBJDIR)/exc.o:         my_math.h
$(OBJDIR)/gauleg.o:      my_math.h
$(OBJDIR)/main.o:        alloc_main.h main.h
$(OBJDIR)/potat.o:       alloc_potat.h free_potat.h

# compiling procedures
$(FOROBJS): $(OBJDIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@ 

$(COBJS): $(OBJDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(INPUTOBJ): $(SOURCES) $(ALLH) Makefile
	@set date = $(DATE) ; \
	echo "Making $@ with current date" ; \
	$(CC) $(CFLAGS) -DCURRENT_DATE=$$date -c input.c -o $@;

.PHONY: clean

clean:
	rm $(OBJS) $(OUT)
