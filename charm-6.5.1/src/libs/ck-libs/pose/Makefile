# Default charmc options 
# - changes to be passed in via simulation makefiles 
# **************************************************
CHARM=../../../..
CHARMINC=$(CHARM)/include
CHARMBIN=$(CHARM)/bin
CHARMLIB=$(CHARM)/lib
CHARMC=$(CHARM)/bin/charmc $(OPTS)
POSEINCDIR = $(CHARMINC)/pose

# POSE object definitions - DO NOT MODIFY
#***********************************
OBJECTS=	eventID.o stats.o srtable.o cancel.o memory_temporal.o \
		mempool.o pvtobj.o lbObject.o ldbal.o gvt.o event.o eqheap.o \
		evq.o sim.o rep.o strat.o seq.o cons.o opt.o opt2.o opt3.o \
		spec.o adapt.o adapt2.o adapt3.o adapt4.o adapt5.o pose.o
SEQOBJECTS=     $(OBJECTS:.o=.seq.o)

LIB = libmodulepose.a
SEQLIB = libmoduleseqpose.a
LIBDEST =  $(CHARMLIB)/$(LIB)
LIBSEQDEST =  $(CHARMLIB)/$(SEQLIB)

# POSE LIBRARY GENERATION RULES -
# DO NOT MODIFY
# *******************************

all:  $(LIBDEST) $(LIBSEQDEST)

.PHONY:  pose seqpose

pose: $(LIBDEST)

seqpose:  $(LIBSEQDEST) 

$(LIBDEST):	$(OBJECTS)
	$(CHARMC) -o $(LIBDEST) $(OBJECTS)
	test ! -d $(POSEINCDIR) && mkdir $(POSEINCDIR) || true
	/bin/cp -f  *.h $(POSEINCDIR)
	/bin/cp libmodulepose.dep $(CHARMLIB)/.
	/bin/cp -f  etrans.pl $(CHARMBIN)

$(LIBSEQDEST):	$(SEQOBJECTS)
	$(CHARMC) -o $(LIBSEQDEST) $(SEQOBJECTS)
	test ! -d $(POSEINCDIR) && mkdir $(POSEINCDIR) || true
	/bin/cp -f  *.h $(POSEINCDIR)
	/bin/cp libmodulepose.dep $(CHARMLIB)/libmoduleseqpose.dep
	/bin/cp -f  etrans.pl $(CHARMBIN)

# *******************************

docs: pose.doxy adapt2.C adapt3.C adapt4.C adapt5.C eventID.C lbObject.C spec.C adapt2.h adapt3.h adapt4.h adapt5.h eventID.h lbObject.h pose.C spec.h pose.ci adapt.C ldbal.C srtable.C adapt.h memory_temporal.C mempool.C ldbal.ci srtable.h memory_temporal.ci mempool.ci pose.h cancel.C stats.C cancel.h ldbal.h pvtobj.C stats.ci memory_temporal.h mempool.h pvtobj.h chpt.h seq.C cons.C evq.C rep.C stats.h seq.h cons.h evq.h opt2.C rep.h opt2.h strat.C eqheap.C gvt.C sim.C strat.h eqheap.h gvt.ci opt3.C sim.ci opt3.h event.C gvt.h opt.C sim.h event.h opt.h
	doxygen pose.doxy

# HOUSE-KEEPING RULES
#********************
clean:
	rm -f charmrun conv-host *.o
	rm -f *.def.h *.decl.h *.stamp
	rm -f gmon.out #*#
	rm -f core *~
	rm -f TAGS *.headers
	rm -f libmodulepose.a
	rm -f libmoduleseqpose.a

CIDEPENDFILE = Make.cidepends
DEPENDFILE = Make.depends

include $(CHARM)/tmp/charm.mk
include $(CIDEPENDFILE)
include $(DEPENDFILE)

depends: cidepends.default
	echo "Building ci files ..."
	for cifile in $(wildcard *.ci); do \
	    $(CHARMC) $$cifile ; \
	done
	echo "Creating " $(DEPENDFILE) " ..."; 	\
	if [ -f $(DEPENDFILE) ]; then \
           /bin/cp -f $(DEPENDFILE) $(DEPENDFILE).old; \
        fi; \
	echo '#generated by make depends' > $(DEPENDFILE); \
        for i in $(OBJECTS) ; do \
	      SRCFILE=`basename $$i .o`.C ; \
              echo "checking dependencies for $$i : $$SRCFILE" ; \
              g++ -MM -Wno-deprecated -I$(CHARMINC) $$SRCFILE >> $(DEPENDFILE); \
              echo '	$$(CHARMC) -o '$$i $$SRCFILE >> $(DEPENDFILE) ; \
        done; 	\
        for i in $(SEQOBJECTS) ; do \
	      SRCFILE=`basename $$i .seq.o`.C ; \
              echo "checking dependencies for $$i : $$SRCFILE" ; \
              g++ -MM -MT $$i -Wno-deprecated -I$(CHARMINC) $$SRCFILE >> $(DEPENDFILE); \
              echo '	$$(CHARMC) -o '$$i -DSEQUENTIAL_POSE=1 $$SRCFILE >> $(DEPENDFILE) ; \
        done; 

