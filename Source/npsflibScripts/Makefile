include ~/config/make.conf
include $(NPSFPATH)/common.mk

NPSFARGS = -I$(NPSFINCPATH) -L$(NPSFLIBPATH) -lnpsf
GSLARGS  = -I$(GSLINCPATH) -L$(GSLLIBPATH) -lgsl -lgslcblas
FFTWARGS = -I$(FFTWINCPATH) -L$(FFTWLIBPATH) -lfftw3
MATARGS  = -I$(MATINCPATH) -L$(MATLIBPATH) $(MATLIBLIST)
LDFLAGS  = $(NPSFARGS) $(GSLARGS) $(FFTWARGS) $(MATARGS) $(ALGLIBARGS)

all: libnpsf.a
libnpsf.a: $(NPSFLIBPATH)/libnpsf.a
	$(MAKE) -C $(NPSFLIBPATH)

eigenstatesBO: eigenstatesBO.cpp $(NPSFLIBPATH)/libnpsf.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

doc:
	$(DOXYGEN)

clean:
	rm -f *.o
	rm -f $(NPSFPATH)/$(NPSFLIBPATH)/*.a
	rm -f $(NPSFPATH)/$(NPSFSRCPATH)/*.[a,o,d]


.PHONY: all doc clean