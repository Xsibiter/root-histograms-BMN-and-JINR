
#OPT           =-02
OPT           =
 
CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC
LD            = g++
LDFLAGS       =
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
#LIBS          = $(ROOTLIBS) $(SYSLIBS) -lEG -lg2c
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lEG -lgfortran
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)



# Fortran:
F77           = gfortran
#F77           = g77
F77FLAGS      = -fPIC
F77OPT        = $(OPT)
CXXOUT        = -o # keep whitespace after "-o"


#----------------------------

FASTMCO        = CFs_Led_PRF_KK.o ltran12.o fsiw.o fsiini.o

FASTMCS        = CFs_Led_PRF_KK.cxx
 
FASTMCF        = ltran12.F fsiw.F fsiini.F      
       
              
TARGET	    = CFs_Led_PRF_KK
#------------------------------------------------------------------------------
$(TARGET):       $(FASTMCO)
		$(LD)  -O2 $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(FASTMCO) $(TARGET)

%.o : %.cxx
	$(CXX) -O2 $(CXXFLAGS) -c $<

%.o: %.F
	$(F77) $(F77OPT) $(F77FLAGS) $(CXXOUT)$@ -c $<
