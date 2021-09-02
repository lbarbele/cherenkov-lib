BINDIR := bin
OBJDIR := obj
SRCDIR := src

DIRECTORIES := $(BINDIR) $(OBJDIR)

#vpath %.cpp $(SRCDIR) # search for .cpp files under ./src
#vpath % $(BINDIR)     # search for other files under ./bin

CXXFLAGS += `root-config --cflags`
ROOTLIBS =  `root-config --libs`

all: lib demo
lib: $(BINDIR)/libcerang.so
demo: $(BINDIR)/demo

$(BINDIR)/demo: $(OBJDIR)/demo.o $(BINDIR)/libcerang.so | $(BINDIR)
	$(CXX) $^ -o $@ $(ROOTLIBS)

$(BINDIR)/libcerang.so: $(OBJDIR)/cerang.o | $(BINDIR)
	$(CXX) -shared $< -o $@ 
	
$(OBJDIR)/demo.o: $(SRCDIR)/demo.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/cerang.o: $(SRCDIR)/cerang.cpp $(SRCDIR)/cerang.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -fpic -c $< -o $@
	
$(DIRECTORIES):
	mkdir -p $@

.PHONY: clean

clean:
	rm -rvf $(DIRECTORIES)
