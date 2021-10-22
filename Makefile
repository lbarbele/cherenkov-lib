SRCDIR := src
OBJDIR := obj
BINDIR := bin

CXXFLAGS += -I$(SRCDIR) -fpic `root-config --cflags`
LDFLAGS += `root-config --libs`

DEPS = $(addprefix $(OBJDIR)/, cherenkov-angular.o)

vpath %.cpp $(SRCDIR) # search for .cpp files under $(SRCDIR) 
vpath %.h $(SRCDIR)   # search for .h files under $(SRCDIR) 
vpath %.o $(OBJDIR)   # search for .o files under $(OBJDIR) 
vpath % $(BINDIR)     # search for anything onde $(BINDIR)

.PHONY: lib demo clean

lib: libcerang.so

clean: 
	rm -rvf *.o bin obj
	
libcerang.so: $(DEPS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -shared $^ -o $(BINDIR)/$@
	
demo: $(OBJDIR)/demo.o libcerang.so
	$(CXX) $(LDFLAGS) $^ -o $(BINDIR)/$@
	
$(OBJDIR)/%.o: %.cpp %.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@	
	
$(OBJDIR):
	mkdir $(OBJDIR)
	
$(BINDIR):
	mkdir $(BINDIR)






#BINDIR := bin
#OBJDIR := obj
#SRCDIR := src

#DIRECTORIES := $(BINDIR) $(OBJDIR)

#vpath %.cpp $(SRCDIR) # search for .cpp files under ./src
#vpath % $(BINDIR)     # search for other files under ./bin

#CXXFLAGS += `root-config --cflags`
#ROOTLIBS =  `root-config --libs`

#all: lib demo
#lib: $(BINDIR)/libcerang.so
#demo: $(BINDIR)/demo
#tests: $(BINDIR)/tests

#$(BINDIR)/tests: $(OBJDIR)/tests.o $(BINDIR)/libcerang.so | $(BINDIR)
#	$(CXX) $^ -o $@ $(ROOTLIBS)

#$(BINDIR)/demo: $(OBJDIR)/demo.o $(BINDIR)/libcerang.so | $(BINDIR)
#	$(CXX) $^ -o $@ $(ROOTLIBS)

#$(BINDIR)/libcerang.so: $(OBJDIR)/cerang.o | $(BINDIR)
#	$(CXX) -shared $< -o $@ 
#	
#$(OBJDIR)/tests.o: $(SRCDIR)/tests.cpp | $(OBJDIR)
#	$(CXX) $(CXXFLAGS) -c $< -o $@
#	
#$(OBJDIR)/demo.o: $(SRCDIR)/demo.cpp | $(OBJDIR)
#	$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJDIR)/cerang.o: $(SRCDIR)/cerang.cpp $(SRCDIR)/cerang.h | $(OBJDIR)
#	$(CXX) $(CXXFLAGS) -fpic -c $< -o $@
#	
#$(DIRECTORIES):
#	mkdir -p $@

#.PHONY: clean

#clean:
#	rm -rvf $(DIRECTORIES)
