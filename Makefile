SRCDIR := src
OBJDIR := obj
BINDIR := bin

CXXFLAGS += -I$(SRCDIR) -fpic `root-config --cflags`
LDFLAGS += `root-config --libs` -lgsl -lcblas -lm

DEPS = $(addprefix $(OBJDIR)/, cherenkov-angular.o cherenkov-longi.o)

vpath %.cpp $(SRCDIR) # search for .cpp files under $(SRCDIR) 
vpath %.h $(SRCDIR)   # search for .h files under $(SRCDIR) 
vpath %.o $(OBJDIR)   # search for .o files under $(OBJDIR) 

.PHONY: lib demo clean

lib: $(BINDIR)/libcerang.so
demo: $(BINDIR)/demo

clean: 
	rm -rvf *.o bin obj
	
$(BINDIR)/libcerang.so: $(DEPS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -shared $^ -o $@
	
$(BINDIR)/%: $(OBJDIR)/%.o $(BINDIR)/libcerang.so
	$(CXX) $(LDFLAGS) $^ -o $@

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
$(OBJDIR):
	mkdir $(OBJDIR)
	
$(BINDIR):
	mkdir $(BINDIR)
