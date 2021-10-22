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

clean: 
	rm -rvf *.o bin obj
	
$(BINDIR)/libcerang.so: $(DEPS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -shared $^ -o $@
	
demo: $(OBJDIR)/demo.o $(BINDIR)/libcerang.so
	$(CXX) $(LDFLAGS) $^ -o $(BINDIR)/$@
	
$(OBJDIR)/%.o: %.cpp %.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@	
	
$(OBJDIR):
	mkdir $(OBJDIR)
	
$(BINDIR):
	mkdir $(BINDIR)
