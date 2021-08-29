BINDIR := bin
OBJDIR := obj
LIBDIR := lib
SRCDIR := src

DIRECTORIES := $(BINDIR) $(OBJDIR) $(LIBDIR)

vpath %.cpp $(SRCDIR) # search for .cpp files under ./src
vpath % $(BINDIR)     # search for other files under ./bin

CXXFLAGS += `root-config --cflags`
ROOTLIBS =  `root-config --libs`

default:
	@echo Starting parallel compilation
	@make -j4 test

test: $(OBJDIR)/test.o $(OBJDIR)/cerang.o | $(BINDIR)
	$(CXX) $^ -o $(BINDIR)/$@ $(ROOTLIBS)

$(OBJDIR)/%.o: %.cpp $(SRCDIR)/*.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DIRECTORIES):
	mkdir -p $@

.PHONY: clean

clean:
	rm -rvf $(DIRECTORIES)