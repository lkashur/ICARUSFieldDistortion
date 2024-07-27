CXXFLAGS += -I. $(shell root-config --cflags) -g
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

PROGRAMS = sce_dx_aggregator

all:	clean $(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cpp -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
