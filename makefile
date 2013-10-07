EXEDIR := scripts
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin

CXX := $(shell root-config --cxx)
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Wshadow $(shell root-config --cflags) -O2 -I $(INCDIR)
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags)
LDLIBS := $(shell root-config --libs) -lMinuit

vpath %.cpp $(SRCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.so $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

# Add new executables to this list
all: make_plots.exe skim_file.exe stack_histos.exe draw_abcd_ratio_plots.exe make_sig_plots.exe fix_skimmed_file.exe calc_abcd.exe count_specific_mass_events.exe draw_npv_plot.exe

# List any object files your executable need to be linked with
$(EXEDIR)/draw_npv_plot.exe: draw_npv_plot.o pu_constants.o
$(EXEDIR)/count_specific_mass_events.exe: count_specific_mass_events.o
$(EXEDIR)/calc_abcd.exe: calc_abcd.o weights.o
$(EXEDIR)/fix_skimmed_file.exe: fix_skimmed_file.o
$(EXEDIR)/make_sig_plots.exe: make_sig_plots.o
$(EXEDIR)/generate_cfa_class.exe: generate_cfa_class.o
$(EXEDIR)/stack_histos.exe: stack_histos.o
$(EXEDIR)/draw_abcd_ratio_plots.exe: draw_abcd_ratio_plots.o
$(EXEDIR)/skim_file.exe: skim_file.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o
$(EXEDIR)/make_plots.exe: make_plots.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))
-include $(MAKEDIR)/cfa.d

$(MAKEDIR)/%.d: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM -MG -MF $@ $< 
	sed -i 's#$*.o#$(OBJDIR)/$*.o $(MAKEDIR)/$*.d#g' $@

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXEDIR)/%.exe:
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# cfa.cpp and cfa.hpp need special treatment.
# Need to figure out how to recompile cfa.* when generate_cfa_class.* changes
# This may also unecessarily run twice (harmless but annoying)
$(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp: $(EXEDIR)/generate_cfa_class.exe example_root_file.root
	./$< example_root_file.root

.PHONY: clean

clean:
	-rm -rf $(EXEDIR)/*.exe $(OBJDIR)/*.o $(MAKEDIR)/*.d $(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp *.exe *.o *.d
	./scripts/remove_backups.sh
