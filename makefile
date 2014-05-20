EXEDIR := scripts
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wswitch-default -Wswitch-enum -Wundef -Wunused -Wvariadic-macros -Wwrite-strings -Wabi -Wctor-dtor-privacy -Wnon-virtual-dtor -Wstrict-null-sentinel -Wsign-promo -Wsign-compare #-Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code 
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Wshadow -Woverloaded-virtual -Wold-style-cast $(EXTRA_WARNINGS) $(shell root-config --cflags) -O2 -I $(INCDIR)
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
all: make_plots.exe skim_file.exe stack_histos.exe draw_abcd_ratio_plots.exe make_sig_plots.exe calc_abcd.exe count_specific_mass_events.exe draw_npv_plot.exe make_cutflow_table.exe calc_abcd_new.exe make_reduced_tree.exe efficiencies.exe full_sim_vs_fast_sim.exe pileup_plots.exe qcd_plots.exe qcd_study.exe n_minus_one.exe abcd_obs_compare.exe abcd_predictions.exe compare_qcd_samples.exe qcd_systematic.exe fix_sample.exe composition_systematic.exe new_composition_systematic.exe

# List any object files your executable oneed to be linked with
$(EXEDIR)/new_composition_systematic.exe: new_composition_systematic.o utils.o math.o
$(EXEDIR)/composition_systematic.exe: composition_systematic.o utils.o math.o
$(EXEDIR)/fix_sample.exe: fix_sample.o utils.o timer.o
$(EXEDIR)/qcd_systematic.exe: qcd_systematic.o utils.o timer.o
$(EXEDIR)/compare_qcd_samples.exe: compare_qcd_samples.o utils.o timer.o
$(EXEDIR)/abcd_predictions.exe: abcd_predictions.o utils.o timer.o math.o
$(EXEDIR)/abcd_obs_compare.exe: abcd_obs_compare.o
$(EXEDIR)/n_minus_one.exe: n_minus_one.o timer.o utils.o plotter.o
$(EXEDIR)/qcd_study.exe: qcd_study.o timer.o utils.o plotter.o
$(EXEDIR)/qcd_plots.exe: qcd_plots.o timer.o
$(EXEDIR)/pileup_plots.exe: pileup_plots.o utils.o plotter.o
$(EXEDIR)/full_sim_vs_fast_sim.exe: full_sim_vs_fast_sim.o timer.o
$(EXEDIR)/efficiencies.exe: efficiencies.o timer.o
$(EXEDIR)/draw_npv_plot.exe: draw_npv_plot.o pu_constants.o
$(EXEDIR)/count_specific_mass_events.exe: count_specific_mass_events.o
$(EXEDIR)/calc_abcd.exe: calc_abcd.o weights.o
$(EXEDIR)/calc_abcd_new.exe: calc_abcd_new.o weights.o math.o abcd_calculator.o abcd_count.o
$(EXEDIR)/make_sig_plots.exe: make_sig_plots.o
$(EXEDIR)/generate_cfa_class.exe: generate_cfa_class.o
$(EXEDIR)/stack_histos.exe: stack_histos.o
$(EXEDIR)/draw_abcd_ratio_plots.exe: draw_abcd_ratio_plots.o
$(EXEDIR)/skim_file.exe: skim_file.o lib_jet_met_objects.so event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o cfa_skimmer.o event_handler.o in_json_2012.o
$(EXEDIR)/make_plots.exe: make_plots.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o in_json_2012.o cfa_plotter.o
$(EXEDIR)/make_cutflow_table.exe: make_cutflow_table.o cutflow.o
$(EXEDIR)/make_reduced_tree.exe: make_reduced_tree.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o in_json_2012.o reduced_tree_maker.o

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))
-include $(MAKEDIR)/cfa.d

$(MAKEDIR)/%.d: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MM -MG -MF $@ $< 
	sed -i'' 's#$*.o#$(OBJDIR)/$*.o $(MAKEDIR)/$*.d#g' $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# This is a bit ugly. Shouldn't need the dependency explicitly.
$(EXEDIR)/%.exe: $(OBJDIR)/%.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# cfa.cpp and cfa.hpp need special treatment. Probably cleaner ways to do this.
$(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp: dummy_cfa.all
.SECONDARY: dummy_cfa.all
dummy_cfa.all: $(EXEDIR)/generate_cfa_class.exe example_root_file.root
	./$< $(word 2,$^)
.PRECIOUS: generate_cfa_class.o

.DELETE_ON_ERROR:

.PHONY: clean

clean:
	-rm -rf $(EXEDIR)/*.exe $(OBJDIR)/*.o $(MAKEDIR)/*.d $(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp *.exe *.o *.d
	./scripts/remove_backups.sh
