

RC     := root-config
ifeq ($(findstring $(MAKECMDTARGET),clean),)
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
RCONFIG := $(wildcard $(shell $(RC) --incdir)/RConfigure.h)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
ifeq ($(ROOTSYS),)
ROOTSYS = ..
endif
include $(ROOTSYS)/etc/Makefile.arch
endif
endif

#------------------------------------------------------------------------------

EVENTO        = Event.$(ObjSuf) EventDict.$(ObjSuf)
EVENTS        = Event.$(SrcSuf) EventDict.$(SrcSuf)
EVENTSO       = libEvent.$(DllSuf)
ifeq ($(PLATFORM),win32)
EVENTLIB      = libEvent.lib
else
EVENTLIB      = $(shell pwd)/$(EVENTSO)
endif


READTREEO    = read_tree.$(ObjSuf)
READTREES    = read_tree.$(SrcSuf)

CONVERTTOMASKO    = convert_to_mask.$(ObjSuf)
CONVERTTOMASKS    = convert_to_mask.$(SrcSuf)

CREATETREEFROMTXTO    = create_tree_from_txt.$(ObjSuf)
CREATETREEFROMTXTS    = create_tree_from_txt.$(SrcSuf)

CREATETREEO    = create_tree.$(ObjSuf)
CREATETREES    = create_tree.$(SrcSuf)

PLOTTREEO    = plot_tree.$(ObjSuf)
PLOTTREES    = plot_tree.$(SrcSuf)

WRITEALLHITSO    = write_all_hits.$(ObjSuf)
WRITEALLHITSS    = write_all_hits.$(SrcSuf)

WRITEALLPMTSO    = write_all_pmts.$(ObjSuf)
WRITEALLPMTSS    = write_all_pmts.$(SrcSuf)

WRITEDETECTORO    = write_detector.$(ObjSuf)
WRITEDETECTORS    = write_detector.$(SrcSuf)

PLOTRADIOACTIVITYHITSO    = plot_radioactivity_hits.$(ObjSuf)
PLOTRADIOACTIVITYHITSS    = plot_radioactivity_hits.$(SrcSuf)

CREATEANGLESO    = create_angles.$(ObjSuf)
CREATEANGLESS    = create_angles.$(SrcSuf)

PLOTANGLESO    = plot_angles.$(ObjSuf)
PLOTANGLESS    = plot_angles.$(SrcSuf)


EXECUTABLEREADTREE         = read_tree$(ExeSuf)
EXECUTABLECONVERTTOMASK         = convert_to_mask$(ExeSuf)
EXECUTABLECREATETREEFROMTXT         = create_tree_from_txt$(ExeSuf)
EXECUTABLECREATETREE         = create_tree$(ExeSuf)
EXECUTABLEPLOTTREE         = plot_tree$(ExeSuf)
EXECUTABLEWRITEALLHITS         = write_all_hits$(ExeSuf)
EXECUTABLEWRITEALLPMTS         = write_all_pmts$(ExeSuf)
EXECUTABLEWRITEDETECTOR         = write_detector$(ExeSuf)
EXECUTABLEPLOTRADIOACTIVITYHITS         = plot_radioactivity_hits$(ExeSuf)
EXECUTABLECREATEANGLES         = create_angles$(ExeSuf)
EXECUTABLEPLOTANGLES         = plot_angles$(ExeSuf)

OBJS          = $(EVENTO) $(READTREEO) $(CONVERTTOMASKO) $(CREATETREEFROMTXTO) $(CREATETREEO) $(PLOTTREEO) $(WRITEALLHITSO) $(WRITEALLPMTSO) $(WRITEDETECTORO) $(PLOTRADIOACTIVITYHITSO) $(CREATEANGLESO) $(PLOTANGLESO)

PROGRAMS      = $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS)  $(EXECUTABLEWRITEALLPMTS)  $(EXECUTABLEWRITEDETECTOR) $(EXECUTABLEPLOTRADIOACTIVITYHITS) $(EXECUTABLECREATEANGLES) $(EXECUTABLEPLOTANGLES)

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------


all:            $(PROGRAMS)

$(EVENTSO):     $(EVENTO)
ifeq ($(ARCH),aix5)
		$(MAKESHARED) $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
		@echo "$@ done"

$(EXECUTABLEREADTREE):       $(EVENTSO) $(READTREEO)
		$(LD) $(LDFLAGS) $(READTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECONVERTTOMASK):       $(EVENTSO) $(CONVERTTOMASKO)
		$(LD) $(LDFLAGS) $(CONVERTTOMASKO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECREATETREEFROMTXT):       $(EVENTSO) $(CREATETREEFROMTXTO)
		$(LD) $(LDFLAGS) $(CREATETREEFROMTXTO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECREATETREE):       $(EVENTSO) $(CREATETREEO)
		$(LD) $(LDFLAGS) -L${WCSIM_BASE_DIR} -I${WCSIM_BASE_DIR} -lWCSimRoot $(CREATETREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTTREE):       $(EVENTSO) $(PLOTTREEO)
		$(LD) $(LDFLAGS) $(PLOTTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


$(EXECUTABLEWRITEALLHITS):       $(EVENTSO) $(WRITEALLHITSO)
		$(LD) $(LDFLAGS) $(WRITEALLHITSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEWRITEALLPMTS):       $(EVENTSO) $(WRITEALLPMTSO)
		$(LD) $(LDFLAGS) $(WRITEALLPMTSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


$(EXECUTABLEWRITEDETECTOR):       $(EVENTSO) $(WRITEDETECTORO)
		$(LD) $(LDFLAGS) $(WRITEDETECTORO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTRADIOACTIVITYHITS):       $(EVENTSO) $(PLOTRADIOACTIVITYHITSO)
		$(LD) $(LDFLAGS) $(PLOTRADIOACTIVITYHITSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECREATEANGLES):       $(EVENTSO) $(CREATEANGLESO)
		$(LD) $(LDFLAGS) $(CREATEANGLESO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTANGLES):       $(EVENTSO) $(PLOTANGLESO)
		$(LD) $(LDFLAGS) $(PLOTANGLESO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"




clean:
		@rm -f $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS) $(EXECUTABLEWRITEALLPMTS) $(EXECUTABLEWRITEDETECTOR) $(EXECUTABLEPLOTRADIOACTIVITYHITS) $(EXECUTABLECREATEANGLES) $(EXECUTABLEPLOTANGLES) $(PROGRAMS)  $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS) $(EXECUTABLEWRITEALLPMTS) $(EXECUTABLEWRITEDETECTOR) $(EXECUTABLEPLOTRADIOACTIVITYHITS) $(EXECUTABLECREATEANGLES) $(EXECUTABLEPLOTANGLES) $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def \
		   *.pcm

.SUFFIXES: .$(SrcSuf)

###
 
Event.$(ObjSuf): Event.h
read_tree.$(ObjSuf): Event.h
convert_to_mask.$(ObjSuf): Event.h
create_tree_from_txt.$(ObjSuf): Event.h
create_tree.$(ObjSuf): Event.h
plot_tree.$(ObjSuf): Event.h
write_all_hits.$(ObjSuf): Event.h
write_all_pmts.$(ObjSuf): Event.h
write_detector.$(ObjSuf): Event.h
plot_radioactivity_hits.$(ObjSuf): Event.h
create_angles.$(ObjSuf): Event.h
plot_angles.$(ObjSuf): Event.h

EventDict.$(SrcSuf): Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
#	$(ROOTCLING) -f $@ -c $^
	rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

