GENPARAMS_TGT=GUD_GENPARAMS.h gud_readgenparams.F
TRAITPARAMS_TGT=GUD_TRAITPARAMS.h gud_readtraitparams.F
TRAITS_TGT=GUD_TRAITS.h gud_init_traits.F gud_readtraits.F
EXF_TGT=GUD_EXF.h GUD_EXF_NAMELIST.h GUD_EXF_DEFAULTS.h gud_exf_init_fixed.F gud_exf_init_varia.F gud_exf_load.F
CHECK_TGT=gud_check.F

TGTS=$(GENPARAMS_TGT) $(TRAITPARAMS_TGT) $(TRAITS_TGT) $(CHECK_TGT) $(DIAGS_TGT) $(EXF_TGT) $(CHECK_TGT)

all: $(TGTS)

$(GENPARAMS_TGT): params.template gentraits
	./gentraits $<

$(TRAITPARAMS_TGT): traitparams.template gentraits
	./gentraits $<

$(TRAITS_TGT): traits.template gentraits
	./gentraits $<

$(EXF_TGT): genexf $(EXF_TGT:=.template)
	./genexf $(EXF_TGT)

$(CHECK_TGT): gencheck $(CHECK_TGT:=.template) GUD_SIZE.h GUD_INDICES.h
	./gencheck $@

CLEAN:
	@make clean

Clean:
	@make clean

clean: 
	-rm -f $(TGTS)

clean_old:
