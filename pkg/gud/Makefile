GENPARAMS_TGT=GUD_GENPARAMS.h gud_readgenparams.F
TRAITPARAMS_TGT=GUD_TRAITPARAMS.h gud_readtraitparams.F
TRAITS_TGT=GUD_TRAITS.h gud_init_traits.F gud_readtraits.F
CHECK_TGT=gud_check.F
DIAGS_TGT=gud_diagnostics_init.F
EXF_TGT=GUD_EXF.h GUD_EXF_NAMELIST.h GUD_EXF_DEFAULTS.h gud_exf_init_fixed.F gud_exf_init_varia.F gud_exf_load.F

TGTS=$(GENPARAMS_TGT) $(TRAITPARAMS_TGT) $(TRAITS_TGT) $(CHECK_TGT) $(DIAGS_TGT) $(EXF_TGT)

all: $(TGTS)

$(GENPARAMS_TGT): params.template gentraits
	./gentraits $<

$(TRAITPARAMS_TGT): traitparams.template gentraits
	./gentraits $<

$(TRAITS_TGT): traits.template gentraits
	./gentraits $<

$(CHECK_TGT): GUD_SIZE.h GUD_INDICES.h mkcheck
	./mkcheck > $@

$(DIAGS_TGT): diags.template mkdiags
	./mkdiags < $< > $@

$(EXF_TGT): exf.template genexf
	./genexf < $<

CLEAN:
	@make clean

Clean:
	@make clean

clean: 
	-rm -f $(TGTS)

clean_old:
