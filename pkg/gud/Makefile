GENPARAMS_TGT=GUD_GENPARAMS.h gud_readgenparams.F
TRAITPARAMS_TGT=GUD_TRAITPARAMS.h gud_readtraitparams.F
TRAITS_TGT=GUD_TRAITS.h gud_init_traits.F gud_readtraits.F

TGTS=$(GENPARAMS_TGT) $(TRAITPARAMS_TGT) $(TRAITS_TGT)

all: $(TGTS) cog

cog: FORCE
	python -m cogapp/__main__ -c -r gud_check.F GUD_EXF.h GUD_EXF_NAMELIST.h GUD_EXF_DEFAULTS.h gud_exf_init_fixed.F gud_exf_init_varia.F gud_exf_load.F

FORCE:

$(GENPARAMS_TGT): params.template gentraits
	./gentraits $<

$(TRAITPARAMS_TGT): traitparams.template gentraits
	./gentraits $<

$(TRAITS_TGT): traits.template gentraits
	./gentraits $<

CLEAN:
	@make clean

Clean:
	@make clean

clean: 
	-rm -f $(TGTS)

clean_old:
