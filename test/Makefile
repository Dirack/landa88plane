# Makefile
# 
# Purpose: Run all tests with 'make'.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programmer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

DEP= ../raytrace.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c ../layer2d.c ../model2d.c
TESTS=test_hubral.x
EXP=test_model3 test_model4
FILES=betas.rsf rnips.rsf m0s.rsf t0s.rsf

define EXPERIMENT
	mkdir -p $(1)/$(2)/cre $(1)/$(2)/cds
	scons -C $(1) nit=200 cds=n
	cd $(1); for f in $(FILES) $(notdir $(wildcard $(1)/otsemb_* $(1)/vspline_* $(1)/zspline*)); do <$$f sfput --out=stdout > $(2)/cre/$$f; done
	scons -C $(1) nit=200 cds=y
	cd $(1); for f in $(FILES) $(notdir $(wildcard $(1)/otsemb_* $(1)/vspline_* $(1)/zspline*)); do <$$f sfput --out=stdout > $(2)/cds/$$f; done
endef

all:	test exp

test:	$(TESTS)

exp:	$(EXP)

test_hubral.x:	test_hubral.c ../grid2.c ../atela.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_model3:	model3/SConstruct
	scons -C model3 nit=200 inverted-original.view

test_model4:	model4/SConstruct
	scons -C model4 nit=200 inverted-original.view

fullexp:	model4
	$(foreach expdir,$^,$(call EXPERIMENT,$(expdir),exp_$(shell date +%d-%m-%Y)))

.PHONY:	$(TESTS)
