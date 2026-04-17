all: figures/00_KFuns.svg figures/01_heatPropagation.svg figures/02_Irrigation01.svg figures/03_SFS.svg figures/04_cryo01.svg 

figures/00_KFuns.svg: 00_Properties.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/00*

figures/01_heatPropagation.svg: 01_HeatingSurface.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/01*

figures/02_Irrigation01.svg: 02_InfiltrationPulse.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/02*

figures/03_SFS.svg: 03_SeasonalFreezeThaw.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/03*

figures/04_cryo01.svg: 04_Cryosuction.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/04*

.PHONY: clean

clean:
	rm -rf figures/*.svg
	rm -rf pars_def.txt const_def.txt

