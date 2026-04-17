all: figures/00_KFuns.png figures/01_heatPropagation.png figures/02_Irrigation01.png figures/03_SFS.png figures/04_cryo01.png 

figures/00_KFuns.png: 00_Properties.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/00*

figures/01_heatPropagation.png: 01_HeatingSurface.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/01*

figures/02_Irrigation01.png: 02_InfiltrationPulse.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/02*

figures/03_SFS.png: 03_SeasonalFreezeThaw.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/03*

figures/04_cryo01.png: 04_Cryosuction.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"
	touch figures/04*

.PHONY: clean

clean:
	rm -rf figures/*.png
	rm -rf pars_def.txt const_def.txt

