all: figures/01_heatPropagation.png figures/02_Irrigation01.png figures/04_cryo01.png figures/05_SFS_01.png figures/06_PFS_01.png

figures/01_heatPropagation.png: 01_HeatingSurface.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"

figures/02_Irrigation01.png: 02_InfiltrationPulse.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"

# figures/03_layered.png: 03_layeredSoil.ipynb
# 	# Do nothing
# 	#jupyter nbconvert --to notebook --execute --inplace "$<"

figures/04_cryo01.png: 04_Cryosuction.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"

figures/05_SFS_01.png: 05_SeasonallyFrozenSoil.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"

figures/06_PFS_01.png: 06_PermafrostSoil.ipynb
	jupyter nbconvert --to notebook --execute --inplace "$<"

.PHONY: clean

clean:
	rm -rf figures/*.png
	rm -rf pars_def.txt const_def.txt

