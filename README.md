# hippdwi_HcECT
Registration and metric extraction for low-bval DWI (requires hippunfold, and prepdwi --no-regT1)

## instructions:
1. clone the repo to a scratch folder on graham
2. update the config.yml file to point to hippunfold/prepdwi folders
3. do a dry-run: `snakemake -np`
4. run the job using an interactive session, or regularSubmit:
```
regularSubmit -j ShortFat snakemake --use-singularity -j32
```
