Each of the folders in this zip file (Data S2) contains iCre1355 in GAMS format with code to enable
running FBA, FVA and single gene deletion analysis under photoautotrophich (Auto), mixotrophic (mixo)
and heterotrophic conditions (hetero).

Each folder contains the following files:
1. iCre1355.gms - model in GAMS format
2. iCre1355_Gene_transcript_map.gms - file providing a mapping between genes and transcripts
3. iCre1355_GPR.gms - the GPR rules from iCre1355
4. Simulation_conditions.gms - specification of simulation conditions for simulations under photoautotrophic,
mixotrophic or heterotrophic growth.
5. FBA_FVA_KO_analysis.gms - code to run FBA, FVA and single deletion analysis (for all genes). The code for
FVA and single gene deletion analysis are commented out. You just need to uncomment them to use them.
