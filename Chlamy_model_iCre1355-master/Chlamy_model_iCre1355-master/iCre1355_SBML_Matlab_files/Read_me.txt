This zip file (Data S2) contains 3 files:

1. iCre1355_auto.xml - for simulation of photoautotrophic growth
2. iCre1355_mixo.xml - for simulation of mixotrophic growth
3. iCre1355_hetero.xml - for simulation of heterotrophic growth

Each SBML file has all the constraints (uptake rates etc) required for prediction of growth under each condition. This can be 
overriden prior to simulation.

The gene associations within all 3 SBML files is based on transcripts (not genes). Thus, to conduct a single gene deletion
analysis, these transcripts have to be mapped to genes. The transcript to gene map is provided in "transcript_gene_map.mat".
"singleGeneKOsim.m" contains matlab code to conduct FBA analysis and single gene deletion analysis under photoautotrophic,
mixotrophic and heterotrophic conditions.
