**Include the simulation conditions for mixotrophic growth
$include Simulation_conditions.gms


Parameter
c(j) used to define the objective function for FBA
t(transcripts) used to capture state of transcripts (on-off)
g(genes) used to capture state of genes (on-off);

Variables
v(j) flux values through reaction in network
Obj  this is the value of the objective function for the FBA solutions;

Equations
massbalance(i) mass balance equations for each metabolite
calcobj calculates the dot product of the c vector the flux vector;

massbalance(i).. sum( j,S(i,j)*v(j) )=e=0;
calcobj.. Obj=e=sum( j,c(j)*v(j) );

Model FBA /massbalance, calcobj/;

v.lo(j)=LowerLimits(j);
v.up(j)=UpperLimits(j);

*make sure all genes and transcripts are turned on initially
g(genes)=1;
t(transcripts)=1;

**Knock out individual genes here
*g('Cre02.g120100')=0;

$include iCre1355_Gene_transcript_map.gms

**Knock out individual transcripts here
*t('Cre10_g442600.t1.2')=0;

$include iCre1355_GPR.gms

c(j)=0;
c('R_Biomass_Chlamy_mixo')=1;

**Run FBA maximizing for biomass
solve FBA using lp maximizing Obj;

$ontext
***Code to carry out all single gene deletions (remove comments to use)
file result /All_single_gene_KO_mixo.txt/;
result.pw=500; result.pc=6; result.ps=130;
put result;
Scalar wt_gr;
c(j)=0;
c('R_Biomass_Chlamy_mixo')=1;
solve FBA using lp maximizing Obj;
wt_gr=Obj.l;
v.lo(j)=LowerLimits(j);
v.up(j)=UpperLimits(j);

put "Gene_ID";put "WT_Biomass";put "Mutant_Biomass";
put /;
alias(subg,genes);
loop(subg,
	c(j)=0;
	g(subg)=0;
$include iCre1355_Gene_transcript_map.gms
$include iCre1355_GPR.gms
	c('R_Biomass_Chlamy_mixo')=1;
	solve FBA using lp maximizing Obj;
	put subg.tl;put wt_gr:8:12;
	if(FBA.modelstat ne 1,
                 put 'Infeasible';
	else put Obj.l:8:12;
	);
	put /;
	v.lo(j)=LowerLimits(j);
	v.up(j)=UpperLimits(j); 
	g(genes)=1;
	t(transcripts)=1;
);
putclose;

$offtext

$ontext

****Code to run Flux variability analysis (FVA) (remove comments to use)
file result /FVA_mixo.txt/;
result.pw=500; result.pc=6; result.ps=130;
put result;
Scalar wt_gr;
c(j)=0;
c('R_Biomass_Chlamy_mixo')=1;
solve FBA using lp maximizing Obj;
wt_gr=Obj.l;
v.lo(j)=LowerLimits(j);
v.up(j)=UpperLimits(j);
alias(j,subj)
loop(subj,
	c(j)=0;
	c(subj)=1;
	put subj.tl;
	v.fx('R_Biomass_Chlamy_mixo')=wt_gr;
	solve FBA using lp maximizing Obj;
	put Obj.l:8:12;
	solve FBA using lp minimizing Obj;
	put Obj.l:8:12;put FBA.modelstat; put /;
	v.lo(j)=LowerLimits(j);
	v.up(j)=UpperLimits(j);	
);

$offtext
