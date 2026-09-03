**Include model file
$include iCre1355.gms

UpperLimits(j)=Vmax;

***Mixotrophic conditions
**Carbon sources
**Acetate
LowerLimits('R_EX_ac_e')=-2;

*NITROGEN SOURCE
**NH4
LowerLimits('R_EX_nh4_e')=-0.5;

**O2
LowerLimits('R_EX_o2_e')=-10;

*Exchange of required nutrients
**H2O
LowerLimits('R_EX_h2o_e')=-10;
**H+
LowerLimits('R_EX_h_e')=-10;
**PO4
LowerLimits('R_EX_pi_e')=-10;
**SO4
LowerLimits('R_EX_so4_e')=-10;
**Fe
LowerLimits('R_EX_fe2_e')=-10;
**Mg2+
LowerLimits('R_EX_mg2_e')=-10;
**Na
LowerLimits('R_EX_na1_e')=-10;


*ATP maintenance
LowerLimits('R_ATPM_NGAM')=0.183;
UpperLimits('R_ATPM_NGAM')=0.183;

***Dark simulations (heterotrophic)
LowerLimits('R_EX_co2_e')=0;
UpperLimits('R_STARCH300DEGRA')=0;
UpperLimits('R_STARCH300DEGRB')=0;
UpperLimits('R_STARCH300DEGR2A')=6.52604E-05;
UpperLimits('R_STARCH300DEGR2B')=6.52604E-05;
UpperLimits('R_GLPThi')=0;
UpperLimits('R_ATPSh')=0;
UpperLimits('R_GAPDH_nadp_hi')=0;
UpperLimits('R_MDH_nadp_hi')=0;
UpperLimits('R_MDHC_nadp_hr')=0;
LowerLimits('R_MDHC_nadp_hr')=0;
UpperLimits('R_PPDKh')=0;
UpperLimits('R_IDPh')=0;
UpperLimits('R_PRUK')=0;
UpperLimits('R_RBPCh')=0;
UpperLimits('R_SBP')=0;
LowerLimits('R_H2Oth')=0;
*************************

**Block other biomass reactions
UpperLimits('R_Biomass_Chlamy_auto')=0;
UpperLimits('R_Biomass_Chlamy_mixo')=0;
