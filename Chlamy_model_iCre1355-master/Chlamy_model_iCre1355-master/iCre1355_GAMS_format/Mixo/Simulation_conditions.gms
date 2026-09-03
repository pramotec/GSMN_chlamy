**Include model file
$include iCre1355.gms

UpperLimits(j)=Vmax;

***Mixotrophic conditions
**Carbon sources
**Acetate
LowerLimits('R_EX_ac_e')=-2;
**CO2
LowerLimits('R_EX_co2_e')=-2;

*NITROGEN SOURCE
**NH4
LowerLimits('R_EX_nh4_e')=-1;

**O2 (required because dummy reaction removes O2 evolved from photosynthesis)
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
*Light
LowerLimits('R_EX_photonVis_e')=-80;


*ATP maintenance
LowerLimits('R_ATPM_NGAM')=2.85;
UpperLimits('R_ATPM_NGAM')=2.85;

******phototrophic constraints
*Light inhibited reactions (Chang et al. 2011)
UpperLimits('R_PCHLDR')=0;
UpperLimits('R_PFKh')=0;
UpperLimits('R_G6PADHh')=0;
UpperLimits('R_G6PBDHh')=0;
UpperLimits('R_H2Oth')=0;
UpperLimits('R_FBAh')=0;
LowerLimits('R_FBAh')=0;
UpperLimits('R_STARCH300DEGRA')=6.52604E-05;
UpperLimits('R_STARCH300DEGRB')=6.52604E-05;
UpperLimits('R_STARCH300DEGR2A')=0;
UpperLimits('R_STARCH300DEGR2B')=0;

***Light spectral reactions (comment out desired light source. Here design_growth spectral reaction used)
UpperLimits('R_PRISM_solar_litho')=0;
UpperLimits('R_PRISM_solar_exo')=0;
UpperLimits('R_PRISM_fluorescent_warm_18W')=0;
UpperLimits('R_PRISM_fluorescent_cool_215W')=0;
UpperLimits('R_PRISM_metal_halide')=0;
UpperLimits('R_PRISM_high_pressure_sodium')=0;
UpperLimits('R_PRISM_growth_room')=0;
UpperLimits('R_PRISM_incandescent_60W')=0;
UpperLimits('R_PRISM_white_LED')=0;
UpperLimits('R_PRISM_red_LED_array_653nm')=0;
UpperLimits('R_PRISM_red_LED_674nm')=0;
*UpperLimits('R_PRISM_design_growth')=51.8841087821545;
********************************************************

**Block other biomass reactions
UpperLimits('R_Biomass_Chlamy_auto')=0;
UpperLimits('R_Biomass_Chlamy_hetero')=0;
