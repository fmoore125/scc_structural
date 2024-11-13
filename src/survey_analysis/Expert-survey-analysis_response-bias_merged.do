* Start log
clear all 
*set more off
   
**************************************************************************
* Analysis corresponds to results reported in SI sections S.2.2.4 and S.2.2.5 of Moore et al. 2024.
* The data for the analysis here cannot be shared to protect the anonymity of respondents. 
* We here solely document the code to analyse it as well as results at an aggregate level.
* Results are insignificant expect where explicitly mentioned with corresponding p-values.
**************************************************************************
   
* Housekeeping 

use ".../merged-SCC-data.dta"

**************************************************************************

**************************************************************************
  
 gen Europe =.
 replace Europe=1 if continent=="Europe"
 replace Europe=0 if continent=="NorthAmerica" | continent=="Asia" | continent=="Australia & Oceania" 
 
  gen NorthAmerica =.
 replace NorthAmerica=1 if continent=="NorthAmerica"
 replace NorthAmerica=0 if continent=="Europe" | continent=="Asia" | continent=="Australia & Oceania" 
 
   gen Asia =.
 replace Asia=1 if continent=="Asia"
 replace Asia=0 if continent=="Europe" | continent=="NorthAmerica" | continent=="Australia & Oceania" 
 
    gen Oceania =.
 replace Oceania=1 if continent=="Australia & Oceania"
 replace Oceania=0 if continent=="Europe" | continent=="NorthAmerica" | continent=="Asia" 
 
 
**************************************************************************
*Compare identified respondents with non-respondents, which include anonymous respondents
**************************************************************************
    
ttest Europe, by(identified_participant)
ttest NorthAmerica, by(identified_participant)
ttest Asia, by(identified_participant)
ttest Oceania, by(identified_participant)

ttest male, by(identified_participant)
tab male identified_participant, chi

ttest year_of_phd, by(identified_participant)
ttest documents , by(identified_participant)
ttest hindex , by(identified_participant)

ttest median_20102030_scc , by(identified_participant)
ttest median_20102030_dr, by(identified_participant)

ttest frameworkexpansion , by(identified_participant)
ranksum frameworkexpansion , by(identified_participant)

sum earth_system_updates, d
ttest earth_system_updates  , by(identified_participant)
ranksum earth_system_updates  , by(identified_participant)

ttest carboncycle  , by(identified_participant)
ttest climatemodel  , by(identified_participant)
ttest  tippingpoints  , by(identified_participant)
ttest  tippingpoints2 , by(identified_participant)
ttest  persistentgrowthdamages , by(identified_participant)
ttest  epsteinzin  , by(identified_participant)
ttest limitedlysubstitutablegoods  , by(identified_participant)
ttest  ambiguitymodeluncertainty  , by(identified_participant)
ttest inequalityaversion  , by(identified_participant)
ttest learning  , by(identified_participant)
ttest alternativeethicalapproachesnotd  , by(identified_participant)

ranksum carboncycle  , by(identified_participant)
ranksum climatemodel  , by(identified_participant)
ranksum  tippingpoints  , by(identified_participant)
ranksum  tippingpoints2 , by(identified_participant)
ranksum  persistentgrowthdamages , by(identified_participant)
ranksum  epsteinzin  , by(identified_participant)
ranksum  ambiguitymodeluncertainty  , by(identified_participant)
ranksum limitedlysubstitutablegoods  , by(identified_participant)
ranksum inequalityaversion  , by(identified_participant)
ranksum learning  , by(identified_participant)
ranksum alternativeethicalapproachesnotd  , by(identified_participant)

ttest empiricalimprovement , by(identified_participant)
*--> higher, p=0.0160
ranksum empiricalimprovement , by(identified_participant)
*--> higher, p=0.0106

**************************************************************************
*Strategic response bias
**************************************************************************

ttest scc_true_central, by(anonym)
ttest scc_true_2p5perc, by(anonym)
ttest scc_true_97p5perc, by(anonym)

ttest scc_lit_central, by(anonym)
*signif, p=0.0151
ttest scc_lit_2p5perc, by(anonym)
ttest scc_lit97p5perc, by(anonym)

gen scc_wedge_25 = scc_true_2p5perc - scc_lit_2p5per
gen scc_wedge_975 = scc_true_97p5perc - scc_lit97p5perc

ttest scc_wedge_value, by(anonym)
ttest scc_wedge_25, by(anonym)
ttest scc_wedge_975, by(anonym)

ttest sccwedge_drivers_other, by(anonym) 
* More weight on other drivers, p=0.0040
ttest sccwedge_drivers_other_number, by(anonym)

gen q2_comment=. if q2_num_comments==.
replace q2_comment=1 if q2_num_comments>0
replace q2_comment=0 if q2_num_comments==0
tab q2_comment anonym, chi

tab sccwedge_drivers_other_mention1 anonym, chi

gen q4_comment=. if next_stepsnumberofmentions==.
replace q4_comment=1 if next_stepsnumberofmentions>0
replace q4_comment=0 if next_stepsnumberofmentions==0
tab q4_comment anonym, chi
 
ttest next_stepsnumberofmentions, by(anonym)
 
ttest sccwedge_drivers_pure_time_pref, by(anonym)
* less weight on pure time discounting, p=0.0902

 
 **************************************************************************

 sum scc_lit_central, d
 
ttest scc_wedge_value==0
*p=0.0000
ttest scc_wedge_value==0 if scc_wedge==1 
*p=0.0000
ttest scc_wedge_value==0 if scc_wedge==2 
*p=0.0975
 
 *** Corr lit strengh and improvement
olog q2_meanscc_earth_system q2_include_earth_system
olog q2_meanscc_tipping_climate_syste q2_include_tipping_climate_syste
*-->p=0.004
olog q2_meanscc_tipping_damages q2_include_tipping_damages 
olog q2_meanscc_persistent_damages q2_include_persistent_damage
olog q2_meanscc_substitutability q2_include_substitutability 
 *-->borderline p=0.106
 olog q2_meanscc_epsteinzin q2_include_epsteinzin
 *-->p=0.000
olog q2_meanscc_model_uncertainty q2_include_model_uncertainty
  *-->p=0.086
olog q2_meanscc_learning q2_include_learning
 olog q2_meanscc_distributional_weight q2_include_distributional_weight
 
 
 ****
 regr sccwedge_drivers_substitutabilit q2_include_substitutability
 *-->p=0.012
 

  *** 
gen SCCwedge= scc_true_central-scc_lit_central

pwcorr SCCwedge sccwedge_drivers_earth_system sccwedge_drivers_tipping_damage sccwedge_drivers_persistent_dama sccwedge_drivers_substitutabilit sccwedge_drivers_epsteinzin sccwedge_drivers_model_uncertain sccwedge_drivers_learning sccwedge_drivers_distributional sccwedge_drivers_endogenous_tech sccwedge_drivers_ethical_approac sccwedge_drivers_adaptation sccwedge_drivers_tipping_climate sccwedge_drivers_other sccwedge_drivers_pure_time_pref sccwedge_drivers_emuc sccwedge_drivers_damage_function, star(0.05)

regr SCCwedge sccwedge_drivers_earth_system sccwedge_drivers_tipping_damage sccwedge_drivers_persistent_dama sccwedge_drivers_substitutabilit sccwedge_drivers_epsteinzin sccwedge_drivers_model_uncertain sccwedge_drivers_learning sccwedge_drivers_distributional sccwedge_drivers_endogenous_tech sccwedge_drivers_ethical_approac sccwedge_drivers_adaptation sccwedge_drivers_tipping_climate sccwedge_drivers_other, r
regr SCCwedge sccwedge_drivers_pure_time_pref sccwedge_drivers_emuc sccwedge_drivers_damage_function, r
   
**************************************************************************
*Compare early versus late respondents
**************************************************************************

drop if scc_lit_central==.
sum interviewnumberongoing if interviewnumberongoing<1000, d
gen early_v_late = .
replace early_v_late=1 if interviewnumberongoing<645
replace early_v_late=0 if interviewnumberongoing>644 & interviewnumberongoing<1000

ttest scc_true_central, by(early_v_late)
ttest scc_true_2p5perc, by(early_v_late)
ttest scc_true_97p5perc, by(early_v_late)

ttest scc_lit_central, by(early_v_late)
ttest scc_lit_2p5perc, by(early_v_late)
ttest scc_lit97p5perc, by(early_v_late)


ttest scc_wedge_positive, by(early_v_late)
tab scc_wedge_positive early_v_late, chi

ttest scc_wedge_value, by(early_v_late)
ttest scc_wedge_25, by(early_v_late)
ttest scc_wedge_975, by(early_v_late)

ttest sccwedge_drivers_other, by(early_v_late) 
ttest sccwedge_drivers_other_number, by(early_v_late)


tab q2_comment early_v_late, chi

tab sccwedge_drivers_other_mention1 early_v_late, chi


gen next_stepsnumberofmention= .
replace next_stepsnumberofmention=. if next_stepsnumberofmentions==.
replace next_stepsnumberofmention=0 if next_stepsnumberofmentions==0
replace next_stepsnumberofmention=1 if next_stepsnumberofmentions>0

tab next_stepsnumberofmention early_v_late, chi
*Chi2 p=0.007
 
ttest next_stepsnumberofmentions, by(early_v_late)
 ttest sccwedge_drivers_pure_time_pref, by(early_v_late)


**************************************************************************
*Additional analyses of the merged dataset
**************************************************************************

regr scc_true_central  median_20102030_scc, r
regr scc_lit_central  median_20102030_scc, r

drop if identified_participant==0
drop if median_20102030_scc>1000

 ttest scc_true_central, by(male)
ttest scc_lit_central, by(male)
ttest scc_wedge_value, by(male)


 **************
 
 gen q2_include_earth_system_sagree = .
 replace q2_include_earth_system_sagree= 1 if q2_include_earth_system==4 | q2_include_earth_system==5
  replace q2_include_earth_system_sagree= 0 if q2_include_earth_system==1 | q2_include_earth_system==2 | q2_include_earth_system==3

 probit q2_include_earth_system_sagree earth_system_updates
 
 
  gen q2_include_tipping_climate_sa = .
 replace q2_include_tipping_climate_sa= 1 if q2_include_tipping_climate_syste==4 | q2_include_tipping_climate_syste==5
  replace q2_include_tipping_climate_sa= 0 if q2_include_tipping_climate_syste==1 | q2_include_tipping_climate_syste==2 | q2_include_tipping_climate_syste==3
 
 probit q2_include_tipping_climate_sa tippingpoints
 
   gen q2_include_tipping_damages_sa = .
 replace q2_include_tipping_damages_sa= 1 if q2_include_tipping_damages==4 | q2_include_tipping_damages==5
  replace q2_include_tipping_damages_sa= 0 if q2_include_tipping_damages==1 | q2_include_tipping_damages==2 | q2_include_tipping_damages==3
 
 probit q2_include_tipping_damages_sa tippingpoints2
 
    gen q2_include_persistent_damage_sa = .
 replace q2_include_persistent_damage_sa= 1 if q2_include_persistent_damage==4 | q2_include_persistent_damage==5
  replace q2_include_persistent_damage_sa= 0 if q2_include_persistent_damage==1 | q2_include_persistent_damage==2 | q2_include_persistent_damage==3
 
 probit q2_include_persistent_damage_sa persistentgrowthdamages
 
     gen q2_include_substitutability_sa = .
 replace q2_include_substitutability_sa= 1 if q2_include_substitutability==4 | q2_include_substitutability==5
  replace q2_include_substitutability_sa= 0 if q2_include_substitutability==1 | q2_include_substitutability==2 | q2_include_substitutability==3
 
 probit q2_include_substitutability_sa limitedlysubstitutablegoods
 
     gen q2_include_epsteinzin_sa = .
 replace q2_include_epsteinzin_sa= 1 if q2_include_epsteinzin==4 | q2_include_epsteinzin==5
  replace q2_include_epsteinzin_sa= 0 if q2_include_epsteinzin==1 | q2_include_epsteinzin==2 | q2_include_epsteinzin==3
  
 probit q2_include_epsteinzin_sa epsteinzin
 
     gen q2_include_model_uncertainty_sa = .
 replace q2_include_model_uncertainty_sa= 1 if q2_include_model_uncertainty==4 | q2_include_model_uncertainty==5
  replace q2_include_model_uncertainty_sa= 0 if q2_include_model_uncertainty==1 | q2_include_model_uncertainty==2 | q2_include_model_uncertainty==3
 
 probit q2_include_model_uncertainty_sa ambiguitymodeluncertainty
 
     gen q2_include_learning_sa = .
 replace q2_include_learning_sa= 1 if q2_include_learning==4 | q2_include_learning==5
  replace q2_include_learning_sa= 0 if q2_include_learning==1 | q2_include_learning==2 | q2_include_learning==3
 
 probit q2_include_learning_sa learning
 
     gen q2_include_distributional_sa = .
 replace q2_include_distributional_sa= 1 if q2_include_distributional_weight==4 | q2_include_distributional_weight==5
  replace q2_include_distributional_sa= 0 if q2_include_distributional_weight==1 | q2_include_distributional_weight==2 | q2_include_distributional_weight==3
  
 olog q2_include_distributional_sa inequalityaversion
 reg q2_include_distributional_weight inequalityaversion
 
probit q2_include_distributional_sa inequalityaversion

 **************
 reg sccwedge_drivers_earth_system earth_system_updates
  reg sccwedge_drivers_tipping_climate tippingpoints
 reg sccwedge_drivers_tipping_damage tippingpoints2

  reg sccwedge_drivers_persistent_dama persistentgrowthdamages
    *sign

 reg sccwedge_drivers_substitutabilit limitedlysubstitutablegoods
 
  reg sccwedge_drivers_epsteinzin epsteinzin
    *sign
 reg sccwedge_drivers_model_uncertain ambiguitymodeluncertainty
 reg sccwedge_drivers_learning learning
 reg sccwedge_drivers_distributional_ inequalityaversion
    *sign
	
  **************
   
  gen sccwedge_drivers_weight_convent =.
   replace sccwedge_drivers_weight_convent= sccwedge_drivers_pure_time_pref + sccwedge_drivers_emuc + sccwedge_drivers_damage_function

    gen sccwedge_drivers_weight_struct =.
   replace sccwedge_drivers_weight_struct=  sccwedge_drivers_earth_system + sccwedge_drivers_tipping_damage  + sccwedge_drivers_persistent_dama  + sccwedge_drivers_substitutabilit  + sccwedge_drivers_epsteinzin  + sccwedge_drivers_model_uncertain  + sccwedge_drivers_learning  + sccwedge_drivers_distributional_  + sccwedge_drivers_endogenous_tech  + sccwedge_drivers_ethical_approac  + sccwedge_drivers_adaptation +  sccwedge_drivers_tipping_climate

     gen sccwedge_drivers_weight_str_meta =.
   replace sccwedge_drivers_weight_str_meta=  sccwedge_drivers_earth_system + sccwedge_drivers_tipping_damage  + sccwedge_drivers_persistent_dama  + sccwedge_drivers_substitutabilit  + sccwedge_drivers_epsteinzin  + sccwedge_drivers_model_uncertain  + sccwedge_drivers_learning  + sccwedge_drivers_distributional_  + sccwedge_drivers_ethical_approac  + sccwedge_drivers_tipping_climate

 **************

