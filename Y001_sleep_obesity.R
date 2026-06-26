
library(starUKB)
library(starr)
library(starx)
library(tidyverse)
library(survival)
library(rms)
library(tableone)
library(mice)
library(mediation)
library(broom)
library(writexl)
library(readxl)
library(ggplot2)

set.seed(20250611)
out_dir = "outputs"
derived_dir = "derived"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(derived_dir,recursive=TRUE,showWarnings=FALSE)
source("ukb_init.R")
source("my_func.R")

p_fmt = function(p) ifelse(is.na(p),NA_character_,ifelse(p<0.001,"<0.001",sprintf("%.3f",p)))
ci_fmt = function(est,lo,hi) sprintf("%.3f (%.3f-%.3f)",est,lo,hi)
exposure_vars = c("sle_group5","sm","zyx1","dh1","btdd1")
exposure_labels = c(sle_group5="Sleep duration",sm="Insomnia",zyx1="Chronotype",dh1="Snoring",btdd1="Daytime sleepiness")
term_labels = c("sle_group51"="Sleep duration <=6 h","sle_group52"="Sleep duration >=10 h","sm2"="Insomnia sometimes","sm3"="Insomnia usually","zyx11"="Evening chronotype","dh11"="Snoring yes","btdd11"="Daytime sleepiness sometimes","btdd12"="Daytime sleepiness often/all of the time")
cov_m1 = c("age","eth")
cov_m2 = c(cov_m1,"edu","alc","smk","met","tdi","dm","hp","food_score")
cov_m3 = c(cov_m2,"num.bb","cc","lcs","byy","hrt")
cov_m4 = c(cov_m3,"bmi")
primary_vars = c(exposure_vars,cov_m3,"td_obr.x","d_obr.x","bmi")

cox_tidy = function(fit,exposure,analysis){tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(str_detect(term,paste0("^",exposure))) %>% mutate(analysis=analysis,exposure=exposure,n=fit$n,events=fit$nevent,term_label=recode(term,!!!term_labels,.default=term),HR_95CI=ci_fmt(estimate,conf.low,conf.high),p_value_fmt=p_fmt(p.value),.before=1)}
run_cox = function(dat,exposure,covars,label,timevar="td_obr.x",eventvar="d_obr.x"){coxph(reformulate(c(exposure,covars),response=paste0("Surv(",timevar,",",eventvar,")")),data=dat) %>% cox_tidy(exposure,label)}
run_cox_set = function(dat,covars,label,timevar="td_obr.x",eventvar="d_obr.x") map_dfr(exposure_vars,~run_cox(dat,.x,covars,label,timevar,eventvar))

# 1. UK Biobank field extraction and variable construction

dislist = tribble(~icd10,~name,~cn,~icd9,~self,"E66[0|8|9]","obr","Obesity",NA,NA,"I10","hp2","Hypertension",NA,NA,"E1[0-4]","dm2","Diabetes",NA,NA)
datalist1 = tribble(
  ~factor,~id_column,~name,~cn,~from,~to,~nato,~id2,
  F,50,"hei_cm","Height",NULL,NULL,NA,0,F,48,"wc_cm","Waist circumference",NULL,NULL,NA,0,F,49,"hc_cm","Hip circumference",NULL,NULL,NA,0,F,21002,"wei","Weight",NULL,NULL,NA,0,
  T,2443,"dm1","Diabetes history",c(-3,-1),c(NA,NA),NA,0,T,2724,"pau","Menopause",c(-3,2,3),c(NA,NA,NA),NA,0,F,2734,"num.bb","Number of live births",NULL,NULL,NA,0,F,2714,"cc","Age at menarche",NULL,NULL,NA,0,
  T,2774,"lcs","Pregnancy loss history",c(-3,-1),c(NA,NA),NA,0,T,2784,"byy","Oral contraceptive use",c(-3,-1),c(NA,NA),NA,0,T,2814,"hrt","HRT use",c(-3,-1),c(NA,NA),NA,0,
  T,3591,"op.ut","Hysterectomy",c(-3,-5),c(NA,NA),NA,0,T,2834,"op.ov","Oophorectomy",c(-3,-5),c(NA,NA),NA,0,
  F,1160,"time","Sleep duration",NULL,NULL,NA,0,T,1180,"zyx","Chronotype",c(-3,-1),c(NA,NA),NA,0,T,1200,"sm","Insomnia",c(-3),c(NA),NA,0,T,1210,"dh","Snoring",c(-3,-1),c(NA,NA),NA,0,T,1220,"btdd","Daytime sleepiness",c(-3,-1),c(NA,NA),NA,0,
  F,1289,"cooked.veg","Cooked vegetables",c(-3,-1,-10),c(NA,NA,NA),NA,0,F,1299,"salad","Raw vegetables/salad",c(-3,-1,-10),c(NA,NA,NA),NA,0,F,1309,"fruit","Fresh fruit",c(-3,-1,-10),c(NA,NA,NA),NA,0,
  T,1329,"oil.fish","Oily fish",c(-3,-1),c(NA,NA),NA,0,T,1339,"fish","Non-oily fish",c(-3,-1),c(NA,NA),NA,0,T,1349,"pro.meet","Processed meat",c(-3,-1),c(NA,NA),NA,0,
  T,1369,"beef","Beef",c(-3,-1),c(NA,NA),NA,0,T,1379,"sheep","Lamb/mutton",c(-3,-1),c(NA,NA),NA,0,T,1389,"pig","Pork",c(-3,-1),c(NA,NA),NA,0,
  F,30710,"s_crp","CRP",NULL,NULL,NA,0,F,30000,"s_num.wbc","White blood cell count",NULL,NULL,NA,0,F,30140,"s_num.ne","Neutrophil count",NULL,NULL,NA,0,F,30120,"s_num.lyn","Lymphocyte count",NULL,NULL,NA,0,F,30130,"s_num.mono","Monocyte count",NULL,NULL,NA,0,
  F,30080,"s_plt","Platelet count",NULL,NULL,NA,0,F,30070,"s_rdw","RDW",NULL,NULL,NA,0,F,30100,"s_mpv","MPV",NULL,NULL,NA,0,F,30090,"s_pct","Platelet crit",NULL,NULL,NA,0,F,30110,"s_pdw","PDW",NULL,NULL,NA,0,
  F,30870,"s_tri","Triglycerides",NULL,NULL,NA,0,F,30690,"s_tc","Total cholesterol",NULL,NULL,NA,0,F,30780,"s_ldl","LDL cholesterol",NULL,NULL,NA,0,F,30760,"s_hdl","HDL cholesterol",NULL,NULL,NA,0,
  F,30630,"s_zza","Apolipoprotein A",NULL,NULL,NA,0,F,30640,"s_zzb","Apolipoprotein B",NULL,NULL,NA,0,F,30790,"s_zdba","Lipoprotein A",NULL,NULL,NA,0,F,30750,"s_hba1c","HbA1c",NULL,NULL,NA,0,F,30740,"s_glu","Glucose",NULL,NULL,NA,0,F,30880,"s_sua","Urate",NULL,NULL,NA,0,F,30800,"s_eatra","Oestradiol",NULL,NULL,NA,0)
datalist2 = tribble(~factor,~id_column,~name,~cn,~from,~to,~nato,~id2,~id3,T,6150,"hp1","Hypertension",c(-7,-3,1,2,3,4),c(NA,NA,0,0,0,1),NA,0,c(0,1,2,3),T,6153,"med1","Medication 1",c(-7,-3,-1,1,2,3,4,5),c(NA,NA,NA,1,2,3,4,5),NA,0,c(0,1,2),T,6177,"med2","Medication 2",c(-7,-3,-1,1,2,3),c(NA,NA,NA,1,2,3),NA,0,c(0,1,2,3))

ddd = UKB_DF[1] %>% ukb_add_cov() %>% ukb_add_icd10(dislist) %>% ukb_add_commondrug(c("dr_atdba","dr_atht","dr_atlp")) %>% ukb_add_data(datalist1) %>% ukb_add_data(datalist2)
ddd = ddd %>% mutate(hp1=as.factor(if_else(rowSums(across(starts_with("hp1_0."),~replace_na(as.numeric(as.character(.x))>=1,FALSE)))>0,1,0)),hp=as.factor(if_else(hp1==1|dr_atht==1,1,0)),dm=as.factor(if_else(dm1==1|dr_atdba==1,1,0)),whr=wc_cm/hc_cm,zyx1=as.factor(cge(zyx,org=1:4,swi=c(0,0,1,1))),dh1=as.factor(cge(dh,org=c(1,2),swi=c(1,0))),btdd1=as.factor(case_when(btdd==0~0,btdd==1~1,btdd %in% c(2,3)~2,TRUE~NA_real_)),sle_group5=as.factor(case_when(time>=7&time<=9~0,time<=6~1,time>=10~2,TRUE~NA_real_)),eth2=as.factor(case_when(eth==1~1,eth %in% 2:6~2,TRUE~NA_real_)),tdi_group=as.factor(ntile(tdi,3)))
r4 = cge(ddd$oil.fish,org=0:5,swi=c(0,0.5,1,3,5.5,7)); r5 = cge(ddd$fish,org=0:5,swi=c(0,0.5,1,3,5.5,7)); r1 = cge(ddd$beef,org=0:5,swi=c(0,0.5,1,3,5.5,7)); r2 = cge(ddd$sheep,org=0:5,swi=c(0,0.5,1,3,5.5,7)); r3 = cge(ddd$pig,org=0:5,swi=c(0,0.5,1,3,5.5,7))
ddd = ddd %>% mutate(fish_score=if_else(r4+r5>=2,1,0),fruit_score=if_else(fruit>=3,1,0),veg_score=if_else(cooked.veg>=4|salad>=4|cooked.veg+salad>=4,1,0),pro.meet_score=if_else(pro.meet %in% 0:2,1,0),red.meet_score=if_else(r1+r2+r3<2,1,0),food_score=rowSums(across(c(fish_score,fruit_score,veg_score,pro.meet_score,red.meet_score)),na.rm=FALSE),food_mod=as.factor(case_when(food_score %in% 4:5~0,food_score %in% 2:3~1,food_score %in% 0:1~2,TRUE~NA_real_)),s_nlr=s_num.ne/s_num.lyn,s_plr=s_plt/s_num.lyn,s_lmr=s_num.lyn/s_num.mono,s_mpvlr=s_mpv/s_num.lyn,s_sii=s_plt*s_num.ne/s_num.lyn)

source_pool = ddd %>% filter(sex==0,age>=40,pau==1,op.ov==0,op.ut==0,between(time,3,20))
dd6 = source_pool %>% drop_na(all_of(c(exposure_vars,cov_m3))) %>% filter(bmi>=18.5,bmi<30)
cohort_flow = tibble(step=c("Total UK Biobank participants","Women","Age >=40 years","Postmenopausal or clear menopause status","No ovarian/uterine resection","Plausible sleep duration","Complete sleep and covariates","Baseline BMI 18.5-<30 kg/m2"),n=c(nrow(ddd),nrow(filter(ddd,sex==0)),nrow(filter(ddd,sex==0,age>=40)),nrow(filter(ddd,sex==0,age>=40,pau==1)),nrow(filter(ddd,sex==0,age>=40,pau==1,op.ov==0,op.ut==0)),nrow(source_pool),nrow(source_pool %>% drop_na(all_of(c(exposure_vars,cov_m3)))),nrow(dd6))) %>% mutate(excluded=lag(n)-n)
followup = dd6 %>% filter(td_obr.x>=0) %>% summarise(n=n(),events=sum(d_obr.x==1,na.rm=TRUE),median_followup_years=survfit(Surv(td_obr.x,1-d_obr.x)~1,data=.) %>% summary() %>% pluck("table") %>% .["median"]/365.25,mean_followup_years=mean(td_obr.x,na.rm=TRUE)/365.25)
write_csv(cohort_flow,file.path(out_dir,"cohort_flow.csv")); write_csv(followup,file.path(out_dir,"followup_summary.csv")); saveRDS(dd6,file.path(derived_dir,"analysis_cohort.rds"))

# 2. Baseline table, main Cox models and RCS

table1_vars = c("age","eth","alc","smk","met","tdi","edu","dm","hp","bmi","whr","food_score","food_mod","num.bb","cc","lcs","byy","hrt","time","sle_group5","sm","zyx1","dh1","btdd1")
table1 = CreateTableOne(vars=table1_vars,strata="d_obr.x",data=dd6,includeNA=TRUE,test=FALSE) %>% print(smd=TRUE,quote=FALSE,noSpaces=TRUE,printToggle=FALSE) %>% as.data.frame() %>% rownames_to_column("Characteristic")
main_cox = bind_rows(run_cox_set(dd6,cov_m1,"Model 1"),run_cox_set(dd6,cov_m2,"Model 2"),run_cox_set(dd6,cov_m3,"Model 3")) %>% mutate(sleep_behaviour=recode(exposure,!!!exposure_labels)) %>% arrange(sleep_behaviour,term,analysis)
write_xlsx(list(Table1_baseline_characteristics_SMD=table1,Table2_main_Cox=main_cox),file.path(out_dir,"Tables_1_2_primary_results.xlsx"))
rcs_dat = dd6 %>% drop_na(time,d_obr.x,td_obr.x,all_of(cov_m3)); dd = datadist(rcs_dat); options(datadist="dd")
rcs_fit = cph(as.formula(paste("Surv(td_obr.x,d_obr.x) ~ rcs(time,4) +",paste(cov_m3,collapse="+"))),data=rcs_dat,x=TRUE,y=TRUE,surv=TRUE)
rcs_pred = Predict(rcs_fit,time,ref.zero=TRUE,fun=exp) %>% as_tibble(); rcs_anova = anova(rcs_fit) %>% as.data.frame() %>% rownames_to_column("term")
write_xlsx(list(rcs_prediction=rcs_pred,rcs_anova=rcs_anova),file.path(out_dir,"Figure_1_RCS_source_data.xlsx"))

# 3. Sensitivity analyses, missing data, PH diagnostics and interaction tests

sensitivity_2yr = run_cox_set(dd6 %>% filter(td_obr.x>=730),cov_m3,"Exclude first 2 years")
sensitivity_trim_sleep = run_cox_set(dd6 %>% filter(time<12),cov_m3,"Trim extreme sleep duration")
sensitivity_whr = run_cox_set(dd6,c(cov_m3,"whr"),"Additional WHR")
sensitivity_bmi = run_cox_set(dd6,cov_m4,"Additional baseline BMI")
write_xlsx(list(Table_S3_exclude_first_2_years=sensitivity_2yr,Table_S4_trimmed_sleep=sensitivity_trim_sleep,Table_S5_additional_WHR=sensitivity_whr,Table_S23_additional_BMI=sensitivity_bmi),file.path(out_dir,"Sensitivity_analyses.xlsx"))
missing_summary = source_pool %>% summarise(across(all_of(primary_vars),~sum(is.na(.x)))) %>% pivot_longer(everything(),names_to="variable",values_to="missing_n") %>% mutate(n_before=nrow(source_pool),missing_pct=100*missing_n/n_before,non_missing_n=n_before-missing_n)
complete_flag = source_pool %>% mutate(complete_case=if_else(if_all(all_of(primary_vars),~!is.na(.x)),"Included complete-case","Excluded due to missingness"))
missing_comparison = CreateTableOne(vars=c("age","eth","edu","tdi","smk","alc","met","food_score","dm","hp","bmi","time","sle_group5","sm","zyx1","dh1","btdd1","num.bb","cc","lcs","byy","hrt"),strata="complete_case",data=complete_flag,includeNA=TRUE,test=FALSE) %>% print(smd=TRUE,quote=FALSE,noSpaces=TRUE,printToggle=FALSE) %>% as.data.frame() %>% rownames_to_column("Characteristic")
write_xlsx(list(Table_S16_missingness=missing_summary,Table_S17_complete_case_comparison=missing_comparison),file.path(out_dir,"Tables_S16_S17_missing_data.xlsx"))
imp_dat = source_pool %>% select(all_of(primary_vars)) %>% mutate(across(where(is.factor),as.factor)); ini = mice(imp_dat,maxit=0,printFlag=FALSE)
imp = mice(imp_dat,m=20,maxit=10,method=ini$method,predictorMatrix=ini$predictorMatrix,seed=20250611,printFlag=FALSE)
mi_results = map_dfr(1:20,~complete(imp,.x) %>% mutate(across(all_of(exposure_vars),as.factor),d_obr.x=as.numeric(d_obr.x)) %>% run_cox_set(cov_m3,paste0("MI dataset ",.x)))
write_xlsx(list(Table_S18_MI_model_outputs=mi_results),file.path(out_dir,"Table_S18_multiple_imputation_sensitivity.xlsx"))
ph_results = map_dfr(exposure_vars,function(x){fit = coxph(reformulate(c(x,cov_m3),response="Surv(td_obr.x,d_obr.x)"),data=dd6); ph = cox.zph(fit); as.data.frame(ph$table) %>% rownames_to_column("term") %>% mutate(exposure=exposure_labels[x],n=fit$n,events=fit$nevent,.before=1)})
write_xlsx(list(Table_S21_PH_assumption=ph_results),file.path(out_dir,"Table_S21_PH_assumption_tests.xlsx"))
dd_int = dd6 %>% mutate(age_group=as.factor(if_else(age<60,0,1)),bmi_group=as.factor(if_else(bmi<25,0,1)),tdi_group=as.factor(ntile(tdi,3)),eth2=as.factor(eth2))
subgroups = tribble(~group,~subgroup,~remove_cov,"age_group","Age group",list("age"),"bmi_group","Baseline BMI group",list("bmi"),"hrt","HRT use",list("hrt"),"eth2","Ethnicity",list("eth"),"dm","Diabetes history",list("dm"),"hp","Hypertension history",list("hp"),"tdi_group","TDI tertile",list("tdi"))
interaction_res = crossing(exposure=exposure_vars,subgroups) %>% pmap_dfr(function(exposure,group,subgroup,remove_cov){covars = cov_m3[!cov_m3 %in% unlist(remove_cov)]; fit0 = coxph(reformulate(c(exposure,group,covars),response="Surv(td_obr.x,d_obr.x)"),data=dd_int); fit1 = coxph(reformulate(c(paste0(exposure,"*",group),covars),response="Surv(td_obr.x,d_obr.x)"),data=dd_int); lr = anova(fit0,fit1,test="LRT"); tibble(sleep_behaviour=exposure_labels[exposure],subgroup_variable=subgroup,n=fit1$n,events=fit1$nevent,df_for_interaction=lr$Df[2],chisq_for_interaction=lr$Chisq[2],p_for_interaction=lr$`Pr(>|Chi|)`[2],p_for_interaction_fmt=p_fmt(p_for_interaction))})
write_xlsx(list(Table_S20_interaction_tests=interaction_res),file.path(out_dir,"Table_S20_formal_interaction_tests.xlsx"))
repeat_bmi = tribble(~factor,~id_column,~name,~cn,~from,~to,~nato,~id2,F,21001,"bmi_repeat","Repeat measured BMI",NULL,NULL,NA,1)
repeat_df = UKB_DF[1] %>% ukb_add_data(repeat_bmi) %>% mutate(eid=as.character(eid)); dd_repeat = dd6 %>% mutate(eid=as.character(eid)) %>% left_join(repeat_df,by="eid") %>% filter(!is.na(bmi_repeat)) %>% mutate(d_bmi_rep=as.integer(bmi_repeat>=30))
repeat_summary = dd_repeat %>% summarise(n_repeat=n(),events=sum(d_bmi_rep==1,na.rm=TRUE),non_events=sum(d_bmi_rep==0,na.rm=TRUE))
repeat_logistic = map_dfr(exposure_vars,function(x){fit = glm(reformulate(c(x,cov_m3),response="d_bmi_rep"),data=dd_repeat,family=binomial()); tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(str_detect(term,paste0("^",x))) %>% mutate(sleep_behaviour=exposure_labels[x],exposure=x,n=nobs(fit),events=sum(dd_repeat$d_bmi_rep==1,na.rm=TRUE),OR_95CI=ci_fmt(estimate,conf.low,conf.high),p_value_fmt=p_fmt(p.value),.before=1)})
write_xlsx(list(summary=repeat_summary,logistic_results=repeat_logistic),file.path(out_dir,"Table_S19_repeat_BMI_logistic_sensitivity.xlsx"))

# 4. Exploratory indirect-association analyses

dd6 = dd6 %>% mutate(sle_short=as.factor(case_when(time>=7&time<=9~0,time<=6~1,TRUE~NA_real_)),sle_long=as.factor(case_when(time>=7&time<=9~0,time>=10~1,TRUE~NA_real_)),insomnia_usual=as.factor(case_when(sm==1~0,sm==3~1,TRUE~NA_real_)))
met_df = UKB_DF[1] %>% ukb_add_met() %>% mutate(eid=as.character(eid)); met_names = names(met_df)[names(met_df)!="eid"] %>% head(249)
dd_met = dd6 %>% mutate(eid=as.character(eid)) %>% left_join(select(met_df,eid,all_of(met_names)),by="eid")
serum_mediators = c("s_crp","s_num.wbc","s_num.ne","s_num.lyn","s_num.mono","s_plt","s_rdw","s_mpv","s_pct","s_pdw","s_tri","s_tc","s_ldl","s_hdl","s_zza","s_zzb","s_zdba","s_hba1c","s_glu","s_sua","s_nlr","s_plr","s_lmr","s_mpvlr","s_sii")
indirect_exposures = tribble(~exposure,~label,"sle_short","Short sleep duration","sle_long","Long sleep duration","insomnia_usual","Usual insomnia","zyx1","Evening chronotype","dh1","Snoring")
fit_indirect = function(dat,exposure,mediator,label){dt = dat %>% select(all_of(c(exposure,"d_obr.x",mediator,cov_m3))) %>% drop_na(); model_m = lm(reformulate(c(exposure,cov_m3),response=mediator),data=dt); model_y = glm(reformulate(c(exposure,mediator,cov_m3),response="d_obr.x"),data=dt,family=binomial()); res = mediate(model_m,model_y,treat=exposure,mediator=mediator,boot=TRUE,sims=1000); tibble(exposure=label,mediator=mediator,n=nrow(dt),events=sum(dt$d_obr.x==1),acme=res$d0,acme_low=res$d0.ci[1],acme_high=res$d0.ci[2],p_acme=res$d0.p,ade=res$z0,ade_low=res$z0.ci[1],ade_high=res$z0.ci[2],p_ade=res$z0.p,total=res$tau.coef,total_low=res$tau.ci[1],total_high=res$tau.ci[2],p_total=res$tau.p,proportion=res$n0)}
serum_results = indirect_exposures %>% pmap_dfr(function(exposure,label) map_dfr(serum_mediators,~fit_indirect(dd6,exposure,.x,label))) %>% group_by(exposure) %>% mutate(fdr_acme=p.adjust(p_acme,"fdr")) %>% ungroup()
met_results = indirect_exposures %>% pmap_dfr(function(exposure,label) map_dfr(met_names,~fit_indirect(dd_met,exposure,.x,label))) %>% group_by(exposure) %>% mutate(fdr_acme=p.adjust(p_acme,"fdr")) %>% ungroup()
met_sample_size = met_results %>% group_by(exposure) %>% summarise(number_of_metabolites_assessed=n_distinct(mediator),participants_min=min(n,na.rm=TRUE),participants_median=median(n,na.rm=TRUE),participants_max=max(n,na.rm=TRUE),events_min=min(events,na.rm=TRUE),events_median=median(events,na.rm=TRUE),events_max=max(events,na.rm=TRUE),.groups="drop")
write_xlsx(list(Table_S13_blood_biomarkers=serum_results,Table_S14_plasma_metabolites=met_results,Table_S24_metabolomics_sample_size=met_sample_size),file.path(out_dir,"Tables_S13_S14_S24_indirect_association.xlsx"))

# 5. Figures and computational environment

fig1 = ggplot(rcs_pred,aes(x=time,y=yhat)) + geom_hline(yintercept=1,linetype="dashed",linewidth=0.4) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.20) + geom_line(linewidth=0.8) + labs(x="Sleep duration, h/day",y="Hazard ratio for ICD-coded incident obesity",caption="Restricted cubic spline with 7 h/day as reference; P-overall = 0.002 and P-nonlinearity = 0.024.") + theme_classic(base_size=11)
ggsave(file.path(out_dir,"Figure_1_RCS_sleep_duration.png"),fig1,width=6.5,height=4.5,dpi=300,bg="white")
fig2_data = bind_rows(serum_results %>% mutate(type="Blood biomarkers"),met_results %>% mutate(type="Plasma metabolites")) %>% filter(p_acme<0.05) %>% mutate(mediator_label=str_replace_all(mediator,"_"," ")) %>% group_by(exposure,type) %>% slice_min(order_by=p_acme,n=20,with_ties=FALSE) %>% ungroup()
fig2 = ggplot(fig2_data,aes(x=acme,y=reorder(mediator_label,acme))) + geom_vline(xintercept=0,linetype="dashed",linewidth=0.3) + geom_errorbar(aes(xmin=acme_low,xmax=acme_high),width=0.2) + geom_point(size=1.8) + facet_grid(type~exposure,scales="free_y",space="free_y") + labs(x="Average indirect effect with 95% CI",y=NULL,caption="Exploratory indirect-association signals; findings should not be interpreted as causal mediation.") + theme_classic(base_size=9) + theme(strip.text=element_text(face="bold"),axis.text.y=element_text(size=7))
ggsave(file.path(out_dir,"Figure_2_exploratory_indirect_association_signals.png"),fig2,width=10,height=8,dpi=300,bg="white")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
packages = c("R","starUKB","starr","starx","tidyverse","survival","rms","mice","mediation","broom","writexl","tableone")
versions = tibble(component=packages,version=map_chr(packages,~if_else(.x=="R",as.character(getRversion()),as.character(packageVersion(.x)))))
write_xlsx(list(Table_S22_computational_environment=versions),file.path(out_dir,"Table_S22_computational_environment.xlsx"))
