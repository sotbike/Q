# Per minute blood loss in RPN/RAPN vs. OPN (Q)

Study title: Comparing robotic and open partial nephrectomy under the prism of surgical precision: A meta-analysis of the average blood loss rate as a novel variable.

Computational procedure: Q variable formulation, meta-analysis (MA), subgroup analysis (SGA) & meta-regression analysis (MRA), sensitivity analysis (SA) at five levels.

Included items: individual study data ("Datasets" folder), derived data ("Derived datasets" folder) & R code ("Code" folder).

Outcome: expected value & 95% confidence interval (CI95%) of the mean difference (MD) of Q when comparing robotic (RPN / RAPN) vs. open partial nephrectomy (OPN). RPN / RAPN is considered the experimental arm, while OPN the control arm.

Meta-analysis: R code is included for the implementation of a random effects model; subgroups: studies published before & after 2018, studies with & without patient matching, single- / multicenter studies, stratification of studies according to ROBINS-I tool class.

Meta-regression analysis: R code is included for the implementation of a random effects model; subgroups: studies published before & after 2018, studies with & without patient matching, single- / multicenter studies, stratification according to ROBINS-I tool, moderators: year of publication, number of quality stars assigned according to the Newcastle - Ottawa Scale (NOS).

Sensitivity analyses: Level 1: Re-analysis in a subset of studies with relatively increased accuracy of reported results based on their CI95% (OSA). Level 2: Re-analysis in the subset of low risk of bias (ROB) studies with patient matching (ML). Level 3: Re-analysis in a subset of studies with relatively large sample sizes based on their patient populations (LS). Level 4: Re-analysis in the subset of multicenter studies with patient matching (MM). Level 5: Exploration of the change in the MD of Q between robotic and open partial nephrectomy for the theoretical range of values obtained by the correlation coefficient (r) between the original variables (EBL & OT).

Instructions for the replication of results: 
   1. Create a desktop folder named "Data" and save the csv files provided in the "Dataseets" folder.
   2. Make sure you provide the appropriate absolute / relative paths for your system.
   3. The code creates various desktop folders with appropriate names to save the derived data and plots.
   4. Make sure to save the Excel files derived as csv files in the "Data" desktop folder.
