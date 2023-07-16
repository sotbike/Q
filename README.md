# Q

Study title: A quantified comparison between robotic and open partial nephrectomy regarding the level of surgical precision through the formulation of a novel variable. A meta-analysis of comparative studies accompanied by sensitivity analysis.

Computational procedure: Q variable formulation, meta-analysis (MA), subgroup analysis (SGA) & meta-regression analysis (MRA), sensitivity analysis (SA) at two levels.

Included items: individual study data, sensitivity analysis data, derived data & R code.

Outcome: expected value & 95% confidence interval (CI95%) of the mean difference (MD) of Q when comparing robotic (RPN / RAPN) vs. open partial nephrectomy (OPN). RPN / RAPN is considered the experimental arm, while OPN the control arm.

Meta-analysis: R code is included for the implementation of a random effects model; subgroups: studies published before & after 2018, studies with & without patient matching, single- / multicenter studies, stratification of studies according to ROBINS-I tool class.

Meta-regression analysis: R code is included for the implementation of a random effects model; subgroups: studies published before & after 2018, studies with & without patient matching, single- / multicenter studies, stratification according to ROBINS-I tool, moderators: year of publication, number of quality stars assigned according to the Newcastle - Ottawa Scale (NOS).

Sensitivity analysis: Level 1: Re-analysis in a subset of studies with relatively increased precison of reported results. Level 2: Observation of the change in the MD of Q between robotic and open partial nephrectomy for the various values obtained by the correlation coefficient (r) between the original variables (EBL & OT).

Instructions for the replication of results: 
   1. Create a desktop folder named "Data" and save the csv files provided.
   2. Make sure you provide the appropriate absolute / relative paths for your system.
   3. The code creates various desktop folders with appropriate names to save the derived data and plots.
