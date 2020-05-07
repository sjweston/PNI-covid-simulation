library(lavaan)
library(tidyverse)

set.seed(04102020) # April 10, 2020

# write general syntax for model describing outside influences on both 
# ses and suceptability to COVID-19

mod.syntax = "
covid ~ param1*ses
inflammation ~ param2*ses + param3*covid
"

# possible effect sizes
possible_r = c(.1, .5)

# possible sample sizes
possible_n = c(25, 75, 125, 200, 300)

#bootstrap
bootn = 200


# have to figure out how to code this
#hard coded version
all_combinations = expand_grid(
  param1 = possible_r,
  param2 = .3,
  param3 = possible_r,
  N = possible_n
)

total_sim = nrow(all_combinations)*bootn

#fixed value
fix_param1 = numeric(length = total_sim)
fix_param2 = numeric(length = total_sim)
fix_param3 = numeric(length = total_sim)
fix_N = numeric(length = total_sim)
#estimated relationships
est_ses_inflam_nocontrol = numeric(length = total_sim)
est_ses_inflam_control_covid = numeric(length = total_sim)
#pvalues
pval_ses_inflam_nocontrol = numeric(length = total_sim)
pval_ses_inflam_control_covid = numeric(length = total_sim)


sim_num = 0
for(i in 1:nrow(all_combinations)){
  
  #replace placeholders with specific parameter estimates
  specific.model = gsub("param1", all_combinations$param1[i], mod.syntax)
  specific.model = gsub("param2", all_combinations$param2[i], specific.model)
  specific.model = gsub("param3", all_combinations$param3[i], specific.model)
  
  for(k in 1:bootn){
    sim_num = sim_num+1
    fix_param1[sim_num] = all_combinations$param1[i]
    fix_param2[sim_num] = all_combinations$param2[i]
    fix_param3[sim_num] = all_combinations$param3[i]
    fix_N[sim_num] = all_combinations$N[i]
    simulated <- simulateData(specific.model, sample.nobs=all_combinations$N[i])
    #no covariates
    test_nocontrol = lm(inflammation~ses, data = simulated)
    est_ses_inflam_nocontrol[sim_num] = test_nocontrol$coefficients["ses"]
    pval_ses_inflam_nocontrol[sim_num] = coef(summary(test_nocontrol))["ses", "Pr(>|t|)"]
    #control covid
    test_covid = lm(inflammation~ses + covid, data = simulated)
    est_ses_inflam_control_covid[sim_num] = test_covid$coefficients["ses"]
    pval_ses_inflam_control_covid[sim_num] = coef(summary(test_covid))["ses", "Pr(>|t|)"]
  }
}

sim_results = data.frame(
  N = fix_N,
  popParam = fix_param2,
  ses2covid = fix_param1,
  covid2inflammation = fix_param3,
  est_nocov = est_ses_inflam_nocontrol,
  pval_nocov = pval_ses_inflam_nocontrol,
  est_covid = est_ses_inflam_control_covid,
  pval_covid = pval_ses_inflam_control_covid
)

save(sim_results, file = "sim_mediation.Rdata")
load("sim_mediation.Rdata")
sim_results = sim_results %>%
  mutate(ID = row_number(),
         indirect = ses2covid*covid2inflammation,
         total = indirect + popParam,
         ses2covid = factor(ses2covid, 
                               levels = c(.1, .5),
                               labels = c("SES weak \n cause of COVID-19", 
                                          "SES strong \n cause of COVID-19")),
         covid2inflammation = factor(covid2inflammation, 
                                  levels = c(.1,  .5),
                                  labels = c("COVID-19 weak \n cause of inflammation", 
                                             "COVID-19 strong \n cause of inflammation")))

man.colors = c("#1A85FF", "#D41159")

sim_results %>%
  gather("key", "value", which(grepl("_", names(.)))) %>%
  separate("key", into = c("stat", "control"), sep = "_") %>%
  spread("stat", "value") %>%
  ggplot(aes(x = as.factor(N), y = est, color = control, fill = control)) +
  geom_hline(aes(yintercept = indirect, linetype = "Indirect Effect"))+
  geom_hline(aes(yintercept = popParam, linetype = "Direct Effect"))+
  geom_hline(aes(yintercept = total, linetype = "Total Effect"))+
  geom_boxplot(alpha = .3) + 
  scale_fill_manual("Covariate(s)", 
                    values = man.colors,
                    labels = c("Controlling for COVID-19",
                               "Not controlling for COVID-19"))+
  labs(x = "Sample Size", 
       y = "Difference between estimate and parameter",
       linetype = "",
       title = "Simulation of regression estimates when COVID-19 is a mediating variable") +
  scale_color_manual(values = man.colors,
                     "Covariate(s)", 
                     labels = c("Controlling for COVID-19",
                                "Not controlling for COVID-19")) +
  facet_grid(ses2covid~covid2inflammation) +
  theme_minimal() + 
  theme(legend.position = "bottom", plot.title.position = "plot")

ggsave("simulation_mediation_estimate.pdf", width = 9, height = 5)

sim_results %>%
  mutate(ID = row_number()) %>%
  gather("key", "value", which(grepl("_", names(.)))) %>%
  separate("key", into = c("stat", "control"), sep = "_") %>%
  spread("stat", "value") %>%
  group_by(N, control, ses2covid, covid2inflammation) %>%
  mutate(sig = ifelse(pval < .05, 1,0)) %>%
  summarize(power = sum(sig)/n()) %>%
  filter(N < 300) %>%
  ggplot(aes(x = N, y = power, color = control)) +
  geom_line()+ 
  labs(x = "Sample Size", 
       y = "Power",
       color = "Covariate(s)") +
  facet_grid(ses2covid~covid2inflammation)
ggsave("simulation_confound_power.pdf", width = 8, height = 5)
