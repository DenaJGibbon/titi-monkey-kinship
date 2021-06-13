
library(dplyr)
library(ggplot2)
library(tidybayes)
library(Hmisc)
library(brms)
library(sjPlot)
library(ggpubr)
library(tidyr)


# Dena's example -----
titi.kinship.mat.in <- 
  read.csv('titi.kinship.mat.csv')

row.names(titi.kinship.mat.in) <- colnames(titi.kinship.mat.in) # add rownames to kinship matrix

# check diagonal of kinship matrix
diag(as.matrix(titi.kinship.mat.in)) # these should all be 1's?

datos.chirp <- read.csv('chirp.features.addbodyweight.Oct2020.csv')

unique(colnames(titi.kinship.mat.in)[!(colnames(titi.kinship.mat.in) %in% datos.chirp$name)])
unique(datos.chirp$name[!(datos.chirp$name %in% colnames(titi.kinship.mat.in))])
# the names in the data and the kinship matrix don't line up

# subset both to a list of the in common names
A <- titi.kinship.mat.in
A2 <- A[row.names(A) %in% unique(datos.chirp$name),]

A2 <- A[row.names(A) %in% unique(datos.chirp$name), colnames(A) %in% unique(datos.chirp$name)]

datos.chirp <- datos.chirp[datos.chirp$name %in% row.names(A), ]

# check histograms to see if normally distributed
hist(datos.chirp$note.rate)
hist(datos.chirp$duration)
hist(log(datos.chirp$low.freq))
hist((datos.chirp$high.freq))
# these actually don't look too bad to me, not sure they necessarily need to be transformed with the exception
# of maybe high.freq



ggpubr::ggscatter(data=datos.chirp,x='sex',y='note.rate',#scales='free',
                                          color=cbbPalette[1],fill=cbbPalette[1],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('Note Rate (Notes/s)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Note Rate (Notes/s)")

# Data summary plots
Noterate.violin.chirp <- ggpubr::ggviolin(data=datos.chirp,x='sex',y='note.rate',#scales='free',
                 color=cbbPalette[1],fill=cbbPalette[1],
                 add = "boxplot", add.params = list(fill = "white"),
                 panel.labs = list(feature = c('Note Rate (Notes/s)')),
                 )+theme(legend.position = "none") +
                  xlab('Sex')+ylab("Note Rate (Notes/s)")

Duration.violin.chirp <- ggpubr::ggviolin(data=datos.chirp,x='sex',y='duration',#scales='free',
                                          color=cbbPalette[2],fill=cbbPalette[2],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('Duration (s)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Duration (s)")


High.freq.violin.chirp <- ggpubr::ggviolin(data=datos.chirp,x='sex',y='high.freq',#scales='free',
                                          color=cbbPalette[3],fill=cbbPalette[3],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('High Frequency (Hz)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("High Frequency (Hz)")

Low.freq.violin.chirp <- ggpubr::ggviolin(data=datos.chirp,x='sex',y='low.freq',#scales='free',
                                           color=cbbPalette[4],fill=cbbPalette[4],
                                           add = "boxplot", add.params = list(fill = "white"),
                                           panel.labs = list(feature = c('Low Frequency (Hz)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Low Frequency (Hz)")


chirp.descriptive <- cowplot::plot_grid(Noterate.violin.chirp,Duration.violin.chirp,
                   High.freq.violin.chirp,Low.freq.violin.chirp,
                   labels = c('A',"B",'C',"D"),label_x = 0.9)

cowplot::plot_grid(pulse.descriptive,NULL,chirp.descriptive,
                   rel_heights = c(.4,.1,.4),
                   labels = c('Pulse','Chirp'),nrow=3,label_x = 0.9)


datos.chirp$high.freq <- log(datos.chirp$high.freq)
datos.chirp$low.freq <- log(datos.chirp$low.freq)

# Build full model for note rate
# NOTE: This model appears to converge and appears to be genetic signature in this feature
# Estimate: 0.43 +/- 0.18 
model_noterate.chirp_full <- brm(note.rate ~ weight_kg + sex + (1|gr(name, cov = A)) + (1|indiv_id), 
  data = datos.chirp, family = gaussian(), data2 = list(A = A),
    prior = c(prior(normal(0,10), "b"),    
              prior(normal(0,50), "Intercept"),
              prior(student_t(3,0,1), "sd"),
              prior(student_t(3,0,1), "sigma")),
  sample_prior = TRUE, chains = 4, cores = 4, iter = 4000,control = list( adapt_delta = .99))

summary(model_noterate.chirp_full)
plot(model_noterate.chirp_full)

brms::conditional_effects(model_noterate.chirp_full)
brms::fixef(model_noterate.chirp_full)

# ICC calculation
hyp_noterate.chirp <- paste("sd_name__Intercept^2 /", 
                            "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")

# sigma-squared is the within-animal
(hyp_noterate.chirp <- hypothesis(model_noterate.chirp_full, hyp_noterate.chirp, class = NULL))

draws_noterate.chirp <- model_noterate.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_noterate.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")

draws_noterate.chirp <- model_noterate.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_noterate.chirp$Genetic))
values <- draws_noterate.chirp$Genetic
Genetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_noterate.chirp$Inter))
values <- draws_noterate.chirp$Inter
Inter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_noterate.chirp$Intra))
values <- draws_noterate.chirp$Intra
Intra <- cbind.data.frame(ICC, values)

ICCvalues.noterate <- rbind.data.frame(Genetic,Inter,Intra)

ggpubr::ggdensity(data=ICCvalues.noterate,x='values',fill='ICC')

# Build full model for duration.chirp
# NOTE: This model does NOT converge
model_duration.chirp_full <- 
  update(model_noterate.chirp_full, formula. = duration ~ ., newdata=datos.chirp)

summary(model_duration.chirp_full)
plot(model_duration.chirp_full)
brms::conditional_effects(model_duration.chirp_full)

hyp_duration.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_duration.chirp <- hypothesis(model_duration.chirp_full, hyp_duration.chirp, class = NULL))

draws_duration.chirp <- model_duration.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_duration.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")

draws_duration.chirp <- model_duration.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_duration.chirp$Genetic))
values <- draws_duration.chirp$Genetic
Genetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_duration.chirp$Inter))
values <- draws_duration.chirp$Inter
Inter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_duration.chirp$Intra))
values <- draws_duration.chirp$Intra
Intra <- cbind.data.frame(ICC, values)

ICCvalues.duration <- rbind.data.frame(Genetic,Inter,Intra)

ggpubr::ggdensity(data=ICCvalues.duration,x='values',fill='ICC')


# High frequency
model_high.freq.chirp_full <- 
  update(model_noterate.chirp_full, formula. = high.freq ~ ., newdata=datos.chirp)

summary(model_high.freq.chirp_full)
plot(model_high.freq.chirp_full)

brms::conditional_effects(model_high.freq.chirp_full)

hyp_high.freq.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_high.freq.chirp <- hypothesis(model_high.freq.chirp_full, hyp_high.freq.chirp, class = NULL))

draws_high.freq.chirp <- model_high.freq.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_high.freq.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")


model_low.freq.chirp_full <- 
  update(model_noterate.chirp_full, formula. = low.freq ~ ., newdata=datos.chirp)

summary(model_low.freq.chirp_full)
plot(model_low.freq.chirp_full)

brms::conditional_effects(model_low.freq.chirp_full)


hyp_low.freq.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_low.freq.chirp <- hypothesis(model_low.freq.chirp_full, hyp_low.freq.chirp, class = NULL))

draws_low.freq.chirp <- model_low.freq.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_low.freq.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")

# High frequency
model_bandwidth.chirp_full <- 
  update(model_noterate.chirp_full, formula. = bandwidth.mean ~ ., newdata=datos.chirp)

summary(model_bandwidth.chirp_full)
plot(model_bandwidth.chirp_full)

brms::conditional_effects(model_bandwidth.chirp_full)

hyp_bandwidth.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_bandwidth.chirp <- hypothesis(model_bandwidth.chirp_full, hyp_bandwidth.chirp, class = NULL))

draws_bandwidth.chirp <- model_bandwidth.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_bandwidth.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")


model_low.freq.chirp_full <- 
  update(model_noterate.chirp_full, formula. = low.freq ~ ., newdata=datos.chirp)

summary(model_low.freq.chirp_full)
plot(model_low.freq.chirp_full)

brms::conditional_effects(model_low.freq.chirp_full)


hyp_low.freq.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_low.freq.chirp <- hypothesis(model_low.freq.chirp_full, hyp_low.freq.chirp, class = NULL))

draws_low.freq.chirp <- model_low.freq.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_low.freq.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")

# Density plots together
draws_noterate.chirp$Feature <- rep('Note Rate',nrow(draws_noterate.chirp))
draws_duration.chirp$Feature <- rep('Duration',nrow(draws_noterate.chirp))
draws_high.freq.chirp$Feature <- rep('High Frequency',nrow(draws_high.freq.chirp))
draws_low.freq.chirp$Feature <- rep('Low Frequency',nrow(draws_low.freq.chirp))

combinedchirps.df <- rbind.data.frame(draws_noterate.chirp,draws_duration.chirp,draws_high.freq.chirp,draws_low.freq.chirp)

data_long <- tidyr::gather(combinedchirps.df, ICC, value, Genetic:Intra, factor_key=TRUE)
data_long

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ChirpPlot <- ggpubr::ggdensity(data=data_long,x='value',fill='ICC',facet.by = 'Feature',
                  palette = colorBlindBlack8[1:3] )+
  ylab('')+xlab('ICC')#+ theme(legend.position = "none")

ChirpPlot <- ChirpPlot+ guides(fill=guide_legend(title="Chirp ICC"))

cowplot::plot_grid(ChirpPlot,PulsePlot,labels = c('A)','B)'),label_x = 0.95,
                   nrow=2)

# Coefficient plots pulse
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Chirp.plot <- sjPlot::plot_models(model_noterate.chirp_full,model_duration.chirp_full,
                                  model_high.freq.chirp_full,model_low.freq.chirp_full,
                                  axis.lim = c(-2,6),
                                  rm.terms='b_Intercept',
                                  colors = cbbPalette[1:4],
                                  axis.labels = c('Sex (M)', 'Weight (kg)'), legend.title = 'Chirp features',
                                  m.labels = c('Note Rate (Notes/s)','Duration (s)','High Frequency (Hz)','Low Frequency (Hz)'))


Chirp.plot 

cowplot::plot_grid(Pulse.plot,Chirp.plot,labels=c("A","B"),label_x = 0.9)

coefplot(model_note.rate.chirp_full)
sjPlot::tab_model(model_note.rate.chirp_full)

bayesplot::mcmc_areas(model_noterate.chirp_full,color='red',
                      pars = c("b_weight_kg","b_sexM"))

get_variables(model_noterate.chirp_full)

combined <- rbind(mcmc_intervals_data(model_noterate.chirp_full,pars =c("b_weight_kg","b_sexM")), 
                  mcmc_intervals_data(model_duration.chirp_full,pars =c("b_weight_kg","b_sexM")))
combined$model <- rep(c("Model 1", "Model 2"), each = 4/2)

# make the plot using ggplot 
library(ggplot2)
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined$model == "Model 2", 0, 0.1))
ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, size=2)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
  geom_point(position = pos, color="black")


# Note rate
noterate.chirp.modeldraws <- model_noterate.chirp_full %>%
  spread_draws(b_weight_kg,b_sexM) 

noterate.chirp.modeldraws$feature.category <- rep('noterate',nrow(noterate.chirp.modeldraws))

noterate.chirp.modeldraws_long <- gather(noterate.chirp.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

#Duration
duration.chirp.modeldraws <- model_duration.chirp_full %>%
  spread_draws(b_weight_kg,b_sexM) 

duration.chirp.modeldraws$feature.category <- rep('duration',nrow(duration.chirp.modeldraws))

duration.chirp.modeldraws_long <- gather(duration.chirp.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

# High frequency
high.freq.chirp.modeldraws <- model_high.freq.chirp_full %>%
  spread_draws(b_weight_kg,b_sexM) 

high.freq.chirp.modeldraws$feature.category <- rep('high.freq',nrow(high.freq.chirp.modeldraws))

high.freq.chirp.modeldraws_long <- gather(high.freq.chirp.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

# Low frequency
low.freq.chirp.modeldraws <- model_low.freq.chirp_full %>%
  spread_draws(b_weight_kg,b_sexM) 

low.freq.chirp.modeldraws$feature.category <- rep('low.freq',nrow(low.freq.chirp.modeldraws))

low.freq.chirp.modeldraws_long <- gather(low.freq.chirp.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)


allchirpmodel_draws <- rbind.data.frame(noterate.chirp.modeldraws_long,
                 duration.chirp.modeldraws_long,
                 high.freq.chirp.modeldraws_long,low.freq.chirp.modeldraws_long)

feature.names <- 
  c('noterate'='Note Rate (Notes/s)',
    'duration'='Duration (s)',
    'high.freq'='High Frequency (Hz)',
    'low.freq' ='Low Frequency (Hz)')

allchirpmodel_draws$feature <- dplyr::recode_factor(allchirpmodel_draws$feature, 
                                                    'b_sexM'='Sex (M)')

allchirpmodel_draws$feature <- dplyr::recode_factor(allchirpmodel_draws$feature, 
                                                    'b_weight_kg'='Weight (kg)')

                     
chirp.coef.plot <- allchirpmodel_draws %>%
  ggplot(aes( x = measurement, y=feature, fill=feature.category)) +
  scale_fill_manual(values=cbbPalette[1:4])+
  stat_halfeye(.width = c(.90, .5), alpha = 0.65)+
  geom_vline(xintercept=0,linetype = "dashed")+
  facet_wrap(~feature.category,scales ='free',labeller = as_labeller(feature.names))+
  theme(legend.position = "none")+xlab('Estimates')+ylab('Predictor')

cowplot::plot_grid(NULL,pulse.coef.plot,NULL,chirp.coef.plot,nrow=4,rel_heights = c(0.05,.4,0.05,.4),
                   label_x = 0,
                   label_y = 1.1,
                   labels = c('','A. Pulse coeffcients','','B. Chirp coeffcients'))




