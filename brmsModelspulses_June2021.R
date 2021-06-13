
library(dplyr)
library(ggplot2)
library(tidybayes)
library(Hmisc)
library(brms)


# Dena's example -----
titi.kinship.mat.in <-read.csv('titi.kinship.mat.csv')
row.names(titi.kinship.mat.in) <- colnames(titi.kinship.mat.in) # add rownames to kinship matrix
# check diagonal of kinship matrix
diag(as.matrix(titi.kinship.mat.in)) # these should all be 1's?

datos.pulse <- read.csv('pulse.features.addbodyweight.Oct2020.csv')

unique(colnames(titi.kinship.mat.in)[!(colnames(titi.kinship.mat.in) %in% datos.pulse$name)])
unique(datos.pulse$name[!(datos.pulse$name %in% colnames(titi.kinship.mat.in))])
# the names in the data and the kinship matrix don't line up

# subset both to a list of the in common names
A <- titi.kinship.mat.in
A2 <- A[row.names(A) %in% unique(datos.pulse$name),]

A2 <- A[row.names(A) %in% unique(datos.pulse$name), colnames(A) %in% unique(datos.pulse$name)]

datos.pulse <- datos.pulse[datos.pulse$name %in% row.names(A), ]

# check histograms to see if normally distributed
hist(datos.pulse$note.rate)
hist(datos.pulse$duration)
hist(datos.pulse$low.freq)
hist(datos.pulse$high.freq)
# these actually don't look too bad to me, not sure they necessarily need to be transformed with the exception
# of maybe high.freq

# Data summary plots
Noterate.violin.pulse <- ggpubr::ggviolin(data=datos.pulse,x='sex',y='note.rate',#scales='free',
                                          color=cbbPalette[1],fill=cbbPalette[1],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('Note Rate (Notes/s)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Note Rate (Notes/s)")

Duration.violin.pulse <- ggpubr::ggviolin(data=datos.pulse,x='sex',y='duration',#scales='free',
                                          color=cbbPalette[2],fill=cbbPalette[2],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('Duration (s)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Duration (s)")


High.freq.violin.pulse <- ggpubr::ggviolin(data=datos.pulse,x='sex',y='high.freq',#scales='free',
                                           color=cbbPalette[3],fill=cbbPalette[3],
                                           add = "boxplot", add.params = list(fill = "white"),
                                           panel.labs = list(feature = c('High Frequency (Hz)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("High Frequency (Hz)")

Low.freq.violin.pulse <- ggpubr::ggviolin(data=datos.pulse,x='sex',y='low.freq',#scales='free',
                                          color=cbbPalette[4],fill=cbbPalette[4],
                                          add = "boxplot", add.params = list(fill = "white"),
                                          panel.labs = list(feature = c('Low Frequency (Hz)')),
)+theme(legend.position = "none") +
  xlab('Sex')+ylab("Low Frequency (Hz)")


pulse.descriptive <- cowplot::plot_grid(Noterate.violin.pulse,Duration.violin.pulse,
                   High.freq.violin.pulse,Low.freq.violin.pulse,
                   labels = c('A',"B",'C',"D"),label_x = 0.9)


datos.pulse$high.freq <- log(datos.pulse$high.freq)
datos.pulse$low.freq <- log(datos.pulse$low.freq)
#datos.pulse[,c('duration','low.freq','high.freq','note.rate')] <- scale(log(datos.pulse[,c('duration','low.freq','high.freq','note.rate')]))

# Build full model for note rate
# NOTE: This model appears to converge and appears to be genetic signature in this feature
# Estimate: 0.43 +/- 0.18 
model_noterate.pulse_full <- brm(note.rate ~ weight_kg + sex + (1|gr(name, cov = A)) + (1|indiv_id), 
                                 data = datos.pulse, family = gaussian(), data2 = list(A = A),
                                 prior = c(prior(normal(0,10), "b"),    
                                           prior(normal(0,50), "Intercept"),
                                           prior(student_t(3,0,1), "sd"),
                                           prior(student_t(3,0,1), "sigma")),
                                 sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 2000, 
                                 control = list( adapt_delta = .99))

summary(model_noterate.pulse_full)
plot(model_noterate.pulse_full)

# ICC calculation
hyp_noterate.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_noterate.pulse <- hypothesis(model_noterate.pulse_full, hyp_noterate.pulse, class = NULL))

draws_noterate.pulse <- model_noterate.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_noterate.pulse, aes(sigma)) + geom_histogram(color = "black", fill = "white")

ggplot(draws_noterate.pulse, aes(sd_name__Intercept, sd_indiv_id__Intercept)) +
  geom_point()

draws_noterate.pulse <- model_noterate.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_noterate.pulse$ICCgenetic))
values <- draws_noterate.pulse$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_noterate.pulse$ICCinter))
values <- draws_noterate.pulse$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_noterate.pulse$ICCintra))
values <- draws_noterate.pulse$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.noterate <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

ggpubr::ggdensity(data=ICCvalues.noterate,x='values',fill='ICC')


# Build full model for duration.pulse
# NOTE: This model does NOT converge
model_duration.pulse_full <- 
  update(model_noterate.pulse_full, formula. = duration ~ ., newdata=datos.pulse)

summary(model_duration.pulse_full)
plot(model_duration.pulse_full)

hyp_duration.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_duration.pulse <- hypothesis(model_duration.pulse_full, hyp_duration.pulse, class = NULL))

draws_duration.pulse <- model_duration.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_duration.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")



# High frequency
model_high.freq.pulse_full <- 
  update(model_noterate.pulse_full, formula. = high.freq ~ ., newdata=datos.pulse)

summary(model_high.freq.pulse_full)
plot(model_high.freq.pulse_full)

brms::conditional_effects(model_high.freq.pulse_full)

hyp_high.freq.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_high.freq.pulse <- hypothesis(model_high.freq.pulse_full, hyp_high.freq.pulse, class = NULL))

draws_high.freq.pulse <- model_high.freq.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_high.freq.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")

# Low frequecy
model_low.freq.pulse_full <- 
  update(model_noterate.pulse_full, formula. = low.freq ~ ., newdata=datos.pulse)

summary(model_low.freq.pulse_full)
plot(model_low.freq.pulse_full)

brms::conditional_effects(model_low.freq.pulse_full)

hyp_low.freq.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_low.freq.pulse <- hypothesis(model_low.freq.pulse_full, hyp_low.freq.pulse, class = NULL))

draws_low.freq.pulse <- model_low.freq.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(Genetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Inter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(Intra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_low.freq.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")


# Density plots together
draws_noterate.pulse$Feature <- rep('Note Rate',nrow(draws_noterate.pulse))
draws_duration.pulse$Feature <- rep('Duration',nrow(draws_noterate.pulse))
draws_high.freq.pulse$Feature <- rep('High Frequency',nrow(draws_high.freq.pulse))
draws_low.freq.pulse$Feature <- rep('Low Frequency',nrow(draws_low.freq.pulse))

combinedpulses.df <- rbind.data.frame(draws_noterate.pulse,draws_duration.pulse,draws_high.freq.pulse,draws_low.freq.pulse)

data_long_pulse <- tidyr::gather(combinedpulses.df, ICC, value, Genetic:Intra, factor_key=TRUE)
data_long_pulse

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



PulsePlot <- ggpubr::ggdensity(data=data_long_pulse,x='value',fill='ICC',facet.by = 'Feature',
                  palette =colorBlindBlack8[1:3]  )+
  ylab('')+xlab('ICC')#+ theme(legend.position = "none")

PulsePlot <- PulsePlot+ guides(fill=guide_legend(title="Pulse ICC"))


PulsePlot

# Coefficient plots pulse


Pulse.plot <- sjPlot::plot_models(model_noterate.pulse_full,model_duration.pulse_full,
                                  model_high.freq.pulse_full,model_low.freq.pulse_full,
                                  axis.lim = c(-2,6),
                                  rm.terms='b_Intercept',colors = cbbPalette[1:4],
                                  axis.labels = c('Sex (M)', 'Weight (kg)'), legend.title = 'Pulse features',
                                  m.labels = c('Note Rate (Notes/s)','Duration (s)','High Frequency (Hz)','Low Frequency (Hz)'))


Pulse.plot 

# Pulse coefficient plots
# Note rate
noterate.pulse.modeldraws <- model_noterate.pulse_full %>%
  spread_draws(b_weight_kg,b_sexM) 

noterate.pulse.modeldraws$feature.category <- rep('noterate',nrow(noterate.pulse.modeldraws))

noterate.pulse.modeldraws_long <- gather(noterate.pulse.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

#Duration
duration.pulse.modeldraws <- model_duration.pulse_full %>%
  spread_draws(b_weight_kg,b_sexM) 

duration.pulse.modeldraws$feature.category <- rep('duration',nrow(duration.pulse.modeldraws))

duration.pulse.modeldraws_long <- gather(duration.pulse.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

# High frequency
high.freq.pulse.modeldraws <- model_high.freq.pulse_full %>%
  spread_draws(b_weight_kg,b_sexM) 

high.freq.pulse.modeldraws$feature.category <- rep('high.freq',nrow(high.freq.pulse.modeldraws))

high.freq.pulse.modeldraws_long <- gather(high.freq.pulse.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)

# Low frequency
low.freq.pulse.modeldraws <- model_low.freq.pulse_full %>%
  spread_draws(b_weight_kg,b_sexM) 

low.freq.pulse.modeldraws$feature.category <- rep('low.freq',nrow(low.freq.pulse.modeldraws))

low.freq.pulse.modeldraws_long <- gather(low.freq.pulse.modeldraws, feature, measurement, b_weight_kg:b_sexM, factor_key=TRUE)


allpulsemodel_draws <- rbind.data.frame(noterate.pulse.modeldraws_long,
                                        duration.pulse.modeldraws_long,
                                        high.freq.pulse.modeldraws_long,low.freq.pulse.modeldraws_long)

feature.names <- 
  c('noterate'='Note Rate (Notes/s)',
    'duration'='Duration (s)',
    'high.freq'='High Frequency (Hz)',
    'low.freq' ='Low Frequency (Hz)')

allpulsemodel_draws$feature <- dplyr::recode_factor(allpulsemodel_draws$feature, 
                                                    'b_sexM'='Sex (M)')

allpulsemodel_draws$feature <- dplyr::recode_factor(allpulsemodel_draws$feature, 
                                                    'b_weight_kg'='Weight (kg)')


allpulsemodel_draws %>%
  ggplot(aes( x = measurement, y=feature, fill=feature.category)) +
  scale_fill_manual(values=cbbPalette[1:4])+
  stat_halfeye(.width = c(.90, .5))+
  geom_vline(xintercept=0,linetype = "dashed")+
  facet_wrap(~feature.category,scales ='free',labeller = as_labeller(feature.names))+
  theme(legend.position = "none")+xlab('Estimates')+ylab('Predictor')


