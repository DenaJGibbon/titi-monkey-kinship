
library(dplyr)
library(ggplot2)
library(tidybayes)
library(Hmisc)
library(brms)


# Dena's example -----
titi.kinship.mat.in <-read.csv('/Users/denasmacbook/titii_acoustic_distance/titi.kinship.mat.csv')
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
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

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
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

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
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

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
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_low.freq.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")


# PCA of all features
pulse.pca <- prcomp(datos.pulse[,3:11])
summary(pulse.pca)
pulse.pca$x[,1:2]

pca1 <- pulse.pca$x[,1]
pca2 <- pulse.pca$x[,2]
datos.pulse.pca <- cbind.data.frame(datos.pulse,pca1,pca2)

model_pca1.pulse_full <- 
  update(model_noterate.pulse_full, formula. = pca1 ~ ., newdata=datos.pulse.pca)

summary(model_pca1.pulse_full)
plot(model_pca1.pulse_full)

brms::conditional_effects(model_pca1.pulse_full)

hyp_pca1.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_pca1.pulse <- hypothesis(model_pca1.pulse_full, hyp_pca1.pulse, class = NULL))

draws_pca1.pulse <- model_pca1.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_pca1.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")


# PCA of all features

model_pca2.pulse_full <- 
  update(model_noterate.pulse_full, formula. = pca2 ~ ., newdata=datos.pulse.pca)

summary(model_pca2.pulse_full)
plot(model_pca2.pulse_full)

brms::conditional_effects(model_pca2.pulse_full)

hyp_pca2.pulse <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_pca2.pulse <- hypothesis(model_pca2.pulse_full, hyp_pca2.pulse, class = NULL))

draws_pca2.pulse <- model_pca2.pulse_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_pca2.pulse, aes(ICC)) + geom_histogram(color = "black", fill = "white")


# Density plots together
draws_noterate.pulse$Feature <- rep('Note Rate',nrow(draws_noterate.pulse))
draws_duration.pulse$Feature <- rep('Duration',nrow(draws_noterate.pulse))
draws_high.freq.pulse$Feature <- rep('High Frequency',nrow(draws_high.freq.pulse))
draws_low.freq.pulse$Feature <- rep('Low Frequency',nrow(draws_low.freq.pulse))

combinedpulses.df <- rbind.data.frame(draws_noterate.pulse,draws_duration.pulse,draws_high.freq.pulse,draws_low.freq.pulse)

ggpubr::ggdensity(data=combinedpulses.df,x='ICC',fill='Feature',facet.by = 'Feature',
                  palette = matlab::jet.colors(length(unique( combinedpulses.df$Feature))) )+
  ylab('')+ theme(legend.position = "none")



