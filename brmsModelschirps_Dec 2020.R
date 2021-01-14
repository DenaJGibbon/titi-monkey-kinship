
library(dplyr)
library(ggplot2)
library(tidybayes)
library(Hmisc)
library(brms)


# Dena's example -----
titi.kinship.mat.in <- 
  read.csv('/Users/denasmacbook/titii_acoustic_distance/titi.kinship.mat.csv')

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
  
datos.chirp$high.freq <- log(datos.chirp$high.freq)
datos.chirp$low.freq <- log(datos.chirp$low.freq)
#datos.chirp[,c('duration','low.freq','high.freq','note.rate')] <- scale(log(datos.chirp[,c('duration','low.freq','high.freq','note.rate')]))

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
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_noterate.chirp$ICCgenetic))
values <- draws_noterate.chirp$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_noterate.chirp$ICCinter))
values <- draws_noterate.chirp$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_noterate.chirp$ICCintra))
values <- draws_noterate.chirp$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.noterate <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

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
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_duration.chirp$ICCgenetic))
values <- draws_duration.chirp$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_duration.chirp$ICCinter))
values <- draws_duration.chirp$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_duration.chirp$ICCintra))
values <- draws_duration.chirp$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.duration <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

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
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

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
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_low.freq.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")



# Number of notes (need poisson)
datos.chirp$n.notes <- log(datos.chirp$n.notes)

model_nnotes.chirp_full<- 
  update(model_noterate.chirp_full, formula. = n.notes ~ .,
         newdata=datos.chirp)


summary(model_nnotes.chirp_full)
plot(model_nnotes.chirp_full)

brms::conditional_effects(model_nnotes.chirp_full)


hyp_nnotes.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_nnotes.chirp <- hypothesis(model_nnotes.chirp_full, hyp_nnotes.chirp, class = NULL))

draws_nnotes.chirp <- model_nnotes.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_nnotes.chirp, aes(ICC)) + geom_density(color = "black", fill = "white")

draws_nnotes.chirp <- model_nnotes.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_nnotes.chirp$ICCgenetic))
values <- draws_nnotes.chirp$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_nnotes.chirp$ICCinter))
values <- draws_nnotes.chirp$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_nnotes.chirp$ICCintra))
values <- draws_nnotes.chirp$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.nnotes <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

ggpubr::ggdensity(data=ICCvalues.nnotes,x='values',fill='ICC')


# PCA of all features
chirp.pca <- 
  prcomp(datos.chirp[,c(4:9,11:12)])

summary(chirp.pca)

chirp.pca$x[,1:2]

pca1 <- chirp.pca$x[,1]
pca2 <- chirp.pca$x[,2]
datos.chirp.pca <- cbind.data.frame(datos.chirp,pca1,pca2)

model_pca1.chirp_full <- 
  update(model_noterate.chirp_full, formula. = pca1 ~ ., newdata=datos.chirp.pca)

summary(model_pca1.chirp_full)

plot(model_pca1.chirp_full)

brms::conditional_effects(model_pca1.chirp_full)

hyp_pca1.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_pca1.chirp <- hypothesis(model_pca1.chirp_full, hyp_pca1.chirp, class = NULL))

draws_pca1.chirp <- model_pca1.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_pca1.chirp$ICCgenetic))
values <- draws_pca1.chirp$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_pca1.chirp$ICCinter))
values <- draws_pca1.chirp$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_pca1.chirp$ICCintra))
values <- draws_pca1.chirp$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.pca1 <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

ggpubr::ggdensity(data=ICCvalues.pca1,x='values',fill='ICC')





# PCA of all features
datos.chirp.pca <- subset(datos.chirp.pca,pca2 < 3000)

model_pca2.chirp_full <- 
  update(model_noterate.chirp_full, formula. = pca2 ~ ., newdata=datos.chirp.pca)

summary(model_pca2.chirp_full)
plot(model_pca2.chirp_full)

brms::conditional_effects(model_pca2.chirp_full)

hyp_pca2.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_pca2.chirp <- hypothesis(model_pca2.chirp_full, hyp_pca2.chirp, class = NULL))

draws_pca2.chirp <- model_pca2.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICCgenetic = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCinter = sd_indiv_id__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))%>%
  mutate(ICCintra = sigma/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ICC <- rep('genetic',length(draws_pca2.chirp$ICCgenetic))
values <- draws_pca2.chirp$ICCgenetic
ICCgenetic <- cbind.data.frame(ICC, values)

ICC <- rep('inter',length(draws_pca2.chirp$ICCinter))
values <- draws_pca2.chirp$ICCinter
ICCinter <- cbind.data.frame(ICC, values)

ICC <- rep('intra',length(draws_pca2.chirp$ICCintra))
values <- draws_pca2.chirp$ICCintra
ICCintra <- cbind.data.frame(ICC, values)

ICCvalues.pca2 <- rbind.data.frame(ICCgenetic,ICCinter,ICCintra)

ggpubr::ggdensity(data=ICCvalues.pca2,x='values',fill='ICC')



# PCA of all features
datos.chirp$bandwidth.mean <- log(datos.chirp$bandwidth.mean)

model_bandwidth.chirp_full <- 
  update(model_noterate.chirp_full, formula. = bandwidth.mean ~ ., newdata=datos.chirp)

summary(model_bandwidth.chirp_full)
plot(model_bandwidth.chirp_full)

brms::conditional_effects(model_bandwidth.chirp_full)

hyp_bandwidth.chirp <- paste("sd_name__Intercept^2 /", "(sd_name__Intercept^2 + sd_indiv_id__Intercept^2 + sigma^2) = 0")
(hyp_bandwidth.chirp <- hypothesis(model_bandwidth.chirp_full, hyp_bandwidth.chirp, class = NULL))

draws_bandwidth.chirp <- model_bandwidth.chirp_full %>%
  spread_draws(b_weight_kg, b_sexM, sd_name__Intercept, sd_indiv_id__Intercept, sigma) %>%
  mutate(ICC = sd_name__Intercept/(sd_name__Intercept + sd_indiv_id__Intercept + sigma))

ggplot(draws_bandwidth.chirp, aes(ICC)) + geom_histogram(color = "black", fill = "white")


# Density plots together
draws_noterate.chirp$Feature <- rep('Note Rate',nrow(draws_noterate.chirp))
draws_duration.chirp$Feature <- rep('Duration',nrow(draws_noterate.chirp))
draws_high.freq.chirp$Feature <- rep('High Frequency',nrow(draws_high.freq.chirp))
draws_low.freq.chirp$Feature <- rep('Low Frequency',nrow(draws_low.freq.chirp))

combinedchirps.df <- rbind.data.frame(draws_noterate.chirp,draws_duration.chirp,draws_high.freq.chirp,draws_low.freq.chirp)

ggpubr::ggdensity(data=combinedchirps.df,x='ICC',fill='Feature',facet.by = 'Feature',
                  palette = matlab::jet.colors(length(unique( combinedchirps.df$Feature))) )+
  ylab('')+ theme(legend.position = "none")

