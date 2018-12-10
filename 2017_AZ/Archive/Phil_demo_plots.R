#############################################
# Plot expectations from each model
# against data
#############################################
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

library(ggplot2)
library(reshape2)

source("Phil_prep_inputs.R")
plot_outfile_IN_h1 = 'Compare_data/KG_INSIGHT_test_h1.png'
plot_outfile_IN_h3 = 'Compare_data/KG_INSIGHT_test_h3.png'
plot_outfile_AZ_h1 = 'Compare_data/KG_AZ_test_h1.png'
plot_outfile_AZ_h3 = 'Compare_data/KG_AZ_test_h3.png'
plot_outfile_IN_h1_ssn = 'Compare_data/KG_INSIGHT_ssn_h1.png'
plot_outfile_IN_h3_ssn = 'Compare_data/KG_INSIGHT_ssn_h3.png'
plot_outfile_IN_h1_country = 'Compare_data/KG_INSIGHT_country_h1.png'
plot_outfile_IN_h3_country = 'Compare_data/KG_INSIGHT_country_h3.png'


############################################
##  Make AZ plots first
############################################


## Load AZ data
source('Inputs_multinomial.R')

flu_season_labels = c('1993-1994',
                      '1994-1995',
                      '2002-2003',
                      '2003-2004',
                      '2004-2005',
                      '2005-2006',
                      '2006-2007',
                      '2007-2008',
                      '2008-2009',
                      '2009-2010',
                      '2010-2011',
                      '2011-2012',
                      '2012-2013',
                      '2013-2014',
                      '2014-2015')



demo_group_h1_exp = demo_group * (rowSums(inc_data$h1_inc))
demo_group_h3_exp = demo_group * (rowSums(inc_data$h3_inc))

demo_group_h1_plot = melt(as.matrix(demo_group_h1_exp), id="Season", variable.name = "Age_Group", value.name="Incidence")
demo_group_h3_plot = melt(as.matrix(demo_group_h3_exp), id="Season", variable.name = "Age_Group", value.name="Incidence")
h1obs_plot = melt(as.matrix(inc_data$h1_inc), id="Season", variable.name = "Age_Group", value.name = "Incidence")
h3obs_plot = melt(as.matrix(inc_data$h3_inc), id="Season", variable.name = "Age_Group", value.name = "Incidence")
names(demo_group_h1_plot) = names(demo_group_h3_plot) = names(h1obs_plot) = names(h3obs_plot) = c('Season', 'Age_Group', 'Incidence')


h1obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> h1obs_plot
h1obs_plot$Incidence = h1obs_plot$Incidence - demo_group_h1_plot$Incidence
h1obs_plot$Season_f = factor(h1obs_plot$Season, levels=flu_season_labels)

h3obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> h3obs_plot
h3obs_plot$Incidence = h3obs_plot$Incidence - demo_group_h3_plot$Incidence
h3obs_plot$Season_f = factor(h3obs_plot$Season, levels=flu_season_labels)


# Plotting parameters
text_sizes = theme(legend.position= "none", legend.title=element_blank(), 
                   plot.title = element_text(size = 18, face = "bold") , legend.text=element_text(size=12),
                   axis.title.x = element_text(face="bold", size=12),
                   axis.title.y = element_text(face="bold", size=12),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

pal = c(H1N1="#1f78b4", H3N2="#fb9a99")
plot_width = 6.6
plot_height = 5.5

# Demography plots
models_h1p <- ggplot(h1obs_plot, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H1N1") #+ ylim(-40, 150)
models_h1plots <- models_h1p + theme_bw()+ 
                  text_sizes +
                  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_AZ_h1, dpi=300, plot = models_h1plots, width = plot_width, height = plot_height )

models_h3p <- ggplot(h3obs_plot, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#fb9a99", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H3N2") #+ ylim(-55, 60)
models_h3plots <- models_h3p + theme_bw()+ 
                  text_sizes + 
                  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_AZ_h3, dpi=300, plot = models_h3plots, width = plot_width, height = plot_height )






#### Repeat for INSIGHT data (Country x Season)
## dem_tab
## H1.summary
## H3.summary
IN_dem_h1_exp = dem_tabH1 *rowSums(H1.summary)
IN_dem_h3_exp = dem_tabH3 *rowSums(H3.summary)
IN_dem_H1_plot = melt(as.matrix(IN_dem_h1_exp), id="Season", variable.name = "Age_Group", value.name="Incidence")
IN_dem_H3_plot = melt(as.matrix(IN_dem_h3_exp), id="Season", variable.name = "Age_Group", value.name="Incidence")


IN_h1obs_plot = melt(as.matrix(H1.summary), id="Season", variable.name = "Age_Group", value.name = "Incidence")
IN_h3obs_plot = melt(as.matrix(H3.summary), id="Season", variable.name = "Age_Group", value.name = "Incidence")
names(IN_dem_H1_plot) = names(IN_dem_H3_plot) = names(IN_h1obs_plot) = names(IN_h3obs_plot) = c('Season', 'Age_Group', 'Incidence')


#IN_h1obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h1obs_plot
IN_h1obs_plot$Incidence = IN_h1obs_plot$Incidence - IN_dem_H1_plot$Incidence
IN_h1obs_plot$Season_f = factor(IN_h1obs_plot$Season)

#IN_h3obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h3obs_plot
IN_h3obs_plot$Incidence = IN_h3obs_plot$Incidence - IN_dem_H3_plot$Incidence
IN_h3obs_plot$Season_f = factor(IN_h3obs_plot$Season)


# Plotting parameters
text_sizes = theme(legend.position= "none", legend.title=element_blank(), 
                   plot.title = element_text(size = 18, face = "bold") , legend.text=element_text(size=12),
                   axis.title.x = element_text(face="bold", size=12),
                   axis.title.y = element_text(face="bold", size=12),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

pal = c(H1N1="#1f78b4", H3N2="#fb9a99")
plot_width = 6.6
plot_height = 5.5

# Demography plots
models_h1p <- ggplot(IN_h1obs_plot, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H1N1") #+ ylim(-40, 150)
models_h1plots <- models_h1p + theme_bw()+ 
  text_sizes +
  facet_wrap(~ Season, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h1, dpi=300, plot = models_h1plots, width = plot_width, height = plot_height )

models_h3p <- ggplot(IN_h3obs_plot, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#fb9a99", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H3N2") #+ ylim(-55, 60)
models_h3plots <- models_h3p + theme_bw()+ 
  text_sizes + 
  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h3, dpi=300, plot = models_h3plots, width = plot_width, height = plot_height )







#### Repeat for INSIGHT data, aggregated by country
# dem_country_H1 is already weighted by case number
# dem_country_H3  ""
IN_dem_H1_country = melt(as.matrix(dem_country_H1), id="Season", variable.name = "Age_Group", value.name="Incidence")
IN_dem_H3_country = melt(as.matrix(dem_country_H3), id="Season", variable.name = "Age_Group", value.name="Incidence")


IN_h1obs_country = melt(as.matrix(H1.country), id="Season", variable.name = "Age_Group", value.name = "Incidence")
IN_h3obs_country = melt(as.matrix(H3.country), id="Season", variable.name = "Age_Group", value.name = "Incidence")
names(IN_dem_H1_country) = names(IN_dem_H3_country) = names(IN_h1obs_country) = names(IN_h3obs_country) = c('Season', 'Age_Group', 'Incidence')


#IN_h1obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h1obs_plot
IN_h1obs_country$Incidence = IN_h1obs_country$Incidence - IN_dem_H1_country$Incidence
IN_h1obs_country$Season_f = factor(IN_h1obs_country$Season)

#IN_h3obs_country %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h3obs_country
IN_h3obs_country$Incidence = IN_h3obs_country$Incidence - IN_dem_H3_country$Incidence
IN_h3obs_country$Season_f = factor(IN_h3obs_country$Season)


# Plotting parameters
text_sizes = theme(legend.position= "none", legend.title=element_blank(), 
                   plot.title = element_text(size = 18, face = "bold") , legend.text=element_text(size=12),
                   axis.title.x = element_text(face="bold", size=12),
                   axis.title.y = element_text(face="bold", size=12),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

pal = c(H1N1="#1f78b4", H3N2="#fb9a99")
plot_width = 6.6
plot_height = 5.5

# Demography plots
models_h1p <- ggplot(IN_h1obs_country, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H1N1") #+ ylim(-40, 150)
models_h1plots <- models_h1p + theme_bw()+ 
  text_sizes +
  facet_wrap(~ Season, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h1_country, dpi=300, plot = models_h1plots, width = plot_width, height = plot_height )

models_h3p <- ggplot(IN_h3obs_country, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#fb9a99", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H3N2") #+ ylim(-55, 60)
models_h3plots <- models_h3p + theme_bw()+ 
  text_sizes + 
  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h3_country, dpi=300, plot = models_h3plots, width = plot_width, height = plot_height )













#### Repeat for INSIGHT data, aggregated by season
# dem_season_H1 is already weighted by case number
# dem_season_H3  ""
IN_dem_H1_season = melt(as.matrix(dem_season_H1), id="Season", variable.name = "Age_Group", value.name="Incidence")
IN_dem_H3_season = melt(as.matrix(dem_season_H3), id="Season", variable.name = "Age_Group", value.name="Incidence")


IN_h1obs_season = melt(as.matrix(H1.season), id="Season", variable.name = "Age_Group", value.name = "Incidence")
IN_h3obs_season = melt(as.matrix(H3.season), id="Season", variable.name = "Age_Group", value.name = "Incidence")
names(IN_dem_H1_season) = names(IN_dem_H3_season) = names(IN_h1obs_season) = names(IN_h3obs_season) = c('Season', 'Age_Group', 'Incidence')


#IN_h1obs_plot %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h1obs_plot
IN_h1obs_season$Incidence = IN_h1obs_season$Incidence - IN_dem_H1_season$Incidence
IN_h1obs_season$Season_f = factor(IN_h1obs_season$Season)

#IN_h3obs_season %>% mutate(Season=ifelse(Season!=2009.5, paste(Season-1, Season, sep='-'), '2009pdm')) -> IN_h3obs_season
IN_h3obs_season$Incidence = IN_h3obs_season$Incidence - IN_dem_H3_season$Incidence
IN_h3obs_season$Season_f = factor(IN_h3obs_season$Season)


# Plotting parameters
text_sizes = theme(legend.position= "none", legend.title=element_blank(), 
                   plot.title = element_text(size = 18, face = "bold") , legend.text=element_text(size=12),
                   axis.title.x = element_text(face="bold", size=12),
                   axis.title.y = element_text(face="bold", size=12),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

pal = c(H1N1="#1f78b4", H3N2="#fb9a99")
plot_width = 6.6
plot_height = 5.5

# Demography plots
models_h1p <- ggplot(IN_h1obs_season, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#1f78b4", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H1N1") #+ ylim(-40, 150)
models_h1plots <- models_h1p + theme_bw()+ 
  text_sizes +
  facet_wrap(~ Season, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h1_ssn, dpi=300, plot = models_h1plots, width = plot_width, height = plot_height )

models_h3p <- ggplot(IN_h3obs_season, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#fb9a99", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H3N2") #+ ylim(-55, 60)
models_h3plots <- models_h3p + theme_bw()+ 
  text_sizes + 
  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_IN_h3_ssn, dpi=300, plot = models_h3plots, width = plot_width, height = plot_height )


