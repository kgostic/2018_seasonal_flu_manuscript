#############################################
# Plot expectations from each model
# against data
#############################################

library(ggplot2)
library(reshape2)

source("prep_inputs.R")
plot_outfile_IN_h1 = 'KG_INSIGHT_test_h1.png'
plot_outfile_IN_h3 = 'KG_INSIGHT_test_h3.png'
plot_outfile_AZ_h1 = 'KG_AZ_test_h1.png'
plot_outfile_AZ_h3 = 'KG_AZ_test_h3.png'

flu_season_labels = c('2007-2008',
                      '2008-2009',
                      '2009pdm',
                      '2009-2010',
                      '2010-2011',
                      '2011-2012',
                      '2012-2013',
                      '2013-2014',
                      '2014-2015',
                      '2015-2016',
                      '2016-2017',
                      '2017-2018')

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
ggsave(plot_outfile_h1, dpi=300, plot = models_h1plots, width = plot_width, height = plot_height )

models_h3p <- ggplot(h3obs_plot, aes(x=Age_Group, y=Incidence)) + 
  geom_bar(stat = "identity", fill = "#fb9a99", color = "black") +
  xlab("Age Group") + ylab('Excess cases') +
  ggtitle("H3N2") #+ ylim(-55, 60)
models_h3plots <- models_h3p + theme_bw()+ 
                  text_sizes + 
                  facet_wrap(~ Season_f, ncol= 3, scales="free_y")
ggsave(plot_outfile_h3, dpi=300, plot = models_h3plots, width = plot_width, height = plot_height )
