library(dplyr)
library(ggplot2)
library(scales)
setwd('/Users/eleno/Dropbox/Medea/03_ddrad/03_ovules')
map <- read.table('flipped_mapping_order.txt',h=T,sep='\t')
colnames(map)
#map$position <- map$position/1000000
p<-ggplot(map, aes(x = position, y = cM)) +
  geom_point(size=.5) +  # Scatter plot for two numeric variables
  facet_wrap(~ lg, nrow = 4, ncol = 4,scales = "free") +  # Create a 2-row, 7-column grid
  scale_x_continuous(
    labels = function(x) x / 1e6,  # Convert y-axis values to megabases
    name = "Physical position (Mb)"  # Update the label to indicate the scale is in megabases
   ) +#scale_x_continuous(
  #   breaks = pretty_breaks(n = 4),  # Limit to a maximum of 4 breaks with whole numbers
  #   labels = scales::label_number(accuracy = 1)  # Ensure whole numbers
  # )+
  labs(
    y = "Genetic distance (cM)"
  )+
  theme_minimal()+ # Use a clean theme
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # Add box around each facet
    axis.title.x = element_text(margin = margin(t = 20)),
    strip.background = element_blank(),  # Optional: remove grey background from facet labels
    strip.text = element_text(face = "bold")  # Optional: bold facet labels
  )
p

ggsave("ovule_linkage_map_plot.jpg", plot = p, width = 6, height = 8, units = "in", dpi = 300)
