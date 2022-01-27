require(WebPower)
require(magrittr)
require(dplyr)
require(ggplot2)
n = p = list()
effect_size = power = seq(0.05,0.95, by = 0.01)
k = 0 

for(i in 1:length(effect_size)){
  for(j in 1:length(power)){
    k = k + 1
    # Estimation of statistical power for Pearson correlation
    # Calculate webpower analysis of sample size applied to correlation for specific effect size and power
    n[[k]] = wp.correlation(n=NULL, 
                            r=effect_size[i], #	Effect size or correlation. 
                            alpha = 0.05/20000,
                            power=power[j], 
                            alternative="two.sided") %>%
      unclass() %>%
      as.data.frame()
    
    # Estimation of statistical power for logistic regression
    # Calculate webpower analysis of sample size applied to logistic regression for specific effect size and power
    res = try(wp.logistic(n=NULL,
                          p0 = effect_size[i],
                          p1 = 1-effect_size[i],
                          # r=effect_size[i], #	Effect size or correlation.
                          alpha = 0.05/20000,
                          power=power[j],
                          alternative="two.sided", 
                          parameter = c(0,1),
                          family = "normal"))
    
    if(class(res)!= "try-error"){
      p[[k]] = res %>%
        unclass() %>%
        as.data.frame()
    }
  }
}
n %<>% dplyr::bind_rows()
p %<>% dplyr::bind_rows()

# Assign a label when the power is over 0.8.
# The function ceiling gets the upper integer value of a float.
#n %<>%
#  mutate(label = ifelse(power >= 0.8, 
#                        ceiling(n), NA)) 

#p %<>%
#  mutate(label = ifelse(power >= 0.8, 
#                        ceiling(n), NA))



# Plot the power in function of the effect size (r) and sample size (n)
output_sample_size_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation.png"
ggplot(n) +
  aes(x = n, y = r, 
      colour = power,
      #label = label,
      group = power) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9))+
  # geom_text() +
  #geom_label_repel() + 
  scale_color_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  labs(title = "Sample Size for correlation", 
       x = "Number of samples",
       y = "Correlation"
       #caption = "To be used in the co-expression networks. \nSize per group."
       ) +
  theme(legend.position = "bottom")
ggsave(
  output_sample_size_plot,
  dpi = 1200
)


# Plot a more specific figure
require(ggrepel)
output_sample_size_reduced_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_reduced.png"
n_reduced = n %>% filter((n <= 1000) & (r < 0.3))
prediction <- data.frame(x=x_axis, y=y_axis)
n_reduced %<>% mutate(label = ifelse(((power == 0.95) & (n < 600) & (n >550)), "Scipher", ifelse(((power == 0.95) & (n < 760) & (n >720)), "GTEx Whole Blood", NA)))
ggplot(n_reduced) +
  aes(x = n, y = r, 
      colour = power,
      #label = label,
      group = power) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  geom_function(fun = fit1, colour = "red")+
  #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3))+
  # geom_text() +
  #geom_label_repel() + 
  scale_color_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  labs(#title = "Sample Size for correlation", 
       x = "Number of samples",
       y = "Correlation"
       #caption = "To be used in the co-expression networks. \nSize per group."
  ) +
  theme(legend.position = "bottom")
ggsave(
  output_sample_size_reduced_plot,
  dpi = 1200
)

n_reduced %>% filter((power == 0.95) & (n < 600) & (n > 550))
n_reduced %>% filter((power == 0.95) & (n < 760) & (n > 720))

# Plot applied to differential expression analysis
require(ggrepel)
ggplot(p) +
  aes(x = p0, 
      y = n, 
      colour = power, 
      group = power, 
      label = label) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     lim = c(0,0.5))+
  # geom_text() +
  geom_label_repel() + 
  scale_color_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  labs(title = "Sample Size for logistic regression", 
       x = "Effect Size", 
       caption = "To be used in the differential expression analysis. \nSize total.")+
  theme(legend.position = "bottom")



