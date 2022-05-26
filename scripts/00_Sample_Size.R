require(WebPower)
require(magrittr)
require(dplyr)
require(ggplot2)
options(bitmapType='cairo')


############
# MANUALLY #
############

N = seq(4, 15000, 1)

total_num_genes = 18884
#total_num_genes = 20000
total_num_edges = (total_num_genes*(total_num_genes-1)/2)
alpha_level = 0.05/ (total_num_edges)
alpha_level = 0.05/ (total_num_genes*total_num_genes)/2
Za = qnorm(alpha_level, lower.tail=F)
Zb = qnorm(0.8, lower.tail=F)

#r = seq(0.1, 1, by = 0.1)
#r = 0.3
#Cr =  0.5 * log ((1 + r)/(1 - r))
#N = ((Za + Zb)/Cr)^2+3

#Za_cal = pnorm( Cr*sqrt(SS - 3) - Zb)
a = (Za + Zb)/(0.5 * sqrt(N-3))
r = (exp(a) - 1) / (exp(a) + 1)
ss_df = data.frame(N=N, r=r)
print(ss_df %>% filter(N == 10198))

#ne * Za_cal
ggplot(ss_df) +
  aes(x = N, y = r) +
  geom_point(shape = "circle", size = 1.5, col="#810f7c", alpha = 0.5) +
  #scale_color_distiller(palette = "BuPu", direction = 1) +
  labs(
    x = "Number of samples",
    y = "Pearson critical correlation",
    title = "Number of samples for critical correlation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17L,
                              face = "bold"),
    axis.title.y = element_text(size = 16L,
                                face = "bold"),
    axis.title.x = element_text(size = 16L,
                                face = "bold"),
    axis.text = element_text(size = 15L)
  )

output_sample_size_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation.png"
ggsave(
  output_sample_size_plot,
  dpi = 1200,
  width = 8750,
  height = 7000,
  units = c("px")
)

# Reduced plot
ss_reduced_df = ss_df %>% filter((N >= 10150) & (N <= 10250))

ggplot(ss_reduced_df) +
  aes(x = N, y = r) +
  geom_point(shape = "circle", size = 3, col="#810f7c", alpha = 0.5) +
  labs(
    x = "Number of samples",
    y = "Pearson critical correlation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17L,
                              face = "bold"),
    axis.title.y = element_text(size = 16L,
                                face = "bold"),
    axis.title.x = element_text(size = 16L,
                                face = "bold"),
    axis.text = element_text(size = 15L)
  )

output_sample_size_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_reduced.png"
ggsave(
  output_sample_size_plot,
  dpi = 1200,
  width = 8750,
  height = 7000,
  units = c("px")
)


r = seq(0.05, 1, by = 0.01)
Cr =  0.5 * log ((1 + r)/(1 - r))
N = ((Za + Zb) / Cr)^2 + 3
ss_df2 = data.frame(N=N, r=r)
plot(N, r)



############
# WEBPOWER #
############

n = p = list()
effect_size = power = seq(0.05,0.95, by = 0.01)
k = 0 
num_genes = 14266

for(i in 1:length(effect_size)){
  for(j in 1:length(power)){
    k = k + 1
    # Estimation of statistical power for Pearson correlation
    # Calculate webpower analysis of sample size applied to correlation for specific effect size and power
    n[[k]] = wp.correlation(n=NULL, 
                            r=effect_size[i], #	Effect size or correlation. 
                            alpha = 0.05/(num_genes*num_genes/2),
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
                          alpha = 0.05/(num_genes*num_genes/2),
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
output_sample_size_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_webpower.png"
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
  theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 13), plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15))

ggsave(
  output_sample_size_plot,
  dpi = 1200,
  width = 8750,
  height = 7000,
  units = c("px")
)


# Plot a more specific figure
require(ggrepel)
output_sample_size_reduced_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_reduced_webpower.png"
n_reduced = n %>% filter((n <= 1000) & (r < 0.4))
#prediction <- data.frame(x=x_axis, y=y_axis)
n_reduced %<>% mutate(label = ifelse(((power == 0.95) & (n < 600) & (n >550)), "Scipher", ifelse(((power == 0.95) & (n < 760) & (n >720)), "GTEx Whole Blood", NA)))
ggplot(n_reduced) +
  aes(x = n, y = r, 
      colour = power,
      #label = label,
      group = power) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
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
  theme(legend.position = "none", plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15))

ggsave(
  output_sample_size_reduced_plot,
  dpi = 1200
)

print(n_reduced %>% filter((power == 0.95) & (n < 600) & (n > 550)))
print(n_reduced %>% filter((power == 0.95) & (n < 760) & (n > 720)))

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




###########
# HICLIMR #
###########

library(HiClimR)
results_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("n", "cor", "p.value"))))
n_list = seq(10, 1000, by = 10)
for(size in n_list){
  result = minSigCor(n = size, alpha = 0.05, r = seq(0, 1, by = 1e-6))
  results_df = rbind(results_df, data.frame(n=size, cor=result$cor, p.value=result$p.value))
}

output_sample_size_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_hiclimr.png"
ggplot(results_df) +
  aes(x = n, y = cor) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5, col="#810f7c") +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9))+
  # geom_text() +
  #geom_label_repel() + 
  #scale_color_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  labs(title = "Sample Size for correlation", 
       x = "Number of samples",
       y = "Correlation"
       #caption = "To be used in the co-expression networks. \nSize per group."
  ) +
  theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 13), plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15))

ggsave(
  output_sample_size_plot,
  dpi = 1200,
  width = 8750,
  height = 7000,
  units = c("px")
)


# Plot a more specific figure
output_sample_size_reduced_plot = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/theoretical_sample_size_vs_correlation_hiclimr_reduced.png"
results_reduced_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("n", "cor", "p.value"))))
n_list = seq(500, 600, by = 1)
for(size in n_list){
  result = minSigCor(n = size, alpha = 0.05, r = seq(0, 1, by = 1e-6))
  results_reduced_df = rbind(results_reduced_df, data.frame(n=size, cor=result$cor, p.value=result$p.value))
}

ggplot(results_reduced_df) +
  aes(x = n, y = cor) +
  geom_line(size = 0.5) +
  geom_point(size = 3, alpha = 0.5, col="#810f7c") +
  #scale_y_continuous(breaks = c(0.1, 0.5, 0.9))+
  # geom_text() +
  #geom_label_repel() + 
  #scale_color_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  labs(title = "Sample Size for correlation", 
       x = "Number of samples",
       y = "Correlation"
       #caption = "To be used in the co-expression networks. \nSize per group."
  ) +
  theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 13), plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15))

ggsave(
  output_sample_size_reduced_plot,
  dpi = 1200
)

print(results_reduced_df %>% filter(n == 584))



