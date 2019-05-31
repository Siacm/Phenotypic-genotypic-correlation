install.packages("openxlsx", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("plotly", dependencies=TRUE)

library(openxlsx)
library(ggplot2)
library(plotly)

#######Using GGPLOTLY. OTHER PLOTLY GRAPH BELOW################

theme_set(
  theme_bw() + 
    theme(legend.position = "right")
)

getSheetNames("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks_cov_notadjusted/updated_SRforbubbleplot_notadjusted.xlsx")
SR_forbubble<-read.xlsx("~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks_cov_notadjusted/updated_SRforbubbleplot_notadjusted.xlsx", sheet = "Short_SR_list" )
colnames(SR_forbubble)

ggplot(SR_forbubble, aes(group=Gene, x = Identity , y = Coverage)) + 
  geom_point(aes(color = Sensitivity, size = Count)) +
  scale_fill_manual(values = c('#C61951', '#1972A4')) +
  xlim(80,105) + scale_y_continuous(breaks = seq(0, 110, 10), limits = c(0, 100)) + scale_size_continuous(breaks = seq(5, 580, 285), limits = c(1, 600), range = c(2, 15))+ labs(title = "Identity vs Coverage In All Detected Samples", x= "Identity (%)", y=("Coverage (%)")) +
  theme(title = element_text(size = 20), legend.title = element_blank(), legend.text=element_text(size=20), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15))
ggsave('~/Documents/PhD/Chapter_2/amr_files/resistance_workbooks/bubble_image.png', width = 16, height = 9, dpi = 100)

bubbleplot<-ggplot(SR_forbubble, aes(group=Gene, x = Identity , y = Coverage)) + 
  geom_point(aes(color = Sensitivity, size = Count)) +
  scale_fill_manual(values = c('#C61951', '#1972A4')) +
  xlim(80,105) + scale_y_continuous(breaks = seq(0, 110, 10), limits = c(0, 100)) + scale_size_continuous(breaks = seq(5, 580, 285), limits = c(1, 600), range = c(2, 15))+ labs(title = "Identity vs Coverage In All Detected Samples", x= "Identity (%)", y=("Coverage (%)")) +
  theme(title = element_text(size = 20), legend.title = element_blank(), legend.text=element_text(size=20), axis.text.x=element_text(size=13), axis.text.y =element_text(size=15))


ggplotly(bubbleplot)

#######Using PLOTLY.################

# We reuse the color scheme from the plotly website
colors <- c('#C61951', '#1972A4')

# Using Plotly
r1 <- plot_ly(
  SR_forbubble, x = ~`Identity`, y = ~`Coverage`,
  color = ~`Sensitivity`, type = "scatter",
  mode="markers", colors=~colors, size=~`Count`,
  marker = list(symbol = 'circle', sizemode = 'diameter',
                line = list(width = 2, color = '#FFFFFF'), opacity=0.4),
  text = ~paste(sep='','Gene:', Gene,
                '<br>Identity', Identity,
                '%', '<br>Coverage:', Coverage,
                '<br>Sensivity:', Sensitivity,'<br>Count:'
                , Count)) %>%
  layout(
    title="% Identity vs % Coverage",
    
    xaxis = list(title = '% Identity',
                 gridcolor = 'rgb(255, 255, 255)',
                 range=c(80,105),
                 zerolinewidth = 1,
                 ticklen = 15,
                 gridwidth = 2),
    yaxis = list(title = '% Coverage',
                 gridcolor = 'rgb(255, 255, 255)',
                 range=c(0,110),
                 zerolinewidth = 1,
                 ticklen = 15,
                 gridwith = 2),
    paper_bgcolor = 'rgb(243, 243, 243)',
    plot_bgcolor = 'rgb(243, 243, 243)'
  )

r1


r2 <- plot_ly() %>% 
  add_markers(x = 1, 
              y = seq_len(length(unique(SR_forbubble$Count))),
              size = sort(unique(SR_forbubble$Count)),
              showlegend = F, 
              marker = list(sizeref=0.03,sizemode="area")) %>% #####change for bubble sizes
  layout(
    annotations = list(
      list(x = 0.2, 
           y = 1, 
           text = "", 
           showarrow = F, 
           xref='paper', 
           yref='paper')),
    xaxis = list(
      zeroline=F,
      showline=F,
      showticklabels=F,
      showgrid=F),
    yaxis=list(
      showgrid=F,
      tickmode = "array",
      tickvals = seq_len(length(unique(SR_forbubble$Count))),
      ticktext = sort(unique(SR_forbubble$Count))))
 r2

 
 #######option two for r2##################
 legend.sizes = seq(1, 600, 250)
 ax = list(zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
 mk = list(sizeref=0.03, sizemode="area", color="black")
 
 r3 = plot_ly() %>%
   add_markers(x = 1, y = legend.sizes, size = legend.sizes, showlegend = F, marker = mk) %>%
   layout(xaxis = ax, yaxis = list(showgrid = FALSE))
 
 r3

### save using EXPORT to HTML