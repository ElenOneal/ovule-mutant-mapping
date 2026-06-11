suppressMessages(library(dplyr))
suppressMessages(library(plotly))
suppressMessages(library(stringr))

plot_lg <- function(order){
  colnames(order) <- c('lg','contig','pos','dist')
  order <- split(order,f=order$lg)
  plotLMorder <- lapply(order, function(x) {
     lgtxt <- paste('',x$lg[1])
    plot_ly(data = x, x= ~dist, y= ~pos,color=~contig,type='scatter',colors = c('red','blue','green','black','purple','cyan','grey'),mode='markers') %>% layout(showlegend=TRUE,title = '', xaxis = list(title = lgtxt),yaxis = list(title = 'Distance'))
  })
}
