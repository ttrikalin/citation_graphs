library(igraph)

mydata <- read.table('../lala.txt')

# fix the names of the columns
names <- rownames(mydata)
colnames(mydata) <- names 

G <- graph.adjacency(mydata, mode='directed')
