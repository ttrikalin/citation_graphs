# write the figure with the descriptive of the main path
pdf(file = paste(project, '.G.spc.pdf', sep=''))
myplot(fix.vertex.colors(G.spc), offset=0)
dev.off(dev.prev())


# this is the colored version of G, the main graph -- log of hub scores weighting 
G.color <- fix.vertex.colors(G)

hs <- hub.score(G.color)
mult <- 5000
min.add <- 0.05 
H <- abs(hs$vector)
H <- (log(H*mult)>0) * log(H*mult)
H <- (H==0)*min.add + H 

pdf(file = paste(project, '.G.color.random.pdf', sep=''))
plot(G.color, vertex.size=H, edge.width=0.1, edge.arrow.size=0.2)
dev.off(dev.prev())
try(dev.off())

# this is the same but ordered by year 
pdf(file = paste(project, '.G.color.year.pdf', sep=''))
myplot(G.color, width=0.1, v.size=H)
dev.off(dev.prev())
try(dev.off())

# this is the graph of the index publications -- only 
pdf(file = paste(project, '.G.index.year.pdf', sep=''))
myplot(G.all.index, offset=0.2)
dev.off(dev.prev())
try(dev.off())


# this is the graph of the index publications -- only 
pdf(file = paste(project, '.G.mesh.design.pdf', sep=''))
G.temp<-G.mesh.primary.data
# temporarily break the colored encoding of rcts and obs
V(G.temp)[V(G.temp)$design=='rct']$color ='red'
V(G.temp)[V(G.temp)$design=='obs']$color ='blue'
myplot(G.temp, offset=0.05)
G.temp <- NULL  # make G.temp vanish!
dev.off(dev.prev())
try(dev.off())
