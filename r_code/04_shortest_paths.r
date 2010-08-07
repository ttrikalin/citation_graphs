# Now let's get the shortest paths between RCTs and OBS in the main graph  
paths.G.main <- get.shortest.path.distribution(G)

fname <- paste(project, '.paths.rct.obs.txt', sep='')
write('Main graph',  file = fname)
write.table(paths.G.main, file = fname , append=TRUE)

paths.G.all.index <- get.shortest.path.distribution(G.all.index)
write('\n\nIndex studies subgraph',  file = fname , append=TRUE)
write.table(paths.G.all.index, file = fname , append=TRUE)

# Now let's get the number of communications between the animal only and the 
# human or animal nodes

animal.only <- V(G.color)[(V(G.color)$mesh.animal==TRUE & V(G.color)$mesh.human==FALSE)]
human <- V(G.color)[V(G.color)$mesh.human==TRUE]

paths.G.human.animal <- get.shortest.path.distribution(G.color, nodes1=human , nodes2=animal.only, names=c('human', 'animal.only'))
write('\n\nMain graph - Human vs animal.only',  file = fname , append=TRUE)
write.table(paths.G.human.animal, file = fname , append=TRUE)


# Now let's get the number of communications between the mesh-RCT studies and the mesh OBS studies

mesh.rct <- V(G.color)[V(G.color)$design=='rct']
mesh.obs <- V(G.color)[V(G.color)$design=='obs']

paths.G.mesh.designs <- get.shortest.path.distribution(G.color, nodes1=mesh.rct , nodes2=mesh.obs, names=c('mesh.rct', 'mesh.obs'))
write('\n\nMain graph - Mesh RCTs vs Observational',  file = fname , append=TRUE)
write.table(paths.G.mesh.designs, file = fname , append=TRUE)

