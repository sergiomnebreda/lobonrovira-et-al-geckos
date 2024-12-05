#########################
# # # # 00.REVIEW # # # #
#########################

rm(list = ls())
library(geomorph)

load("all.dt.bin")

# Species info ###
spec <- read.csv("info20Dec2022.csv", sep = ",")
str(spec)
rownames(spec) <- spec$Species

# Shape data ###
# Ear
ear <- dt.fixed$ear

all(dimnames(ear)[[3]]%in%spec$Species)
dimnames(ear)[[3]][which((dimnames(ear)[[3]]%in%spec$Species)==F)]
all(spec$Species%in%dimnames(ear)[[3]])
spec$Species[which((spec$Species%in%dimnames(ear)[[3]])==F)]
spec.ear <- spec[dimnames(ear)[[3]],]
all(spec.ear$Species==dimnames(ear)[[3]])

# Jaw
jaw <- dt.fixed$jaw

all(dimnames(jaw)[[3]]%in%spec$Species)
all(spec$Species%in%dimnames(jaw)[[3]])
spec$Species[which((spec$Species%in%dimnames(jaw)[[3]])==F)]
spec.jaw <- spec[dimnames(jaw)[[3]],]
all(spec.jaw$Species==dimnames(jaw)[[3]])

# Mandibles
mand <- dt.fixed$mand

all(dimnames(mand)[[3]]%in%spec$Species)
all(spec$Species%in%dimnames(mand)[[3]])
spec$Species[which((spec$Species%in%dimnames(mand)[[3]])==F)]
spec.mand <- spec[dimnames(mand)[[3]],]
all(spec.mand$Species==dimnames(mand)[[3]])

# SUPERIMPOSITION ###
# Ear
ear.semi <- read.table("ear_semi.txt", sep = " ")
gpa.ear <- gpagen(ear, curves = ear.semi)
plotOutliers(gpa.ear$coords, inspect.outliers = T) 
# 
# graniticolus - lms 7 and 8 switched?
# sp_AngDRC - idem
# bivittis - check lm 1

# Jaw
gpa.jaw <- gpagen(jaw)
plotOutliers(gpa.jaw$coords, inspect.outliers = T) # Nothing special

# Mandibles
gpa.mand <- gpagen(mand)
plotOutliers(gpa.mand$coords, inspect.outliers = T) # Idem

# SAVE TO EXTERNAL FILE ###
all.gpa <- list(ear = list(gpa = gpa.ear, spec = spec.ear),
                jaw = list(gpa = gpa.jaw, spec = spec.jaw),
                mand = list(gpa = gpa.mand, spec = spec.mand))

save(all.gpa, file = "all.gpa.bin")



###############################
# # # # 00.REVIEW_PHYLO # # # #
###############################

rm(list = ls())

library(ape)

spec <- read.csv("info20Dec2022.csv", sep = ",")
str(spec)

# Read phylo ###
tr <- read.nexus("Lygo5_MCC.nex")

spec$Species[which(spec$Species%in%tr$tip.label==F)]
#[1] "Lygodactylus_angularis_paurospilus" "Lygodactylus_broadleyi"            
#[3] "Lygodactylus_decaryi"               "Lygodactylus_depressus"            
#[5] "Lygodactylus_fischeri"              "Lygodactylus_grandisonae"          
#[7] "Lygodactylus_inexpectatus"          "Lygodactylus_insularis"            
#[9] "Lygodactylus_manni"                 "Lygodactylus_nyaneka"              
#[11] "Lygodactylus_scheffleri"            "Lygodactylus_scorteccii"           
#[13] "Lygodactylus_somalicus_battersbyi" 

spec.phylo <- spec[which(spec$Species%in%tr$tip.label==T),]
sp.drop <- tr$tip.label[which(tr$tip.label%in%spec.phylo$Species==F)]
tr.spec <- drop.tip(tr, sp.drop)

phy.spec <- list(tr = tr.spec, spec = spec.phylo)
save(phy.spec, file = "phy.spec.MATCHED.bin")



#######################
# # # # 01.PCAs # # # #
#######################

rm(list = ls())
library(geomorph)

# Import superimposed data ###
# And split into datasets
load("all.gpa.bin")

ear.sh <- all.gpa$ear$gpa$coords
ear.cs <- all.gpa$ear$gpa$Csize
ear.spec <- all.gpa$ear$spec

jaw.sh <- all.gpa$jaw$gpa$coords
jaw.cs <- all.gpa$jaw$gpa$Csize
jaw.spec <- all.gpa$jaw$spec

mand.sh <- all.gpa$mand$gpa$coords
mand.logcs <- log(all.gpa$mand$gpa$Csize)
write.csv(mand.logcs, 'skull.cs.inter.csv')
mand.spec <- all.gpa$mand$spec
write.csv(mand.spec, 'skull.clasif.inter.csv')


# PCAs ###
pdf("PCAs.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "EAR", font = 2)
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Continent))
legend("topleft", legend = levels(as.factor(ear.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green"), bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "red"), bty = "n")

# Jaw
pca.jaw <- gm.prcomp(jaw.sh)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Continent))
legend("topleft", legend = levels(as.factor(jaw.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

# Mandible
pca.mand <- gm.prcomp(mand.sh)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand)
text(par()$usr[1]*0.7, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Continent))
legend("topleft", legend = levels(as.factor(mand.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "blue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

dev.off()


###
# Visualize shape change using deformation grids ###
###

# Ear
# Javi, aqui pasa algo raro, hay que revisar los datos del oido...
msh.ear <- mshape(ear.sh)
# minPC1 vs. mean shape
plotRefToTarget(pca.ear$shapes$shapes.comp1$min, msh.ear, method = "vector") 
# maxPC1 vs. mean shape
plotRefToTarget(pca.ear$shapes$shapes.comp1$max, msh.ear, method = "vector") 

# Jaw
msh.jaw <- mshape(jaw.sh)
# minPC1 vs. mean shape
plotRefToTarget(pca.jaw$shapes$shapes.comp1$min, msh.jaw, method = "vector") 
# maxPC1 vs. mean shape
plotRefToTarget(pca.jaw$shapes$shapes.comp1$max, msh.jaw, method = "vector") 

# Mandible
msh.mand <- mshape(mand.sh)
# minPC1 vs. mean shape
plotRefToTarget(pca.mand$shapes$shapes.comp1$min, msh.mand, method = "vector") 
# maxPC1 vs. mean shape
plotRefToTarget(pca.mand$shapes$shapes.comp1$max, msh.mand, method = "vector") 



###################################
# # # # 01b.phyPCA_versions # # # #
###################################

rm(list = ls())
library(geomorph)
library(phytools)

# Import superimposed data ###
load("all.gpa.bin")

# Import phylogeny 
load("phy.spec.MATCHED.bin")
tr <- phy.spec$tr

# Split into datasets and match with phylo
tr.drop <- tr$tip.label[which(tr$tip.label%in%dimnames(all.gpa$ear$gpa$coords)[[3]]==F)]
ear.tr <- drop.tip(tr, tr.drop)
ear.sh <- all.gpa$ear$gpa$coords[,,ear.tr$tip.label]
ear.cs <- all.gpa$ear$gpa$Csize[ear.tr$tip.label]
ear.spec <- all.gpa$ear$spec[ear.tr$tip.label,]

jaw.sh <- all.gpa$jaw$gpa$coords[,,tr$tip.label]
jaw.cs <- all.gpa$jaw$gpa$Csize[tr$tip.label]
jaw.spec <- all.gpa$jaw$spec[tr$tip.label,]

mand.sh <- all.gpa$mand$gpa$coords[,,tr$tip.label]
mand.cs <- all.gpa$mand$gpa$Csize[tr$tip.label]
mand.spec <- all.gpa$mand$spec[tr$tip.label,]

###
# PhyloMorphospace ###
# Traditional PCA (as in script 01a, but with reduced data to match phylo)
# Projection of phylogeny on morphospace
# Allows to visually review to what extent the data align with the phylogeny
# This will be formally tested later through phylogenetic signal analysis
pdf("phylomorphospace.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear, phylo = T, xlim = c(-0.15, 0.12),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.9, par()$usr[4]*0.95, labels = "EAR", font = 2)

# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw,  phylo = T, xlim = c(-0.075, 0.15),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand,  phylo = T, xlim = c(-0.1, 0.12),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.8, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)

dev.off()

###
# phyloPCA ###
# Rotates the data such as to reduce any effect of phylogeny
pdf("phyPCAs.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "EAR", font = 2)
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Continent))
legend("topleft", legend = levels(as.factor(ear.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green"), bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "red"), bty = "n")

# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Continent))
legend("topleft", legend = levels(as.factor(jaw.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand)
text(par()$usr[1]*-0.99, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Continent))
legend("topleft", legend = levels(as.factor(mand.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "blue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

dev.off()

###
# PaCA ###
# Aligns the data to phylogenetic signal, rather than maxVAR
# In a way, it allows to visualize the phylogenetic component of shape variation

pdf("PaCA.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear, phylo = T, xlim = c(-0.15, 0.12),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.9, par()$usr[4]*0.95, labels = "EAR", font = 2)

# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw,  phylo = T, xlim = c(-0.025, 0.05),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand,  phylo = T, xlim = c(-0.075, 0.12),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.8, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)

dev.off()



# A DIFFERENT WAY TO PLOT PCAs #

phyloPCA <- gm.prcomp(mand.sh, phy = tr)   
summary(phyloPCA)
plot(phyloPCA, axis1 = 1, axis2 = 2)
plot(phyloPCA, phylo = TRUE, main = "phyloPCA", loadings = TRUE , phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
arrows(0, 0, phyloPCA$rotation[,1]*100, phyloPCA$rotation[,2]*100, length = 0.2)



# Extract PC axes for plotting
PCAvalues <- data.frame(phyloPCA$x)

df<-data.frame(cs = log(all.gpa$mand$gpa$Csize[tr$tip.label]), 
               earcs = log(all.gpa$ear$gpa$Csize[tr$tip.label]),
               groups = as.factor(all.gpa$mand$spec[tr$tip.label, "Group"]),
               cont = as.factor(all.gpa$mand$spec[tr$tip.label, "Continent"]),
               macro = as.factor(all.gpa$mand$spec[tr$tip.label, "Macrohabitat"]),
               micro = as.factor(all.gpa$mand$spec[tr$tip.label, "Microhabitat"]),
               svl = log(all.gpa$mand$spec[tr$tip.label, "SVL"]),
               PC1 = phyloPCA$x[,1], PC2 = phyloPCA$x[,2], PC3 = phyloPCA$x[,3], 
               row.names = as.factor(all.gpa$mand$spec[tr$tip.label, "Species"]))

# Plot
Plot1 <- ggplot(df, aes(PC2, PC3, label = row.names(df),  colour = cont)) + geom_point(aes(size = svl))
Plot2 <- print(Plot1 
               + geom_text(color = "black", size = 0.5, check_overlap = F))

phylomorph.PC23<-phylomorphospace(tr, PCAvalues[,2:3], colors = "grey", label = "horizontal", lwd = 1, xlab = "PC2", ylab = "PC3")



###################################
# # # # 01c.phyPCA_versions # # # #
###################################

rm(list = ls())
library(geomorph)
library(phytools)

# Import superimposed data ###
load("all.gpa.bin")

# Import phylogeny 
load("phy.spec.MATCHED.bin")
tr <- phy.spec$tr

# Split into datasets and match with phylo
tr.drop <- tr$tip.label[which(tr$tip.label%in%dimnames(all.gpa$ear$gpa$coords)[[3]]==F)]
ear.tr <- drop.tip(tr, tr.drop)
ear.sh <- all.gpa$ear$gpa$coords[,,ear.tr$tip.label]
ear.cs <- all.gpa$ear$gpa$Csize[ear.tr$tip.label]
ear.spec <- all.gpa$ear$spec[ear.tr$tip.label,]

jaw.sh <- all.gpa$jaw$gpa$coords[,,tr$tip.label]
jaw.cs <- all.gpa$jaw$gpa$Csize[tr$tip.label]
jaw.spec <- all.gpa$jaw$spec[tr$tip.label,]

mand.sh <- all.gpa$mand$gpa$coords[,,tr$tip.label]
mand.cs <- all.gpa$mand$gpa$Csize[tr$tip.label]
mand.spec <- all.gpa$mand$spec[tr$tip.label,]

###
# PhyloMorphospace ###
# Traditional PCA (as in script 01a, but with reduced data to match phylo)
# Projection of phylogeny on morphospace
# Allows to visually review to what extent the data align with the phylogeny
# This will be formally tested later through phylogenetic signal analysis
pdf("phylomorphospace.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear, phylo = T, xlim = c(-0.15, 0.12), pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.9, par()$usr[4]*0.95, labels = "EAR", font = 2)






# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw,  phylo = T, xlim = c(-0.075, 0.15),pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand,  phylo = T, xlim = c(-0.1, 0.12),pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.8, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)

dev.off()




pcs.pca.ear <- pca.ear$x
pcs.pca.skull <- pca.mand$x
pcs.pca.jaw <- pca.jaw$x

ear.df <- data.frame(PC1 = pcs.pca.ear[,1], PC2 = pcs.pca.ear[,2], PC3 = pcs.pca.ear[,3],
                     cs = log(all.gpa$ear$gpa$Csize[ear.tr$tip.label]),
                     skullcs = log(all.gpa$mand$gpa$Csize[ear.tr$tip.label]),
                     groups = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Group"]),
                     cont = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Continent"]),
                     macro = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Macrohabitat"]),
                     micro = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Microhabitat"]),
                     svl = log(all.gpa$ear$spec[ear.tr$tip.label, "SVL"]),
                     row.names = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Species"]))

skull.df <- data.frame(PC1 = pcs.pca.skull[,1], PC2 = pcs.pca.skull[,2], PC3 = pcs.pca.skull[,3],
                     cs = log(all.gpa$mand$gpa$Csize[tr$tip.label]),
                     groups = as.factor(all.gpa$mand$spec[tr$tip.label, "Group"]),
                     cont = as.factor(all.gpa$mand$spec[tr$tip.label, "Continent"]),
                     macro = as.factor(all.gpa$mand$spec[tr$tip.label, "Macrohabitat"]),
                     micro = as.factor(all.gpa$mand$spec[tr$tip.label, "Microhabitat"]),
                     svl = log(all.gpa$mand$spec[tr$tip.label, "SVL"]),
                     row.names = as.factor(all.gpa$mand$spec[tr$tip.label, "Species"]))

jaw.df <- data.frame(PC1 = pcs.pca.jaw[,1], PC2 = pcs.pca.jaw[,2], PC3 = pcs.pca.jaw[,3],
                       cs = log(all.gpa$jaw$gpa$Csize[tr$tip.label]),
                       groups = as.factor(all.gpa$jaw$spec[tr$tip.label, "Group"]),
                       cont = as.factor(all.gpa$jaw$spec[tr$tip.label, "Continent"]),
                       macro = as.factor(all.gpa$jaw$spec[tr$tip.label, "Macrohabitat"]),
                       micro = as.factor(all.gpa$jaw$spec[tr$tip.label, "Microhabitat"]),
                       svl = log(all.gpa$jaw$spec[tr$tip.label, "SVL"]),
                       row.names = as.factor(all.gpa$jaw$spec[tr$tip.label, "Species"]))

Plot1 <- ggplot(jaw.df, aes(PC1, PC2, label = row.names(jaw.df), colour = micro, size = cs*2)) + geom_point(alpha = 1)
Plot2 <- print(Plot1 
               + geom_text(color = "black", size = 2, check_overlap = T))




###
# phyloPCA ###
# Rotates the data such as to reduce any effect of phylogeny
pdf("phyPCAs.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "EAR", font = 2)
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Continent))
legend("topleft", legend = levels(as.factor(ear.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green"), bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.ear, pch = 21, bg = as.factor(ear.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(ear.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "red"), bty = "n")

# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw)
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Continent))
legend("topleft", legend = levels(as.factor(jaw.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "cornflowerblue"), 
       bty = "n")
plot(pca.jaw, pch = 21, bg = as.factor(jaw.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(jaw.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr, GLS = T)
layout(matrix(1:4, nrow = 2))
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand)
text(par()$usr[1]*-0.99, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Continent))
legend("topleft", legend = levels(as.factor(mand.spec$Continent)), pch = 21, 
       pt.bg = c("black", "red", "green", "cornflowerblue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Macrohabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Macrohabitat)), 
       pch = 21, pt.bg = c("black", "red", "green", "blue"), bty = "n")
plot(pca.mand, pch = 21, bg = as.factor(mand.spec$Microhabitat))
legend("topleft", legend = levels(as.factor(mand.spec$Microhabitat)), 
       pch = 21, pt.bg = c("black", "red", "green"), bty = "n")

dev.off()

###
# PaCA ###
# Aligns the data to phylogenetic signal, rather than maxVAR
# In a way, it allows to visualize the phylogenetic component of shape variation

pdf("PaCA.pdf")
# Ear
pca.ear <- gm.prcomp(ear.sh, phy = ear.tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.ear, phylo = T, xlim = c(-0.15, 0.12), pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.9, par()$usr[4]*0.95, labels = "EAR", font = 2)

# Jaw
pca.jaw <- gm.prcomp(jaw.sh, phy = tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.jaw,  phylo = T, xlim = c(-0.025, 0.05),pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.85, par()$usr[4]*0.95, labels = "JAW", font = 2)

# Mandible
pca.mand <- gm.prcomp(mand.sh, phy = tr, align.to.phy = T)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.25, 0.35, 0),
    cex.axis = 0.7, font.lab = 2)
plot(pca.mand,  phylo = T, xlim = c(-0.075, 0.12),pch = 21, bg = as.factor(mand.spec$Group),
     phylo.par = list(anc.states = F,
                      edge.color = "grey 35",
                      edge.width = 1.5,
                      tip.txt.cex = 0.7,
                      node.txt.col = "n"))
text(par()$usr[1]*0.8, par()$usr[4]*0.95, labels = "MANDIBLE", font = 2)

dev.off()




nodes.temp <- unique( tr$edge[,1] ) #cambiar arbol de origen

tr$node.label <- nodes.temp #le damos nombre a los nodos para que sepa reconocerlos al hacer los analisis

eigendirections.pca.paca <- function( data , phy ) {
  
  coords.matrix <- two.d.array(data) 
  nodes.temp <- unique( phy$edge[,1] )
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy )
  names( descendants.temp ) <- nodes.temp
  discarded.descendants <- as.vector(names(which(lapply(descendants.temp, length)<5)))
  descendants.temp.mod <- descendants.temp[names(descendants.temp) %in% discarded.descendants == FALSE]
  node.data.list <- lapply( descendants.temp.mod , FUN = subset.rows , data = coords.matrix )
  
  
  vectors.matrix1 <- matrix( ncol = length(node.data.list) , nrow = 60) #cambiar nrow -> skull = 363, jaw = 60, ear = 492
  
  for (i in 1:length(node.data.list)) {
    node1.data <- node.data.list[[i]]
    GM.node1<-arrayspecs(as.matrix(node1.data), dim(as.matrix(node1.data))[2]/3, 3)
    pc.node1.vector<-gm.prcomp(GM.node1)$rotation[,1] #cambiar [,x] para definir con qué pc hacerlo
    name.node1 <- names(node.data.list)[i]
    
    vectors.matrix1[,i] <- pc.node1.vector
  }
  
  colnames(vectors.matrix1) <- names(node.data.list)
  
  
  vectors.matrix2 <- matrix( ncol = length(node.data.list) , nrow = 60) #cambiar nrow -> skull = 363, jaw = 60, ear = 492
  
  for (i in 1:length(node.data.list)) {
    node2.data <- node.data.list[[i]]
    GM.node2<-arrayspecs(as.matrix(node2.data), dim(as.matrix(node2.data))[2]/3, 3)
    data.names <- dimnames(GM.node2)[[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    phy.prunned <- drop.tip( phy , drop.taxa )  
    paca.node2.vector<-gm.prcomp(GM.node2, phy = phy.prunned, align.to.phy = TRUE, GLS = TRUE)$rotation[,1] #cambiar [,x] para definir con qué pc hacerlo
    name.node2 <- names(node.data.list)[i]
    
    vectors.matrix2[,i] <- paca.node2.vector      
  }
  
  colnames(vectors.matrix2) <- names(node.data.list)
  
  angles.vectors <- matrix(ncol = 1 , nrow = length(node.data.list))  
  
  for (x in 1:length(node.data.list)) {
    angles.vectors[x,] <- angle( vectors.matrix1[,x] , vectors.matrix2[,x], degree = T)
  }  
  
  rownames(angles.vectors) <- names(node.data.list)
  
  #this block makes angles all angles below 90 degrees: because PLS vectors are arbitrary in direction the maximum angle is orthogonal	
  pca.paca.angles.mod <- as.numeric(angles.vectors)
  names(pca.paca.angles.mod) <- row.names(angles.vectors)
  pca.paca.angles.mod.wrong <-  pca.paca.angles.mod [which(  pca.paca.angles.mod > 90 )]
  pca.paca.angles.mod.corrected <- (pca.paca.angles.mod.wrong -180) * -1
  pca.paca.angles.mod[names(pca.paca.angles.mod.corrected)] <- pca.paca.angles.mod.corrected
  angles.vectors <- pca.paca.angles.mod
  
  angles.vectors 
  
  
}  #modificar algunos parametros dentro de la funcion  #### FUNCIONAAAAAAAAA VAAAAMOOOOOS

  angles.1.pca.paca.skull <- eigendirections.pca.paca( data = all.gpa$mand$gpa$coords[,,tr$tip.label] , phy = tr)

  angles.1.pca.paca.ear <- eigendirections.pca.paca( data = all.gpa$ear$gpa$coords[,,ear.tr$tip.label] , phy = ear.tr)
  
  angles.1.pca.paca.jaw <- eigendirections.pca.paca( data = all.gpa$jaw$gpa$coords[,,tr$tip.label] , phy = tr)


# generates colours for angles 

edge.values <- angles.1.pca.paca.jaw
colours.rgb <- c("#ff6200","#fff367","#e8e8e8") #naranja rojizo-amarillo-gris ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.6) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.angles<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.angles[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

tree.temp <- tr
variable = angles.1.pca.paca.jaw

#tree with pca-paca.angles

plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T)
nodelabels( node = as.numeric(names(angles.1.pca.paca.jaw)), frame = "n" , cex = 3 , pch = 21 , bg = BG.angles )
angles.1.pca.paca.jaw








########################
# # # # 04.pGLSs # # # #
########################


                    
ear.gdf.pgls <- geomorph.data.frame(shape = all.gpa$ear$gpa$coords[,,ear.tr$tip.label], 
                                   cs = log(all.gpa$ear$gpa$Csize[ear.tr$tip.label]),
                                   skullcs = log(all.gpa$mand$gpa$Csize[ear.tr$tip.label]),
                                   groups = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Group"]),
                                   cont = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Continent"]),
                                   macro = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Macrohabitat"]),
                                   micro = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Microhabitat"]),
                                   svl = log(all.gpa$ear$spec[ear.tr$tip.label, "SVL"]),
                                   phy = ear.tr,
                                   row.names = as.factor(all.gpa$ear$spec[ear.tr$tip.label, "Species"]))
                    
jaw.gdf.pgls <- geomorph.data.frame(shape = all.gpa$jaw$gpa$coords[,,tr$tip.label], 
                                   cs = log(all.gpa$jaw$gpa$Csize[tr$tip.label]), 
                                   skullcs = log(all.gpa$mand$gpa$Csize[tr$tip.label]),
                                   groups = as.factor(all.gpa$jaw$spec[tr$tip.label, "Group"]),
                                   cont = as.factor(all.gpa$jaw$spec[tr$tip.label, "Continent"]),
                                   macro = as.factor(all.gpa$jaw$spec[tr$tip.label, "Macrohabitat"]),
                                   micro = as.factor(all.gpa$jaw$spec[tr$tip.label, "Microhabitat"]),
                                   svl = log(all.gpa$jaw$spec[tr$tip.label, "SVL"]),
                                   phy = tr,
                                   row.names = as.factor(all.gpa$jaw$spec[tr$tip.label, "Species"]))
                    
mand.gdf.pgls <- geomorph.data.frame(shape = all.gpa$mand$gpa$coords[,,tr$tip.label], 
                                    cs = log(all.gpa$mand$gpa$Csize[tr$tip.label]), 
                                    earcs = log(all.gpa$ear$gpa$Csize[tr$tip.label]),
                                    groups = as.factor(all.gpa$mand$spec[tr$tip.label, "Group"]),
                                    cont = as.factor(all.gpa$mand$spec[tr$tip.label, "Continent"]),
                                    macro = as.factor(all.gpa$mand$spec[tr$tip.label, "Macrohabitat"]),
                                    micro = as.factor(all.gpa$mand$spec[tr$tip.label, "Microhabitat"]),
                                    svl = log(all.gpa$mand$spec[tr$tip.label, "SVL"]),
                                    phy = tr,
                                    row.names = as.factor(all.gpa$mand$spec[tr$tip.label, "Species"]))

### pGLSs
# allometry
pgls.size.cs <- procD.pgls(shape ~ svl, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.size.cs)
    plot.pgls.size.cs <- plot(pgls.size.cs, type = "regression", reg.type = "RegScore", predictor = mand.gdf.pgls$cs)
    RegScores.skullCS<-plot.pgls.size.cs$RegScore
    #comparisons.allometry<-pairwise(pgls.skullCS.sinpremax.unique, groups = gdf.pgls.sinpmx$clades.b, covariate = gdf.pgls.sinpmx$log.brainconvhull)
    #summ.comparisons.allom<-summary(comparisons.allometry, test.type = "dist") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
    #summ.comparisons.allom
    pred.pgls.skullCS<-plot(pgls.size.cs, type="regression", reg.type="PredLine", predictor = mand.gdf.pgls$cs)
    PredLine.skullCS <- pred.pgls.skullCS$PredLine  
    reg.lines.x <- as.vector(pred.pgls.skullCS$plot.args$x)
    names(reg.lines.x) <- mand.gdf.pgls$row.names
    reg.lines.y <- pred.pgls.skullCS$plot.args$y
    names(reg.lines.y) <- mand.gdf.pgls$row.names
    df.skullCS.plots<-data.frame(RegScores.skullCS = RegScores.skullCS, 
                                 log.skullCS = mand.gdf.pgls$cs,
                                 groups = mand.gdf.pgls$groups,
                                 reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
    rownames(df.skullCS.plots) <- row.names(mand.gdf.pgls)
    Plot1 <- ggplot(df.skullCS.plots, aes(log.skullCS, RegScores.skullCS, label = rownames(df.skullCS.plots), colour = groups))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
    Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))
    
pgls.size.skullcs <- procD.pgls(shape ~ earcs, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.size.skullcs) 
    plot.pgls.size.earcs <- plot(pgls.size.earcs, type = "regression", reg.type = "RegScore", predictor = mand.gdf.pgls$earcs)
pgls.size.svl <- procD.pgls(shape ~ svl, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.size.svl)
  
    preds.allometry<-shape.predictor(pgls.size.cs$GM$pgls.fitted, x = plot.pgls.size.cs$RegScore, Intercept = F, 
                                     predmin = min(plot.pgls.size.cs$RegScore), predmax = max(plot.pgls.size.cs$RegScore))
      
      mesh.original.skull <- read.ply("capensis_Skull_clean.ply", ShowSpecimen = F, addNormals = F ) #asegurarse de que no tiene el binary encoding marcado en MeshLab
      coords.original.skull <- all.gpa$mand$gpa$coords[,,"Lygodactylus_capensis"]
      M.allometry <- pgls.size.cs$GM$pgls.mean
      
      library(hot.dots)
      library(tibble)
      
        per_lm_variance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
        {
        variances <- rowSums(apply(shape.data, c(1, 2), var))
        cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
        cols <- cols1(100)
        x = (log10(variances))
        xlims <- NULL
        tol <- 1e-06
        xlims <- range(x) + c(-tol, tol)
        nbin = 100
        breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
        whichColor <- function(p, cols, breaks) {
          i <- 1
          while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
              1
          cols[i]
        }
        variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
        variance_table <- tibble(Per_Lm_Variance = variances, Log_Variance = x,Variance_Colors = variancecolors)
        return(variance_table)
      }
        my.variances.allometry <- per_lm_variance(shape.data = pgls.size.cs$GM$pgls.fitted)
  
      mesh.mean <- tps3d(mesh.original.skull, refmat = coords.original.skull, tarmat = M.allometry)
      shade3d(mesh.mean, col= 8, alpha = 0.3)
      
      open3d()
      par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
      Sys.sleep(1)
      mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
      
      mesh.min <- tps3d(mesh.original.skull, refmat = coords.original.skull, tarmat = preds.allometry$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
      shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
      spheres3d(preds.allometry$predmin, col = my.variances.allometry$Variance_Colors, radius = 0.01) #radius cambia el tamaño del LM
      
      next3d()
      mesh.max <- tps3d(mesh.original.skull, refmat = coords.original.skull, tarmat = preds.allometry$predmax, lambda = 0.5)
      shade3d(mesh.max, col= 8, alpha = 0.3) 
      spheres3d(preds.allometry$predmax, col = my.variances.allometry$Variance_Colors, radius = 0.01)
      
      
  

pgls.sizes.model1 <- procD.pgls(shape ~ cs + earcs, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.sizes.model1)
pgls.sizes.model2 <- procD.pgls(shape ~ cs + svl, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.sizes.model2)
pgls.sizes.model3 <- procD.pgls(shape ~ earcs + svl, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.sizes.model3) 
pgls.sizes.model4 <- procD.pgls(shape ~ cs + earcs + svl, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.sizes.model4) 
  
# ecology & geography
pgls.ecolgeogr.cont <- procD.pgls(shape ~ cont, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.cont)
pgls.ecolgeogr.macro <- procD.pgls(shape ~ macro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.macro) 
pgls.ecolgeogr.micro <- procD.pgls(shape ~ micro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.micro)  
  pgls.ecolgeogr.microgroups <- procD.pgls(shape ~ micro*groups, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
    summary(pgls.ecolgeogr.microgroups)
      microgroups<-interaction(mand.gdf.pgls$micro, mand.gdf.pgls$groups)
        PW1<-pairwise(pgls.ecolgeogr.microgroups, groups = microgroups)
          summary(PW1,test.type = "dist", confidence = 0.95, stat.table = TRUE)
  pgls.ecolgeogr.micromacro <- procD.pgls(shape ~ micro*macro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
    summary(pgls.ecolgeogr.micromacro)
      micromacro<-interaction(mand.gdf.pgls$micro, mand.gdf.pgls$macro)
        PW2<-pairwise(pgls.ecolgeogr.micromacro, groups = micromacro)
          summary(PW2,test.type = "dist", confidence = 0.95, stat.table = TRUE)
  pgls.ecolgeogr.microcont <- procD.pgls(shape ~ micro*cont, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
    summary(pgls.ecolgeogr.microcont)
      microcont<-interaction(mand.gdf.pgls$micro, mand.gdf.pgls$cont)
        PW3<-pairwise(pgls.ecolgeogr.microcont, groups = microcont)
          summary(PW3,test.type = "dist", confidence = 0.95, stat.table = TRUE)
  
pgls.ecolgeogr.model1 <- procD.pgls(shape ~ cont * macro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.model1)
pgls.ecolgeogr.model2 <- procD.pgls(shape ~ cont * micro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.model2) 
pgls.ecolgeogr.model3 <- procD.pgls(shape ~ macro * micro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.model3)
pgls.ecolgeogr.model4 <- procD.pgls(shape ~ cont + macro + micro, phy = phy, data = ear.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolgeogr.model4)   

# combined  
pgls.ecolbysize.model1 <- procD.pgls(shape ~ micro + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model1)                    
pgls.ecolbysize.model2 <- procD.pgls(shape ~ macro + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model2)
pgls.ecolbysize.model3 <- procD.pgls(shape ~ cont + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model3)
pgls.ecolbysize.model3.5 <- procD.pgls(shape ~ groups + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model3.5)
  
pgls.ecolbysize.model4 <- procD.pgls(shape ~ micro + earcs:cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model4)
pgls.ecolbysize.model5 <- procD.pgls(shape ~ macro + micro + earcs*cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model5)
pgls.ecolbysize.model6 <- procD.pgls(shape ~ macro + micro + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model6)  
pgls.ecolbysize.model7 <- procD.pgls(shape ~ micro:earcs + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model7)
pgls.ecolbysize.model8 <- procD.pgls(shape ~ micro:macro + cs, phy = phy, data = mand.gdf.pgls, iter = 9999, SS.type = "II")
  summary(pgls.ecolbysize.model8) 
  
  
###########################################################################
# # # # 05.DISPARITY, DTT, PHYLOGENETIC SIGNAL & EVOLUTIONARY RATES # # # #
###########################################################################

# DISPARITY
  
  disparity.skull.groups <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ groups, groups = ~ groups, data = mand.gdf.pgls, iter=999)
    summary(disparity.skull.groups)
  disparity.skull.cont <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ cont, groups = ~ cont, data = mand.gdf.pgls, iter=999)
    summary(disparity.skull.cont)  
  disparity.skull.macro <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ macro, groups = ~ macro, data = mand.gdf.pgls, iter=999)
    summary(disparity.skull.macro)
  disparity.skull.micro <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ micro, groups = ~ micro, data = mand.gdf.pgls, iter=999)
    summary(disparity.skull.micro)
    disparity.skull.microgroups <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ micro*groups, groups = ~ micro*groups, data = mand.gdf.pgls, iter=9999, stat.table=TRUE)
    summary(disparity.skull.microgroups, formula = TRUE)
    disparity.skull.micromacro <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ micro*macro, groups = ~ micro*macro, data = mand.gdf.pgls, iter=9999, stat.table=TRUE)
    summary(disparity.skull.micromacro)
    disparity.skull.microcont <- morphol.disparity(all.gpa$mand$gpa$coords[,,tr$tip.label] ~ micro*cont, groups = ~ micro*cont, data = mand.gdf.pgls, iter=9999, stat.table=TRUE)
    summary(disparity.skull.microcont)
  
  disparity.jaw.groups <- morphol.disparity(all.gpa$jaw$gpa$coords[,,tr$tip.label] ~ groups, groups = ~ groups, data = jaw.gdf.pgls, iter=999)
    summary(disparity.jaw.groups)
  disparity.jaw.cont <- morphol.disparity(all.gpa$jaw$gpa$coords[,,tr$tip.label] ~ cont, groups = ~ cont, data = jaw.gdf.pgls, iter=999)
    summary(disparity.jaw.cont)  
  disparity.jaw.macro <- morphol.disparity(all.gpa$jaw$gpa$coords[,,tr$tip.label] ~ macro, groups = ~ macro, data = jaw.gdf.pgls, iter=999)
    summary(disparity.jaw.macro)
  disparity.jaw.micro <- morphol.disparity(all.gpa$jaw$gpa$coords[,,tr$tip.label] ~ micro, groups = ~ micro, data = jaw.gdf.pgls, iter=999)
    summary(disparity.jaw.micro)
    
  disparity.ear.groups <- morphol.disparity(all.gpa$ear$gpa$coords[,,ear.tr$tip.label] ~ groups, groups = ~ groups, data = ear.gdf.pgls, iter=999)
    summary(disparity.ear.groups)
  disparity.ear.cont <- morphol.disparity(all.gpa$ear$gpa$coords[,,ear.tr$tip.label] ~ cont, groups = ~ cont, data = ear.gdf.pgls, iter=999)
    summary(disparity.ear.cont)  
  disparity.ear.macro <- morphol.disparity(all.gpa$ear$gpa$coords[,,ear.tr$tip.label] ~ macro, groups = ~ macro, data = ear.gdf.pgls, iter=999)
    summary(disparity.ear.macro)
  disparity.ear.micro <- morphol.disparity(all.gpa$ear$gpa$coords[,,ear.tr$tip.label] ~ micro, groups = ~ micro, data = ear.gdf.pgls, iter=999)
    summary(disparity.ear.micro)
    
  
# DTT  (Disparity Through Time)
  matrix.coords.skull <- two.d.array(all.gpa$mand$gpa$coords[,,tr$tip.label])           
DTT.skull <- dtt(phy = tr, data = matrix.coords.skull, index = c("avg.sq"), mdi.range = c(0,30), nsim = 1000, CI = 0.95, plot = TRUE)      
  
  matrix.coords.ear <- two.d.array(all.gpa$ear$gpa$coords[,,ear.tr$tip.label])           
DTT.ear <- dtt(phy = ear.tr, data = matrix.coords.ear, index = c("avg.sq"), mdi.range = c(0,30), nsim = 1000, CI = 0.95, plot = TRUE) 

  matrix.coords.jaw <- two.d.array(all.gpa$jaw$gpa$coords[,,tr$tip.label])           
DTT.jaw <- dtt(phy = tr, data = matrix.coords.jaw, index = c("avg.sq"), mdi.range = c(0,30), nsim = 1000, CI = 0.95, plot = TRUE)


# K multivariate
physignal(all.gpa$mand$gpa$coords[,,tr$tip.label], phy = tr, iter = 9999)
physignal(all.gpa$ear$gpa$coords[,,ear.tr$tip.label], phy = ear.tr, iter = 9999)
physignal(all.gpa$jaw$gpa$coords[,,tr$tip.label], phy = tr, iter = 9999)

# K univariate
phylosig(all.gpa$mand$gpa$Csize[tr$tip.label], tree = tr, method ="K", test = T, nsim = 9999)
phylosig(all.gpa$ear$gpa$Csize[ear.tr$tip.label], tree = ear.tr, method ="K", test = T, nsim = 9999)
phylosig(all.gpa$jaw$gpa$Csize[tr$tip.label], tree = tr, method ="K", test = T, nsim = 9999)


# Rates of shape evolution
  # skull
  data.spec.tr <- all.gpa$mand$spec[tr$tip.label,]
gp.skull.group <- factor(data.spec.tr$Group)
  names(gp.skull.group) <- data.spec.tr$Species
gp.skull.cont <- factor(data.spec.tr$Continent)
  names(gp.skull.cont) <- data.spec.tr$Species
gp.skull.macro <- factor(data.spec.tr$Macrohabitat)
  names(gp.skull.macro) <- data.spec.tr$Species
gp.skull.micro <- factor(data.spec.tr$Microhabitat)
  names(gp.skull.micro) <- data.spec.tr$Species
  
ER.skull.groups <- compare.evol.rates(A = all.gpa$mand$gpa$coords , phy = tr , method = "simulation" , gp = gp.skull.group , iter = 999)
  summary(ER.skull.groups)
    plot(ER.skull.groups)
      ER.skullCS.groups <- compare.evol.rates(A = all.gpa$mand$gpa$Csize , phy = tr , method = "simulation" , gp = gp.skull.group , iter = 999)
        summary(ER.skullCS.groups)
ER.skull.cont <- compare.evol.rates(A = all.gpa$mand$gpa$coords , phy = tr , method = "simulation" , gp = gp.skull.cont , iter = 999)
  summary(ER.skull.cont)
    plot(ER.skull.cont)
      ER.skullCS.cont <- compare.evol.rates(A = all.gpa$mand$gpa$Csize , phy = tr , method = "simulation" , gp = gp.skull.cont , iter = 999)
        summary(ER.skullCS.cont)
ER.skull.macro <- compare.evol.rates(A = all.gpa$mand$gpa$coords , phy = tr , method = "simulation" , gp = gp.skull.macro , iter = 999)
  summary(ER.skull.macro)
    plot(ER.skull.macro)
      ER.skullCS.macro <- compare.evol.rates(A = all.gpa$mand$gpa$Csize , phy = tr , method = "simulation" , gp = gp.skull.macro , iter = 999)
        summary(ER.skullCS.macro)
ER.skull.micro <- compare.evol.rates(A = all.gpa$mand$gpa$coords , phy = tr , method = "simulation" , gp = gp.skull.micro , iter = 999)
  summary(ER.skull.micro)
    plot(ER.skull.micro)
      ER.skullCS.micro <- compare.evol.rates(A = all.gpa$mand$gpa$Csize , phy = tr , method = "simulation" , gp = gp.skull.micro , iter = 999)
        summary(ER.skullCS.micro)
 
  # labyrinth
  data.spec.tr <- all.gpa$ear$spec[ear.tr$tip.label,]
gp.ear.group <- factor(data.spec.tr$Group)
  names(gp.ear.group) <- data.spec.tr$Species
gp.ear.cont <- factor(data.spec.tr$Continent)
  names(gp.ear.cont) <- data.spec.tr$Species
gp.ear.macro <- factor(data.spec.tr$Macrohabitat)
  names(gp.ear.macro) <- data.spec.tr$Species
gp.ear.micro <- factor(data.spec.tr$Microhabitat)
  names(gp.ear.micro) <- data.spec.tr$Species
    
ER.ear.groups <- compare.evol.rates(A = all.gpa$ear$gpa$coords , phy = ear.tr , method = "simulation" , gp = gp.ear.group , iter = 999)
  summary(ER.ear.groups)
    plot(ER.ear.groups)
    ER.earCS.groups <- compare.evol.rates(A = all.gpa$ear$gpa$Csize , phy = ear.tr , method = "simulation" , gp = gp.ear.group , iter = 999)
    summary(ER.earCS.groups)
ER.ear.cont <- compare.evol.rates(A = all.gpa$ear$gpa$coords , phy = ear.tr , method = "simulation" , gp = gp.ear.cont , iter = 999)
  summary(ER.ear.cont)
    plot(ER.ear.cont)
    ER.earCS.cont <- compare.evol.rates(A = all.gpa$ear$gpa$Csize , phy = ear.tr , method = "simulation" , gp = gp.ear.cont , iter = 999)
    summary(ER.earCS.cont)
ER.ear.macro <- compare.evol.rates(A = all.gpa$ear$gpa$coords , phy = ear.tr , method = "simulation" , gp = gp.ear.macro , iter = 999)
  summary(ER.ear.macro)
    plot(ER.ear.macro)
    ER.earCS.macro <- compare.evol.rates(A = all.gpa$ear$gpa$Csize , phy = ear.tr , method = "simulation" , gp = gp.ear.macro , iter = 999)
    summary(ER.earCS.macro)
ER.ear.micro <- compare.evol.rates(A = all.gpa$ear$gpa$coords , phy = ear.tr , method = "simulation" , gp = gp.ear.micro , iter = 999)
  summary(ER.ear.micro)
    plot(ER.ear.micro)
    ER.earCS.micro <- compare.evol.rates(A = all.gpa$ear$gpa$Csize , phy = ear.tr , method = "simulation" , gp = gp.ear.micro , iter = 999)
    summary(ER.earCS.micro)

# Mandible
  data.spec.tr <- all.gpa$jaw$spec[tr$tip.label,]
gp.jaw.group <- factor(data.spec.tr$Group)
  names(gp.jaw.group) <- data.spec.tr$Species
gp.jaw.cont <- factor(data.spec.tr$Continent)
  names(gp.jaw.cont) <- data.spec.tr$Species
gp.jaw.macro <- factor(data.spec.tr$Macrohabitat)
  names(gp.jaw.macro) <- data.spec.tr$Species
gp.jaw.micro <- factor(data.spec.tr$Microhabitat)
  names(gp.jaw.micro) <- data.spec.tr$Species

ER.jaw.groups <- compare.evol.rates(A = all.gpa$jaw$gpa$coords , phy = tr , method = "simulation" , gp = gp.jaw.group , iter = 999)
  summary(ER.jaw.groups)
    plot(ER.jaw.groups)
    ER.jawCS.groups <- compare.evol.rates(A = all.gpa$jaw$gpa$Csize , phy = tr , method = "simulation" , gp = gp.jaw.group , iter = 999)
    summary(ER.jawCS.groups)
ER.jaw.cont <- compare.evol.rates(A = all.gpa$jaw$gpa$coords , phy = tr , method = "simulation" , gp = gp.jaw.cont , iter = 999)
  summary(ER.jaw.cont)
    plot(ER.jaw.cont)
ER.jaw.macro <- compare.evol.rates(A = all.gpa$jaw$gpa$coords , phy = tr , method = "simulation" , gp = gp.jaw.macro , iter = 999)
  summary(ER.jaw.macro)
    plot(ER.jaw.macro)
ER.jaw.micro <- compare.evol.rates(A = all.gpa$jaw$gpa$coords , phy = tr , method = "simulation" , gp = gp.jaw.micro , iter = 999)
  summary(ER.jaw.micro)
    plot(ER.jaw.micro)
    

###########################################################
# # # # 06.INTRASPECIFIC VS INTERSPECIFIC VARIATION # # # #
###########################################################

coords.inter<-two.d.array(all.gpa$mand$gpa$coords)
      write.csv(coords.inter, 'coords.inter.csv')
coords.intra<-two.d.array(all.gpa$mand$gpa$coords)
      write.csv(coords.intra, 'coords.intra.csv')
        
load("intra.gpa.bin")
load("intra.dt.bin")
mand.spec <- all.gpa$mand$spec
write.csv(mand.spec, 'skull.clasif.intra.csv')


inter.intra.together<-read.csv("coords.all.csv",header=T,sep=";",row.names=1)
coords.all<-arrayspecs(inter.intra.together[,7:369],121,3)
gpa.interintra<-gpagen(coords.all)
coords.inter<-gpa.interintra$coords[,,1:79]
coords.intra<-gpa.interintra$coords[,,80:90]

mand.gdf.interintra <- geomorph.data.frame(sh = gpa.interintra$coords,
                                           cs = inter.intra.together[,6],
                                           groups = as.factor(inter.intra.together[,1]),
                                           cont = as.factor(inter.intra.together[,2]),
                                           macro = as.factor(inter.intra.together[,3]),
                                           micro = as.factor(inter.intra.together[,4]),
                                           interintra = as.factor(inter.intra.together[,5]))

allom.interintra.lm <- procD.lm(sh ~ cs*interintra, SS.type = "I", data = mand.gdf.interintra, iter = 9999)
  anova(allom.interintra.lm)
   pw.allom.lm.groups <- pairwise(allom.interintra.lm, groups = mand.gdf.interintra$interintra, covariate = mand.gdf.interintra$cs)
    summ.pw.allom<-summary(pw.allom.lm.groups, confidence = 0.95, test.type = "VC", angle.type = "deg") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
      summ.pw.allom


mand.gdf.inter <- geomorph.data.frame(sh.inter = coords.inter[,,tr$tip.label],
                                     cs.inter = inter.intra.together[tr$tip.label,5],
                                     groups.inter = as.factor(inter.intra.together[tr$tip.label,1]),
                                     cont.inter = as.factor(inter.intra.together[tr$tip.label,2]),
                                     macro.inter = as.factor(inter.intra.together[tr$tip.label,3]),
                                     micro.inter = as.factor(inter.intra.together[tr$tip.label,4]),
                                     phy = tr)

pgls.size.csinter <- procD.pgls(sh.inter ~ cs.inter, phy = phy, data = mand.gdf.inter, iter = 9999, SS.type = "II")
summary(pgls.size.csinter)


mand.gdf.intra <- geomorph.data.frame(sh.intra = coords.intra,
                                      cs.intra = inter.intra.together[80:90,5],
                                      groups.intra = as.factor(inter.intra.together[80:90,1]),
                                      cont.intra = as.factor(inter.intra.together[80:90,2]),
                                      macro.intra = as.factor(inter.intra.together[80:90,3]),
                                      micro.intra = as.factor(inter.intra.together[80:90,4]))

allom.intra.lm.commgpa <- procD.lm(sh.intra ~ cs.intra, SS.type = "I", data = mand.gdf.intra)
anova(allom.intra.lm.commgpa)

M.gpainterintra <-rbind(coef.inter <- pgls.size.csinter$LM$gls.coefficients[2,],
                        coef.intra <- allom.intra.lm.commgpa$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(M.gpainterintra))*180/pi




rdf <- rrpp.data.frame(cs = log(all.gpa$mand$gpa$Csize), shape = all.gpa$mand$gpa$coords,
                       macro = as.factor(all.gpa$mand$spec$Macrohabitat), micro = as.factor(all.gpa$mand$spec$Microhabitat), 
                       cont = as.factor(all.gpa$mand$spec$Continent), row.names = as.factor(all.gpa$mand$spec$Species))

mand.df <- geomorph.data.frame(sh = all.gpa$mand$gpa$coords, 
                               cs = log(all.gpa$mand$gpa$Csize), 
                               cont = as.factor(all.gpa$mand$spec$Continent),
                               macro = as.factor(all.gpa$mand$spec$Macrohabitat),
                               micro = as.factor(all.gpa$mand$spec$Microhabitat),
                               svl = log(all.gpa$mand$spec$SVL))


# 1: Evolutionary Allometry

allom.intra <- lm.rrpp(shape~cs, data = rdf)
anova(allom.ind)

allom.intra.lm <- procD.lm(sh ~ cs, SS.type = "I", data = mand.df)
anova(allom.intra.lm)


M <-rbind(coef.inter <- pgls.size.cs$LM$gls.coefficients[2,],
          coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(M))*180/pi

Mmicro <-rbind(coef.inter <- pgls.ecolbysize.model1$LM$gls.coefficients[2,],
          coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(Mmicro))*180/pi

Mmacro <-rbind(coef.inter <- pgls.ecolbysize.model2$LM$gls.coefficients[2,],
               coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(Mmacro))*180/pi

Mcont <-rbind(coef.inter <- pgls.ecolbysize.model3$LM$gls.coefficients[2,],
               coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(Mcont))*180/pi

Mspp <-rbind(coef.inter <- pgls.ecolbysize.model3.5$LM$gls.coefficients[2,],
              coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(Mspp))*180/pi



M.gpainterintra <-rbind(coef.inter <- pgls.ecolbysize.model3.5$LM$gls.coefficients[2,],
             coef.intra <- allom.intra.lm$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(Mspp))*180/pi
