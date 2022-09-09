require(treeio)
require(ggtree)
require(dplyr)
source("clean_functions.R") # contains no_()

beast218 <- read.beast("./BEAST2_analysis/beast_218_tree_asterophryinae.nex")
iqtree218 <- read.iqtree("./IQTREE2_analysis/iqtree_218_timetree_asterophryinae.contree")
d <- read.csv("Table1.csv")
if(dir.exists("out")!=TRUE) dir.create("out") # check if output directory out exists, if false create

tree <- full_join(beast218, d, by="label") # 218 sp tree with metadata

################################################################
#### GeoDistance Between Sites using GPS
####
#### Calculate pairwise distances in km between sampling points of individual specimens, 
#### group all frogs within 1.1 km of each other into "sites"
#### INPUT treemetadata.csv
#### OUTPUT grouped by site to sppbysite.csv
#### By Diana Gao, Ethan Hill, Marguerite Butler
################################################################

require("geodist")
require("dplyr")

dd <- d[d$terrane != "notpng",]  # d is treemetadata in Table1.csv

dist <- geodist(dd[c("longitude", "latitude")], measure="geodesic")/1000     # pairwise distance matrix in km
dimnames(dist) <- list(dd$id, dd$id)           

ll <- apply( dist, 1, function(x) x[ x <1.1 ]  )							# list with all values within < 1.1 km of each other
llnames <- unique(lapply(ll, function(x) sort(names(x))))		# ids of all those frogs grouped by 1.1km distance
names(llnames) <- 1:length(llnames)

sites <- sapply(llnames, function(y)       			       # list of sites with the following data
				t(
					sapply(y, 
						function(x) dd[dd$id==x, c("id", "genus", "species", "site", "locality", "latitude", "longitude", "elevation")]
					)
				)
			)		# dropped elevation from this function
						
sites <- lapply(seq_along(sites),    # add site number to each dataframe, add blank line for padding
			function(i) { 
				rbind(
					cbind(site=i, sites[[i]] ), 
					""
					)
				}
		)						
				
sites_all <- do.call(rbind, sites)		# bind all the dataframes together into giant dataframe
write.csv(sites_all, file="out/spbysite.csv", row.names=F)			

######################## 
##     Process metadata for trees
##	  	 Create tip labels for support tree Fig 2a 
##         Add site for multisite species
##         Color code "oddballs"
##         Generate alternate tiplabels
##         d2 is a dataframe with tree metadata, used to join with tree for plotting fig 2 
########################

# alternative tip labels that could be used for plotting dataframe d2 = 218 spp tree:
#		gensp = genus species
# 		tiplabid = genus species id
#		tiplab2	= genus species-site if at multiple sites
# 		tiplab2id = tiplab2 plus id

d$gensp <- paste(d$genus, d$species, sep=" ")
d$tiplabid <- paste(d$gensp, d$id, sep=" ")

# d2 = 218 species metadata
d2 <- merge(data.frame("label"=get_taxa_name(ggtree(tree))), d, by="label")  # 218 species
od <- allduplicates(d2$gensp)

# paste -site if the species occurs at multiple sites
d2$tiplab2 <- paste(d2$gensp, ifelse(od, paste("-", d2$site, sep=""), ""), sep="")
d2$tiplab2id <- paste(d2$tiplab2, d2$id, sep=" ")

# oddball taxa that fall in weird places, odds is pointer to them
oddballs <- c("Paedophryne dekot", "Aphantophryne", "Copiula tyleri")
odds <- unlist(sapply(oddballs, grep, d2$tiplab2))
d2$tipcol <- "black"
d2$tipcol[odds] <- "blue"

genera <- unique(d2$genus) %w/o% c("Oxydactyla", "Sphenophryne", "Aphantophryne A", "Aphantophryne B")   # 17 recognized genera

################################################################
#### Gather BI/ML nodal Support in dataframe (for Figure 2a)
####
#### Consensus topology with both BI and ML support values, 218 species tree
#### read support values from BEAST and IQTREEs, combine into bi218phylo object
#### code support dots dataframes: bidots, mldots and genus
################################################################

require(phytools)
require(treeio)
require(tidytree)
require(ggtree)
require(ape)
require(ggplot2)

bi218phylo <- beast218 %>% beast2phylo  # beast2phylo is in clean_functions.R
bi218phylo <- ladderize(bi218phylo, right=F)
iqtree218phylo <- iqtree218 %>% as.phylo

### aligning the tree to make sure tips, nodes are in same order
obj <- cophylo(bi218phylo, iqtree218phylo)     
#cbind(obj$trees[[1]]$tip.label, obj$trees[[2]]$tip.label)   # one mismatch, row 133-134, use BEAST labels
#cbind(obj$trees[[1]]$node.label, obj$trees[[2]]$node.label)

########## Generate Support Value data.frame for plotting ########################
## supp = data.frame with columns: sup, bi, ml, bicol, mlcol
## BI in probability 0-1, ML in % 0-100, sup in bi/ml, bicol, mlcol by level 0<.7, .7<.8, .8<.9, .9-1  

bis <- obj$trees[[1]]$node.label
mls <- obj$trees[[2]]$node.label
supp <- data.frame(
					node=(Ntip(bi218phylo)+1):(Nedge(bi218phylo)+1), 
					supp=paste(bis, mls, sep = "/"), # BI/ML values
					bi=bis,  # Bayesian posterior probability
					ml=mls,  # Maximum Likelihood bootstraps    
													# Set colors according to support breaks
					bicol = as.character(
											cut(as.numeric(bis), 
											include.lowest=T, 
											breaks=c(.7,.8,.90,1), 
											labels=c("low", "med", "high"))),
					mlcol = as.character(
											cut(as.numeric(mls), 
											include.lowest=T, breaks=c(70,80,90,100), 
											labels=c("low", "med", "high")))
				)

write.csv(supp, file="out/tree_support_dataframe.csv", row.names=F)  # write support dataframe

bi218phylo <- obj$trees[[1]]    ## use BEAST tree - extracted from cophylo obj
bi218phylo$node.label <- supp$supp   ## paste combined support labels 
#write.tree(bi218phylo, file = "data/fig2.tree_biml_multisite.nwk")     ## write combined, cleaned tree

# find basal node for each genus to use in highlighting clades
# first make tipvec, dropping oddball species (probably misindentifications)
tipvec <- grep(paste(c(nospaces(oddballs),"Copiula_sp.7"), collapse="|"), bi218phylo$tip.label, value=T, invert=T) 
# gennode = vector containing node numbers at root of each genus (exclude oddballs)
gennode <- sapply(nospaces(genera), function(x)   # for each set of tips per genus, find MRCA
                            MRCA(bi218phylo,   
                              grep(x, 
								tipvec,
								value=T)
							)
		         )

##################################################
##	Join tree with metadata
##################################################

tree <- bi218phylo  # combined topology 218 tip tree

td <- full_join(bi218phylo, supp, by="node")
td <- full_join(td, d2[c("label", "gensp", "tiplabid", "tiplab2", "tiplab2id", "tipcol", "site", "lifestyle")], by="label")
tdat <- td %>% as_tibble %>% as.data.frame

##################################################
##	Plot Support tree with dumbbells (Figure 2a) 
##################################################

gc <- read.csv("gencolorABC.csv")    ## color code highlight boxes by genus
gcol <- gc$col
names(gcol) <- gc$gen         ## gcol holds the colors=values for the genera=names  

# plot no ID tree: (use label=tiplab2id for labels with IDs)
#		gensp = genus species
# 		tiplabid = genus species id
#		tiplab2	= genus species-site if at multiple sites
# 		tiplab2id = tiplab2 plus id

p <- ggtree(td) + 
	 geom_tiplab(aes(label=tiplab2, color=tipcol), 
	 					size=1.2, 
	 					offset=.05, 
	 					fontface='italic' ) + 
	 scale_color_manual(values = c("black", "blue")) +
     scale_x_continuous(limits = c(0, 30), 
     					breaks=c(-0.71, 4.71, 9.71, 14.71, 19.71),  
     					labels = c(20,15,10,5,0)) +
     coord_cartesian(clip="off") + 
	 theme_tree2(legend.position="none") 
	 
# test tree if you want to show node numbers
#test <- p + geom_text2(aes(label=node), hjust=-.2, size=1.5) +   
#    geom_text2(aes(label=sup), hjust=-1, size=1.5)                 

## adds genus highlight boxes  
q <- p + geom_hilight(node=gennode, fill=gcol[no_(names(gennode))], alpha=0.25)
## adds support dots

yl=10  # vertical placement adjustment for legend
sl=0  # size adjustment for legend
# overlay support dots and support legend
r <-  q + geom_nodepoint( aes(fill=bicol, subset= (as.numeric(bi)>=0.7)), 
 						size=1.27,
						shape=21, 				
						position=position_nudge(x=-.35, y=.8), 
						alpha=1, 
						na.rm=T,
						stroke=.4) +
	geom_nodepoint( aes(fill=mlcol, subset=(as.numeric(ml)>=70)), 
						size=1.27,
						shape=21, 
						position=position_nudge(x=-.35, y=-.8), 
						alpha=1, 
						na.rm=T,
						stroke=.4) +
	scale_fill_manual(values = c("high"="black", "low"="azure", "med"="skyblue")) +
    annotate("text", 0, 215+yl, hjust=.1, size=4+sl, label="Posterior Probability/" ) +
	annotate("text", 0, 211+yl, hjust=.09, size=4+sl, label="Bootstrap Support" ) +
	annotate("point", .1, 206+yl, shape=21, size=3, fill="black" ) +
	annotate("text", 1, 206+yl, hjust=0.01, size=3.5+sl, label="PP>0.9 / BS>90" ) +
	annotate("point", .1, 200+yl,  shape=21, size=3, fill="skyblue" ) +
	annotate("text", 1, 201.5+yl, hjust=0.01, size=3.5+sl, label="0.8<=PP<0.9") +
	annotate("text", 1, 198.5+yl, hjust=0.01, size=3.5+sl, label="80<=BS<90" ) +
	annotate("point", .1, 193+yl,  shape=21, size=3, fill="azure" ) +
	annotate("text", 1, 194+yl, hjust=0.01, size=3.5+sl, label="0.7<=PP<0.8") +
	annotate("text", 1, 191+yl, hjust=0.01, size=3.5+sl, label="70<=BS<80" ) 


pdf(file="out/Fig2a_BIMLsupporthilight10x4_noid.pdf", height=10, width=4)
   print(r)
dev.off()


##################################################
####
#### Lifestyle Analysis (Fig 2b)
####
#### lifestyle data is in td under "ecomorph"
##################################################

require(geiger)
require(ggimage)
lifestyle <- tdat$lifestyle     # data vector should be named w/ tip labels
names(lifestyle) <- tdat$label
lifestyle <- lifestyle[!is.na(lifestyle)]   # data for just tips, named by tiplabels

## The models below will take 10 min or so to run
##
# geigerER <- fitDiscrete(tree, ecomorph, model = "ER")
# geigerER 
# geigerER$opt$aicc  # AICc=320

# geigerARD <- fitDiscrete(tree, ecomorph, model = "ARD", ncores = 4, control = list(method = c("subplex", "L-BFGS-B"), niter = 100))
# geigerARD
# geigerARD$opt$aicc  # AICc=311

geigerSYM <- fitDiscrete(tree, lifestyle, model = "SYM", ncores = 4, control = list(method = c("subplex", "L-BFGS-B"), niter = 100))
# geigerSYM$opt$aicc   # AICc=310
			# SYM.simmap.summary$ace # has prob anc states
ecorate <- geiger:::.Qmatrix.from.gfit(geigerSYM) # extracts rate matrix from model
ecorate <- round(ecorate, 5)

## 2000 stochastic maps generated from geiger rate matrix
strees <- make.simmap(tree, lifestyle, model="SYM", nsim=2000, Q = ecorate)
SYM.simmap.summary <- summary(strees,plot=FALSE)

# print Geiger output to file
sink(file="out/lifestyle_analysis_geiger_output.txt")
  print("SIMMAP summary:")
  print(SYM.simmap.summary)
  print("Geiger best fit model SYM")
  print(geigerSYM)
  print("rate matrix =")
  print(ecorate)
sink()

# # Save Geiger outputs
# write.csv(ecorate, file="out/ecomorph_transition_matrix.csv", row.names=F)   # save ecomorph transition matrix
# write.csv(eco, file="out/ecomorph_nodes.csv", row.names=F)                   # save ecomorph - nodes
# save(geigerSYM, SYM.simmap.summary, ecorate, file="out/ecomorph_anc_geiger.Rdata")      # raw Geiger output

## Get ecomorph states for each node (including tips) for plotting

states <- apply(SYM.simmap.summary$ace, 1, function(x) names(x[x==max(x)]))   # majority rule for node states
       # the $ace object is a dataframe with one row per node (internal first, then tip nodes)

ntips <- Ntip(tree)
nodes <- tdat[c("node", "label")]
nnodes <- dim(nodes)[1]  # same as Ntip() + Nnode() internal nodes
nodes$label[(ntips+1):(nnodes)] <- nodes$node[(ntips+1):(nnodes)]
eco <- c(states[ntips:nnodes], states[1:(nnodes-ntips)]) # tree order: tips -> internal nodes

## make ancestral state pies
## plot phylogenies painted by ecomorphs and with pies

# "Arboreal", "Fossorial", "Scansorial", "Semi-Aquatic", "Terrestrial"
cols <- c("#00ff00", "#ffa500", "#ff1493", "#00bfff", "#0000ff")
names(cols) <- sort(unique(lifestyle))

ancstates <- data.frame(node=rownames(SYM.simmap.summary$ace),  SYM.simmap.summary$ace)
pies <- nodepie(ancstates, cols=2:6, alpha=.7)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

# Create pies only for nodes where the maximum ancestral state prob is less than 70%
pies70 <- ancstates[apply(ancstates, 1, function(x) max(x[-1])<.7),]

##
##	Plot Ecomorph tree with pies (fig 2b)
##
##  Use theme_get() to see theme options
e <- ggtree(tree, aes(color=eco), size=.5) +
     theme(legend.position=c(0.9,.92), 					# customize the legend placement, text, title size
     		legend.text=element_text(size=11), 
     		legend.title=element_text(size=14)
     		) + 
  	 scale_x_reverse() +
     guides(color = guide_legend(override.aes=list(size=2), title = "Lifestyles")) + # customize legend line size
 	 scale_colour_manual(values = cols) +
  	 geom_inset(pies[pies70$node], width = .25, height = .25, hjust = -.06, vjust = .4, reverse_x=T)   # show pies < .70


pdf(file="out/Fig2b_ecomorph10x2.pdf", height=10, width=2)
  print(e)
dev.off()

##################################################
####
#### Genus Circle Plot (fig 3a)
####
##################################################

# gcol contains a color palette for genera for tree plotting (from support tree above)

# keep one representative per genus in the tree
tokeep <- nospaces(c("Hylophorbus_proekes", "Mantophryne_menziesi", "Callulops sp.1", "Asterophrys sp.", "Xenorhina sp.1", "Genyophryne thomsoni", "Liophryne aff. dentata.1", "Copiula sp.1", "Austrochaperina A parkeri", "Barygenys sp.", "Austrochaperina B macrorhyncha", "Austrochaperina C aff. palmipes.1", "Choerophryne sp.1", "Oreophryne B ezra", "Cophixalus sp.1", "Paedophryne titan", "Oreophryne A sp.1"))
keep <- sapply( tokeep, grep, tree$tip.label, value=T)
names(keep) <- NULL
ctree <- ladderize(keep.tip(tree, keep), right=T)  # genus-level tree

# reduce tip labels to genus only
# (ID) (genus or genus A/B/C) (species) keep only \\2
gtips <- sub("(^[A-Za-z0-9]*-)([A-Za-z]*|[A-Za-z]*_[ABC])(_.*$)", "\\2", ctree$tip.label)
gtips <- gsub("_", " ", gtips)
 
pdf(file="out/Fig3a_fan.pdf")
  plot.phylo(ctree, "fan", edge.width=1, label.offset = 1, cex=1.2, show.tip.label=F)
  tiplabels(pch = 15, cex=2, col = gcol[no_(gtips)])
dev.off()

##################################################
####
#### LTT analysis and plot (fig 3b)
####
##################################################

lttdat<-ltt(tree,plot=FALSE)
pdf(file = "out/Fig3b_lttplot.pdf", height=10, width = 10)
  par(mgp = c(3,1.2,0), mar=c(5,6,4,1)+.1)
  plot(lttdat,log.lineages=T,col="black",lwd=2, show.tree=F, xlab = "Time in MYA", ylab = "ln(Number of Lineages)", yaxt='n', xaxt='n', transparency=.6, cex.lab = 2.5)
  axis(2, cex.axis = 2, lwd.ticks = 1.5, cex.lab = 2)
  axis(1, at=c(0,5,10,15,20), labels=c(20, 15, 10, 5, 0), cex.axis = 2, lwd.ticks = 1.5, cex.lab = 2, )
  text(x = .7, y = 5, pos=4, paste("Gamma =", round(lttdat$gamma, 2)), cex = 2)
  text(x = .7, y = 4.5, pos=4, paste("P-value =", signif(lttdat$p, 3)), cex = 2)
dev.off()
     
##################################################
####
#### Community line plots (fig 5a,b,c)
####
##################################################

#### custom functions ################################
paint.clade <- function( tree, x, name="clade1", startreg = NULL, ... ) {

	## get list of nodes matching each x
	if (is.null(x)) {
		return( reg <- rep("anc", times=tree@nnodes) )
    } else {
	if (class(x)=="character") node = lapply( x, function(x) grep( x, tree@nodelabels, ... )) else
		if (class(x)=="integer" | class(x)=="numeric") node = x else
		return("node must be character or numeric")

	lin <- tree@lineages[unlist(node)]   ## list of lineages
	shared <- sharedlin(lin)
	st <- unique(unlist(sapply( lin, function(x) x %w/o% shared )))  # shared history

	if (is.null(startreg)) {
		reg <- rep("anc", times=tree@nnodes)
		names(reg) <- tree@nodes
		reg[st] = name
	} else {
		reg <- startreg
		reg[st] = name     # paint regimes
	}

    return( reg )
  }
}

sharedlin <- function( lin ){
	l <- length(lin)
	if (l==1) return( lin[[1]][-1] )
	  else 	return( intersect( lin[[ l ]], sharedlin( lin[ -l ] )  ) )
}

mrca <- function( tt, x, ... ) {

	if (class(x)=="character") node = sapply( x, function(x) grep( x, tt@nodelabels, ... )) else
		if (class(x)=="integer" | class(x)=="numeric") node = x else
		return("x must be character or numeric")

	lin <- tt@lineages[node]

    return( sharedlin( lin )[1] )
}

# Function to drop the id number 
dropid <- function(x) sub("^[A-Za-z0-9]*-", "", x)		
# grab just genus species names  (labels are ID-Genus species-site)
getbinom <- function(x) no_(sub("-.*$", "", dropid(x)))

### Community overplotting function
community_overplot <- function(otree, gensp, community="Cliffside", name="Cliffside Camp Kamiali"){
	
  ss <- grep(community, otree@nodelabels, ignore.case=T, )
  tnc <- gensp
  tnc[as.numeric(otree@nodes %w/o% ss)] <- ""

  par(mar=c(2,0,0,0))
  plot(otree, text_opts = list(cex=.5, offset=0.1), margin=c(.3), regimes=paint.clade(otree, community, name=name, ignore.case=T), labels=tnc, palette=mypal2, xaxt="n", legend=F)
  axis(1, at=ticks-adj, labels=labs, pos=-0.01, padj=-1.8, cex.axis=.75, lwd.ticks=.5, tcl=-.3)
  mtext("Time MY Ago", cex=.75, side=1, line=.6, at=10)
  title(name, cex.main=.75, line=-1)
}
#### custom functions end ################################

require(ape)
require(ouch)
require(dplyr)

# Convert to phylo then to ouch format 
tibb <- as_tibble(td)  ## paste site onto tip labels, save as phylo
tibb$label <- paste(tibb$label, tibb$site, sep="-")
tibb$label[219:length(tibb$label)] <- NA
tree <- as.phylo(tibb)
otree <- ape2ouch(tree, scale=F)    # 218 species - ape version is tree
# get tip labels, then get genus species from the labels
gensp <- getbinom(no_(otree@nodelabels))


# settings for the plots
mypal2 <- function(x) {return(c("grey70", "blue"))}
theight <- otree@depth   ## 19.71611
adj <- 20-theight
ticks <- c(0,5,10,15,20)
labs <- c(20,15,10,5,0)

### Make community lineplots 
pdf(paste("out/Fig5a_Cliffside_community.pdf"), width=4, height=7)
  community_overplot(otree, gensp, "Cliffside", "Cliffside Camp Kamiali")
dev.off()

pdf(paste("out/Fig5b_Mwatebu_community.pdf"), width=4, height=7)
  community_overplot(otree, gensp, "Mwatebu", "Mwatebu, Normanby Island")
dev.off()

pdf(paste("out/Fig5c_Gerebu_community.pdf"), width=4, height=7)
  community_overplot(otree, gensp, "Maru Ruama", "Maru Ruama, Mt. Gerebu")
dev.off()


##################################################
####
#### Supplementary Figures - Anc State Rec
####
##################################################
tree <- bi218phylo
tree <- ladderize(tree, right = F)
symcol <- setNames(c("#00ff00", "#ffa500", "#ff1493", "#00bfff", "#0000ff"),
                   c("arboreal", "fossorial", "scansorial", "semi-Aquatic", "terrestrial"))

geigerSYM <- fitDiscrete(tree, lifestyle, model = "SYM", ncores = 4, control = list(method = c("subplex", "L-BFGS-B"), niter = 100))
# geigerSYM$opt$aicc   # AICc=310
# SYM.simmap.summary$ace # has prob anc states
ecorate <- geiger:::.Qmatrix.from.gfit(geigerSYM) # extracts rate matrix from model
ecorate <- round(ecorate, 5)

## 2000 stochastic maps generated from geiger rate matrix
strees <- make.simmap(tree, lifestyle, model="SYM", nsim=2000, Q = ecorate)
SYM.simmap.summary <- summary(strees,plot=FALSE)

pdf("Figure5.pdf", height = 11, width = 8.5)
  plot(SYM.simmap.summary, colors = symcol, fsize = 0.3, ftype = "i")
  add.simmap.legend(colors = symcol, x = 1, y = 220, prompt = F)
dev.off()

##################################################
####
#### Supplementary Figures - tree plots
####
##################################################
require(ape)
require(ggtree)
require(ggplot2)
require(treeio)

# Read in the trees
tree <- read.beast("./BEAST2_analysis/beast_tree_asterophryinae.nex")
iqtree <- read.iqtree("./IQTREE2_analysis/iqtree_timetree_asterophryinae.contree") 
nuc <- read.beast("./BEAST2_analysis/beast_tree_asterophryinae_nuclear.nex")
mito <- read.beast("./BEAST2_analysis/beast_tree_asterophryinae_mitochondrial.nex")

# Function to plot the trees
plot_SI_tree <- function(tree, ttitle="Bayesian Inference Phylogeny (BEAST2)", tipmargin=5) {
	
  # get tree height
  height <- max(node.depth.edgelength(as.phylo(tree))) 
  # get support values
  support <- if ("posterior" %in% get.fields(tree)) round(as_tibble(tree)$posterior, 2) else as_tibble(tree)$UFboot  

  # generate tree with ggtree
  s <- ggtree(tree, size=.3) + 
	 geom_tiplab(size=1.2, offset=.05, fontface='italic' ) + 
     geom_text2(aes(label=support), size = 1.75, vjust=-.3, hjust = 1.1, nudge_x=-.05) +
     scale_x_continuous(limits = c(0, height+tipmargin), breaks=(height - (round(height/5):0)*5),  labels = (round(height/5):0)*5 ) +
     coord_cartesian(clip="off") + 
	 theme_tree2(legend.position="none") +
     ggtitle(ttitle) +
     xlab("Time (MYA)") +
     theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(vjust = 0.5, size = 12, colour = 'black'), 
        axis.line.x = element_line(colour = 'black', size = .5, linetype = 1),
        axis.ticks.x = element_line(size = .5, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.ticks.length=unit(.1, "cm"))
  print(s)    
}

# Plot trees to pdf
pdf(file="out/SupplementaryFigures.pdf", width=8.5, height=11)
  plot_SI_tree(tree, "Bayesian Inference Phylogeny (BEAST2)")
  plot_SI_tree(iqtree, "Maximum Likelihood Phylogeny (IQTREE)")
  plot_SI_tree(nuc, "Bayesian Inference Phylogeny (BEAST2)\nNuclear Loci Only")
  plot_SI_tree(mito, "Bayesian Inference Phylogeny (BEAST2)\nMitochondrial Loci Only", tipmargin=3)
dev.off()

# Plotting figures for DiB
pdf(file="Figure1.pdf", width=8.5, height=11)
  plot_SI_tree(tree, "Bayesian Inference Phylogeny (BEAST2)")
dev.off()

pdf(file="Figure2.pdf", width=8.5, height=11)
plot_SI_tree(iqtree, "Maximum Likelihood Phylogeny (IQTREE)")
dev.off()

pdf(file="Figure3.pdf", width=8.5, height=11)
plot_SI_tree(nuc, "Bayesian Inference Phylogeny (BEAST2)\nNuclear Loci Only")
dev.off()

pdf(file="Figure4.pdf", width=8.5, height=11)
plot_SI_tree(mito, "Bayesian Inference Phylogeny (BEAST2)\nMitochondrial Loci Only", tipmargin=3)
dev.off()


