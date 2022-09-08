# Functions for cleaning metadata and trees

######################## DATA CLEANING FUNCTIONS #####################
"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
nospaces <- function(x) gsub(" ", "_", x)  # replace spaces with "_"
no_ <- function(x) gsub("_", " ", x)  # replace spaces with "_"
delspaces <- function(x) gsub(" ", "", x)  # deletes spaces
### Better Tiplabel functions
getid <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\1", no_(x) )  # split tiplabel by - grab id
getsp <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\2", no_(x) )
getsite <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\3", no_(x) )
getlocality <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\4", no_(x) )
getidsp <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\1-\\2", no_(x) )
getspsite <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\2-\\3", no_(x) ) 
getspsitelocality <- function(x) gsub("(.*)\\-(.*)\\-(.*)\\-(.*)","\\2-\\3-\\4", no_(x) ) 


clean_elevation <- function( dat ) {
  loc <- dat$Locality.Selector

  pattern <- "([[:digit:]]+[ ]*)-([ ]*[[:digit:]]+) m|([[:digit:]]+) m|([[:digit:]]+[ ]*)-([ ]*[[:digit:]]+) M|([[:digit:]]+) M"   # search for ʻ0-100 mʻ or ʻ100 mʻ
  lout <- regmatches(loc, gregexpr(pattern, loc))
  lout[lout=="character(0)"] <- ""
  elevation <- sapply(lout, function(x) return(x[length(x)]) )   # elevation data
  dat$elevationrange <- gsub(" [Mm]$","", elevation)
  dat$elevation <- sapply(strsplit(dat$elevationrange, "-"), function(x) mean(as.numeric(unlist(x))))
  return(dat)
}

clean_species <- function( dat ){
  species <- trimws(dat$species)
  species <- sub("(Genus of Microhylidae)", "sp.", species, fixed=T)  # replace with sp.
  species <- gsub("[ \t\r\n]+Kraus 2013", "", species)
  species <- gsub("Albericus", "Choerophryne", species)
  species <- gsub("Metamagnusia", "Asterophrys", species)
  species <- gsub("Pseudocallulops", "Asterophrys", species)
  species <- gsub("Oninia", "Asterophrys", species)
  species <- gsub("sp[0-9]+$", "sp.", species)
  
  dat$gensp <- trimws(species)
  dat$genus <- trimws(gsub(" .*$", "", dat$gensp))
  dat$species <- trimws(sub(".*? ", "", dat$gensp))
    
  return(dat)	
}

clean_gps <- function(dat){

  localnotes <- dat$ORIGLOCAL
  dat$sitelocality <- paste(trimws(dat$site), trimws(dat$locality), sep="-")  # paste site, locality

# Cleans the gps data
  gps <- dat$gps
  gps <- gsub("[°º]", "", gps)            # remove degree symbols
  gps <- gsub("^(.*)S(.* *)E", "-\\1\\2", gps)  # if S it should be - , if E should be +
  gps <- gsub("^(.*)N(.* *)E", "\\1\\2", gps)   # if N it should be + , if E should be +
  gps <- gsub("^(.*)S(.* *)W", "-\\1-\\2", gps)  # if S it should be - , if W should be -
  gps <- gsub("^(.*)N(.* *)W", "\\1-\\2", gps)  # if N it should be + , if W should be -
  gps <- gsub("^([-0-9.0-9]*)[ ,]*([-0-9.0-9]*)[ ]*", "\\1, \\2", gps)  # tidy up white spaces ,
  gps_clean <- gps

  dat$lat <- as.numeric(sub(",.*", "", gps))
  dat$lon <- as.numeric(sub(".*,", "", gps))
  dat$gps_raw <- dat$gps
  dat$gps <- gps
  return(dat)
}

pasteprefix <- function (id, prefix) {
	id <- as.character(id)
	id[!is.na(id)] <- paste(prefix, id[!is.na(id)], sep="")
	return(id)
}

# converts beast to phylo, adding posterior probability to node labels
beast2phylo <- function(beast){
	tibb <- as_tibble(beast)
	istip <- !is.na(tibb$label)    # is it a tip? 
    phylo <- ladderize(as.phylo(tibb), right=F)
    phylo$node.label <- round(as.numeric(tibb$posterior[!istip]),2)  # support values at nodes	
	return(phylo)
}

## drops any columns of a dataframe that are lists
flatten <- function(dat){
  oo <- sapply(dat, is.list)           					# need to flatten listd by dropping list elements
  dat <- dat[!oo]
  return(dat)	
}

## drops internal node rows from tree dataframes, id should have NAs in internal nodes 
tipsonly <- function( dat, id ) {
  istip <- !is.na(id)
  dat <- dat[istip,]
  id <- id[istip]
  rownames(dat) <- id
  return(dat)
}

allduplicates <- function(x) {
	y <- duplicated(x)
	z <- duplicated(x, fromLast = TRUE)
    return(y + z >0)	
}

make.transparent <- function(color="gold", transparency=25) {
  co <- col2rgb(color)
  trco <- rgb(co["red",], co["green",], co["blue",], max=255, alpha=transparency)
  return(trco)
}
###################################################################
