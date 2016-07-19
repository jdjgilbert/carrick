
setwd('~/Dropbox/_work/R analyses/Carrick/')

tax <- read.table('~/Downloads/ott/taxonomy1.tsv', sep='\t', comment.char='')
tree <- read.tree('~/Downloads/opentree5.0/output/grafted_solution/grafted_solution.tre')
tree1 <- tree
tree1$tip.label <- as.character(tax$V5[match(sub('ott','',tree$tip.label), tax$V1)])
tree1$node.label <- as.character(tax$V5[match(sub('ott','',tree$node.label), tax$V1)])
grep('Hexapoda', tree1$node.label) # 27829

hextree <- extract.clade(tree1, grep('Hexapoda', tree1$node.label)+Ntip(tree1))

Ntip(hextree) #3098

head(hextree$tip.label)
[1] "Diplura"       "Monocondylia"  "Odonata"       "Ephemeroptera" "Phasmatodea"   "Blaberoidea"  


head(hextree$node.label)
[1] "Hexapoda"   NA           "Insecta"    "Dicondylia" NA           "Neoptera"  



write.nexus(hextree, file='hexapoda_opentol.nex')


cat_spp <- scan(what='character') ## read in cat spp list

cat_spp[grep('^[A-Z]\\.', cat_spp)] <- c("Eotetranychus_willametei","Clavigralla_tomentosicollis","Nephotettix_nigropictus")


intersect(cat_spp, sub(' ','_',hextree$tip.label))  ## 0
length(intersect(sub('_.*', '', cat_spp), sub(' .*','',hextree$tip.label)))  ## 17

cat_fam <- scan(what='character', sep='\n')

length(unique(cat_fam)) # 100
length(intersect(cat_fam, hextree$node.label)) #19


### how many are left after we account for the families and genera in the phylogeny


check <- rep(0, length(cat_spp))
check[cat_spp %in% c(hextree$tip.label, hextree$node.label)] <- 1
check[sub('_.*', '', cat_spp) %in% sub(' .*','', c(hextree$tip.label, hextree$node.label))] <- 1
check[cat_fam %in% c(hextree$node.label,hextree$tip.label)] <- 1
table(check)


### Use taxonomy to place other spp/families?

length(intersect(cat_fam, tax$V5))  #79 out of 100
length(intersect(sub(' .*','',cat_spp), tax$V5))  #0

length(intersect(intersect(cat_fam, hextree$node.label), tax$V5))  ### all 19

which(tax$V5 %in% intersect(cat_fam, hextree$node.label))
o

known_parent <- as.numeric(as.character(unique(tax$V3[tax$V5 %in% intersect(cat_fam, hextree$node.label)])))


known_parent_id <- as.character(unique(tax$V5[tax$V3 %in% known_parent]))

length(intersect(cat_fam, known_parent_id)) # 28

check[cat_fam %in% intersect(cat_fam, known_parent_id)] <- 1

table(check)
check
  0   1 
242 104 

parent_id <- as.character(unique(tax$V5[tax$V1 %in% known_parent]))

known_parent2 <- as.numeric(as.character(unique(tax$V3[tax$V1 %in% parent_id])))
known_parent2_id <- as.character(unique(tax$V5[tax$V3 %in% known_parent2]))
