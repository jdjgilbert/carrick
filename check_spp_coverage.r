
## Checking Cat's list of species against mine

## First load up my list of species with associated taxonomy

require(ape)
require(digest)

## Get phylogeny
Ntip(phy <- read.nexus('Data/wphy.wheeler.20130726.nex'))
# 1935
digest(phy$tip.label, algo='md5', serialize=F)  ## md5 summary of phylogeny tips
#[1] "2970395c9154c946689fcb7deb964267"

digest(phy$node.label, algo='md5', serialize=F)  ## md5 summary of phylogeny nodes  -- interestingly this often changes - why?
# [1] "947de11a14b6d4476362008e527548d6"

## Get taxonomy
nrow(tax <- read.csv('data/hierarchical_taxonomy_2011_01.csv'))
# 1740
head(tax)
#  id      sbp spc     cla        sbc       ifc      dvn           sbd        ord         sbo              ifo            spf
#  1 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera Endopterygota Coleoptera    Adephaga             <NA>           <NA>
#    2 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera  Paraneoptera  Hemiptera Heteroptera       Nepomorpha       Nepoidea
#  3 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera Endopterygota Coleoptera   Polyphaga Staphyliniformia Staphylinoidea
#  4 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera          <NA> Dermaptera Forficulina             <NA>  Forficuloidea
#  5 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera Endopterygota Coleoptera   Polyphaga     Cucujiformia Chrysomeloidea
#  6 NA Hexapoda  NA Insecta Dicondylia Pterygota Neoptera  Paraneoptera  Hemiptera Heteroptera  Pentatomomorpha  Pentatomoidea
#  fam              sbf  tri             gen        sp
#  1        Carabidae             <NA> <NA>            Abax       spp
#  2   Belostomatidae             <NA> <NA>          Abedus  herberti
#  3        Silphidae        Silphinae <NA>      Ablattaria laevigata
#  4     Forficulidae             <NA> <NA>   Acanthocordax  papuanus
#  5        Bruchidae             <NA> <NA> Acanthoscelides  obtectus
#  6 Acanthosomatidae Acanthosomatinae <NA>     Acanthosoma   griseum

digest(as.character(unlist(tax)), algo='md5', serialize=F)  ## md5 summary of phylogeny tips
#[1] "d4cd0dabcf4caa22ad92fab40844c786"

## separate genera from species in spp list
length(phy_genera <- do.call(rbind, strsplit(phy$tip.label, '_'))[,1])
head(phy_genera)
#[1] "Neomachilellus" "Neomachilellus" "Neomachilellus" "Meinertellus"   "Lepisma"        "Thermobia"     

## How many genera are there on the tree
length(unique(phy_genera))
## 1194
## How many genera on the tree are in the taxonomy and can therefore be assigned to a family
length(intersect(phy_genera, as.character(tax$gen)))
# 802

### Downloading ITIS genus-level database so that I can cross reference using the taxonomy

head(itisgen <- read.table('Data/2016-01-19_ITIS_arthropoda_genus_level.csv',  sep='|', nrows=37007, na.strings=""))
 
names(itisgen) <- c('field1','tsn','unit_ind1','unit_name1','unit_ind2','unit_name2','unit_ind3','unit_name3','unit_ind4','unit_name4','unnamed_taxon_ind','usage','unaccept_reason','credibility_rtng','completeness_rtng','currency_rating','phylo_sort_sequence','initial_time_stamp','parent_tsn','taxon_author_id','hybrid_author_id','kingdom_id','rank_id','update_date','uncertain_prnt_ind','name_usage','complete_name')

head(itisgen)


## are all "parent" TSNs in the list of TSNs?
with(itisgen, all(parent_tsn %in% tsn))
## F

with(itisgen, length(which(!parent_tsn %in% tsn)))
## 4397 parent taxa do not have entries in their own right [out of 37007]

## Try to match parents to offspring taxa
head(itisgen$parent_tax <- itisgen$unit_name1[match(itisgen$parent_tsn, itisgen$tsn)])
# [1] Hypogastruridae <NA>            Hypogastruridae Hypogastruridae <NA>            Hypogastruridae
# 23211 Levels: Aaata Aaroniella Aaroniellinae Aata Abablemma Abacaria Abacetini Abachrysa ... Zyzzyx

## Check Schaefferia is in the Hypogastruridae
## YES

length(intersect(phy_genera, itisgen$unit_name1))
## 544 [out of 1256]

## How many genera are in either of the taxonomy or the ITIS taxonomy
length(intersect(phy_genera, c(as.character(itisgen$unit_name1), as.character(tax$gen))))
# 942


### Import Cat's species list
nrow(cat <- read.csv('Data/catherine_carrick_2016-01-19.csv'))
# 344

head(cat$genus <- sub(' .*','',cat$Species2))
#[1] "T."              "Eotetranychus"   "E."              "Eutetranychus"   "Metatetranychis" "Panonychus"     

head(catsp <- strsplit(as.character(cat$Species2), ' '))
head(cat$sp <- sapply(catsp, function(x) {x[2]} ))
cat$binomial <- paste(cat$genus, cat$sp, sep='_')


## How many genera are they from?
length(unique(cat$genus))
# [1] 263

## can we assign Cat's genera to families using ITIS?
length(intersect(cat$genus, itisgen$unit_name1))
# 152
## Need Cat to do the rest

head(cat)


  
## Write out handy reference dataset for placing children taxa into parent taxa
itistax <- subset(itisgen, select=c('unit_name1','parent_tax'))
names(itistax)<-c('taxon','parent_taxon')
 write.csv(itistax, file="Data/ITIS_parent_list.csv")

length(unplaced_genera <- unique(setdiff(cat$genus, unlist(c(phy_genera, itistax$taxon)))))
# 219

length(intersect(cat$Family, itisgen$parent_tax))


library(RNeXML)

### Use RNeXML to extract family-level tree structure from tree of life 
### 
       





### ALTERNATIVELY: Use simple taxonomy based tree for now.

## Use Wheeler order-level tree as base
## Below order level, construct hierarchical taxonomy table from ITIS
## Below family level, create star phylogenies
  
head(a <- subset(itisgen, select=c('unit_name1', 'parent_tax')))
names(a) <- c('t1','t2')
head(ca <- subset(a, t1 %in% cat$Family))

ca$t3 <- a$t2[match(ca$t2, a$t1)]
ca$t4 <- a$t2[match(ca$t3, a$t1)]
ca$t5 <- a$t2[match(ca$t4, a$t1)]
ca$t6 <- a$t2[match(ca$t5, a$t1)]
ca$t7 <- a$t2[match(ca$t6, a$t1)]
ca$t8 <- a$t2[match(ca$t7, a$t1)]
ca$t9 <- a$t2[match(ca$t8, a$t1)]
ca$t10 <- a$t2[match(ca$t9, a$t1)]
ca$t11 <- a$t2[match(ca$t10, a$t1)]
ca$t12 <- a$t2[match(ca$t11, a$t1)]
ca$t13 <- a$t2[match(ca$t12, a$t1)]
ca$t14 <- a$t2[match(ca$t13, a$t1)]
ca$t15 <- a$t2[match(ca$t14, a$t1)]
ca$t16 <- a$t2[match(ca$t15, a$t1)]
ca$t17 <- a$t2[match(ca$t16, a$t1)]
all(is.na(ca$t17))
#[1] TRUE

write.csv(ca, file='cat_family_taxonomy_1.csv')
## Now use Excel to align the columns!

cattax <- subset(read.csv(file='cat_family_taxonomy.csv', colClasses='character'), select=c('t1','t2','t3','t4','t5','t6','t7','t8'))

cattax$t2[which(cattax$t2=='')]<-cattax$t1[which(cattax$t2=='')]
cattax$t3[which(cattax$t3=='')]<-cattax$t2[which(cattax$t3=='')]
cattax$t4[which(cattax$t4=='')]<-cattax$t3[which(cattax$t4=='')]
cattax$t5[which(cattax$t5=='')]<-cattax$t4[which(cattax$t5=='')]
cattax$t6[which(cattax$t6=='')]<-cattax$t5[which(cattax$t6=='')]
cattax$t7[which(cattax$t7=='')]<-cattax$t6[which(cattax$t7=='')]
cattax$t8[which(cattax$t8=='')]<-cattax$t7[which(cattax$t8=='')]
 cattax$t9<-rep('Animalia', nrow(cattax))

for (i in 1:9) { cattax[,i]<-as.factor(cattax[,i]) }

#plot(as.phylo.formula(~t8/t7/t6/t5/t4/t3/t2/t1, cattax))

wheelertree <- read.caic('Hexapoda.card')
wheelertree$tip.label <- sub('Neuropteroidea','Neuroptera', wheelertree$tip.label)
wheelertree$tip.label <- sub('Phasmatodea','Phasmida', wheelertree$tip.label)
wheelertree$tip.label <- sub('Myriapoda','Arachnida', wheelertree$tip.label)

for (i in unique(cattax$t5))  {
    xx <- subset(cattax, t5 == i) 
      f4 <- ~t4/t3/t2/t1
      f3 <- ~t3/t2/t1
      f2 <- ~t2/t1
      f1 <- ~t1
    form <- f1
    if(length(unique(xx$t2))>1) form <- f2
    if(length(unique(xx$t3))>1) form <- f3
    if(length(unique(xx$t4))>1) form <- f4
    assign(as.character(i), as.phylo.formula(form, xx))
}

phy <- wheelertree
 for (i in intersect(unique(cattax$t5), phy$tip.label)) {
    phy <- bind.tree(phy, get(i), where=grep(i, phy$tip.label))
}

Ntip(phy)
# 111

#code from treelego.r
makestarclade <- function(splist) {
  collapse.singles(read.tree(text=paste('(',paste(splist, collapse=','), ');', sep='')))
}

for (i in unique(cat$Family)) {
  phy <- bind.tree(phy, makestarclade(cat$binomial[which(cat$Family==i)]), where=grep(paste('^',i,'$',sep=''), phy$tip.label))  
}

Ntip(phy)
#[1] 365

all(cat$binomial %in% phy$tip.label)
## [1] TRUE

phy <- compute.brlen(phy, power=0.8)

## some labels are duplicated
phy$tip.label[which(duplicated(phy$tip.label))] <- paste(phy$tip.label[which(duplicated(phy$tip.label))], '_1', sep='')

write.nexus(phy, file='phylogeny_from_wheeler_plus_taxonomy_20160713.nex')
