library(xQTLStats)
library(qs)
library(qs2)
#yuck
library(tidyverse)
library(AlphaSimR)
library(Rfast)
library(ggpubr)
library(susieR)

#optional 
#library(vcfR)

#pull from github
xQTLSims.dir = '/home/jbloom/Dropbox/code/QuantGenSandbox/'

source.dir=paste0(xQTLSims.dir, 'R/')

#function to simulate crosses and additional helper functions
source(paste0(source.dir, 'helperFxs.R'))
source(paste0(source.dir, 'makeCountTables.R'))


#unique chromosomes yeast
uchr=paste0('chr', as.character(as.roman(1:16)))

#yeast genetic maps
gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
gmap=gmaps[['A']]

vcf.dir='/data1/yeast/reference/pop_vcfs/'
#preprocessed genotype data from Joseph's yeast collection (filtered for only biallelic SNPs)
mega_filtered_processed=(paste0(vcf.dir, '1011/mega_filtered_processed.qs2'))
#careful, this thing blows up to 10GB in R
gt=qs2::qs_read(mega_filtered_processed)
attr(gt, 'fix')$ID=paste0(data.table::tstrsplit(rownames(gt) ,'_')[[1]],'_', data.table::tstrsplit(rownames(gt), '_')[[2]])


# GWAS-related sims -------------------------------------------------------------------------
    #build data structure for GWAS.variants 
    mnames=data.frame('chrom'=attr(gt,'fix')[,'CHROM'], 'mname'=rownames(gt), 'pos'=attr(gt,'fix')[,'POS'], stringsAsFactors=F)
    #this has to be structured with strains as rows and markers as columns 
    GWAS.variants=list(g=t(gt), mnames=mnames)
    source(paste0(source.dir, 'GWAS_preprocessing.R'))
    source(paste0(source.dir, 'GWAS_dofastLMM_mvnpermute.R'))

    #this takes a minute
    GWAS.variants=filterGWASvariants(GWAS.variants, af.cutoff=.05, na.frq=.05)
    #calc relatedness matrix A, for the sake of the tutorial aspect #impute missing genotype values to make the code and math less annoying, and hard call them
    GWAS.variants=calcA(GWAS.variants, return.imp=T)
    GWAS.variants$g=ifelse(GWAS.variants$g<.5,0,1)
    svdA=fastLMMsvd(GWAS.variants$A)

    Gin=t(GWAS.variants$g)
    attr(Gin, 'fix')=attr(gt,'fix')[rownames(attr(gt, 'fix')) %in% rownames(Gin),]

    FB=createFounderPop(Gin,colnames(Gin),gmap, X.only=F,X.drop=F, filterNAs=T)
    SP=SimParam$new(FB)
    genMap=getGenMap(FB)
    genMap$maf=attr(Gin, 'fix')$MAF #[match(genMap$id, rownames(Gin))]

    rc.vecs=simRareCommon(genMap, nQTL.r=20, nQTL.c=5)
    QTL.sims=make.QTL.sims(rc.vecs, o.h2.norm=T, o.h2=.4)
    #f designates QTL effects for hermaphrodites 
   
    simFR=simPheno(FB, genMapMarkers=genMap$id, QTL.sims=QTL.sims, returnG=F)


    test=dofastLMM_mvnpermute(y=simFR$simy,X=NULL,GWAS.variants, svdA=svdA, nperm=500, REML=T) 
    plot(test$REML.nlp, col=test$REMLsig+1)
    abline(v=match(rc.vecs$o.add.qtl.ind, genMap$id))
# -------------------------------------------------------------------------------------------------------



# biparental cross sims ------------------------------------------------------------------------------------
    #consider setting to 10,000 for individual segregant mapping and 100,000 for xQTL 
    max.per.gen=1e4
    #meta.results=list()
    #setup sims using mega pop then mini sims per biparental
    #we'll do BY and RM to start, for example 
    p.names=c('S288C', 'AAA')
    FB=createFounderPop(gt,p.names,gmap)
    SP=SimParam$new(FB)
    genMaps=getGenMap(FB)
    nid=paste0(data.table::tstrsplit(genMaps$id, '_')[[1]],'_', data.table::tstrsplit(genMaps$id, '_')[[2]])
    genMaps$maf=attr(gt, 'fix')$MAF[match(nid, attr(gt, 'fix')$ID)]

    FB=newPop(FB, simParam=SP)
    f1=makeCross(FB, matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
    f2=makeDH(f1, nDH=max.per.gen, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) ,
    FR=f2

    rc.vecs=simRareCommon(genMaps, nQTL.r=20, nQTL.c=5, seed=20)
    QTL.sims=make.QTL.sims(rc.vecs, o.h2.norm=T, o.h2=.4)
    # for QTL sims set return G = T 
    simFR=simPheno(FR, genMapMarkers=genMaps$id, QTL.sims=QTL.sims, returnG=T)
    #-
    #genotypes are simFR$G, phenotypes are simFR$simy
    #do QTL mapping
# ----------------------------------------------------------------------------------------


# for xQTL, can resuse the above bit and set max.per.gen to 100,000  -----------------------
    depth=50
    sel.frac=.1
    max.per.gen=1e5
    #for xQTL sims, set returnG = F to help with memory management
    #over write simFR
    f2=makeDH(f1, nDH=max.per.gen, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) ,
    FR=f2
    simFR=simPheno(FR, genMapMarkers=genMaps$id, QTL.sims=QTL.sims, returnG=F)
   
    countdf.h=simSequencingRefAlt(simFR$simy,FR, genMaps$id, depth=depth, sel.frac=sel.frac, lower.tail=F)
        #if you set sel.frac=1 sample the population of existing genotypes without QTL effects 
        #countdf.l=simSequencingRefAlt(y=NULL, FR,genMaps$id, depth=depth, sel.frac=1 , lower.tail=F)
    ds.ind=sort(sample(max.per.gen, max.per.gen*sel.frac))
    countdf.l=simSequencingRefAlt(y=NULL, FR[ds.ind],genMaps$id, depth=depth, sel.frac=1 , lower.tail=F)

    countdf.h=phaseBiparental(countdf.h, p.names[1], FB, genMaps)
    countdf.l=phaseBiparental(countdf.l, p.names[1], FB, genMaps)
    #  plot(countdf.h$p1/(countdf.h$p1+countdf.h$p2))
    #  points(countdf.h$expected, col='red') #alt/(countdf$al

    #-----------------------------
    test  = calcAFD(countdf.h, experiment.name='high1',sample.size=1e4, bin.width=500, sel.strength=.1, uchr=unique(genMaps$chr) ) 
    test2 = calcAFD(countdf.l, experiment.name='unsel1',sample.size=1e4, bin.width=500, sel.strength=1, uchr=unique(genMaps$chr) )
    results=calcContrastStats(results=list(test, test2), L='_high1', R='_unsel1')
    plotSummary(results) 
#-----------------------------------------------------------------------------------------------
