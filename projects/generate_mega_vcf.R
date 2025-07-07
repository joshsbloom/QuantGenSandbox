library(vcfR)
library(qs)
library(qs2)

xQTLSims.dir = '/home/jbloom/Dropbox/code/QuantGenSandbox/'

source.dir=paste0(xQTLSims.dir, 'R/')
#data.dir=paste0(xQTLSims.dir, 'data/')
#project.dir=paste0(xQTLSims.dir, 'projects/032924/')

#function to simulate crosses, treatement of X is incomplete/broken
#and additional helper functions
#source(paste0(source.dir, 'simWormCrosses.R'))
source(paste0(source.dir, 'helperFxs.R'))
#source(paste0(source.dir, 'makeCountTables.R'))


#unique chromosomes 
#yeast
uchr=paste0('chr', as.character(as.roman(1:16)))
#uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

#elegans
#gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')
#gmap=restructureGeneticMap(gmap.file)
gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
gmap=gmaps[['A']]


#------------------------------------------------------------------------------------------------------------------------------------------------------
#read new vcfs 
vcf.dir='/data1/yeast/reference/pop_vcfs/'

    #was used to retain the set of 1011 samples to keep, saved in pop_vcfs
    #samples.to.keep=c(1,which(nchar(colnames(vcf@gt))==3))
    #filter variants in the 1011 strain collection and only retain SNPs
    setwd(vcf.dir)
    #bcftools view -S samples.txt -v snps -m2 -M2 chrI.norm.vcf.gz |  bcftools filter -e 'AC==0 || AC==AN'  | bcftools view -O z -o test.vcf.gz
    for(u in uchr) {
        print(u)
        system(paste0('bcftools view -S samples.txt -v snps -m2 -M2 ', u, ".norm.vcf.gz |  bcftools filter -e 'AC==0 || AC==AN'  | bcftools view -O z -o", ' 1011/', u, '.vcf.gz'))
    }


    #now read in each vcf and save as R objects for faster access
    #also add BY (S288C) given it is missing in the collection
    #vcfS=list()
    for( u in uchr) {
    #test.vcf.gz
    #}
        print(u)
        yeast.ref.vcf=paste0(vcf.dir, '1011/', u, '.vcf.gz') #'/data1/yeast/reference/pop_vcfs/test.vcf.gz') #', u, '.norm.vcf.gz')
        vcf=vcfR::read.vcfR(yeast.ref.vcf) #, cols=samples.to.keep)
        gt=vcfR::extract.gt(vcf, as.numeric=T)

        print('# of variants')
        print(nrow(gt))

        gcnt=rowSums(!is.na(gt))
        acnt=rowSums(gt, na.rm=T)
        af=acnt/gcnt
        monomorphic= (af==0 | af==1) # ((acnt/gcnt)==0)
        
        vcf=vcf[!monomorphic,]
        gt=gt[!monomorphic,]
        print(nrow(gt))

        gcnt=gcnt[-which(monomorphic)] #rowSums(!is.na(gt))
        acnt=acnt[-which(monomorphic)] #rowSums(gt, na.rm=T)
        af=acnt/gcnt

        #simplest to just prune multi-allelic sites
        pos=getPOS(vcf)
        multiallelic.pos=unique(pos[which(duplicated(pos))])
        biallelic=!(pos %in% multiallelic.pos)
        
        vcf=vcf[biallelic,]
        gt=gt[biallelic,]

        print(nrow(gt))

        ref.hack="0/0:0,0:0:0:.:.:0,0,0:." 
        S288C=rep(ref.hack, nrow(gt))
        vcf@gt=cbind(vcf@gt, S288C)
        S288C=rep(0, nrow(gt))
        gt=cbind(gt, S288C)

        qsave(list(gt=gt, vcf=vcf), file=paste0(vcf.dir, '1011/', u, '.gs'))
        #vcfS[[u]]=vcf
    }


    #manually filter and retain the freaking snps 
    v1=qread('/data1/yeast/reference/pop_vcfs/1011/chrI.gs')
    vcf=v1$vcf
    #manually norm the supid vcf so GATK doesn't choke 
    vcf@fix[,'REF']= substr(vcf@fix[,'REF'],1,1)
    vcf@fix[,'ALT']= substr(vcf@fix[,'ALT'],1,1)
    #snpL = nchar(vcf@fix[,'REF'])==1  & nchar(vcf@fix[,'ALT'])==1 
    gt=v1$gt
    #vcf=vcf #[snpL,]
    #gt=gt #[snpL,]
    for(u in uchr[-1]) { 
        print(u)
        vtemp=qread(paste0('/data1/yeast/reference/pop_vcfs/1011/', u, '.gs'))
        vtemp$vcf@fix[,'REF']= substr(vtemp$vcf@fix[,'REF'],1,1)
        vtemp$vcf@fix[,'ALT']= substr(vtemp$vcf@fix[,'ALT'],1,1)
        #snpL = nchar(vtemp$vcf@fix[,'REF'])==1  & nchar(vtemp$vcf@fix[,'ALT'])==1 
        vcf@fix=rbind(vcf@fix, vtemp$vcf@fix) #[snpL,])
        vcf@gt=rbind(vcf@gt, vtemp$vcf@gt) #[snpL,])
        gt=rbind(gt, vtemp$gt) #[snpL,])
    }

    # save vcf and gt for future processing and simulations 
    qsave(list(gt=gt, vcf=vcf), file=paste0(vcf.dir, '1011/mega_filtered.qs'))


    #hack a vcf with het calls at each variant site for use with gatk ASEReadCounter
    vcfout=vcf
    vcfout@gt=vcfout@gt[,c(1,2)]
    vcfout@gt[,1]='GT:DP'
    vcfout@gt[,2]='0/1:.'
    write.vcf(vcfout,'/data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf')
    #then on the command line, bgzip and gatk IndexFeatureFile


    # figure out reticulate behavior
    #library(reticulate)
    #use_condaenv('genomics')
    #maybe switch to absolute python path 
    #system(paste0('
    bcftools view -I /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf -O z -o /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz 
    #'), intern=T)

    #bcftools norm -a -m -snps -f /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta mega_filtered.vcf.gz > mega_filtered_snps.vcf
    #use_condaenv('genomics')
    #system(paste0("
    gatk IndexFeatureFile -I  /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz
    #"), intern=T)

    #------------------------------------------------------------------------------
    #bcftools index mega_filtered.vcf.gz 
    #tabix mega_filtered.vcf.gz
    gatk ASEReadCounter --min-mapping-quality 20 -I /media/hoffman2/lcrisman/Dmagicmarker_042024/bam_files/G1_KANa_S62.bam -V /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz -R /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta -O /data0/xQTLSims/projects/061024/data/G1_KANa_S62.txt
    gatk ASEReadCounter --min-mapping-quality 20 -I /media/hoffman2/lcrisman/Dmagicmarker_042024/bam_files/G2_NATalpha_S63.bam -V /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz -R /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta -O /data0/xQTLSims/projects/061024/data/G2_NATalpha_S63.txt



    #manually removed structural variants, and sorted the thing with bcftools sort 
    #yeast.ref.vcf='/media/hoffman2/jsbloom/reference/rr_parents_no_svar.vcf.gz'
    #vcf=vcfR::read.vcfR(yeast.ref.vcf)
    #vcf=vcf[vcfR::is.biallelic(vcf),]
    ##vcf=vcf[-which(duplicated(paste0(getCHROM(vcf),':', getPOS(vcf)))),]
    #vcf=vcf[!getCHROM(vcf)=='chrM',]

    #gt=vcfR::extract.gt(vcf, as.numeric=T)

#-----------------------------------------------------------------------------------------------------------------------
library(qs)
library(qs2)
library(vcfR)
vcf.dir='/data1/yeast/reference/pop_vcfs/'
mega_filtered=qread(paste0(vcf.dir, '1011/mega_filtered.qs'))
mega_filtered_processed=(paste0(vcf.dir, '1011/mega_filtered_processed.qs2'))
vcf=mega_filtered$vcf
vcf.fix=vcf@fix
gt=mega_filtered$gt
attr(gt, 'fix')=vcf@fix
fix.df=data.frame(attr(gt, 'fix'), stringsAsFactors=F)
fix.df[,'POS']=as.numeric(fix.df[,'POS'])
fix.df[,'QUAL']=as.numeric(fix.df[,'QUAL'])

#calc minor allele frequencies --------------------------------
altcnt=rowSums(gt, na.rm=T)
nna=rowSums(!is.na(gt))
af=altcnt/nna
af[af>.5]=1-af[af>.5]
maf=af
hist(maf)

fix.df$MAF=maf
rownames(fix.df)= rownames(gt)

rownames(attr(gt, 'fix'))=rownames(gt)
attr(gt, 'fix')=fix.df
qs_save(gt, mega_filtered_processed, compress_level=16)

