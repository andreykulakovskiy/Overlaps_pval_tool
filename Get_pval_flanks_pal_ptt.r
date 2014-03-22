require(IRanges)
require(parallel)
require (data.table)
require(plyr)

feature <- "transposase"
work.env <- "bacteria_S6-15_L0-10_M1"
cat(paste(Sys.time(), "This R script has started!", sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)

pal.path <- paste0("D:/DNAPUNCTUATION/ftp/", work.env, "/")
ptt.path <- "D:/DNAPUNCTUATION/all.ptt/"

# pal.path <- paste0("D:/DNAPUNCTUATION/ftp/", work.env, "/")
# ptt.path <- "D:/DNAPUNCTUATION/all.ptt/"

pal.files <- list.files(path=pal.path, pattern="cleaned$", recursive=T, full.names=F)
ptt.files <- list.files(path=ptt.path, pattern="ptt$", recursive=T, full.names=F)

genlen.all <- read.table(file="D:/smpl_files/genome.length.txt", header=F, sep="\t")

#extract NC ID's from .ptt files
m <- regexpr("\x2F(.+?).ptt", x<-ptt.files, perl=T)
ptt.nc.id <-as.vector(regmatches(x, m))
ptt.nc.id <- substring(ptt.nc.id, first=2, last=10)

#extract NC ID's from .pal.cleaned files
m <- regexpr("\x2F(.+?).fna", x<-pal.files, perl=T)
pal.nc.id <-as.vector(regmatches(x, m))
pal.nc.id <- substring(pal.nc.id, first=2, last=10)

#prepare data.frames with .ptt and .pal.cleaned files
ptt.df <- data.frame(ptt.file = paste0(ptt.path,ptt.files), nc.id = ptt.nc.id)
pal.df <- data.frame(pal.file = paste0(pal.path,pal.files), nc.id = pal.nc.id)

#data.frame with both .ptt and pal.cleaned files, will be passed to our parallel env.
all.df <- merge(ptt.df, pal.df,by = intersect("nc.id", "nc.id"))
all.df$nc.id <-as.character(all.df$nc.id)
all.df$ptt.file <-as.character(all.df$ptt.file)
all.df$pal.file <-as.character(all.df$pal.file)

#add new one feature
coverage.data <- read.table('s615.f50.trp.aroundone.txt', header=T, sep='\t')

all.df <- subset(all.df, subset = (all.df$nc.id %in% coverage.data$nc))


GetGenLen <- function (nc.id=NULL) {
    genlen.subset <- subset(genlen.all, subset=(genlen.all$V1==nc.id));
    genlen<-genlen.subset$V2
    return(genlen)
}

GetPTTFile <- function (ptt.file = NULL, product = NULL){
    ptts <- try(fread(input=ptt.file))
    if ( any(class(ptts) == "try-error") ) {
        ptts <- NULL
        cat(paste(Sys.time(), "error in PTT reading!: ",ptt.file, sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)  
    }else{
        ptt.product.pattern <- paste0("*",product)
        ptts<-subset(x=ptts, grepl(ptt.product.pattern, ptts$Product))
        
        list.tmp <- strsplit(as.character(ptts$Location), "\\..")
        StartEnd <- ldply(list.tmp)
        if(nrow(StartEnd)>0) colnames(StartEnd) <- c("Start", "End")
        
        ptts$Start <- as.numeric(StartEnd$Start)
        ptts$End <- as.numeric(StartEnd$End)
        ptts$Size <- as.numeric(ptts$End-ptts$Start)
        
        ptts <- subset(ptts, ptts$Size>0)
        class(ptts)<- "data.frame"
        
    }
    return(ptts)    
}

pval.find <-function (nc.id=NULL, ptt.file=NULL, pal.file=NULL){
    ###logging###
    cat(paste(Sys.time(), "start: ",nc.id, ptt.file, pal.file, sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)
    ###/logging###
    PTTs <- GetPTTFile(ptt.file,feature)
    genlen <- GetGenLen(nc.id)
    PALs <- try(fread(input = pal.file, sep = "\t"))
    
    if ( any(class(PALs)!="try-error" ) & any(class(PTTs) != "NULL") & nrow(PTTs)>0){
        PALs.IR <- reduce(IRanges(start=PALs$V1, end=PALs$V2))
        
        ##MAybe there's an ERROR!
        leftside <- IRanges(start=PTTs$Start-50, end=PTTs$Start)
        rightside <-IRanges(start=PTTs$End, end = PTTs$End+50)
        PTTs.IR <-reduce(c(leftside, rightside)) 
        
        
        hits.real <- sum(countOverlaps(query = PALs.IR, subject = PTTs.IR))
        pval <- 0
        
        maxcoord<-genlen - (max(PTTs.IR@width)+50)
        reps <- length(PTTs.IR)
        if (reps >0 & maxcoord >0){
            for (i in 1:1000){
                start.rnd <- sort(sample(maxcoord,reps))
                end.rnd <-start.rnd+PTTs$Size
                
                #PTTs.IR УЖЕ содержат фланки, следовательно, сейчас мы сгенерим в 2 раза больше отрезков
                leftside <- IRanges(start=start.rnd-50, end=start.rnd)
                rightside <-IRanges(start=end.rnd, end = end.rnd+50)
                PTTs.rand.IR <-reduce(c(leftside, rightside))
                
                hits.rand <-sum(countOverlaps(query=PALs.IR, subject=PTTs.rand.IR))
                if (hits.rand >= hits.real) pval=pval+1
            }
            
        }else{
            cat(paste(Sys.time(), "error in reps or maxcoord: ",nc.id, ptt.file, pal.file, sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)    
            pval = -1
        }   
    }else{
        cat(paste(Sys.time(), "PTT or PAL err: ",nc.id, ptt.file, pal.file, sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)
        pval= -2
    }
    cat(paste(Sys.time(), "result: ",nc.id, pval, sep="\t"), file=paste0(work.env, "_", feature,"_flank50_PALvsPTT.results.txt"), sep="\n", append=T)
    return(pval)
}

cat(paste(Sys.time(), "Making cluster!", sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)
cl.local<-makeCluster(getOption("cl.cores", 3))
clusterExport(cl=cl.local, varlist=c("all.df", "genlen.all", "pval.find", "GetGenLen","GetPTTFile", "work.env", "feature"))
clusterEvalQ(cl.local,require(IRanges))
clusterEvalQ(cl.local,require(data.table))
clusterEvalQ(cl.local,require(plyr))
all.df$rslts<-parApply(cl.local, all.df[,c('nc.id','ptt.file','pal.file')], 1, function(x) pval.find(x[1],x[2],x[3]))
stopCluster(cl.local)
cat(paste(Sys.time(), "Stopping Cluster", sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)

write.table(x<-all.df, file=paste0(work.env, "_", feature,"_flank50_PALvsPTT.results2.txt"), sep="\t")

cat(paste(Sys.time(), "DONE!", sep="\t"), file=paste0(work.env, "_", feature, "_flank50_PALvsPTT.log"), sep="\n", append=T)
