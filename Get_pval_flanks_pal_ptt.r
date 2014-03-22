worktime.start = Sys.time() #Get time of script's start...
require(IRanges)
require(parallel)
require (data.table)
require(plyr)
require(RCurl)

#====== List of parameters: ======
pal.path <- "D:/DNAPUNCTUATION/ftp/bacteria_S6-15_L0-10_M1" #folder with .pal.cleaned files
ptt.path <- "D:/DNAPUNCTUATION/all.ptt/" #folder with .ptt files
feature = 'transposase' #feature in genes markdown(PTT files)
genlenfile = 'D:/smpl_files/genome.length.txt'    
flank_length = 50  #size of flank
tag = 'TEST' #testing
#=================================

#Now let's create names of log.file, results and final summary
palindrome.type = basename(pal.path) # Example: /ftp/bacteria_S6-15_L0-10_M1 ===> bacteria_S6-15_L0-10_M1  
main.dir = getwd()# Script executing dir
results.path =file.path(main.dir, 'pal_ptt_results')
dir.create(results.path, showWarnings=F) #path to results folder for current work
cur.res.name = paste(tag, palindrome.type, feature, 'flank', flank_length,'pal_ptt', sep = '_') #common name for current work

log.file = file.path(results.path, paste0(cur.res.name,'.log'))
out.file = file.path(results.path, paste0(cur.res.name,'.out'))
all.file = file.path(results.path, paste0(cur.res.name,'.summary'))

#get all needed data
pal.files <- list.files(path=pal.path, pattern="cleaned$", recursive=T, full.names=F)
ptt.files <- list.files(path=ptt.path, pattern="ptt$", recursive=T, full.names=F)
genlen.all <- read.table(file=genlenfile, header=F, sep="\t")

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
#i'm not sure that 3 lines below needed...
all.df$nc.id <-as.character(all.df$nc.id)
all.df$ptt.file <-as.character(all.df$ptt.file)
all.df$pal.file <-as.character(all.df$pal.file)

#function returns length of genome with id=nc.id
GetGenLen <- function (nc.id=NULL) {
    genlen.subset <- subset(genlen.all, subset=(genlen.all$V1==nc.id));
    genlen<-genlen.subset$V2
    return(genlen)
}

#functon returns coordinates of feature from given .ptt file
GetPTTFile <- function (ptt.file = NULL, product = NULL){
    ptts <- try(fread(input=ptt.file))
    if ( any(class(ptts) == "try-error") ) {
        ptts <- NULL
        cat(paste(Sys.time(), "ERROR! in PTT reading!: ",ptt.file, sep="\t"), file=log.file, sep="\n", append=T)  
    }else{
        ptt.product.pattern <- paste0("*",product)
        ptts<-subset(x=ptts, grepl(ptt.product.pattern, ptts$Product))
        
        #Separate string starts..ends to starts and ends
        list.tmp <- strsplit(as.character(ptts$Location), "\\..")
        StartEnd <- ldply(list.tmp)
        
        if(nrow(StartEnd)>0) colnames(StartEnd) <- c("Start", "End")
        
        ptts$Start <- as.numeric(StartEnd$Start)
        ptts$End <- as.numeric(StartEnd$End)
        ptts$Size <- as.numeric(ptts$End-ptts$Start)
        
        ptts <- subset(ptts, ptts$Size>0) #filter 0-sized features
        class(ptts)<- "data.frame"
        
    }
    return(ptts)    
}

cat(paste('Time', 'nc.id','pval','hits.real','reps','genlen','maxcoord' sep="\t"), file=out.file, sep="\n", append=T)
pval.find <-function (nc.id=NULL, ptt.file=NULL, pal.file=NULL){
    cat(paste(Sys.time(), "[START:]: ",nc.id, ptt.file, pal.file, sep="\t"), file=log.file, sep="\n", append=T)

    PTTs <- GetPTTFile(ptt.file,feature)
    genlen <- GetGenLen(nc.id)
    PALs <- try(fread(input = pal.file, sep = "\t"))
    
    if ( any(class(PALs)!="try-error" ) & any(class(PTTs) != "NULL") & nrow(PTTs)>0){
        pals.ir <- reduce(IRanges(start=PALs$V1, end=PALs$V2))
        
        #Get flanks of feature from .ptt file
        leftside <- IRanges(start=PTTs$Start-flank_length, end=PTTs$Start)
        rightside <-IRanges(start=PTTs$End, end = PTTs$End+flank_length)
        flanks.ir <-reduce(c(leftside, rightside)) 
        
        
        hits.real <- sum(countOverlaps(query = pals.ir, subject = flanks.IR))
        pval <- 0
        
        maxcoord<-genlen - (max(PTTs$Size)+flank_length)
        reps <- nrow(PTTs)
        if (reps >0 & maxcoord >0){
            for (i in 1:1000){
                start.rnd <- sort(sample(maxcoord,reps))
                end.rnd <-start.rnd+PTTs$Size
                
                
                leftside <- IRanges(start=start.rnd-flank_length, end=start.rnd)
                rightside <-IRanges(start=end.rnd, end = end.rnd+flank_length)
                flanks.rand.ir <-reduce(c(leftside, rightside))
                
                hits.rand <-sum(countOverlaps(query=pals.ir, subject=flanks.rand.ir))
                if (hits.rand > hits.real) pval=pval+1
            }
            
        }else{
            cat(paste(Sys.time(), "ERROR! in reps or maxcoord: ",nc.id, ptt.file, pal.file, sep="\t"), file=log.file, sep="\n", append=T)    
            pval = -1
        }   
    }else{
        cat(paste(Sys.time(), "ERROR! in reading .pal or empty PTT: ",nc.id, ptt.file, pal.file, sep="\t"), file=log.file, sep="\n", append=T)
        pval= -2
    }
    cat(paste(Sys.time(), nc.id, pval, hits.real, reps, genlen, maxcoord, sep="\t"), file=out.file, sep="\n", append=T) #'Time', 'nc.id','pval','hits.real','reps','genlen','maxcoord'
    return(pval)
}

send_sms <- function(time.start, time.end){
    nodename = Sys.info()['nodename']
    filename = cur.res.name
    status = 'DONE'
    time.total = paste0(ceiling(difftime(time.end,time.start, units='mins')[[1]]), " mins")
    msg.text = paste(filename,nodename,status,time.total, sep = ' | ')

    api_id = 'dfdbb9ca-08dd-6954-cd3b-ee8556c9a9b5'
    phone = '79269326848'
    
    sms <- postForm(uri='http://sms.ru/sms/send', api_id = api_id, to = phone, text=msg.text, translit = 1)
    return(sms)
}


cat(paste(Sys.time(), "Making cluster!", sep="\t"), file=log.file, sep="\n", append=T)
cl.local<-makeCluster(getOption("cl.cores", 3))
clusterExport(cl=cl.local, varlist=c("all.df", "genlen.all", "pval.find", "GetGenLen","GetPTTFile", 'log.file','out.file', "feature"))
clusterEvalQ(cl.local,require(IRanges))
clusterEvalQ(cl.local,require(data.table))
clusterEvalQ(cl.local,require(plyr))
all.df$rslts<-parApply(cl.local, all.df[,c('nc.id','ptt.file','pal.file')], 1, function(x) pval.find(x[1],x[2],x[3]))
stopCluster(cl.local)
cat(paste(Sys.time(), "Stopping Cluster", sep="\t"), file=log.file, sep="\n", append=T)

write.table(x<-all.df, file=all.file, sep="\t")

worktime.finish = Sys.time() #Time when spript finished


cat(paste(Sys.time(), "DONE!", sep="\t"), file=log.file, sep="\n", append=T)
send_sms(worktime.start, worktime.finish)



