#! /usr/bin/Rscript 

rm(list=ls()) 

############################
# PART 1: prepare analyses #
############################

# collect arguments
args = commandArgs(TRUE)

# parse arguments (we expect the form --arg=value)
parseArgs = function(x) strsplit(sub("^--", "", x), "=")
argsL = as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) = as.data.frame(do.call("rbind", parseArgs(args)))$V1
args = argsL
rm(argsL)

# default setting when no all arguments passed or help needed
if("--help" %in% args | is.null(args$sync) | is.null(args$names) | is.null(args$o)) {
  cat("
      R SCRIPT FOR CONVERTING A SYNC FILE INTO A VCF FILE
      This script converts a sync file into a VCF file with a format that is compatible for SNP annotation with snpEff
      This script requires the following R packages to be installed: data.table, tidyr, stringr, YAPSA, rio

      Usage: Rscript ./sync2VCF.R --sync=dummy.sync --names=names.txt --o=test
  
      Mandatory arguments:
      --sync=File            - sync file (e.g., output PoPoolation2), for example: test.sync
      --names=File           - txt file with the names of the pools on the first line, separated by tabs
      --o=String             - name of the output file, for example: test
      --help                 - print this text

      ## Input sync file
          Sync file of the format given as output of PoPoolation2's mpileup2sync.jar script
          Order of columns: contig/chromosome, site position, reference allele, allele count for each pool in the format A:T:C:G:N:Del
          Pools are given in the same order as the maximum coverage .txt file
          Code assumes sync file only contains biallelic positions
          Format example:
          
              Contig1	77	A	0:0:0:0:0:0	16:0:1:0:0:0	0:0:0:0:0:0	3:0:0:2:0:0	0:0:0:0:0:0	17:0:0:0:0:0	0:30:0:0:0:0	2:0:0:0:0:0
              Contig1	93	C	0:0:0:0:0:0	0:0:16:4:0:0	0:0:0:0:0:0	0:0:5:0:0:0	0:1:0:0:0:0	0:0:17:0:0:0	0:0:0:0:0:0	0:0:2:0:0:0
              Contig1	124	G	2:0:0:122:0:0	1:0:0:46:0:0	3:0:0:19:0:0	9:0:0:68:0:0	2:0:0:20:0:0	5:0:0:60:0:0	9:0:0:56:0:0	7:0:0:60:0:0
              Contig1	131	G	2:0:0:1220:0:0	100:0:0:0:0:0	355:0:0:0:0:0	90:0:0:680:0:0	201:0:0:20:0:0	555:0:0:0:0:0	90:0:0:56:0:0	71:0:0:0:0:0
              
      ## Input txt file
          Text file with the names of the pools, separated by tabs
          Pool1   Pool2   Pool3   Pool4   Pool5   Pool6   Pool7   Pool8
              
      ## Notes on the output VCF files
          - REF = reference allele (based on the reference genome)
          - ALT = alternative alleles (can be one or more alleles)
          - INFO field contains the allele counts in sync format
          - The output consists of a VCF file of all SNPs (REF = major allele across all pools, ALT = minor allele across all pools),
            and separate VCF files for each pool (REF = major allele in the pool, ALT =  minor allele in the pool)

      ## Source
          - This script is written by Eveline Pinseel in November 2022, available from https://github.com/evelinepinseel
          - Some parts of the script were adapted from the PoolParty script 'r_frequency.R'
          - Source PoolParty script: https://github.com/StevenMicheletti/poolparty/tree/master/rscripts, downloaded in May 2022")

  q(save="no")
}

## print arguments used for analyses
cat("Sync file that will be converted to VCF: ", args$sync,"\n",sep="")
cat("Text file with pool names: ", args$names,"\n",sep="")
cat("Name of the output files: ", args$o,"\n",sep="")

# load packages
suppressMessages(require("data.table"))
suppressMessages(require("tidyr"))
suppressMessages(require("stringr"))
suppressMessages(require("YAPSA"))
suppressMessages(require("rio"))

# to run in Rstudio:
#sync = as.data.frame(fread(file="test-folder-filtered_all.sync", stringsAsFactors=FALSE, showProgress=FALSE))
#names = read.table("names.txt", header=TRUE)
#alertname = test

# get present working directory
outdir = getwd()

# import sync file
sync = as.data.frame(fread(file=args$sync, stringsAsFactors=FALSE, showProgress=FALSE))

# define output file
alertname = args$o

# get names of the pools
names = read.table(args$names, header=TRUE)
names = colnames(names)

# create header for the sync file
colnames(sync) = c('Chr', 'Pos', 'Ref', names)

paste0("STEP 1: done. The sync file has been succesfully imported. Analysis started.")

#################################
# PART 2: VCF for whole dataset #
#################################

# restructure the .sync file
## break into genotypes and genomic positions heading
counts = sync[, 4:(ncol(sync))]
gpops = sync[, 1:3] 
## convert to matrix for faster write
counts = as.matrix(counts)
## get the working directory, make temporary file
activedir = getwd()
tfile= paste0(outdir,"/","tempRf.tem")
## write temporary file of allele counts, this allows to break up colon-separated alleles
write.table(counts, file=tfile, sep=':', row.names=FALSE, col.names=FALSE, quote=FALSE)
## read counts back in separate columns
counts = fread(tfile, stringsAsFactors=FALSE, sep=':', showProgress=FALSE)
counts = as.data.frame(counts)
## remove temporary file
invisible(file.remove(tfile))
## rename row names
counts = as.data.frame(counts)
rownames = paste0(sync$Chr,':',sync$Pos,':',sync$Ref)
rownames(counts) = rownames
## get reference allele
REF = sync$Ref

# create .sync list for each population separately
## calculate number of populations
npops = ncol(counts)/6
## get summary info from file
totcol = ncol(counts)
thetseq = seq(6,totcol,6)
thebseq = seq(1,totcol,6)
## break up genotypes for each individual into separate dataframes
Var = list()
for (i in (1:npops)) {
  Var[[i]] = counts[, c(thebseq[i]:thetseq[i], ncol(counts))]
  gc()
}
## rename column names
colnames = c("A", "T", "C", "G", "N", "Del")
Var = lapply(Var, setNames, colnames)

# get a list of the major allele for each position
## add in column with dummy data
Var2 = Var
for (i in (1:npops)) {
  V = as.data.frame(Var2[[i]])
  V$E = rep(0,nrow(as.data.frame(Var2[i])))
  Var2[[i]] = V
}
## reshuffle columns
Var2 = Map(function(V) {
  V = subset(V, select=c(7,1:6))
  V
}, Var2)
## get the major allele at each position
MaxC = list()
for (i in (1:npops)) {
  MaxC[[i]] = colnames(Var2[[i]])[max.col(Var2[[i]], ties.method="first")]
  gc()
}
## reset all 'E' values (corresponds with the extra column that was added in), to zero
MaxC = Map(function(V) {
  V = ifelse(V == 'E', NA, V)
  V
}, MaxC)
## make nucleotide list into data frame, printing the primary allele for each population
### when no coverage in a sample, print NA
### when two alleles are as common: take first
MaxC = as.data.frame(MaxC)
colnames(MaxC) = c(1:npops)
MaxC = as.matrix(MaxC)
## turn into numerical matrix
MaxN = chartr("ATCGND", "123456", MaxC)
MaxN = matrix(as.numeric(MaxN), ncol = ncol(MaxN))
## create function
modefunc = function(x){
  tabresult = tabulate(x)
  themode = which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1) themode = themode[[1]]
  return(themode)
}
## find the major allele
MaxN = apply(MaxN, 1, modefunc)
MaxN = as.matrix(MaxN)
MAX = chartr("123456", "ATCGND", MaxN)
rm(Var2) ; invisible(gc())

# get a list of the minor allele for each position (= allele that is not the major allele)
## reset major alleles to zero
VarS = sum_over_list_of_df(Var)
cols = as.vector(MaxN)
rows = seq_len(nrow(VarS))
VarS[cbind(rows, cols)] = 0
## get the minor allele
MIN = colnames(VarS)[max.col(VarS, ties.method="first")]
rm(VarS) ; invisible(gc())

# get IDs from the sync file
CHROM = sync[,1]
POS = sync[,2]

# create vector with alternative allele(s)
## create dataframe with the reference, major, and minor alleles
df = as.data.frame(cbind(REF, MAX, MIN))
## reset major or minor allele to NA if major/minor allele = reference allele
MAX.r = ifelse(df$MAX == df$REF, NA, df$MAX)
MIN.r = ifelse(df$MIN == df$REF, NA, df$MIN)
ALT = paste0(MAX.r, ', ', MIN.r)
## remove NAs from the vector with alternative alleles
ALT = gsub("NA, ", "", ALT) 
ALT = gsub(", NA", "", ALT) 

# create ID and QUAL columns
ID = rep('.', length=nrow(sync)) 
QUAL = rep('.', length=nrow(sync)) 

# combine sync columns into one column
sync2 = sync %>% unite("INFO", c(4:length(sync)), sep= " ", remove = FALSE)
INFO = paste0('SYNC=', sync2[,c(4)])

# create VCF
VCF = as.data.frame(cbind(CHROM, POS, ID, REF, ALT, QUAL, INFO))
colnames(VCF) = c('#CHROM', 'ID', 'POS', 'REF', 'ALT', 'QUAL', 'INFO')

# export VCF files for each pool
export(VCF, file = paste0(alertname, "_all.txt"))

########################################
# PART 3: VCF for each individual pool #
########################################

# remove positions that are not SNPs (pool-by-pool)
Var = Map(function(V) {
  V$count = apply(V, 1, function(x) length(which(x>0)))
  rm = rownames(V[V$count == 0 | V$count == 1, ])
  V = V[! rownames(V) %in% rm, ]
  V = V[,c(1:6)]
  V
}, Var)

# get reference allele for each SNP (pool-by-pool)
REF = Map(function(V) {
  V = as.data.frame(str_split_fixed(rownames(V), ':', 3))
  V = V[,3]
  V
}, Var)

# get a list of the major allele for each pool
MAX = list()
for (i in (1:npops)) {
  MAX[[i]] = colnames(Var[[i]])[max.col(Var[[i]], ties.method="first")]
  gc()
}

# get a list of the minor allele for each pool
## reset all zero values to Inf
Var2 = Var
for (i in (1:npops)) {
  V = as.data.frame(Var2[[i]])
  V[V == 0] = Inf
  Var2[[i]] = V
}
## get the minor allele
min.col = function(m, ...) max.col(-m, ...)
MIN = list()
for (i in (1:npops)) {
  MIN[[i]] = colnames(Var2[[i]])[min.col(Var2[[i]], ties.method="last")]
  gc()
}
rm(Var2) ; invisible(gc())

# create list of vectors with alternative allele(s)
ALT = Map(function(REF, MAX, MIN) {
  df = as.data.frame(cbind(REF, MAX, MIN))
  MAX.r = ifelse(df$MAX == df$REF, NA, df$MAX)
  MIN.r = ifelse(df$MIN == df$REF, NA, df$MIN)
  ALT = paste0(MAX.r, ', ', MIN.r)
  ALT = gsub("NA, ", "", ALT) 
  ALT = gsub(", NA", "", ALT) 
}, REF, MAX, MIN)

# create simple VCF format for each pool
VCF = Map(function(V, REF, ALT) {
  A = V[,1]
  T = V[,2]
  C = V[,3]
  G = V[,4]
  N = V[,5]
  Del = V[,6]
  INFO = paste0('SYNC=', A,':',T,':',C,':',G,':',N,':',Del)
  POS = rownames(V)
  POS = as.data.frame(str_split_fixed(POS, ':', 3))
  POS = POS[,c(1:2)]
  ID = rep('.', length=nrow(V)) 
  VCF = cbind(POS, ID, REF, ALT, ID, INFO)
  colnames(VCF) = c('#CHROM', 'ID', 'POS', 'REF', 'ALT', 'QUAL', 'INFO')
  VCF
}, Var, REF, ALT)

# give names to list of dataframes
names(VCF) = names

# export VCF files for each pool
export_list(VCF, file = paste0(alertname, "_%s.txt"))

paste0("All files have been created.")
