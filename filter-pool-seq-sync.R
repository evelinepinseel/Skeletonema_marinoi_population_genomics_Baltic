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
if("--help" %in% args | is.null(args$sync) | is.null(args$minCOV) | is.null(args$maxCOV) |  is.null(args$MAC) | is.null(args$MAF) | is.null(args$o)) {
  cat("
      R SCRIPT FOR FILTERING POOL-SEQ SYNC FILE
      This script filters a sync file using the specified parameters, and calculates allele frequencies
      This script requires the following R packages to be installed: tidyr, data.table, stringr, matrixStats, YAPSA, gtools

      Usage: Rscript ./filter-pool-seq-sync.R --sync=dummy.sync --minCOV=20 --maxCOV=dummy-depth.txt --MAC=4 --MAF=0.001 --o=/test-folder
  
      Mandatory arguments:
      --sync=File            - sync file (e.g., output PoPoolation2), for example: test.sync
      --minCOV=Integer       - minimum coverage, for example: 10
      --maxCOV=File          - txt file specifying maximum coverage filter per contig per pool, for example: test.txt
      --MAC=Integer          - minimum allele count (MAC), for example: 4
      --MAF=Proportion       - minimum allele frequency (MAF) as an integer (e.g., 0.01 = 1%), for example: 0.05
      --o=String             - name of the output folder (as a subfolder of the working directory), for example: /test-folder
      --help                 - print this text

      ## Input sync file
          Sync file of the format given as output of PoPoolation2's mpileup2sync.jar script
          Order of columns: contig/chromosome, site position, reference allele, allele count for each pool in the format A:T:C:G:N:Del
          Pools are given in the same order as the maximum coverage .txt file
          Format example:
      
              Contig1	77	A	0:0:0:0:0:0	16:0:1:0:0:0	0:0:0:0:0:0	3:0:0:2:0:0	0:0:0:0:0:0	17:0:0:0:0:0	0:30:0:0:0:0	2:0:0:0:0:0
              Contig1	93	C	0:0:0:0:0:0	0:0:16:4:0:0	0:0:0:0:0:0	0:0:5:0:0:0	0:1:0:0:0:0	0:0:17:0:0:0	0:0:0:0:0:0	0:0:2:0:0:0
              Contig1	124	G	2:0:0:122:0:0	1:0:0:46:0:0	3:0:0:19:0:0	9:0:0:68:0:0	2:0:0:20:0:0	5:0:0:60:0:0	9:0:0:56:0:0	7:0:0:60:0:0
              Contig1	131	G	2:0:0:1220:0:0	100:0:0:0:0:0	355:0:0:0:0:0	90:0:0:680:0:0	201:0:0:20:0:0	555:0:0:0:0:0	90:0:0:56:0:0	71:0:0:0:0:0

      ## Minimum coverage filter
          The minimum coverage filter is applied to each pool separately
          For example, when for a given site pool A does not pass the minimum coverage filter, but pool B does: the site is removed (reset to NA) for pool A, but not pool B
          For example, when minCOV > 20, only pools that have at least 20 reads for a site will be retained 
          When for a given site, none of the pools pass the minimum coverage filter, the entire site is removed
          
      ## Maximum coverage filter
          The maximum coverage filter is applied to each site across all pools
          As soon as one pool does not pass the maximum coverage filter, the entire site is removed
          Maximum coverage filtering is done using a custom threshold for each contig in each pool (e.g., twice the average sequencing depth across a contig for each pool)
          This is achieved by giving a .txt file with the desired maximum coverage threshold to this script of the format:
      
                          P1      P2      P3      P4      P5      P6     P7       P8  
              Contig1     256	  184	  160	  181	  193	  133	  148	  165	
              Contig2	  220	  147	  138	  140	  155	  109	  118	  136	
              Contig3	  322	  213	  213	  224	  192	  125	  146	  171	
              
          In above example, a site of contig 1 will only pass the maximum coverage filter if the number of reads is 
          lower than or equal to 256 in P1, 184 in P2, 160 in P3, 181 in P4, 193 in P5, 133 in P6, 148 in P7, and 165 in P8
          If only one maximum coverage threshold is desired for the entire dataset, the input .txt file can be adjusted. 
          For example, for a maximum coverage = 250:
        
                          P1      P2      P3      P4      P5      P6      P7    P8
              Contig1     250	  250	  250	  250	  250	  250	  250   250
              Contig2	  250	  250	  250	  250	  250	  250	  250   250	
              Contig3	  250	  250	  250	  250	  250	  250	  250   250	
        
      ## Minimum allele count (MAC)
          The MAC filter is applied across all libraries combined
          For example, if MAC = 4, a site will only be considered as polymorphic if the cumulative count of the 
          minor allele across all pools in the input .sync file equals at least 4
          In case of multiallelic (> 2 alleles) sites where one allele is below the MAC threshold: 
          the allele counts of the allele that did not pass the MAC will be removed, transforming the site into a biallelic site
          
      ## Minimum allele frequency (MAF)
          This script creates an output file with sites that did not pass the MAF filter
          Sites that did not pass the MAF filter are also removed from the dataset
          
      ## Disclaimer on filtering conditions
          - Filtering is done in this order: minimum coverage, maximum coverage, MAC, multiallelic positions, removal of positions that are no SNPs
          - Multiallelic sites (> 2 alleles) are removed from the dataset
      
      ## Output files
          - Two .sync files: 
                  '-filtered_all.sync' contains all sites that passed filtering
                  '-filtered_complete.sync' contains all complete sites that passed filtering
                  
          - Two .fz files with allele frequencies: 
                  '-filtered_all.fz' contains all sites that passed filtering
                  '-filtered_complete.fz' contains all complete sites that passed filtering
            Format of the .fz files is: (Chr = contig/chromosome, Pos = site/position, Ref = reference allele, 
            MaA =  major allele, MiA = minor allele, followed by the allele frequency of the major allele across all pools): 
              
              Chr       Pos   Ref   MaA   MiA   P1    P2    P3    P4    P5    P6    P7    P8
              Contig1   124   G     G     A     0.984 0.979 0.864 0.883 0.909 0.923 0.862 0.896 
              Contig1   132   G     G     T     1     1     1     0.944 1     1     1     1
              
          - .txt files with positions that did not pass the filters:
                -'_min_coverage-fail.txt': sites that contained pools that did not pass the minimum coverage filter
                -'_max_coverage-fail.txt': sites that contained pools that did not pass the maximum coverage filter
                -'_MAC-fail.txt': sites that contained allales that did not pass the MAC filter
                -'_multiallelic-fail.txt': multiallelic (> 2 alleles) sites
                -'_other-sites-fail.txt': sites that were not SNPs after all filtering steps
                -'_MAF-fail.txt': sites that failed the MAF filter

      ## Source
          - This script is written by Eveline Pinseel in October 2022, available from https://github.com/evelinepinseel
          - Some parts of the script were adapted from the PoolParty script 'r_frequency.R'
          - Source PoolParty script: https://github.com/StevenMicheletti/poolparty/tree/master/rscripts, downloaded in May 2022")

  q(save="no")
}

## print arguments used for analyses
cat("Sync file that will be filtered: ", args$sync,"\n",sep="")
cat("Applied minimum coverage threshold: ", args$minCOV,"\n",sep="")
cat("Applied maximum coverage threshold taken from file: ", args$maxCOV,"\n",sep="")
cat("Applied MAC threshold: ", args$MAC,"\n",sep="")
cat("Applied MAF threshold: ", args$MAF,"\n",sep="")
cat("Data will be stored in the output folder: ", args$o,"\n",sep="")

# load packages
suppressMessages(require("tidyr"))
suppressMessages(require("data.table"))
suppressMessages(require("stringr"))
suppressMessages(require("matrixStats"))
suppressMessages(require("YAPSA"))
suppressMessages(require("gtools"))

# to run in Rstudio:
#sync = as.data.frame(fread(file="dummy.sync", stringsAsFactors=FALSE, showProgress=FALSE))
#minCOV = 40
#maxCOV = read.table(file="dummy-depth.txt", header = TRUE)
#MAC = 2
#MAF = 0.05
#alertname = "/test-folder2"

# get present working directory
outdir = getwd()

# import sync file
sync = as.data.frame(fread(file=args$sync, stringsAsFactors=FALSE, showProgress=FALSE))
colnames(sync) = c("Chr", "Pos", "Ref", names)
rownames = paste0(sync$Chr,':',sync$Pos,':',sync$Ref)
rownames(sync) = rownames

# define filtering specifications
minCOV = as.numeric(args$minCOV)
maxCOV = read.table(file=args$maxCOV, header = TRUE)
MAC = as.numeric(args$MAC)
MAF = as.numeric(args$MAF)

# define output folder name
alertname = args$o

# get names of the pools from the maximum coverage file
names = colnames(maxCOV)

# get total number of SNPs
tot.snp = nrow(sync)

paste0("STEP 1: done. ", tot.snp, " SNPs have been succesfully imported. Analysis started.")

################################
# PART 2: filter the sync file #
################################

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
write.table(counts, file=tfile, sep=':',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
## read counts back in separate columns
counts = fread(tfile, stringsAsFactors=FALSE, sep=':', showProgress=FALSE)
## remove temporary file
invisible(file.remove(tfile))
## rename row names
counts = as.data.frame(counts)
rownames = paste0(sync$Chr,':',sync$Pos,':',sync$Ref)
rownames(counts) = rownames

# prepare maxCOV dataframe for max coverage filtering
## get list of all contig names
contig.list = rownames(maxCOV)
contig.length = list()
for (i in contig.list){
  tmp = nrow(sync[grepl(i, sync$Chr), ])
  contig.length[[i]] = tmp
}
contig.length = as.vector(unlist(contig.length))
## multiply rows in max coverage dataframe with contig.length
maxcov.m = data.frame()
for (i in c(1:length(contig.length))) {
  tmp = maxCOV[rep(i, contig.length[i]),]
  maxcov.m = rbind(maxcov.m, tmp)
}

# remove redundant objects
rm(sync) ; invisible(gc())
rm(contig.length) ; invisible(gc())

# create .sync list for each population separately
## calculate number of populations
npops = ncol(counts)/6
## get summary info from file
totcol = ncol(counts)
thetseq = seq(6,totcol,6)
thebseq = seq(1,totcol,6)
L = nrow(counts) #total number of SNPs in the file
## add column with unique numeric position ID (needed for downstream coverage filtering)
counts$Pos = c(1:nrow(counts))
## break up genotypes for each individual into separate dataframes
Var = list()
for (i in (1:npops)) {
  Var[[i]] = counts[, c(thebseq[i]:thetseq[i], ncol(counts))]
  gc()
}
## rename column names
colnames = c("A", "T", "C", "G", "N", "Del", "Pos")
Var = lapply(Var, setNames, colnames)

# calculate total coverage for each position in each population
for (i in c(1:npops)) {
  V = as.data.frame(Var[i])
  sum = rowSums(V[,c(1:6)])
  Var[[i]]$cov = sum }

# add information on max allowed coverage to Var
Var = Map(function(V, pop) {
  V$max.cov = maxcov.m[,pop]
  V
}, Var, c(1:npops))
rm(maxcov.m) ; invisible(gc())

# remove positions with too low coverage within each pool separately
## get positions with too low coverage in each individual pool
pos.min = list()
for (i in c(1:npops)) {
  df = as.data.frame(Var[i])
  tmp = rownames(df[df$cov < minCOV & df$cov > 0, ])
  tmp2 = df[rownames(df) %in% tmp, ]
  pos.min[[i]] = tmp2$Pos  
}
## remove allele counts of positions from pools with too low coverage
Var = Map(function(V, mic) {
  V$A = ifelse(as.numeric(V$Pos) %in% mic, 0, V$A)
  V$T = ifelse(as.numeric(V$Pos) %in% mic, 0, V$T)
  V$C = ifelse(as.numeric(V$Pos) %in% mic, 0, V$C)
  V$G = ifelse(as.numeric(V$Pos) %in% mic, 0, V$G)
  V$N = ifelse(as.numeric(V$Pos) %in% mic, 0, V$N) 
  V$Del = ifelse(as.numeric(V$Pos) %in% mic, 0, V$Del)
  V
}, Var, pos.min)
## get list of low coverage SNP IDs
pos.min.all = unique(sort(unlist(pos.min)))
## print info on minimum coverage filter
tmp = as.data.frame(Var[[1]])
min.list = rownames(tmp[tmp$Pos %in% pos.min.all, ])
pro.min.cov = round((length(min.list) / tot.snp * 100),2)
paste0("Filtering details: ", length(min.list), " SNPs (", pro.min.cov, "%) are from pools that do not pass the min coverage filter ", minCOV, ", and allele counts for these pools have been removed.")

# remove positions with too high coverage within each population, separately for each contig and pool
## get positions with too high coverage in each individual pool
pos.max = list()
for (i in c(1:npops)) {
  df = as.data.frame(Var[i])
  tmp = rownames(df[df$cov > df$max.cov, ])  
  tmp2 = df[rownames(df) %in% tmp, ]
  pos.max[[i]] = tmp2$Pos
}
## remove allele counts of positions from pools with too low coverage
Var = Map(function(V, mxc) {
  V$A = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$A)
  V$T = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$T)
  V$C = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$C)
  V$G = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$G)
  V$N = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$N) 
  V$Del = ifelse(as.numeric(V$Pos) %in% mxc, 0, V$Del)
  V
}, Var, pos.max)
## get list of low coverage SNP IDs
pos.max.all = unique(sort(unlist(pos.max)))
## print info on minimum coverage filter
tmp = as.data.frame(Var[[1]])
max.list = rownames(tmp[tmp$Pos %in% pos.max.all, ])
pro.max.cov = round((length(max.list) / tot.snp * 100),2)
paste0("Filtering details: ", length(max.list), " SNPs (", pro.max.cov, "%) are from pools that do not pass the max coverage filter, and allele counts for these pools have been removed")

# remove coverage columns
for (i in c(1:npops)) {
  Var[[i]]$cov = NULL
  Var[[i]]$max.cov = NULL}

# remove alleles with Minimum Allele Count (MAC) below threshold
## calculate number of alleles below the MAC
MACF = seq.int(1, MAC - 1, by = 1)
## combine all dataframes into one
Var2 = do.call('rbind', Var)
## take the sum for each position across all populations
Var2 = aggregate(Var2[,1:6], by=list(Var2$Pos), sum)
colnames(Var2)[1] = 'Pos'
Var2$Pos = as.numeric(Var2$Pos)
Var2 = Var2[order(Var2$Pos),]
## select positions with alleles that do not pass the MAC filter
MAC_fail_A = Var2[Var2$A %in% MACF, ]$Pos
MAC_fail_T = Var2[Var2$T %in% MACF, ]$Pos
MAC_fail_C = Var2[Var2$C %in% MACF, ]$Pos
MAC_fail_G = Var2[Var2$G %in% MACF, ]$Pos
MAC_fail_N = Var2[Var2$N %in% MACF, ]$Pos
MAC_fail_Del = Var2[Var2$Del %in% MACF, ]$Pos
## transform vectors into lists
A_list = rep(list(MAC_fail_A), npops)
T_list = rep(list(MAC_fail_T), npops)
C_list = rep(list(MAC_fail_C), npops)
G_list = rep(list(MAC_fail_G), npops)
N_list = rep(list(MAC_fail_N), npops) 
Del_list = rep(list(MAC_fail_Del), npops)
rm(Var2) ; invisible(gc())
## loop over the .sync list to reset positions that failed the MAC threshold to zero
Var = Map(function(V, AF, TF, CF, GF, NF, DelF) {
  V$A = ifelse(as.numeric(V$Pos) %in% AF, 0, V$A)
  V$T = ifelse(as.numeric(V$Pos) %in% TF, 0, V$T)
  V$C = ifelse(as.numeric(V$Pos) %in% CF, 0, V$C)
  V$G = ifelse(as.numeric(V$Pos) %in% GF, 0, V$G)
  V$N = ifelse(as.numeric(V$Pos) %in% NF, 0, V$N) 
  V$Del = ifelse(as.numeric(V$Pos) %in% DelF, 0, V$Del)
  V
}, Var, A_list, T_list, C_list, G_list, N_list, Del_list)
## get a full list of all positions that contained alleles that failed the MAC
pos.mac = unique(sort(c(unlist(A_list), unlist(T_list), unlist(C_list), 
                        unlist(G_list), unlist(N_list), unlist(Del_list))))
tmp = Var[[1]]
mac.list = rownames(tmp[tmp$Pos %in% pos.mac, ])
## print info on MAC filter
pro.mac = round((length(mac.list) / tot.snp * 100),2)
paste0("Filtering details: ", length(mac.list), " SNPs (", pro.mac, "%) contained alleles that did not pass the MAC filter. These alleles were removed.")

# remove column with position numbers
for (i in c(1:npops)) {
  Var[[i]]$Pos = NULL}

# reduce dataset to biallelic alleles only
#alertname = sub('.*/', '', alertname)
## determine polymorphic sites (>2 alleles)
SumA = sum_over_list_of_df(Var)
## get positions with >2 alleles
SumA$count = apply(SumA, 1, function(x) length(which(x>0)))
PolY = rownames(SumA[SumA$count > 2, ])  
## remove positions with >2 alleles
Var  = Map(function(V) {
  Var[[i]] = V[! rownames(V) %in% PolY, ]
}, Var)
## print information on SNPs that were not biallelic
multi.other = rownames(SumA[rownames(SumA) %in% PolY, ])
pro.multi = round((length(PolY) / tot.snp * 100),2)
paste0("Filtering details: ", length(PolY), " SNPs (", pro.multi, "%) represented multiallelic positions, and have been removed.")

# remove positions that are no longer SNPs post filtering
## summarize all populations: total coverage of each position
SumA = sum_over_list_of_df(Var)
## get positions with no or only one allele
SumA$count = apply(SumA, 1, function(x) length(which(x>0)))
pos.zero.one = rownames(SumA[SumA$count == 0 | SumA$count == 1, ])  
## remove positions with no or only one allele
Var  = Map(function(V) {
  Var[[i]] = V[! rownames(V) %in% pos.zero.one, ]
}, Var)
## above code does not remove positions where only one pool has coverage, but has two alleles: remove these also
### identify positions that still need to be removed
SumB = list()
SumB = Map(function(V) {
  SumB[[i]] = rowSums(V[,1:6])
}, Var)
### create list of positions that still need to be removed
SumB = as.data.frame(do.call(cbind, SumB)) 
SumB$count = apply(SumB, 1, function(x) length(which(x>0)))
pos.one = rownames(SumB[SumB$count == 0 | SumB$count == 1, ])  
### remove positions
Var  = Map(function(V) {
  Var[[i]] = V[! rownames(V) %in% pos.one, ]
}, Var)
## print information on positions removed in this step
tmp = sort(c(pos.zero.one, pos.one))
sub.other = rownames(SumA[rownames(SumA) %in% tmp, ])
pro.other = round((length(sub.other) / tot.snp * 100),2)
paste0("Filtering details: ", length(sub.other), " SNPs (", pro.other, "%) were removed after all other filtering steps due to not being a SNP post-filtering.")

paste0("STEP 2: done. SNPs have been succesfully filtered.")

########################################
# PART 3: calculate allele frequencies #
########################################

# get list of reference alleles
Ref.list = list()
Ref.list = Map(function(V) {
  tmp = as.data.frame(V)
  tmp = tibble::rownames_to_column(tmp, "Position")
  tmp = str_split_fixed(tmp$Position, ':', 3)
  Ref.list[[i]] = as.data.frame(cbind(tmp, V))
}, Var)
colnames = c("Chr", "Pos", "Ref", "A", "T", "C", "G", "N", "Del")
Ref.list = lapply(Ref.list, setNames, colnames)

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
ComPos = as.data.frame(chartr("123456", "ATCGND", MaxN))
rm(Var2) ; invisible(gc())

paste0("STEP 3: list of major alleles created.")

# get a list of the minor allele for each position (= allele that is not the major allele)
## reset major alleles to zero
VarS = sum_over_list_of_df(Var)
cols = as.vector(MaxN)
rows = seq_len(nrow(VarS))
VarS[cbind(rows, cols)] = 0
## get the minor allele
RarPos = as.data.frame(colnames(VarS)[max.col(VarS, ties.method="first")])
rm(VarS) ; invisible(gc())

paste0("STEP 3: list of minor alleles created.")

# calculate major allele frequencies
## get total coverage at positions 
Sumz = lapply(Var, function(x) rowSums(x))
## make into data matrix
for (i in (1:npops)) {
  Var[[i]] = as.matrix(Var[[i]]) 
  }
## calculate frequency of major allele
Maxz = lapply(Var, function(x) rowMaxs(x))
Freqz = mapply("/",Maxz,Sumz,SIMPLIFY = FALSE)  
## merge into frequency table of major allele
Freqz = as.data.frame(Freqz)
colnames(Freqz) = names
Freqz = as.matrix(Freqz)

paste0("STEP 3: allele frequencies of major allele calculated.")

# using the major allele, create comparisons index based on matching alleles
## make index, TRUE if the reference major allele is different in population x
## index based on most common allele
MaxC[is.na(MaxC)] = 'X'
idx = (!apply(MaxC, 2, function(x) x == ComPos))

# calculate final allele frequencies
## if the major isn't same between populations, make it the minor allele
Freqz[idx] = 1 - Freqz[idx]
Freqz = round(Freqz, 3)

paste0("STEP 3: final allele frequencies calculated.")

# minor allele calculations
## get minor allele frequency failing positions
mAx = rowMaxs(Freqz, na.rm = TRUE)
mIn = rowMins(Freqz, na.rm = TRUE)
rAnge = as.vector(mAx - mIn)
rAnge = as.matrix(rAnge)
mIn = as.matrix(mIn)
mAx = as.matrix(mAx)
rownames(rAnge) = c()
colnames(rAnge) = "MAFT"
rownames(mAx) = c()
colnames(mAx) = "MAX"
rownames(mIn) = c()
colnames(mIn) = "MIN"
## obtain positions that do not pass the MAF
IDs = as.data.frame(Ref.list[1])[,c(1:3)]
MAFFAIL = cbind (IDs, Freqz, rAnge, mAx, mIn)
MAFFAIL = as.data.frame(MAFFAIL)
MAFFAIL = subset(MAFFAIL[,1:3], ((MAFFAIL$MAX >= (1-MAF) | MAFFAIL$MIN <= MAF) & (MAFFAIL$MAFT <= MAF)) )

paste0("STEP 3: MAF filter has been calculated.")

# create a table with allele frequencies
Freqz.F = cbind(IDs, ComPos, RarPos, Freqz)
colnames(Freqz.F) = c("Chr", "Pos", "Ref", "Maj", "Min", names)

# remove positions that did not pass MAF filter
Freqz.F = Freqz.F[! rownames(Freqz.F) %in% rownames(MAFFAIL), ]

## print information on MAF filter
pro.maf = round((nrow(MAFFAIL) / tot.snp * 100),2)
paste0("Output: ", nrow(MAFFAIL), " SNPS (", pro.maf, "%) do not pass additional population MAF threshold of ", MAF, ", and have been removed.")

paste0("STEP 3: done. Allele frequencies have been calculated.")

##################################
# PART 4: clean and export files #
##################################

# make a frequency table with the filtered data, including positions with missing data
outname = paste0(paste0(outdir, alertname, '/'), alertname, "-filtered_all.fz")
write.table(Freqz.F, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)

# make a frequency table with the filtered data, only including fully resolved positions
Freqz.C = na.omit(Freqz.F)
outname = paste0(paste0(outdir, alertname, '/'), alertname, "-filtered_complete.fz")
write.table(Freqz.C, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)

# rearrange and export fully filtered .sync file
## combine all data frames into one dataframe in .sync format
df = NULL
for (i in c(1:npops)) {
  A = Var[[i]][,1]
  T = Var[[i]][,2]
  C = Var[[i]][,3]
  G = Var[[i]][,4]
  N = Var[[i]][,5]
  Del = Var[[i]][,6]
  x = paste0(A,':',T,':',C,':',G,':',N,':',Del)
  df = cbind(df, x) }
##add column identifiers
colnames(df) = names
## add row identifiers
Fsync = cbind(IDs, df)
## remove sites that did not pass the MAF filter
Fsync = Fsync[! rownames(Fsync) %in% rownames(MAFFAIL), ]
##export file of all positions that passed filtering (including positions with missing data)
outname = paste0(paste0(outdir, alertname, '/'), alertname, "-filtered_all.sync")
write.table(Fsync, file=outname,
            row.names=FALSE, col.names=FALSE, quote=FALSE)
##export file of fully resolved positions that passed filtering
Fsync.C = Fsync[rownames(Fsync) %in% rownames(Freqz.C), ]
outname = paste0(paste0(outdir, alertname, '/'), alertname, "-filtered_complete.sync")
write.table(Fsync.C, file=outname,
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# export files with sites that did not pass filtering
## minimum coverage filter
## get list of low coverage SNP IDs
tmp = as.data.frame(min.list)
sub.min = as.data.frame(str_split_fixed(tmp[,1], ':', 3))
colnames(sub.min) = c('Chr', 'Pos', 'Ref')
outname = paste0(paste0(outdir, alertname, '/'), alertname, "_min-coverage-fail.txt")
write.table(sub.min, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)

## maximum coverage filter
tmp = as.data.frame(max.list)
sub.max = as.data.frame(str_split_fixed(tmp[,1], ':', 3))
colnames(sub.max) = c('Chr', 'Pos', 'Ref')
outname = paste0(paste0(outdir, alertname, '/'), alertname, "_max-coverage-fail.txt")
write.table(sub.max, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)
## MAC
tmp = as.data.frame(mac.list)
sub.mac = as.data.frame(str_split_fixed(tmp[,1], ':', 3))
colnames(sub.mac) = c('Chr', 'Pos', 'Ref')
outname = paste0(paste0(outdir, alertname, '/'), alertname, "_MAC-fail.txt")
write.table(sub.mac, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)

## MAF
outname = paste0(paste0(outdir, alertname, '/'), alertname, "_MAF-fail.txt")
write.table(MAFFAIL, file=outname,
            row.names=FALSE, col.names=FALSE, quote=FALSE)

## multiallelic positions
tmp = as.data.frame(PolY)
PolY = as.data.frame(str_split_fixed(tmp[,1], ':', 3))
colnames(PolY) = c("Chr", "Pos", "Ref")
outname = paste0( paste0(outdir, alertname, '/'), alertname, "_multiallelic.txt") 
write.table(PolY, file=outname, quote=FALSE, row.names=FALSE)

## other removed sites
tmp = as.data.frame(sub.other)
sub.other = as.data.frame(str_split_fixed(tmp[,1], ':', 3))
colnames(sub.other) = c('Chr', 'Pos', 'Ref')
outname = paste0(paste0(outdir, alertname, '/'), alertname, "_other-sites-fail.txt")
write.table(sub.other, file=outname,
            row.names=FALSE, col.names=TRUE, quote=FALSE)

paste0("STEP 4: done. All output files have been created.")

# print the total number of filtered SNPs
pro.totF = round((nrow(Freqz.F) / tot.snp * 100), 2)
pro.totC = round((nrow(Freqz.C) / tot.snp * 100), 2)
paste0(nrow(Freqz.F), " SNPS (", pro.totF, "%) passed all filtering steps.")
paste0(nrow(Freqz.C), " SNPS (", pro.totC, "%) passed all filtering steps and have no missing data.")

paste0("Analysis finished")
