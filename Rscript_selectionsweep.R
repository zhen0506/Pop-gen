library(tidyverse)
library(qqman)
temp = list.files(pattern="*",path = "xpclr/")
myfiles = lapply(paste0("xpclr/",temp), function(x) read.table(x, header = T,sep = "\t") )
temp
xpclr = bind_rows(myfiles, .id = "a")
pops = data.frame(file=temp,a = as.character(1:length(temp))) %>% mutate(pop = str_extract(string =file,pattern = "\\w{2}(?=_eco)"))
chrs= data.frame(chrom = paste0("chr_",rep(1:7,each =4),rep(c("A","B","C",'D'), 7)),chr = 1:28)
xpclr = xpclr %>% left_join(chrs, by="chrom") %>% left_join(pops, by= "a")
xpclr_UF = xpclr %>% filter(pop == "UF")
xpclr_UC = xpclr %>% filter(pop == "UC")


#manhattan plot
max(xpclr$xpclr,na.rm = T)
UF_q95 = quantile(xpclr_UF$xpclr,na.rm = T,probs = 0.95)
UF_q99 = quantile(xpclr_UF$xpclr,na.rm = T,probs = 0.99)
manhattan(xpclr_UF , chr="chr", bp="start", snp="id", p="xpclr", logp=F, ylim = c(0,285),cex=0.75,
          genomewideline = UF_q99, suggestiveline = UF_q95, chrlabs =chrs$chrom , ylab = "XP-CLR",cex.axis=0.75)

UC_q95 = quantile(xpclr_UC$xpclr,na.rm = T,probs = 0.95)
UC_q99 = quantile(xpclr_UC$xpclr,na.rm = T,probs = 0.99)
manhattan(xpclr_UC , chr="chr", bp="start", snp="id", p="xpclr", logp=F, ylim = c(0,285),cex=0.75,
          genomewideline = UC_q99, suggestiveline = UC_q95, chrlabs =chrs$chrom , ylab = "XP-CLR",cex.axis=0.75)

##################################################################
##                   Identify selection sweep                   ##
##                          with xpclr                          ##
##################################################################

sweepID = function(xpclr_UF, prob=0.95,limit = 400000,genome_size = 780000000) {
  ###Identify selection sweep
  ###Three inputs are xpclr output, quantile probility and max gap size to merge a sweep region
  ###Output boundaries of selection sweeps and the region with the highest xpclr value
  UF_q95 = quantile(xpclr_UF$xpclr,na.rm = T,probs = prob)
  xpclr_UF_95 = xpclr_UF %>% filter(xpclr>UF_q95) %>% mutate(dis_f = pos_start-lag(pos_start), dis_n=lead(pos_start)-pos_start)
  #plot(hist(xpclr_UF_95 %>% filter(dis_f>0,dis_f < 1000000) %>% pull(dis_f),na.rm = T,breaks = 50))
  sweep_uf_95 = xpclr_UF_95 %>% mutate(sweep = cumsum(ifelse(dis_f>limit | dis_f < 0 | is.na(dis_f) , 1,0)))
  Ssweep_UF = sweep_uf_95 %>% group_by(chrom, sweep) %>% 
    summarise(region = n(), start=min(pos_start),end=max(pos_stop))  %>% 
    filter(region>2) %>% ungroup()
  df_max_UF = sweep_uf_95 %>% group_by(sweep) %>%  filter(xpclr==max(xpclr)) %>% 
    mutate(max_xplcr=xpclr,max_posS = pos_start,max_posE = pos_stop) %>% select(max_xplcr,max_posS,max_posE )
  Ssweep_UF = Ssweep_UF %>% left_join(df_max_UF, by = "sweep") %>% mutate(size = end-start) %>% 
    rownames_to_column(var = "sweep_ID") %>% select(-sweep)
  Ssweep_UF = Ssweep_UF %>% rename(chr_n=chrom) %>%  rowwise() %>% 
    mutate(sel_coef_m = xpclr_UF %>% filter(chrom==chr_n, pos_start>=start,pos_stop<=end) %>% 
             summarise(sel_coef = mean(sel_coef, na.rm = T)) %>% pull(sel_coef)) 
  ##proportion of genome region under selection
  print(paste("Proportion of genome under selection is",sum(Ssweep_UF$size)/genome_size))
  return(Ssweep_UF)
  }
Ssweep_UC = sweepID(xpclr_UC)
Ssweep_UF = sweepID(xpclr_UF)
mean(Ssweep_UC$sel_coef_m)
mean(Ssweep_UF$sel_coef_m)


##################################################################
##      Intersect/substract sweeps between two populations      ##
##################################################################

library("bedr")

##################################################################
##                   Identify selection sweep                   ##
##                          with RAiSD                          ##
##################################################################

#RAisd data loading
load_raisd = function(Raisd = file1, Chr_start = file2) {
  RAisd = read.table(Raisd, header = F,sep = "\t",comment.char = '/')
  UF_chr = read.table(Chr_start, header = F,sep = " ")
  UF_chr = UF_chr %>% mutate(r_num = as.numeric(str_remove(V1,"://")), row_num = row_number())%>% 
    mutate(start = r_num - row_num + 1, end = lead(start)-1 )
  UF_chr[28,6] = nrow(RAisd)
  UF_chr = UF_chr %>% mutate(len=end-start+1)
  chr = c()
  for (i in 1:28) {
    l = UF_chr$len[i]
    c = UF_chr$V2[i]
    chr = c(chr,rep(c,l))
  }
  RAisd[,"chr"]= chr
  chrs= data.frame(chr = paste0("chr_",rep(1:7,each =4),rep(c("A","B","C",'D'), 7)),id = 1:28)
  RAisd = RAisd %>% left_join(chrs, by=c("chr")) 
  return(RAisd)
}
raisd = load_raisd("RAiSD_Report.UF","UF.chr_start")
raisd_UF_q99 = quantile(raisd$V7,na.rm = T,probs = 0.99)
raisd_UF_q95 = quantile(raisd$V7,na.rm = T,probs = 0.95)
max(raisd$V7,na.rm = T)
head(raisd)

manhattan(raisd, chr="id", bp="V2", snp="V2", p="V7", logp=F, ylim = c(0,1415),cex=0.75,
          genomewideline = raisd_UF_q99, suggestiveline = raisd_UF_q95, chrlabs =chrs$chr , ylab = "RAisd (UF)",cex.axis=0.75)

raisd_uc = load_raisd("RAiSD_Report.UF","UF.chr_start")