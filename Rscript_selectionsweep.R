library(tidyverse)
library(qqman)
library(karyoploteR)
##################################################################
##                   Load xpclr files                          ##
##################################################################
## Load file function
load_xpclr = function(path1 = "./", pattern1 = "*", is.strawberry=T ) {
  temp = list.files(pattern=pattern1,path = path1)
  myfiles = lapply(paste0(path1,temp), function(x) read.table(x, header = T,sep = "\t") )
  xpclr = bind_rows(myfiles, .id = "a")
  if (is.strawberry ) {
    chrs= data.frame(chrom = paste0("chr_",rep(1:7,each =4),rep(c("A","B","C",'D'), 7)),chr = 1:28)
    xpclr = xpclr %>% left_join(chrs, by="chrom")
  }
  return(xpclr)
}
UF_xp = load_xpclr(path1 = "./", pattern1 = "UF_ecoH" )
UC_xp = load_xpclr(path1 = "./", pattern1 = "UC_ecoH" )
heir_xp = load_xpclr(path1 = "./", pattern1 = "heir_ecoH" )

## Manhattan plot function
manhattan_plot = function(xpclr, ylab1= "XP-CLR"){
  max_l = max(xpclr$xpclr,na.rm = T) + 20
  UF_q95 = quantile(xpclr$xpclr,na.rm = T,probs = 0.95)
  UF_q99 = quantile(xpclr$xpclr,na.rm = T,probs = 0.99)
  chr_lablel = unique(str_remove(xpclr$chrom, "chr_"))
  p = manhattan(xpclr , chr="chr", bp="start", snp="id", p="xpclr", logp=F, ylim = c(0,max_l),cex=0.75,
                genomewideline = UF_q99, suggestiveline = UF_q95, chrlabs =chr_lablel, ylab = ylab1,cex.axis=0.75)
  return(p)
}

manhattan_plot(heir_xp)
manhattan_plot(UC_xp)
manhattan_plot(UF_xp)

##################################################################
##                   Identify selective sweep                   ##
##                    using a merge distance                    ##
##################################################################
##function to get genomic intervals under selection 
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
  Ssweep_UF = Ssweep_UF %>% dplyr::rename(chr_n="chrom") %>%  rowwise() %>% 
    mutate(sel_coef_m = xpclr_UF %>% filter(chrom==chr_n, pos_start>=start,pos_stop<=end) %>% 
             summarise(sel_coef = mean(sel_coef, na.rm = T)) %>% pull(sel_coef)) 
  ##proportion of genome region under selection
  print(paste("Proportion of genome under selection is",sum(Ssweep_UF$size)/genome_size))
  return(Ssweep_UF)
}
Ssweep_UC = sweepID(UC_xp, prob = 0.95, limit = 400000)
Ssweep_UF = sweepID(UF_xp, prob = 0.95, limit = 400000)
Ssweep_heir = sweepID(heir_xp, prob = 0.95, limit = 400000)
mean(Ssweep_UC$sel_coef_m)
mean(Ssweep_heir$sel_coef_m)
mean(Ssweep_UF$sel_coef_m)

#write_csv(Ssweep_heir,file="Ssweep_heir.csv")
#write_csv(Ssweep_UF,file="Ssweep_UF.csv")
#write_csv(Ssweep_UC,file="Ssweep_UC.csv")

##################################################################
##                   Identify selective sweep                   ##
##                    using smooth function                     ##
##################################################################

library(npreg)
library(pracma)
##main function to identify selective sweeps using smooth and findpeaks functions
##load the df output from load_xpclr function
sweepfinder = function(UF_UC_xp = df, q = q)
  { UF_UC_xp = UF_UC_xp %>% rownames_to_column(var = "order")
  UF_UC_xp = UF_UC_xp  %>% group_by(chrom) %>% 
    mutate(xpclr_sm = ss(x=order,y=xpclr %>% replace_na(0),m=1,lambda = 0.1, spar = 0.1)$y) %>% 
    ungroup()  
  q95 = quantile(UF_UC_xp$xpclr_sm,q)
  test = findpeaks(UF_UC_xp$xpclr_sm, threshold = q95 )
  xp_df = UF_UC_xp %>% select(chrom,start,stop,order) %>% 
    mutate(order = as.numeric(order))
  test = as.data.frame(test) %>% select(V3,V4) %>% 
    rename("V3"="LB", "V4"="RB") %>% 
    left_join(xp_df, by =c("LB" = "order")) %>% 
    left_join(xp_df, by =c("RB" = "order")) %>% 
    filter(chrom.x == chrom.y) %>% 
    select(chrom.x, start.x,stop.y) 
  print(paste0("With threshold of ",q,", ",
               nrow(test)," peaks were selected, totaling ",
               sum(test$stop.y-test$start.x), 'Mb.'))
  return(list(UF_UC_xp,test))
  }


## data processing for real data 
UF_UC_xp1 = load_xpclr(path1 = "./New_select_083123/", pattern1 = "UF_UC" )
vir_chi_xp1 = load_xpclr(path1 = "./New_select_083123/", pattern1 = "vir_chi" )
pac_chi_xp1 = load_xpclr(path1 = "./New_select_083123/", pattern1 = "pac_chi" )

UF_s95 = sweepfinder(UF_UC_xp = UF_UC_xp1, q = 0.95)
UF_s99 = sweepfinder(UF_UC_xp = UF_UC_xp1, q = 0.99)
vir_s95 = sweepfinder(UF_UC_xp = vir_chi_xp1, q = 0.95)
vir_s99 = sweepfinder(UF_UC_xp = vir_chi_xp1, q = 0.99)
pac_s95 = sweepfinder(UF_UC_xp = pac_chi_xp1, q = 0.95)
pac_s99 = sweepfinder(UF_UC_xp = pac_chi_xp1, q = 0.99)
##for both UF and vir, 0.95 was used, for pac 0.99 was used
##convergent peak for UF and between two species
human_nat.int = bt.merge(bt.intersect(UF_s95[[2]], vir_s99[[2]]),)
human_chi.int = bt.merge(bt.intersect(UF_s95[[2]], pac_s99[[2]]))
div.swp = bind_rows(human_nat.int,human_chi.int,.id = "class")
div.swp = div.swp %>% left_join(UF_s95[[1]] 
                                %>% select(chrom, start, stop,order), by =c("V1"="chrom", "V2"="start")) %>% 
  left_join(UF_s95[[1]] 
            %>% select(chrom, start, stop,order), by =c("V1"="chrom", "V3"="stop"))


sum(div.swp$V3-div.swp$V2)

#build plots for selection sweeps
xp_df = UF_s95[[1]] %>% mutate(order = as.numeric(order)) %>%  
  filter(order %% 10 ==1) %>% 
  left_join(suffix = c("UF","vir"),vir_s99[[1]], by = c("chrom","start")) %>% 
  left_join(suffix = c("","pac"),pac_s99[[1]], by = c("chrom","start")) 
label = xp_df %>% mutate(chrom = str_remove(chrom,"chr_")) %>% group_by(chrom) %>% summarise(start = min(start),order=max(orderUF))
q95 = quantile(UF_s95[[1]]$xpclr_sm,0.95) 
xp_uf = xp_df %>% select(orderUF,chrom,start,stop,xpclr_smUF,xpclr_smvir,xpclr_sm) %>% 
  pivot_longer(cols = c(xpclr_smUF,xpclr_smvir,xpclr_sm), names_to = "x") %>% 
  filter(x == "xpclr_smUF") 
div.swp %>% mutate(xmin=as.numeric(order.x),
                   xmax = as.numeric(order.y))
ggplot() +
  geom_line(data = xp_uf, aes(y=value,x=orderUF),alpha = 1, color = "grey")+
  geom_rect(data = div.swp %>% mutate(xmin=as.numeric(order.x)-500,
                                      xmax = as.numeric(order.y)+500) %>% 
              filter(class == 1), 
            aes(xmin=xmin,xmax = xmax,ymin = -Inf, ymax=Inf),
            fill = "green", alpha = 0.5) +
  geom_rect(data = div.swp %>% mutate(xmin=as.numeric(order.x)-500,
                                      xmax = as.numeric(order.y)+500) %>% 
              filter(class == 2), 
            aes(xmin=xmin,xmax = xmax,ymin = -Inf, ymax=Inf),
            fill = "blue", alpha = 0.5) +
  scale_x_continuous(breaks = label$order, labels = label$chrom ) +
  geom_hline(yintercept =q95 )+
  labs(x="", y = "xpclr_smooth")+
  theme_bw(base_size = 15)
q99 = quantile(vir_s99[[1]]$xpclr_sm,0.99)
xp_df %>% select(orderUF,chrom,start,stop,xpclr_smUF,xpclr_smvir,xpclr_sm) %>% 
  pivot_longer(cols = c(xpclr_smUF,xpclr_smvir,xpclr_sm), names_to = "x") %>% 
  filter(x == "xpclr_smvir") %>% 
  ggplot(aes(y=value,x=orderUF)) +
  geom_line(alpha = 0.5, color = "green") +
  scale_x_continuous(breaks = label$order, labels = label$chrom )+
  geom_hline(yintercept =q99 )+
  labs(x="", y = "xpclr_smooth")+
  theme_bw(base_size = 15)
plot(x)
q99 = quantile(pac_s99[[1]]$xpclr_sm,0.99)
xp_df %>% select(orderUF,chrom,start,stop,xpclr_smUF,xpclr_smvir,xpclr_sm) %>% 
  pivot_longer(cols = c(xpclr_smUF,xpclr_smvir,xpclr_sm), names_to = "x") %>% 
  filter(x == "xpclr_sm") %>% 
  ggplot(aes(y=value,x=orderUF)) +
  geom_line(alpha = 0.5, color = "blue") +
  scale_x_continuous(breaks = label$order, labels = label$chrom )+
  geom_hline(yintercept =q99 )+
  labs(x="", y = "xpclr_smooth")+
  theme_bw(base_size = 15)


##################################################################
##      Intersect/substract sweeps between two populations 
##      8-4-2022 update
##      Farm CPU results updated
##      Overlay AWT and TMY GWAS results
##################################################################

#devtools::install_github("PhanstielLab/bedtoolsr")
options(bedtools.path = "/apps/bedtools/2.30.0/bin")
library("bedtoolsr")

## read sweeps data
Ssweep_heir = read_csv(file="Ssweep_heir.csv")
Ssweep_UF = read_csv(file="Ssweep_UF.csv")
Ssweep_UC = read_csv(file="Ssweep_UC.csv")

## identify early selective sweeps
UF.bed = Ssweep_UF %>% select(chr_n,start,end,sel_coef_m)
UC.bed = Ssweep_UC %>% select(chr_n,start,end,sel_coef_m)
heir.bed = Ssweep_heir %>% select(chr_n,start,end,sel_coef_m)
UF.int = bt.intersect(UF.bed, UC.bed,wa = T)
UFUC.int = bt.intersect(UF.bed, UC.bed)
sum(UFUC.int$V2-UFUC.int$V3)/780000000
#UFUC.int = UFUC.int %>% mutate(size = V3-V2)
heir.int = bt.intersect(UFUC.int, heir.bed)
sum(heir.int$V2-heir.int$V3)/780000000

## parallel selection sweeps
para.bed = bt.intersect(a=UFUC.int, b=heir.bed,v = T)
sum(para.bed$V2-para.bed$V3)/780000000
## divergent selection 
UF.di = bt.intersect(a = UF.bed, b = UFUC.int,v = T)
UC.di = bt.intersect(a = UC.bed, b = UFUC.int,v = T)
sum(UF.di$V2-UF.di$V3)/780000000
sum(UC.di$V2-UC.di$V3)/780000000

###find selection coefficient for each sweep class
library(ggpubr)
library(rstatix)
UF.con = bt.intersect(UF.bed,heir.int,wa=T)
UC.con = bt.intersect(UC.bed,heir.int,wa=T)
heir.con = bt.intersect(heir.bed,heir.int,wa=T)
df.con = bind_rows(UF.con,UC.con,heir.con, .id="Pop") %>% 
  mutate(Pop = recode(Pop,"1"="UF", "2" = "UC", "3" = "Early")) 
stat.test = df.con %>% t_test(V4 ~ Pop)
stat.test <- stat.test %>% add_xy_position(x = "Pop")
p1 = df.con %>%
  ggplot(aes(x=Pop, y= V4)) +
  geom_boxplot(aes(fill=Pop)) + 
  geom_jitter() 
p1 + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(y="Selection coefficient",x="Sweep class") +
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",legend.background = element_rect(colour="black"),legend.spacing.x = unit(0.1,units = "cm"))+ guides(colour = guide_legend(override.aes = list(size=4)))

#################################################################
##             Identify GWAS peaks under selection             ##
#################################################################
##GWAS files
GWAS = read_csv(file = "yielddata/TMY.FarmCPU.TMY.GWAS.Results.csv")
GWAS_awt = read_csv(file = "yielddata/AWT.FarmCPU.predicted.value.GWAS.Results.csv")

#update few positions
ovlap_GWAS = function(GWAS=df,pvalue = 0.05,window = 400000){
  pos = readxl::read_xlsx(path = "yielddata/Map_Faar1.xlsx")
  chr_name = read_table(file = "Farr1_chr_naming.txt")
  names(pos)
  GWAS = GWAS %>% mutate(SNP = str_replace(SNP,"AX","AX-")) %>% left_join(pos[,c(1,2,5)],by = c("SNP"="SNP_ID"))
  GWAS = GWAS %>% left_join(chr_name[,c("Number","Farr1")], by = c("Chromosome" = "Number"))
  GWAS1 = GWAS %>% 
    filter(`FDR_Adjusted_P-values`<pvalue) %>% mutate(Faar1_positions = as.numeric(Faar1_positions)) %>% 
    mutate(start = Faar1_positions-window, end = Faar1_positions+window) %>% 
    mutate(start = if_else(start<0,0,start)) %>% 
    select(Farr1,start,end) %>% filter(!is.na(start))
  GWAS.int = bt.intersect(UFUC.int, GWAS1,wb=T) %>% distinct() 
  GWAS.UC = bt.intersect(bt.intersect(UC.bed, GWAS1,wb=T),GWAS.int,v=T)
  GWAS.UF = bt.intersect(bt.intersect(UF.bed, GWAS1,wb=T),GWAS.int,v=T)
  GWAS_highlight = bind_rows(GWAS.int,GWAS.UC,GWAS.UF,.id = "color")  %>% 
    mutate(color=recode(color,"1"="black","2"="red","3"="green"))  
  return(list(GWAS,GWAS_highlight))
  }

list_TMY = ovlap_GWAS(GWAS,pvalue=0.05)
list_awt = ovlap_GWAS(GWAS_awt,pvalue=0.01)

list_TMY[[2]]
list_awt[[2]]

##manhatton plots for GWAS, color based on selection sweeps categories
#build a custom genome 
chr = read.table(file = "farr1.chr.length", header = F, as.is = T)
chr = chr[1:28,] %>% mutate(start = rep(1,28)) %>% dplyr::rename(chr = V1, end=V2) %>% 
  select(chr, start,end) %>%   mutate(chr = str_remove(chr,"chr_"))
custom.genome <- toGRanges(chr)
#build SNP dataset 
thres = 0.05/5534
GWAS = list_TMY[[1]] %>% filter(!is.na(Faar1_positions)) %>%  mutate(Farr1 = str_remove(Farr1,"chr_"), Faar1_positions = as.numeric(Faar1_positions))
GWAS_bed = GWAS[,c("Farr1", "Faar1_positions","Faar1_positions","P.value")] 
names(GWAS_bed) = c("V1","V2","V3","P.value")

GWAS_bed1 = GWAS_bed %>% left_join(list_TMY[[2]] %>% select(V5,V6,V6,color) %>% mutate(V1=str_remove(V5,"chr_"),V2=V6+400000, V3=V6+400000) %>%
                                           select(V1,V2,V3,color), by=c("V1","V2","V3")) %>% mutate(color = replace_na(color,"grey")) %>% 
  filter(!is.na(V2))

GWAS.data <- toGRanges(as.data.frame(GWAS_bed1))
#names(GWAS.data) <- GWAS2$SNP
GWAS.data$pval = GWAS_bed1$P.value

#genome-wide GWAS view
kp <- plotKaryotype(plot.type=4, genome = custom.genome,cex=1.2)
#kpAddBaseNumbers(kp, add.units = TRUE, cex=1.2)
kp <- kpPlotManhattan(kp, data=GWAS.data, points.col = GWAS.data$color,ymax = 10,genomewideline = -log10(0.0001),
                      suggestiveline = 0, r0=0.55,r1 = 1)
kpAxis(kp, ymin=6, ymax=10,cex = 1.2, r0=0.55,r1 = 1,labels =c(0,5,10))

##add AWT GWAS

GWAS_awt1 = list_awt[[1]] %>% filter(!is.na(Faar1_positions)) %>%  mutate(Farr1 = str_remove(Farr1,"chr_"), Faar1_positions = as.numeric(Faar1_positions))
GWAS_awt1_bed = GWAS_awt1[,c("Farr1", "Faar1_positions","Faar1_positions","P.value")] 
names(GWAS_awt1_bed) = c("V1","V2","V3","P.value")
awt_highlight = list_awt[[2]] %>% arrange(V5, V6) %>% distinct(V6,.keep_all = T)
GWAS_bed_2 = GWAS_awt1_bed %>% left_join(awt_highlight %>% select(V5,V6,V6,color) %>% mutate(V1=str_remove(V5,"chr_"),V2=V6+400000, V3=V6+400000) %>%
                                     select(V1,V2,V3,color), by=c("V1","V2","V3")) %>% mutate(color = replace_na(color,"grey")) %>% 
  filter(!is.na(V2))

GWAS.data2 <- toGRanges(as.data.frame(GWAS_bed_2))
#names(GWAS.data) <- GWAS2$SNP
GWAS.data2$pval = GWAS_bed_2$P.value

#kpAddBaseNumbers(kp, add.units = TRUE, cex=1.2)
kp <- kpPlotManhattan(kp, data=GWAS.data2, points.col = GWAS.data2$color,ymax = 20,genomewideline = -log10(0.00005),
                      suggestiveline = 0, r0=0,r1 = 0.5)
kpAxis(kp, ymin=0, ymax=20,cex = 1.2, r0=0,r1 = 0.5)
write_csv(bind_rows(TMY = GWAS_bed1 %>% filter(P.value<0.0001), AWT = GWAS_bed_2 %>% filter(P.value<0.00005),.id = "trait"),file = "QTLunderselection.csv")

##simulation test for signicance of overlap between GWAS and selection sweeps
row.names(chr) = chr$chr
c=0
for (i in 1:1000) {
  chr_s = sample(chr$chr,1)
  chr_l = chr[chr_s,"end"]
  pos = sample(1:chr_l,1)
  bed_s = data.frame(V1=paste0("chr_",chr_s),V2=ifelse(pos-500000<0,1,pos-500000), V3=ifelse(pos+500000>chr_l,chr_l,pos+500000))
  c=c+nrow(bt.intersect(bed_s,UC.bed))
  
}
dbinom(4, size=11, prob=0.7) 
dbinom(10,size = 26,prob = 0.7)


##################################################################
##                   Overlap selection sweep                   
##                    with chromosomal ancestry    
##                          update 4/9/23
##################################################################

library(ggpmisc)
AD =read.csv(file = "chr_group_virPro.csv")
AD_w = AD %>% select(-X) %>% pivot_wider(names_from = Group,values_from = Fvv.m) %>% 
  mutate(dif = abs(UCD-UF))
##selection sweeps are the main source for ancestory diffience
di.sw = bind_rows(UF.di,UC.di)
di.sw %>% group_by(V1) %>% summarise(size = sum(V3-V2)) %>%
  left_join(AD_w,by=c("V1"="Chr")) %>%
  ggplot(aes(x=dif,y=size)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point()+ labs(x="Ancestry proportion difference",y = "Size of selection sweep") +
  theme_bw(base_size = 15)
write.csv(di.sw %>% group_by(V1) %>% summarise(size = sum(V3-V2)) %>%
            left_join(AD_w,by=c("V1"="Chr")), file = "SS_Anc_linearregression.csv",row.names = F)



