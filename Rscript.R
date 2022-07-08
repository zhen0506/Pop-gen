
##################################################################
##                   kmer based IBD detection                   ##
##                           Zhen Fan                           ##
##                            6-13-22                           ##
##################################################################

library(npreg)
library(tidyverse)
bed = read.table(file = "coverage/FL_12_105_54_filter.bedgraph",header = F)

#signal detection function using spine smooth function 
bimodal_detect <- function(y,wind = 40,lap = 20,thres_adjust = 0) {
  ##spine smoothing 
  bed_sm = ss(x=1:length(y),y=y, m = 1)
  smooth_y = bed_sm$y
  den = density(smooth_y)
  q_90 = quantile(smooth_y,probs = 0.90)
  q_10 = quantile(smooth_y,probs = 0.1)
  ##detect junction of two distributions
  den_min = data.frame(x = den$x, y =den$y) %>% filter( x < q_90,x> q_10) %>% 
    filter(y==min(y))
  nibd_sd = sd(smooth_y[smooth_y<den_min$x])
  nibd_m = mean(smooth_y[smooth_y<den_min$x])
  nsd = floor((den_min$x-nibd_m)/nibd_sd)
  if (nsd>=4) {
    thre =nibd_m + 3*nibd_sd
    
  }else {
    thre = nibd_m + (nsd+thres_adjust)*nibd_sd
    }
  ##peak detection
  df <- data.frame("signal" = rep(0,ceiling(length(smooth_y)/(wind-lap))),
                   "window" = 0,"wstart"=0,"wend"=0,"mean_y" =0)
  for (i in 1:ceiling(length(smooth_y)/(wind-lap))) {
    wstart=(i-1)*wind-(i-1)*lap + 1
    if (i<ceiling(length(smooth_y)/(wind-lap))-1) {
      wend=i*wind-(i-1)*lap
    }else {
      wend = length(smooth_y)
      }
    yw = median(smooth_y[wstart:wend])
    sig_y = ifelse(yw>=thre,1,0)
    df[i,] = c(sig_y,i,wstart,wend,yw)
  }
  return(df)
  }

df1 = bimodal_detect(bed$V4,thres_adjust = 0)

#plot coverage, smooth function and peak detection results
overlay_plot <- function(coverage, df) {
  pos = 1:length(coverage)
  df_s = df %>% filter(signal == 1) %>% mutate(signal=signal-1)
  data.frame("cov" = coverage, "pos" = pos) %>%  
    ggplot(aes(x=pos)) +
    geom_point(aes(y=coverage,color="grey"),alpha = 0.5,show.legend = T) +
    geom_point(aes(y=mean_y,x=wstart,color="red") , data = df,show.legend = T) +
    geom_point(aes(x=wstart,y=signal,color="blue"), data = df_s,show.legend = T) +
    scale_y_continuous(name = "Coverage") +
    scale_colour_manual(name = '', 
                        values =c('grey'='grey','red'='red','blue'='blue'), labels = c('Coverage','Smooth_C',"Signal")) +
    labs(x="Marker number") +
    theme_bw(base_size =15) +
    theme(legend.position = c(0.1,0.85))
}
overlay_plot(bed$V4,df1)

##plot IBD each generation in one plot
path = "./coverage/"
files = list.files(path = path,pattern="*.bedgraph")
genotype = str_extract(files,pattern = "\\w+(?=_filter)")
PI612498_BC = c("FL_12_105_54", "FL_15_80_74", "FL_18_48_46")
PI612498_l = c("BC1","BC2","BC3")
PI551736_BC = c("FL_11_124_34", "FL_12_107_28","FL_18_46_54")
PI551736_l= c("F1","BC2","BC3")
wsize=5000
list_plot = PI551736_BC
overlayP=T
legend_names=PI551736_l


bed_input_main <- function(path, files,genotype, list_plot = genotype, legend_names=genotype, wsize=5000, overlayP=T,wind = 40,lap = 20,thres_adjust = 0) {
  myfiles = lapply(paste0(path,files), function(x) read.table(x,header = F) %>% filter(str_detect(V1,"chr")))
  names(myfiles) = genotype
  chr_l = myfiles %>% bind_rows() %>% group_by(V1) %>% summarise(max_end=max(V3)) %>% mutate(end = ceiling(max_end/5000))
  start=c();end=c();chr=c()
  chr_end = chr_l %>% pull(max_end)
  chr_n = chr_l %>% pull(V1)
  chr_b = chr_l %>% pull(end)
  for (i in 1:length(chr_end)) {
    chr = c(chr,rep(chr_n[i],chr_b[i]))
    start=c(start,c(0:(chr_b[i]-1))*wsize)
    end=c(end,c((1:(chr_b[i]-1))*wsize,chr_end[i]))
  } 
  df.st = data.frame(chr,start,end)
  myfiles = lapply(myfiles, function(x) x = df.st %>% left_join(x,c("chr" = "V1","start"= "V2","end"="V3")) %>% 
                     mutate(V4=replace_na(V4,0)))
  names(myfiles) = genotype
  IBD_list = lapply(myfiles, function(x)  bimodal_detect(x$V4, wind = wind,lap = lap,thres_adjust = thres_adjust))
  names(IBD_list) = genotype
  if (overlayP==T) {
    for (i in list_plot) {
      p1 = overlay_plot(myfiles[[i]]$V4,IBD_list[[i]])
      ggsave(filename = paste0(i,"_overlayplot.png"),units="in",width=8,height=4)
    }
  }  
  ##build a plot with base layers
  IBD_df = bind_rows(IBD_list,.id = "Geno")
  
  #IBD_df = IBD_df %>% mutate(signal=ifelse(Geno=="FL_18_46_54" & window >= 2391 & window <= 2404 , 0, signal)) %>% 
  #  mutate(signal=ifelse(Geno=="FL_18_46_54" & window >= 5718 & window <= 5887 , 0, signal)) %>% 
  #  mutate(signal=ifelse(Geno=="FL_18_46_54" & window <= 526 , 0, signal))
  
  coef_df = data.frame("Geno"=list_plot,"coef"=1:length(list_plot))
  IBD_df %>% filter(Geno %in% list_plot) %>% left_join(coef_df,by = "Geno") %>% mutate(signal = signal * coef)%>% 
    ggplot(aes(x=wstart,y=signal)) +
    geom_point(aes(color = Geno)) + 
    scale_x_continuous(labels = str_remove(chr_n,"chr_"),breaks =cumsum(chr_b) )+
    scale_color_discrete(labels=legend_names, name="chiloensis") +
    coord_cartesian(ylim=c(0.7,length(list_plot)+0.2),xlim = c(0,max(IBD_df$wend)))+
    theme_bw(base_size = 15) +
    labs(y="IBD",x="Chromosome")+
    theme(axis.text.x = element_text(hjust = 1),panel.grid.major.y  = element_blank(),panel.grid.minor.y = element_blank() ,panel.grid.minor.x = element_blank(),
          axis.text.y = element_blank())


  #summary stats table, needs to be returned
  IBD_df = IBD_df %>% left_join(df.st %>% rownames_to_column(var = "id") %>% mutate(id=as.numeric(id)), by=c("wstart"="id")) %>% select(-c(start,end))
  IBD_sum = IBD_df %>%  group_by(Geno,chr) %>% summarise(IBD_p = sum(signal)/n()) %>% 
    pivot_wider(names_from = Geno,values_from = IBD_p ) %>% mutate(chr_length = chr_end)
  write.csv(IBD_sum,file = "IBD_proportion_chr.csv",row.names = F)
  return(IBD_df)
}


B3_over = IBD_df %>%select(-mean_y) %>%  pivot_wider(names_from = Geno,values_from = signal) %>% filter(FL_18_46_54==1 & FL_18_48_46==1)
B3_over %>% group_by(chr) %>% summarise(st=min(wstart),end=max(wend))

#compare retained IBD regions and chi pro
chi = read.csv(file = "chiloensis_pro.csv",header = T)
library(ggpmisc)
IBD_sum %>%  mutate(chr = str_remove(chr,"chr_")) %>% left_join(chi %>% filter(Group=="UF") %>% group_by(chr) %>% summarise(chi_p = median(chi_p)),by = "chr") %>% 
  ggplot(aes(x=chi_p,y=FL_18_46_54)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() + 
  labs(x="Chiloensis Proportion", y="IBD from Chiloensis") +
  theme_bw(base_size = 15)

#yield change, line plot
library(readxl)
path = "./Field_raw/"
files = list.files(path = path,pattern="*.xlsx")
myfiles = lapply(paste0(path,files), function(x) read_xlsx(path = x,col_names = T))
phen_r = bind_rows(myfiles,.id = "Week")
df_name = read_xlsx(path = "Pedigree.xlsx", sheet = "Name")
phen_df = phen_r %>% inner_join(df_name, by = c("Genotype"))
#plot yield
pal = palette("Dark2")
phen_df %>% group_by(ID,Group) %>% summarise("Early_Yield" = sum(MktWt)/4) %>% ungroup() %>% 
  mutate(Group = factor(Group, levels = c("W","F1","BC1", "BC2","BC3","Radiance","Brilliance"))) %>% 
  ggplot(aes(x=Group,y=Early_Yield)) +
  geom_col(aes(fill = Group)) +
  scale_fill_manual(values = c("black","black",pal))+
  labs(x="Generation", y="Early Yield (g/plant)")+
  theme_bw(base_size = 15) +
  theme(legend.position = "")
#plot fruit size
phen_df %>% mutate(size=MktWt/MktNo) %>% 
  mutate(Group = factor(Group, levels = c("W","F1","BC1", "BC2","BC3","Radiance","Brilliance"))) %>% 
  ggplot(aes(x=Group,y=size)) +
  geom_boxplot(aes(fill = Group)) +
  scale_fill_brewer(palette = "Dark2")+
  labs(x="Generation", y="Fruit size (g/fruit)")+
  theme_bw(base_size = 15) +
  theme(legend.position = "")
#Brix
phen_df %>% filter(Brix>0) %>% 
  mutate(Group = factor(Group, levels = c("BC1", "BC2","BC3","Radiance","Brilliance"))) %>% 
  ggplot(aes(x=Group,y=Brix)) +
  geom_boxplot(aes(fill = Group)) +
  scale_fill_brewer(palette = "Dark2")+
  labs(x="Generation", y="Brix(%)")+
  theme_bw(base_size = 15) +
  theme(legend.position = "")
