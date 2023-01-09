library(qiime2R)
library(tidyverse)
library(vegan)
#load file and correct zero values
getwd()
table_240<- read_qza("Data/merge_table_240_noplant_filtered_nous.qza")$data %>% 
  as.data.frame()

taxonomy_240<- read_qza("Data/taxonomy_blast_240_0.97.qza")$data

#remove non-fungi taxa

taxonomy_filter<- taxonomy_240 %>% filter(
  !str_detect(Taxon, "ozoa"))%>% filter(
    !str_detect(Taxon, "helida")) %>% filter(
      !str_detect(Taxon, "ophyta")) %>% filter(
        !str_detect(Taxon, "Ciliophora")) %>% filter(
          !str_detect(Taxon, "Nucleariidae_and_Fonticula_group")) %>% filter(
            !str_detect(Taxon, "Arthrioida")) %>% filter(
              !str_detect(Taxon, "Labyrinthulomycetes"))  %>% filter(
                !str_detect(Taxon, "Apicomplexa")) %>% filter(
                  !str_detect(Taxon, "Bicosoecida")) %>% filter(
                    !str_detect(Taxon, "Breviatea")) %>% filter(
                      !str_detect(Taxon, "Aphelidea"))  %>% filter(
                        !str_detect(Taxon, "Arthropoda")) 

table_filter<-  table_240[match(
  taxonomy_filter$Feature.ID, rownames(table_240)),] %>% drop_na(.)

#load metadata and taxonomy file
meta<-read_tsv("Data/FINALMAP18S") %>% 
  rename(SampleID=`#SampleID`) %>%
  filter(SampleID!="#q2:types")

meta$Compartments<- factor(meta$Type_of_soil,
                           levels = c( "Non-rizospheric", "Rizospheric", 
                                       "Roots", "Uncultivated"),
                           labels = c("Bulk soil", "Rhizosphere",
                                      "Roots", "Uncultivated"))
meta$Watering_regimes<- factor(meta$Treatment,
                               levels = c( "0", "1", "2", "3"),
                               labels = c("Initial","Wet", "Dry", "Extremely dry"))
library(tidyverse)
meta<- meta %>% unite("treats", c("Compartments", "Watering_regimes"), sep = "__", remove = F)

d.pro.0<- table_filter %>% dplyr::select_at(vars(!contains("US")))
d.pro <- t(cmultRepl(t(d.pro.0), method="CZM", output="p-counts"))
d.clr.abund <- CoDaSeq::codaSeq.clr(t(d.pro) )
d.pro.dis<- vegdist(d.clr.abund, method="euclidean")
#fung.dist<-vegdist(decostand(t(table_filter),"hellinger"),"bray")

meta_just<- data.frame(t(d.pro), check.names = F) %>%
    rownames_to_column( var = "SampleID") %>% inner_join(meta)

meta_just  %>% group_by(Compartments) %>% count()
meta_just  %>% group_by(Watering_regimes) %>% count()

library(cowplot)
library(ggpubr)
#treatments
bd<-betadisper(d.pro.dis,factor(meta_just$Watering_regimes),  type = "centroid")
df1<-data.frame(permutest(bd, permutations = 999)$tab)  %>% 
  rename("p-value"="Pr..F.", "Sum \n Squares"= "Sum.Sq","Means \n Squares"=  "Mean.Sq",
         "N. \n Perms"="N.Perm"  ) %>% 
  rownames_to_column( var="Vars") %>%
  mutate_at(3:4, funs(round(., 0)))  %>%
  mutate_at(5:6, as.numeric) %>% 
  mutate_at(5, funs(round(., 2))) %>% 
  mutate_at(5:7, ~replace(., is.na(.), "-"))%>% 
  ggtexttable(., rows = NULL,  theme = ttheme(
    colnames.style = colnames_style(color = "black",
                                    fill = "white",hjust=0, x=0.1, size = 10),
    tbody.style = tbody.style)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
   tab_add_hline(at.row = c(3), row.side = "bottom", 
                linewidth = 3, linetype = 1) 
p1 <- ~{
  par(
    mar = c(3, 3, 1, 1), 
    mgp = c(2, 1, 0)
  )
  boxplot(bd, col=c("#479330","#FFFF00", "#FF0000"), xlab="Watering regimes")}
#plot(bd,label.cex = 0.8 , seg.col=c("#479330","#FFFF00", "#FF0000"),
 #    col=c("#479330","#DAD231", "#FF0000"),
  #   main="", hull=F)
#Compartments
bd1<-betadisper(d.pro.dis,factor(meta_just$Compartments),  type = "centroid")
df2<-data.frame(permutest(bd1, permutations = 999)$tab)  %>% 
  rename("p-value"="Pr..F.", "Sum \n Squares"= "Sum.Sq","Means \n Squares"=  "Mean.Sq",
         "N. \n Perms"="N.Perm"  ) %>% 
  rownames_to_column( var="Vars") %>%
  mutate_at(3:4, funs(round(., 0)))  %>%
  mutate_at(5:6, as.numeric) %>% 
  mutate_at(5, funs(round(., 2))) %>% 
  mutate_at(5:7, ~replace(., is.na(.), "-"))%>% 
  ggtexttable(., rows = NULL,  theme = ttheme(
    colnames.style = colnames_style(color = "black",
                                    fill = "white",hjust=0, x=0.1, size = 10),
    tbody.style = tbody.style)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom", 
                linewidth = 3, linetype = 1) 
p2 <- ~{
  par(
    mar = c(3, 3, 1, 1), 
    mgp = c(2, 1, 0)
  )
  boxplot(bd1, col=c("#800000", "#808000", "#008000"), xlab="Compartments")
}

all<-plot_grid(df1, p1, df2, p2, labels = c("A)", "B)", "C)", "D)"), label_size = 20)
#treatments by compartment
#ro
d.pro1<- table_filter %>% dplyr::select_at(vars(contains("RO"))) %>%filter(rowSums(across(where(is.numeric)))!=0)
d.pro2 <- t(cmultRepl(t(d.pro1), method="CZM", output="p-counts"))
d.clr.abund1 <- CoDaSeq::codaSeq.clr(t(d.pro2) )
d.pro.dis1<- vegdist(d.clr.abund1, method="euclidean")
#fung.dist<-vegdist(decostand(t(table_filter),"hellinger"),"bray")

meta_just1<- data.frame(t(d.pro1), check.names = F) %>%
  rownames_to_column( var = "SampleID") %>% inner_join(meta)
bda<-betadisper(d.pro.dis1,factor(meta_just1$Watering_regimes),  type = "centroid")
dfa<-data.frame(permutest(bda, permutations = 999)$tab)  %>% 
  rename("p-value"="Pr..F.", "Sum \n Squares"= "Sum.Sq","Means \n Squares"=  "Mean.Sq",
         "N. \n Perms"="N.Perm"  ) %>% 
  rownames_to_column( var="Vars") %>%
  mutate_at(3:4, funs(round(., 0)))  %>%
  mutate_at(5:6, as.numeric) %>% 
  mutate_at(5, funs(round(., 2))) %>% 
  mutate_at(5:7, ~replace(., is.na(.), "-"))%>% 
  ggtexttable(., rows = NULL,  theme = ttheme(
    colnames.style = colnames_style(color = "black",
                                    fill = "white",hjust=0, x=0.1, size = 10),
    tbody.style = tbody.style)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom", 
                linewidth = 3, linetype = 1) 
pa <- ~{
  par(
    mar = c(3, 3, 1, 1), 
    mgp = c(2, 1, 0)
  )
  boxplot(bda, col=c("#479330","#FFFF00", "#FF0000"), xlab="Watering regimes")
}


#ri
d.pro1<- table_filter %>% dplyr::select_at(vars(contains("RI"))) %>%filter(rowSums(across(where(is.numeric)))!=0)
d.pro2 <- t(cmultRepl(t(d.pro1), method="CZM", output="p-counts"))
d.clr.abund1 <- CoDaSeq::codaSeq.clr(t(d.pro2) )
d.pro.dis1<- vegdist(d.clr.abund1, method="euclidean")
#fung.dist<-vegdist(decostand(t(table_filter),"hellinger"),"bray")

meta_just1<- data.frame(t(d.pro1), check.names = F) %>%
  rownames_to_column( var = "SampleID") %>% inner_join(meta)
bdb<-betadisper(d.pro.dis1,factor(meta_just1$Watering_regimes),  type = "centroid")
dfb<-data.frame(permutest(bdb, permutations = 999)$tab)  %>% 
  rename("p-value"="Pr..F.", "Sum \n Squares"= "Sum.Sq","Means \n Squares"=  "Mean.Sq",
         "N. \n Perms"="N.Perm"  ) %>% 
  rownames_to_column( var="Vars") %>%
  mutate_at(3:4, funs(round(., 0)))  %>%
  mutate_at(5:6, as.numeric) %>% 
  mutate_at(5, funs(round(., 2))) %>% 
  mutate_at(5:7, ~replace(., is.na(.), "-"))%>% 
  ggtexttable(., rows = NULL,  theme = ttheme(
    colnames.style = colnames_style(color = "black",
                                    fill = "white",hjust=0, x=0.1, size = 10),
    tbody.style = tbody.style)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom", 
                linewidth = 3, linetype = 1) 
pb <- ~{
  par(
    mar = c(3, 3, 1, 1), 
    mgp = c(2, 1, 0)
  )
  boxplot(bdb, col=c("#479330","#FFFF00", "#FF0000"), xlab="Watering regimes")
}
#nr
d.pro1<- table_filter %>% dplyr::select_at(vars(contains("NR"))) %>%filter(rowSums(across(where(is.numeric)))!=0)
d.pro2 <- t(cmultRepl(t(d.pro1), method="CZM", output="p-counts"))
d.clr.abund1 <- CoDaSeq::codaSeq.clr(t(d.pro2) )
d.pro.dis1<- vegdist(d.clr.abund1, method="euclidean")
#fung.dist<-vegdist(decostand(t(table_filter),"hellinger"),"bray")

meta_just1<- data.frame(t(d.pro1), check.names = F) %>%
  rownames_to_column( var = "SampleID") %>% inner_join(meta)
bdc<-betadisper(d.pro.dis1,factor(meta_just1$Watering_regimes),  type = "centroid")
dfc<-data.frame(permutest(bdc, permutations = 999)$tab)  %>% 
  rename("p-value"="Pr..F.", "Sum \n Squares"= "Sum.Sq","Means \n Squares"=  "Mean.Sq",
         "N. \n Perms"="N.Perm"  ) %>% 
  rownames_to_column( var="Vars") %>%
  mutate_at(3:4, funs(round(., 0)))  %>%
  mutate_at(5:6, as.numeric) %>% 
  mutate_at(5, funs(round(., 2))) %>% 
  mutate_at(5:7, ~replace(., is.na(.), "-"))%>% 
  ggtexttable(., rows = NULL,  theme = ttheme(
    colnames.style = colnames_style(color = "black",
                                    fill = "white",hjust=0, x=0.1, size = 10),
    tbody.style = tbody.style)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom", 
                linewidth = 3, linetype = 1) 

pc <- ~{
  par(
    mar = c(3, 3, 1, 1), 
    mgp = c(2, 1, 0)
  )
  boxplot(bdc, col=c("#479330","#FFFF00", "#FF0000"), xlab="Watering regimes")
}
all1<-plot_grid(dfa, pa, dfb, pb, dfc, pc, ncol = 2, nrow = 3,
          labels = c("A) ROOTS","",  "B) RHIZOSPHERE","", "C) BULK SOIL"))
ggsave('Figures/Fig_betadisper1.eps',
       width = 10, height = 6, dpi = 300, plot =all)
ggsave('Figures/Fig_betadisper2.eps',
       width = 10, height = 6, dpi = 300, plot =all1)
