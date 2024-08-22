library(tidyverse)
library(see)
library(magrittr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(tidyverse)
library(patchwork)
library(MetBrewer)

dat=read.csv("data/expert_survey/data_SCC-expert-survey_final_anonymous.csv")

fig3dat=dat%>%
  filter(!is.na(SCCwedge_Drivers_earth_system))%>% #only show subset of respondents that provide breakdown for SCC wedge
  select(contains(c("SCC_Lit","SCC_True_","Interview.number")))%>%
  pivot_longer(!"Interview.number..ongoing.")%>%
  mutate(type=substr(name,5,7))%>%
  mutate(quantile=as.factor(substr(name,9,11)))%>%
  mutate(quantile=fct_collapse(quantile,lower=c("2p5","_2p"),upper=c("7p5","_97"),central=c("Cen","_Ce")))%>%
  filter(!is.na(value))%>%
  pivot_wider(id_cols=c("Interview.number..ongoing.",type),values_from=value,names_from=quantile)
colnames(fig3dat)[1]="id"

#create jittered x variable
fig3dat$xj=jitter(as.numeric(as.factor(fig3dat$type)),amount=0.1)

b=ggplot(fig3dat)+
  geom_violinhalf(data=fig3dat%>%filter(type=="Tru"),aes(x=type,y=central),position=position_nudge(0.17),fill="#FF495C")+
  geom_violinhalf(data=fig3dat%>%filter(type=="Lit"),aes(x=type,y=central),position=position_nudge(-0.17),fill="#E5E059",flip=TRUE)+
  theme_classic()+scale_x_discrete(labels=c("Lit"="Literature","Tru"="Comprehensive"))+theme(text=element_text(size=12))+
  labs(x="",y="2020 SCC ($ per ton CO2)")+scale_fill_manual(values=c("#E5E059","#FF495C"),guide="none")+
  geom_segment(data=pivot_wider(fig3dat,id_cols=id,names_from=type,values_from = c(central,xj)),aes(x=xj_Lit,xend=xj_Tru,y=central_Lit,yend=central_Tru),col="grey50",lwd=0.4)+
  scale_y_continuous(limits=c(-10,450))+expand_limits(x=c(0.25,2.7))+
  geom_point(aes(x=xj,y=central,fill=type),size=2,pch=21,col="black")+
  geom_point(data=fig3dat%>%group_by(type)%>%dplyr::summarise(mean=mean(central)),aes(x=type,y=mean,shape=type),fill="black",size=3,pch=c(15,17))

#Figure 3b - barplot of the SCC wedge
data_barplots <- read_csv2("data/expert_survey/data_SCCwedge.csv")
data_barplots%<>% filter(complete.cases(.))
data_barplots=data_barplots[,-c(19,20)]
data_barplots %<>% filter(rowSums(is.na(.)) != ncol(.))

data_barplots_3b =  data_barplots %>% 
  filter(SCC_wedge_value !=0)%>%
  mutate(across(contains("Driver"), ~(. * abs(SCC_wedge_value) /100)))

# Weight by percent change to SCC estimate
data_barplots_3b_temp <- data_barplots %>%
  filter(SCC_wedge_value == 0) %>%
  mutate(across(contains("Driver"), ~(. * SCC_Lit_Central/100)))

data_barplots_3b %<>% bind_rows(data_barplots_3b_temp)
rm(data_barplots_3b_temp)

# Name Variables
names(data_barplots_3b) <- c("SCC Wedge Value",
                             "Earth System",
                             "Climate\nTipping Points",
                             "Damages\nTipping Points",
                             "Limited\nSubstitutability",
                             "Growth Damages",
                             "Distributional\nWeights",
                             "Epstein-Zin",
                             "Model Uncertainty",
                             "Learning",
                             "Alternative Ethics",
                             "Adaptation",
                             "Tech Progress",
                             "PRTP",
                             "EMUC",
                             "Damage Function",
                             "Other",
                             "Literature SCC")


data_barplots_3b$Subject <- 1:nrow(data_barplots_3b)
data_barplots_3b$"True SCC" <- data_barplots_3b$`SCC Wedge Value` + data_barplots_3b$"Literature SCC"

data_mean_pos <- data_barplots_3b  %>%
  colMeans() %>%
  round(0) %>% #change to 1 for version b
  data.frame() %>%
  t() %>%
  as.data.frame()

data_mean_all_long <- pivot_longer(data_mean_pos, cols = -Subject, names_to = "name", values_to = "value")

# Make grouping variable
data_mean_all_long %<>% 
  mutate(group = case_when(name == "Earth System" | name == "Climate\nTipping Points" ~ "Earth System Structural Changes",
                           name == "Damages\nTipping Points" | name == "Growth Damages" ~ "Climate Damage Structural Changes",
                           name == "Limited\nSubstitutability" | name == "Distributional\nWeights" | name == "Epstein-Zin" | name == "Model Uncertainty" ~ "Utility Function Structural Changes",
                           name == "Alternative Ethics" | name == "Learning" | name == "Endogenous Adaptation" ~ "Other Structural Changes",
                           name == "Pure Time Preference" | name == "Elasticity of Marginal Utility" ~ "Discounting Parametric Changes",
                           name == "Damage Function Parameters" ~ "Climate Damages Parametric Changes",
                           name == "Alternative Ethics" | name == "Learning" | name == "Adaptation" | name == "Tech Progress" ~ "Other Structural Changes",
                           name == "PRTP" | name == "EMUC" ~ "Discounting Parametric Changes",
                           name == "Damage Function" ~ "Climate Damages Parametric Changes",
                           name == "Other" ~ "Other",
                           name == "Literature SCC" ~ "Literature SCC"))

# make factors
data_mean_all_long$Subject %<>% factor()
data_mean_all_long$name %<>% factor()
data_mean_all_long$name %<>% factor(levels = c(names(data_barplots_3b)[c(1:4, 6, 5, 7:18)], "True SCC"))
data_mean_all_long$group %<>% factor()
data_mean_all_long$group %<>% factor(levels = levels(data_mean_all_long$group)[c(4, 1, 8, 7, 3, 2, 6, 5)])


data_mean_all_long_pos <- data_mean_all_long %>%
  filter(!is.na(group)) %>%
  group_by(Subject) %>%
  filter(value > 0 | name == "Literature SCC") %>%
  arrange(name) %>%
  arrange(Subject) %>%
  mutate(position = cumsum(rev(value))) %>%
  mutate(position = rev(position)) %>%
  mutate(position = position - value) %>%
  mutate(position = position + 0.5*value) %>%
  select(c("Subject", "name", "position"))


data_mean_all_long_neg <- data_mean_all_long %>%
  filter(!is.na(group)) %>%
  group_by(Subject) %>%
  filter(value < 0 | name == "Literature SCC") %>%
  arrange(name) %>%
  arrange(Subject) %>%
  mutate(position = cumsum(rev(value))) %>%
  mutate(position = rev(position)) %>%
  mutate(position = position - 0.5*value) %>%
  select(c("Subject", "name", "position"))

data_mean_all_long %<>% 
  left_join(data_mean_all_long_pos, by = c("Subject", "name")) %>%
  left_join(data_mean_all_long_neg, by = c("Subject", "name"))

data_mean_all_long %<>%
  mutate(position = coalesce(position.x, position.y)) 

# Set color palette
cols <- c(brewer.pal(5, "BrBG")[c(1,2)],
          brewer.pal(6, "YlOrRd")[c(2,3)],
          brewer.pal(6, "Greens")[c(3:6)],
          brewer.pal(6, "Blues")[c(3:6)],
          brewer.pal(5, "Reds")[c(4,5)],
          brewer.pal(5, "Purples")[5],
          "lightgray",
          "black")

barplot_mean_all <- data_mean_all_long %>%
  filter(!is.na(group)) %>%
  filter(name != "Literature SCC") %>%
  ggplot(aes(x = Subject)) +
  theme_classic() +
  labs(y =  element_blank(), x = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10),axis.line.x = element_blank()) +
  geom_tile(aes(y = position, height = value, fill = name), width = 0.4) +
  geom_segment(data = filter(data_mean_all_long, name == "Literature SCC"), 
                aes(y = value, yend = value, x =-Inf, xend = 1.17), linetype = "dashed", color = "gray") +
  geom_segment(data = filter(data_mean_all_long, name == "True SCC"), 
               aes(y = value, yend = value, x =-Inf, xend = 1.17), linetype = "dashed", color = "gray") +
  geom_point(data = filter(data_mean_all_long, name == "Literature SCC"), aes(y = value, x = -Inf), shape = 15, fill = "black", size = 4) +
  geom_point(data = filter(data_mean_all_long, name == "True SCC"), aes(y = value, x = -Inf), shape = 17,  fill = "black", size = 4) +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  # geom_segment(data = filter(Lit_True_df_all, Subject == 24), aes(x = Subject, xend = Subject, y = value.x, yend = value.y)) +
  geom_text(aes(x = 1.3, y = position, label = ifelse(value<0,paste0(name,"(-",substr(value,2,2),")"),paste0(name,"(+",value,")"))),hjust = 0)+ #force_pull = 1, force = 0.0003) +
  #geom_label_repel(aes(x = Subject, y = position, label = value), label.size = 0, force_pull = 1, force = 0.0003) +
  #geom_text_repel(data = filter(data_mean_all_long, name == "Literature SCC"), aes(y = value, x = 1.17, label = "Central literature SCC estimate"), 
                  #fontface = "italic", nudge_y = 5,  nudge_x = .13, hjust = 0) + 
  #geom_text_repel(data = filter(data_mean_all_long, name == "True SCC"), aes(y = value, x = 1.17, label = "Central \"true\" SCC estimate"), 
                  #fontface = "italic", nudge_x = .13, nudge_y = -5, hjust = 0) +
  #scale_y_continuous(breaks = c(50, filter(data_mean_all_long, name == "Literature SCC")$value, 75, 100, 125, 150, filter(data_mean_all_long, name == "True SCC")$value),
  #                   minor_breaks = seq(from = 50, to = 165, by = 12.5)) +
  # annotate("text", x = .465, y = filter(data_mean_all_long, name == "Literature SCC")$value, 
  #          label = filter(data_mean_all_long, name == "Literature SCC")$value, hjust = 1, size = 3.5, fontface= "italic", color = "#4d4d4d") +
  # annotate("text", x = .465, y = filter(data_mean_all_long, name == "True SCC")$value, 
  #          label = filter(data_mean_all_long, name == "True SCC")$value, hjust = 1, size = 3.5, fontface= "italic", color = "#4d4d4d") +
  coord_cartesian(xlim = c(1.1, 2.1), clip = "off")

#----------Add expert assessment of discount rate and model structure compared to literature ----------
#show expert assessment of discount rate and model structure compared to literature

#discount rate
#Expert discount rates from Drupp and Hansel 2018 AEJ paper
discountsurvey=read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")

#get literature discount rates
source("src/data_cleaining_scripts/cleaning_master.R")

fig3cdat=data.frame(dr=discountsurvey$SDR,type="Expert")
fig3cdat=rbind(fig3cdat,data.frame(dr=dat$discountrate,type="Literature"))

fig3c=ggplot(fig3cdat)+geom_density(aes(x=dr,group=type,lty=type,col=type),lwd=0.75)+
  theme_classic()+labs(x="Discount Rate",y="")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,12.5))+scale_color_manual(values=c("deepskyblue4","firebrick2"))+
  annotate("text",x=c(6.8,5.5),y=c(0.18,0.29),col=c("firebrick2","deepskyblue4"),label=c("Literature","Expert Assessment"))

#this is output from bayesian meta-analysis over joint probability of structural change inclusion (in bayes.R)
bayespost=read.csv("data/expert_survey/meta-analysis-distribution.csv")

#get fraction of dataset containing structural changes
dat$'Earth System'=ifelse(dat$`Climate Model`==1|dat$`Carbon Cycle`==1,1,NA)

strucfrac=dat%>%
  #dplyr::select("Tipping Points":"Alternative ethical approaches (not Discounted Utilitarianism)","Earth System")%>%
  dplyr::summarise(across(c("Tipping Points":"Learning","Earth System"), ~ sum(.x==1,na.rm=T)/n()))%>%
  pivot_longer(everything())

strucfrac$name[c(1:2,6)]=c("Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability")
colnames(strucfrac)[1]="question"

bayespost$question=factor(bayespost$question,levels=c("Earth System","Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Inequality Aversion","Epstein-Zin","Learning","Ambiguity/Model Uncertainty"))

fig3d=ggplot(bayespost)+geom_density(aes(x=prob,group=question),col="deepskyblue4")+facet_wrap(~as.factor(question),ncol=3)+
  theme_bw()+geom_vline(data=strucfrac,aes(xintercept=value),lty=2,lwd=1,color="firebrick2")+labs(x="Probability of Inclusion",y="")+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.background=element_rect(fill="white"))+xlim(-0.02,1)


# #add variable importance plot
# vip=read.csv("outputs/randomforest_plots/variable_importance_plot_Jan2024.csv")
# vip$column=fct_reorder(vip$column,vip$import)
# 
# fig4c=ggplot(vip,aes(x=0,xend=import,y=column,yend=column,col=Type))+geom_segment(linewidth=4)+
#   scale_color_manual(values=met.brewer(name="Hokusai1", n=4, type="discrete"),labels=c("Damage\nFunction","Other","Parametric\nUncertainty","Model Structure"))+
#   theme_bw()+labs(x="Varaible Importance",y="",color="")

layout="
##B
AAB
AAB
##B
CCB
CCB
CCD
CCD
CCD
"
patchwork=b+barplot_mean_all+fig3d+fig3c+plot_layout(design=layout)
ggsave("figures/Science Revision/figure3.pdf",plot=patchwork,width=10,height=8,units="in")
