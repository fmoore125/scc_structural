#Sep 2020

#read in searches from WebofScience, SCOPUS, and EconLit and discard duplicate DOIs

econlit=read.csv("data/abstract_search/EconLit_20200930.csv")
scopus=read.csv("data/abstract_search/Scopus_20200930.csv")
wos=read.csv("data/abstract_search/WebofScience_20200930.csv")

#scopus is by far the largest, so identify DOIs in other two not in scopus

econlit_diff=econlit[-which(trimws(econlit$digitalObjectIdentifier)%in%tolower(scopus$DOI)),]
#just one paper in econ lit not in scopus - "Federal Coal Program Reform, the Clean Power Plan, and the Interaction of Upstream and Downstream Climate Policies" - Gerarden et al.

econlit_diff=econlit[which(econlit$Missing.DOI.Duplicate.in.SCOPUS==0),] #110 papers in econlit missing DOIs not reproduced in SCOPUS

wos_diff=wos[-which(tolower(wos$DOI)%in%tolower(scopus$DOI)),]
wos_diff=rbind(wos_diff,wos[which(wos$MissingDOI_EconLitDup==0&wos$MissingDOI_ScopusDup==0),]) #add in unique papers with missing DOIs
if(length(which(tolower(wos_diff$Book.DOI)%in%tolower(scopus$DOI)&wos_diff$Book.DOI!=""))>0) wos_diff=wos_diff[-which((wos_diff$Book.DOI%in%tolower(scopus$DOI)&wos_diff$Book.DOI!="")),]
#28 papers in web of science - also including the one non-repeat paper from econlit

compiled=scopus
wos_diff_df=data.frame(Authors=wos_diff$Authors,Author.s..ID=NA,Title=wos_diff$Article.Title,Year=wos_diff$Publication.Year,Source.title=wos_diff$Source.Title,Volume=wos_diff$Volume,Issue=wos_diff$Issue,Art..No.=wos_diff$Article.Number,Page.start=wos_diff$Start.Page,Page.end=wos_diff$End.Page,Page.count=NA,Cited.by=wos_diff$Times.Cited..All.Databases,DOI=wos_diff$DOI,Link=NA,Document.Type=wos_diff$Publication.Type,Publication.Stage=NA,Access.Type=wos_diff$Open.Access.Designations,Source="WebofScience",EID=NA,X=NA)
econlit_diff_df=data.frame(Authors=econlit_diff$Authors,Author.s..ID=NA,Title=econlit_diff$Title,Year=econlit_diff$year,Source.title=econlit_diff$pubtitle,Volume=econlit_diff$volume,Issue=econlit_diff$issue,Art..No.=NA,Page.start=econlit_diff$startPage,Page.end=NA,Page.count=NA,Cited.by=NA,DOI=econlit_diff$digitalObjectIdentifier,Link=econlit_diff$DocumentURL,Document.Type=econlit_diff$ArticleType,Publication.Stage=NA,Access.Type=NA,Source="EconLit",EID=NA,X=NA)

compiled=rbind(compiled,wos_diff_df)
compiled=rbind(compiled,econlit_diff_df)

#subset only to published papers
compiled_papersonly=compiled[which(compiled$Document.Type%in%c("Article","Letter","J","Scholarly Journals")),]

dir.create("outputs/abstract_search", showWarnings = FALSE)
write.csv(compiled_papersonly,"outputs/abstract_search/compiledpapers_20200930_finalreview.csv")
