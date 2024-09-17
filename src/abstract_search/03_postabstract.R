#post abstract screening

#identify additional abstracts for screening based on second review

abstracts=read.csv("data/abstract_search/compiledpapers_20200930_finalreview.csv")

included=read.csv("data/abstract_search/currentpapers.csv")

match=which(abstracts$ID_number%in%included$ID_Number)

#take out abstracts already included in review
abstracts=abstracts[-match,]

#find abstracts not currently included that have a 1 in first, second, or third RA review
toreview=which(abstracts$Abstract.Filter==1|abstracts$Second.NA.Filter==1|abstracts$Third.Abstract.Filter==1)

newabstracts=abstracts[toreview,]

#double check titles don't match ones already being reviewed
check=which(newabstracts$Title%in%included$Title)
newabstracts=newabstracts[-check,]
write.csv(newabstracts,file="outputs/abstract_search/newabstracts.csv")
