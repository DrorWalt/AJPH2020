###############################################################################################
## Russian Twitter Accounts and the Partisan Polarization of Vaccine Discourse, 2015-2017:   ##
##    Supplementary code                                                                     ##
##                                                                                           ##
## Author: Dror Walter                                                                       ##
##                                                                                           ##
## To Cite:                                                                                  ##
## Walter D., Ophir Y. & Hall Jamieson, K. (2020)                                            ##
## Russian Twitter Accounts and the Partisan Polarization of Vaccine Discourse, 2015-2017.   ##
## American Journal of Public Health. http://dx.doi.org/10.2105/AJPH.2019.305564             ##
## First published Online (19 March 2020)                                                    ##
##                                                                                           ##
###############################################################################################


################################ Importing Libraries

options(stringsAsFactors = F)
library(stringi) 
library(stringr)
library(qdap)
library(tm)
library(ggplot2)
library(lubridate)
library(irr)
library(quanteda)
library(ldatuning)
library(topicmodels)
library(xlsx)
library(textcat)
library(parallel)
library(RSQLite)
library(doParallel)
library(scales)
library(lsa)
library(igraph)
library(cld2) 
library(tidyverse)
library(dplyr)
library(rgexf)



############################## Importing and Cleaning the Data

# Read twitter data (available at Twitter Election Integritiy Database)
text.df <- read.csv("raw_data/ira_tweets_csv_hashed.csv", stringsAsFactors = F)
# limit to English tweets - based on twitter metadata
text.en <- text.df[text.df$tweet_language == "en",]
# Convect to lower case
text.en$tweet_text2 <- tolower(text.en$tweet_text)
# limit to English tweets as identified by google cld
text.en$textcatcld2<- cld2::detect_language(text.en$tweet_text)
text.en<-text.en[which(text.en$textcatcld2=='en'),]
# correct encoding
text.en$tweet_text3<-iconv(text.en$tweet_text2, "latin1", "ASCII", sub="")
text.en$tweet_text3<-iconv(text.en$tweet_text3, "UTF-8", "ASCII", sub="")
text.df<-text.en
rm(text.en)

# Cleaning Twitter artifacts (links, images etc.)
text.df$tweet_text3<-gsub("[h][t][t][p][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[h][t][t][p][s][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[p][i][c][.][^[:space:]]+","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[.][c][o][m][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[w][w][w][.][^[:space:]]+","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[.][l][y][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[y][o][u][t][u][.][b][e][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[.][c][o][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[.][c][m][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[.][o][r][g][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[^[:space:]]*[w][a][p][t][o][.][s][t][^[:space:]]*","",text.df$tweet_text3)
text.df$tweet_text3<-gsub("[&][a][m][p]"," ",text.df$tweet_text3)
# arranging the data and renaming columns
data <- text.df
colnames(data)[13]<-"orig_text"
colnames(data)[34]<-"text"
# setting data as date format
data$date2 <- as.Date(data$tweet_time, '%Y-%m-%d')
# adding index column
data$index<-seq(1,nrow(data))
# limiting dates to relevant
data<-data[data$date2<"2018-01-01",]
data<-data[data$date2>"2015-01-01",]


#########################################################################
###########                                                       #######
########### LDA PIPELINE- corpus and hyperparameter grid search   #######
###########                                                       #######
#########################################################################

# removing extremely short text data
removed_short<-subset(data,nchar(as.character(data$text))<6) 
data2<-subset(data,!nchar(as.character(data$text))<6)
# removing duplicate tweets for data3 and saving removed in DF. 
# will be added after model estimation tp calculate sailence of all texts
removed_df<-data2[duplicated(data2$text),] 
data3 <- data2[!duplicated(data2$text),]
# sampling 10% for K and Hyperparameter optimization
data3_10perc_numbers<-sample((1:nrow(data3)),(nrow(data3)/10),replace=FALSE)
data3_10perc<-data3[data3_10perc_numbers,]
# Building the corpus data - Notice the additiomal stopwords and the removal of extremly common/rare words
mycorpus <- corpus(data3_10perc)
stopwords_and_single<-c(stopwords("english"),LETTERS,letters, "t.co", "http", "https", "rt", "p", "amp", "via")
dfm_counts <- dfm(mycorpus,tolower = TRUE, remove_punct = TRUE,remove_numbers=TRUE, 
                  remove = stopwords_and_single,stem = FALSE,
                  remove_separators=TRUE) 
docnames(dfm_counts)<-dfm_counts@docvars$index
# removing extemely common or rare tokens
dfm_counts2<-dfm_trim(dfm_counts, max_docfreq = 0.95, min_docfreq=0.0001,docfreq_type="prop")
# convert to LDA-ready object
dtm_lda <- convert(dfm_counts2, to = "topicmodels",docvars = dfm_counts2@docvars)
full_data<-dtm_lda
# count numer of documents for crossvalidation
n <- nrow(full_data)
# clean temp data
rm(text.df)
rm(dfm_counts)
rm(dfm_counts2)

# Run the crossvalidation loop
print(Sys.time())
# create container for results
MainresultDF<-data.frame(k=c(1),perplexity=c(1),myalpha=c("x"))
MainresultDF<-MainresultDF[-1,]
# set possible alpha and k values
candidate_alpha<- c(0.01, 0.05, 0.1, 0.2, 0.5) # candidates for alpha values
candidate_k <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120) # candidates for how many topics
# run the 10-fold cross validation
for (eachalpha in candidate_alpha) { 
  print ("now running ALPHA:")
  print (eachalpha)
  print(Sys.time())
  cluster <- makeCluster(detectCores(logical = TRUE) - 1) # We are leaving one Core spare. If number of corse on pc is 1, then -1 in this line should be removed.
  registerDoParallel(cluster)
    clusterEvalQ(cluster, {
    library(topicmodels)
  })
  folds <- 10
  splitfolds <- sample(1:folds, n, replace = TRUE)
  clusterExport(cluster, c("full_data", "splitfolds", "folds", "candidate_k"))
  system.time({
    results <- foreach(j = 1:length(candidate_k), .combine = rbind) %dopar%{
      k <- candidate_k[j]
      print(k)
      results_1k <- matrix(0, nrow = folds, ncol = 2)
      colnames(results_1k) <- c("k", "perplexity")
      for(i in 1:folds){
        train_set <- full_data[splitfolds != i , ]
        valid_set <- full_data[splitfolds == i, ]
        
        fitted <- LDA(train_set, k = k, method = "Gibbs",
                      control = list(alpha=eachalpha) )
        
        results_1k[i,] <- c(k, perplexity(fitted, newdata = valid_set))
      }
      return(results_1k)
    }
  })
  stopCluster(cluster)
  results_df <- as.data.frame(results)
  results_df$myalpha<-as.character(eachalpha)
  MainresultDF<-rbind(MainresultDF,results_df)
}
print ("DONE!")
print(Sys.time())

# arrange and examine results
MainresultDF$kalpha=paste0(as.character(MainresultDF$k),MainresultDF$myalpha) 
ggplot(MainresultDF) +geom_boxplot(aes(x=k, y=perplexity, group=kalpha,color=myalpha))

# run additional k and alpha values
candidate_alpha<- c(0.01,0.05)
candidate_k <- c(130,140,150,160) # candidates for how many topics
for (eachalpha in candidate_alpha) { 
  print ("now running ALPHA:")
  print (eachalpha)
  print(Sys.time())
  cluster <- makeCluster(detectCores(logical = TRUE) - 1) # leave one CPU spare...
  registerDoParallel(cluster)
  clusterEvalQ(cluster, {
    library(topicmodels)
  })
  folds <- 10
  splitfolds <- sample(1:folds, n, replace = TRUE)
  clusterExport(cluster, c("full_data", "splitfolds", "folds", "candidate_k"))
  system.time({
    results <- foreach(j = 1:length(candidate_k), .combine = rbind) %dopar%{
      k <- candidate_k[j]
      print(k)
      results_1k <- matrix(0, nrow = folds, ncol = 2)
      colnames(results_1k) <- c("k", "perplexity")
      for(i in 1:folds){
        train_set <- full_data[splitfolds != i , ]
        valid_set <- full_data[splitfolds == i, ]
        
        fitted <- LDA(train_set, k = k, method = "Gibbs",
                      control = list(alpha=eachalpha) )
        
        results_1k[i,] <- c(k, perplexity(fitted, newdata = valid_set))
      }
      return(results_1k)
    }
  })
  stopCluster(cluster)
  NEWresults_df <- as.data.frame(results)
  NEWresults_df$myalpha<-as.character(eachalpha)
  MainresultDF$kalpha<-paste0(as.character(MainresultDF$k),MainresultDF$myalpha)  
  NEWresults_df$kalpha<-paste0(as.character(NEWresults_df$k),NEWresults_df$myalpha) 
  MainresultDF<-rbind(MainresultDF,NEWresults_df)
}
print ("DONE!")
print(Sys.time())
# examine results
ggplot(MainresultDF) +
  geom_boxplot(aes(x=k, y=perplexity, group=kalpha,color=myalpha))+
  geom_smooth(se = TRUE, aes(x=k, y=perplexity,color=myalpha))

# run last k values
candidate_alpha<- c(0.01)
candidate_k <- c(170,180,190,200) # candidates for how many topics
for (eachalpha in candidate_alpha) { 
  print ("now running ALPHA:")
  print (eachalpha)
  print(Sys.time())
  cluster <- makeCluster(detectCores(logical = TRUE) - 1) # leave one CPU spare...
  registerDoParallel(cluster)
  clusterEvalQ(cluster, {
    library(topicmodels)
  })
  folds <- 10
  splitfolds <- sample(1:folds, n, replace = TRUE)
  clusterExport(cluster, c("full_data", "splitfolds", "folds", "candidate_k"))
  system.time({
    results <- foreach(j = 1:length(candidate_k), .combine = rbind) %dopar%{
      k <- candidate_k[j]
      print(k)
      results_1k <- matrix(0, nrow = folds, ncol = 2)
      colnames(results_1k) <- c("k", "perplexity")
      for(i in 1:folds){
        train_set <- full_data[splitfolds != i , ]
        valid_set <- full_data[splitfolds == i, ]
        
        fitted <- LDA(train_set, k = k, method = "Gibbs",
                      control = list(alpha=eachalpha) )
        
        results_1k[i,] <- c(k, perplexity(fitted, newdata = valid_set))
      }
      return(results_1k)
    }
  })
  stopCluster(cluster)
  NEWresults_df <- as.data.frame(results)
  NEWresults_df$myalpha<-as.character(eachalpha)
  MainresultDF$kalpha<-paste0(as.character(MainresultDF$k),MainresultDF$myalpha)  
  NEWresults_df$kalpha<-paste0(as.character(NEWresults_df$k),NEWresults_df$myalpha) 
  MainresultDF<-rbind(MainresultDF,NEWresults_df)
}
print ("DONE!")
print(Sys.time())

# examine full k and alpha resutls
ggplot(MainresultDF) +
  geom_boxplot(aes(x=k, y=perplexity, group=kalpha,color=myalpha))+
  scale_color_discrete(name = "Alpha Levels")+
  xlab("K (Number of Topics)")+
  ylab("Perplexity")
# Identify 2nd derivative max point on perplexity  
MainresultDF_MYALPHA<-MainresultDF[MainresultDF$myalpha==0.01,]
cars.spl <- with(MainresultDF_MYALPHA, smooth.spline(k, perplexity, df = 3))
plot(with(cars, predict(cars.spl, x = MainresultDF_MYALPHA$k, deriv = 2)), type = "l")
abline(v=60)

# Run LDA with optimal values of alpha 0.01 and k=60 on full data
mycorpus <- corpus(data3)
stopwords_and_single<-c(stopwords("english"),LETTERS,letters, "t.co", "http", "https", "rt", "p", "amp", "via")
dfm_counts <- dfm(mycorpus,tolower = TRUE, remove_punct = TRUE,remove_numbers=TRUE, 
                  remove = stopwords_and_single,stem = FALSE,
                  remove_separators=TRUE) 
docnames(dfm_counts)<-dfm_counts@docvars$index
dfm_counts2<-dfm_trim(dfm_counts, max_docfreq = 0.95, min_docfreq=0.0001,docfreq_type="prop")
dtm_lda <- convert(dfm_counts2, to = "topicmodels",docvars = dfm_counts2@docvars)

LDA.60 <- LDA(dtm_lda, k = 60, method = "Gibbs",control = list(alpha=0.01,seed=125231)) 
LDAfit<-LDA.60

#########################################################################
###########                                                       #######
###########            Analyzing the Topics                       #######
###########                                                       #######
#########################################################################

# Printing main files for analysis (WORDS/FREX/TEXTS)
## setting the text column
datacolnum=13

## Funciton to print Beta, Frex and Theta
extract_topic_xls<-function (eachLDA) {
  LDAfit<-eachLDA
  
  mybeta<-data.frame(LDAfit@beta)
  colnames(mybeta)<-LDAfit@terms
  mybeta<-t(mybeta)
  colnames(mybeta)<-seq(1:ncol(mybeta))
  mybeta=exp(mybeta)
  
  #### Now we cycle and print top words for each topic
  nwords=50
  topwords <- mybeta[1:nwords,]
  for (i in 1:LDAfit@k) {
    tempframe <- mybeta[order(-mybeta[,i]),]
    tempframe <- tempframe[1:nwords,]
    tempvec<-as.vector(rownames(tempframe))
    topwords[,i]<-tempvec
  }
  rownames(topwords)<-c(1:nwords)
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(topwords, paste0(kalpha,"_ALLBOTS_Topwords.xlsx"))
  
  #### Get Frex (unique) words
  #### get the beta
  mybeta<-data.frame(LDAfit@beta)
  colnames(mybeta)<-LDAfit@terms
  mybeta<-t(mybeta)
  colnames(mybeta)<-seq(1:ncol(mybeta))
  mybeta=exp(mybeta)
  #### apply FREX formula below
  myw=0.3
  word_beta_sums<-rowSums(mybeta)
  my_beta_for_frex<-mybeta
  for (m in 1:ncol(my_beta_for_frex)) {
    for (n in 1:nrow(my_beta_for_frex)) {
      my_beta_for_frex[n,m]<-1/(myw/(my_beta_for_frex[n,m]/word_beta_sums[n])+((1-myw)/my_beta_for_frex[n,m]))
    }
    print (m)
  }
  #### print top 50 frex:
  nwords=50
  topwords <- my_beta_for_frex[1:nwords,]
  for (i in 1:LDAfit@k) {
    tempframe <- my_beta_for_frex[order(-my_beta_for_frex[,i]),]
    tempframe <- tempframe[1:nwords,]
    tempvec<-as.vector(rownames(tempframe))
    topwords[,i]<-tempvec
  }
  rownames(topwords)<-c(1:nwords)
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(topwords,paste0(kalpha,"_ALLBOTS_TopFREX.xlsx"))
  
  #### TOP TEXTS --->
  data33<-data3
  data33$index<-as.character(data33$index)
  deleted_lda_texts<-(setdiff(data33$index, LDAfit@documents))
  '%!in%' <- function(x,y)!('%in%'(x,y))
  data33<-data33[data33$index %!in% deleted_lda_texts,]
  metadf<-data33
  meta_theta_df<-cbind(metadf[datacolnum],LDAfit@gamma)
  ntext=50
  toptexts <- mybeta[1:ntext,]
  for (i in 1:LDAfit@k) {
    print(i)
    tempframe <- meta_theta_df[order(-meta_theta_df[,i+1]),]
    tempframe <- tempframe[1:ntext,]
    tempvec<-as.vector(tempframe[,1])
    toptexts[,i]<-tempvec
  }
  rownames(toptexts)<-c(1:ntext)
  kalpha<-paste0(as.character(LDAfit@k),"_",gsub("\\.","",as.character(LDAfit@alpha)))
  write.xlsx(toptexts, paste0(kalpha,"_ALLBOTS_TopTexts.xlsx"))
}

## Apply function to model
extract_topic_xls(LDAfit)

# returning the duplicate texts deleted earlier
data33<-data3
data33$index<-as.character(data33$index)
deleted_lda_texts<-(setdiff(data33$index, LDAfit@documents))
#deleted_lda_texts2<-(setdiff(as.character(LDAfit@documents),as.character(data3$doc_id)))
#deleted_lda_texts<-unique(c(deleted_lda_texts1,deleted_lda_texts2))
'%!in%' <- function(x,y)!('%in%'(x,y))
#data33<-data3
data33<-data33[data33$index %!in% deleted_lda_texts,]
metadf<-data33
meta_theta_df<-cbind(metadf,LDAfit@gamma)
removed_df2<-inner_join(removed_df,meta_theta_df,by="text")
removed_df2<-removed_df2[,-c(37:73)]
colnames(removed_df2)<-gsub("\\.x","",colnames(removed_df2))
removed_df2$index<-as.character(removed_df2$index)
meta_theta_df2<-bind_rows(meta_theta_df,removed_df2)
meta_theta_df<-meta_theta_df2
rm(meta_theta_df2)

# marking the users and tweets which are vaccine related
vaccine_words <-  c("vaccin")
vactweets <- meta_theta_df[grep(paste0(vaccine_words, collapse = "|"), meta_theta_df$tweet_text, value=FALSE),]
vacc_users<-unique(vactweets$userid)
meta_theta_df2<-meta_theta_df
meta_theta_df2$vacc_user<-"no_vacc"
meta_theta_df2[which(meta_theta_df2$userid %in% vacc_users),"vacc_user"]<-"yes_vac"


# Running ANTMN function - 
# cite: Walter and Ophir (2019) News Frame Analysis. Communication Methods and Measures 13(4), 248-266  
network_from_LDA<-function(LDAobject,deleted_topics=c(),topic_names=c(),save_filename="",topic_size=c()) {
  # Importing needed packages
  require(lsa) # for cosine similarity calculation
  require(dplyr) # general utility
  require(igraph) # for graph/network managment and output
  
  print("Importing model")
  
  # first extract the theta matrix form the topicmodel object
  theta<-LDAobject@gamma
  # adding names for culumns based on k
  colnames(theta)<-c(1:LDAobject@k)
  
  # claculate the adjacency matrix using cosine similarity on the theta matrix
  mycosine<-cosine(as.matrix(theta))
  colnames(mycosine)<-colnames(theta)
  rownames(mycosine)<-colnames(theta)
  
  # Convert to network - undirected, weighted, no diagonal
  
  print("Creating graph")
  
  topmodnet<-graph.adjacency(mycosine,mode="undirected",weighted=T,diag=F,add.colnames="label") # Assign colnames
  # add topicnames as name attribute of node - importend from prepare meta data in previous lines
  if (length(topic_names)>0) {
    print("Topic names added")
    V(topmodnet)$name<-topic_names
  } 
  # add sizes if passed to funciton
  if (length(topic_size)>0) {
    print("Topic sizes added")
    V(topmodnet)$topic_size<-topic_size
  }
  newg<-topmodnet
  
  # delete 'garbage' topics
  if (length(deleted_topics)>0) {
    print("Deleting requested topics")
    
    newg<-delete_vertices(topmodnet, deleted_topics)
  }
  
  # run community detection and attach as node attribute
  print("Calculating communities")
  
  mylouvain<-(cluster_louvain(newg)) 
  mywalktrap<-(cluster_walktrap(newg)) 
  myspinglass<-(cluster_spinglass(newg)) 
  myfastgreed<-(cluster_fast_greedy(newg)) 
  myeigen<-(cluster_leading_eigen(newg)) 
  
  V(newg)$louvain<-mylouvain$membership 
  V(newg)$walktrap<-mywalktrap$membership 
  V(newg)$spinglass<-myspinglass$membership 
  V(newg)$fastgreed<-myfastgreed$membership 
  V(newg)$eigen<-myeigen$membership 
  
  # if filename is passsed - saving object to graphml object. Can be opened with Gephi.
  if (nchar(save_filename)>0) {
    print("Writing graph")
    write.graph(newg,paste0(save_filename,".graphml"),format="graphml")
  }
  
  # graph is returned as object
  return(newg)
}

costum_topic_labels<-c("Crime", "Disasters and disruptions", "Aphorisms", "Municipalities", "BlackLivesMatter", "Arrests", "Clinton emails", "Victims and deaths", "Leisure activities", "Race whites and blacks", "North Korea", "education", "Sports", "Police brutality", "Extremist groups", "Finances", "Legal issues", "Black history", "Terrorist attacks", "Watch listen", "Criminals", "Presidential elections", "Mixed and science", "Islamist violence", "Obama and Trump", "Laws ", "Sports", "Blessings holidays congrats", "Celebs", "Imitating black southern language", "Trump activities and nominations", "RT of a user writing about loans", "Companies and business", "Diets", "Crime", "Sports", "Immigration", "Pro Trump Anti his enemies", "GOP primaries", "Aphorisms", "Weather", "Anti media and Dems", "Mixed mundane", "Taxes budgets money", "Bills", "Supporting Trump", "Primaries", "Attacks on fake media", "Ukranian nuclear crisis", "Workout", "Sports", "Mixed", "Health and science", "Sports", "Accidents", "Mueller invistigation", "Foreign affairs", "Mixed PP and personal attacks", "Sexual misconduct", "Conservative tweets from pjnet and tcot")

mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           topic_names=mynames,
                           save_filename="TOPIC_ALL_BOT_NET")


#########################################################################
###########                                                       #######
###########        Thematic Communities Analysis                  #######
###########                                                       #######
#########################################################################

# prepare data on user-level aggregation
## average topic loading per user
meta_theta_by_user<-aggregate(x = meta_theta_df2[,c(37:(ncol(meta_theta_df2)-1))], 
                                   by = list(meta_theta_df2$userid), FUN = "mean")
## retweet sum per user
meta_theta_by_user_retweets<-aggregate(x = meta_theta_df2[,c(26)], 
                              by = list(meta_theta_df2$userid), FUN = "sum")
## user engagemnt with vaccines
meta_theta_by_user_ynvac<-aggregate(x = meta_theta_df2[,97], 
                              by = list(meta_theta_df2$userid), FUN = "unique")
## user lifetime volume of tweets
meta_theta_df3<-meta_theta_df2
meta_theta_df3$forsum<-1
meta_theta_by_user_volume<-aggregate(x = meta_theta_df3[,"forsum"], 
                                    by = list(meta_theta_df3$userid), FUN = "sum")
rm(meta_theta_df3)
## combine to user level data
rownames(meta_theta_by_user)<-meta_theta_by_user$Group.1
meta_theta_by_user<-meta_theta_by_user[,-1]
meta_theta_by_user2<-t(meta_theta_by_user)

# correalte users by their topic loading
mybeta_cor<-cosine(as.matrix(meta_theta_by_user2))

# create the network
sem_net_weighted<-graph.adjacency(mybeta_cor,mode="undirected",weighted=T,diag=F,add.colnames="label") # Assign colnames

# add node level metadata from previous aggregation
V(sem_net_weighted)$name<-V(sem_net_weighted)$label
V(sem_net_weighted)$ynvac<-meta_theta_by_user_ynvac$x
V(sem_net_weighted)$size<-meta_theta_by_user_volume$x
V(sem_net_weighted)$retw<-meta_theta_by_user_retweets$x

# Keep significant edges p<0.05
set.seed(763423)
g<-disparity_filter(g=sem_net_weighted,alpha=0.05)

# run community detection
set.seed(433547)
V(g)$louvain<-(cluster_louvain(g)) $membership 

# export the graph
saveAsGEXF(g,"outputTEMP/TROLL_USER_NET_bbone_new.gexf")
print("done!")

# Explore the Communities 
## Find users most central for each community
### get sum of edgeds of each node and inside edges - calculate for each node
nodelist<-list()
for (node in 1:length(V(g))) {
  print(node)
  outside<-strength(g, vids = V(g)[node])
  tempg<-induced_subgraph(g,V(g)$louvain==V(g)$louvain[node])
  inside<-strength(tempg, vids = V(tempg)$label==V(g)[node]$label)
  nodelist[[node]]<-data.frame(
    node=node,label=V(g)[node]$label,inside,comm=V(g)$louvain[node],between=outside,within=inside,commstr=inside/outside)
}
user_comm_df<-do.call(rbind,nodelist)

### grab for each community the top 20 users
top_user_com_df<-data.frame(matrix(NA, nrow = 20, ncol = length(unique(user_comm_df$comm))))
for (i in 1:max(user_comm_df$comm)) {
  print (i)
  temp_df<-user_comm_df[user_comm_df$comm==i,]
  temp_df<-temp_df[order(temp_df$commstr,decreasing = TRUE),]
  towrite<-temp_df$label[1:20]
  top_user_com_df[,i]<-towrite
}

## print top tweets
### go to data - filter by top 20 users
### grab 200 random tweets - print to dataframe and XL
comm_tweets_list<-list()
for (i in 1:max(user_comm_df$comm)) {
    print(i)
    temp_meta_theta_df<-meta_theta_df2[meta_theta_df2$userid %in% top_user_com_df[,i],]
    temp_meta_theta_df<- temp_meta_theta_df[sample(nrow(temp_meta_theta_df), 200), ]
    comm_tweets_list[[i]]<-c(temp_meta_theta_df)
    write.xlsx(temp_meta_theta_df,paste0(as.character(i),"_COMM_200_tweets.xlsx"))
}

## for each comm get vaccine related tweets that appeared in it
vaccine_words <-  c("vaccin")
for (i in 1:max(user_comm_df$comm)) {
  print(i)
  temp_df<-user_comm_df[user_comm_df$comm==i,]
  temp_meta_theta_df<-meta_theta_df2[meta_theta_df2$userid %in% temp_df$label,]
  temp_meta_theta_df <- temp_meta_theta_df[grep(paste0(vaccine_words, collapse = "|"), 
                                                temp_meta_theta_df$tweet_text, value=FALSE),]
  if (nrow(temp_meta_theta_df)>0) {
    write.xlsx(temp_meta_theta_df,paste0(as.character(i),"_COMM_VAC_tweets.xlsx"))
  }
}

# calculate the salience of vaccine related content in each community
## by number of tweets
vacc_comm_ratio<-list()
for (i in 1:max(user_comm_df$comm)) {
  print(i)
  temp_df<-user_comm_df[user_comm_df$comm==i,]
  temp_meta_theta_df<-meta_theta_df2[meta_theta_df2$userid %in% temp_df$label,]
  nrow1<-nrow(temp_meta_theta_df)
  temp_meta_theta_df <- temp_meta_theta_df[grep(paste0(vaccine_words, collapse = "|"), 
                                                temp_meta_theta_df$tweet_text, value=FALSE),]
  nrow2<-nrow(temp_meta_theta_df)
  vacc_comm_ratio[[i]]<-c(i,(nrow2/nrow1))
}
toprint<-t(data.frame(vacc_comm_ratio))
rownames(toprint)<-c()
colnames(toprint)<-c("comm","ratio")
toprint<-as.data.frame(toprint)
toprint[order(toprint$ratio),]

## by number of users
vacc_user_ratio<-list()
for (i in 1:max(user_comm_df$comm)) {
  print(i)
  temp_df<-user_comm_df[user_comm_df$comm==i,]
  temp_meta_theta_df<-meta_theta_df2[meta_theta_df2$userid %in% temp_df$label,]
  nrow1<-length(unique(temp_meta_theta_df$userid))
  temp_meta_theta_df <- temp_meta_theta_df[grep(paste0(vaccine_words, collapse = "|"), 
                                                temp_meta_theta_df$tweet_text, value=FALSE),]
  nrow2<-length(unique(temp_meta_theta_df$userid))
  vacc_user_ratio[[i]]<-c(i,(nrow2/nrow1))
}

toprint<-t(data.frame(vacc_user_ratio))
rownames(toprint)<-c()
colnames(toprint)<-c("comm","ratio")
toprint<-as.data.frame(toprint)
toprint[order(toprint$ratio),]


## print the Top Words for each Community
commlabels<-c("1: Hard News",
              "2: Anti-Trump",
              "3: Pro-Trump",
              "4: Youth Talk and Celebs",
              "5: African American and BLM",
              "6: Mixed International",
              "7: Ukraine",
              "8: Soft News",
              "9: Retweets and Hashtag Games")

comm_as_string<-comm_as_string[-1,]

comm_words <- comm_as_string %>%
  unnest_tokens(word, text) %>%
  count(comm, word, sort = TRUE)

pd<-comm_words %>%
  filter(!word %in% c(stopwords_and_single,"just","will")) %>%  
  arrange(desc(n)) %>%
  mutate(word = factor(word, levels = rev(unique(word)))) %>% 
  group_by(comm) %>% 
  top_n(15) %>% 
  ungroup() %>%
  arrange(comm, n) %>%
  mutate(order = row_number())

ggplot(pd, aes(order, n)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ comm, scales = "free", labeller=labeller(comm = commlabels)) +
  xlab("") +
  ylab("Freq") +
  theme_bw() +
  # Add categories to axis
  scale_x_continuous(
    breaks = pd$order,
    labels = pd$word,
    expand = c(0,0)
  ) +
  ylim(0,120000)+
  coord_flip()



