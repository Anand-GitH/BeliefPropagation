#######################################################################
#Demonstrate belief propogation and d-separation and predictions
#on the DAG-236 observations on fourteen variables from the Danish Heart Clinic
#
#Modified by : Anand
#Modified on: 05/05/2020
######################################################################

rm(list = ls())
graphics.off()

library(gRain)
library(Rgraphviz)
library(gRbase)
library(ggm)
library(gRim)
library(bnlearn)
library(igraph)

##########Function to find d separation in the graph#################################
ffinddsepingraph<-function(cdagmat){
  gnodes<-row.names(cdagmat)
  visitnodes<-vector()
  for(i in 1:length(gnodes)){
    for(j in 1:length(gnodes)){
      if(j!=i & !(gnodes[j] %in% visitnodes)){
        #cat("Node i ",gnodes[i],"Node j ",gnodes[j],"\n")
        if(dSep(cdagmat,gnodes[i],gnodes[j],cond=NULL)){
          cat("Nodes ",gnodes[i]," and ",gnodes[j]," are d-separated given  NULL","\n")
        }
        
        joincond<-vector()
        
        for (k in 1:length(gnodes)){
          if(k!=i & k!=j){
            #cat("Node i ",gnodes[i],"Node j ",gnodes[j]," given Node k ",gnodes[k],"\n")
            if(dSep(cdagmat,first=gnodes[i],second=gnodes[j],cond=gnodes[k])){
              cat("Nodes ",gnodes[i]," and ",gnodes[j]," are d-separated given ",gnodes[k],"\n")
            }
            
            joincond<-append(joincond,gnodes[k],length(joincond))
            
            if(length(joincond)>1){
              #cat("Node i ",gnodes[i],"Node j ",gnodes[j]," given joint cond ",joincond,"\n")
              if(dSep(cdagmat,first=gnodes[i],second=gnodes[j],cond=joincond)){
                cat("Nodes ",gnodes[i]," and ",gnodes[j]," are d-separated given ",joincond,"\n")
              }
            }
          }
        }
      }
    }
    visitnodes<-append(visitnodes,gnodes[i],length(visitnodes))
  }
}

######################################################################################################

#load data cad1 Coronary artery disease data
data(cad1)
?cad1

#236 records and 14 variables
dim(cad1)
names(cad1)

#1-  "Sex"
#8-  "SuffHeartF"
#10- "Hyperchol"   
#11- "Smoker"     
#12- "Inherit"    
#14- "CAD" 

cad2<-cad1[,c(1,8,10,11,12,14)]
names(cad2)
head(cad2)
dim(cad2)

###########Lets construct the network using the few variables mentione##
#Sex        - Female or Male
#Smoker     - Yes or No
#SuffHeartF - Yes or No
#Inherit    - Yes or No
#Hyperchol  - Yes or No
#CAD        - Yes or No 
############With the graph given we see that the Coronary artery disease is derived from other variables########

###############Creating a DAG###################
g <- list(~Sex, ~Smoker|Sex, ~SuffHeartF, ~Inherit|Smoker, ~Hyperchol|Smoker:SuffHeartF, ~CAD|Inherit:Hyperchol)


cdag<-dagList(g)
plot(cdag)
##################Inquire for d-separations################
cdagmat<-as(cdag,"matrix")
ffinddsepingraph(cdagmat)
###########################################################
#To Get conditional probability tables
yn<-c("Yes","No")
mf<-c("Male","Female")
sexcpt<- cptable(~Sex,values=c(0.8,0.2),levels = mf)
smok.sex<- cptable(~Smoker|Sex,values =c(0.6,0.4,0.99,0.01),levels = yn)
suffh<-cptable(~SuffHeartF,values = c(0.6,0.4),levels = yn)
inh.smok<-cptable(~Inherit|Smoker,values = c(0.8,0.2,0.6,0.4),levels = yn)
hyp.smok.suff<-cptable(~Hyperchol|Smoker:SuffHeartF,values = c(0.5,0.5,0.7,0.3,0.9,0.1,0.4,0.6),levels = yn)
cad.inh.hyp<-cptable(~CAD|Inherit:Hyperchol,values = c(0.6,0.4,0.5,0.5,0.8,0.2,0.4,0.6),levels = yn)

##################################################################################
##Build the network
#################################################################################
plist<- compileCPT(list(sexcpt, smok.sex, suffh, inh.smok, hyp.smok.suff, cad.inh.hyp))
grn1<- grain(plist)
summary(grn1)

graphics.off()
plot(grn1)

cptdata<-extractCPT(cad2,grn1$dag)
cptdata

##Compile the network
grn1c<- compile(grn1)
summary(grn1c)
plot(grn1c)

###Propogate the network
grn1c<- propagate(grn1c)
summary(grn1c)

####################Test for different evidences#######################################
##Female with hypercholestrol and check for heart failure and CAD probabilities differ
summary(grn1c)
grn1c.ev<- setFinding(grn1c,nodes=c("Sex","Hyperchol"), states=c("Female","Yes"))

#Probability of finding
getFinding(grn1c.ev)

#Not absorbed  - Before evidence
querygrain(grn1c,nodes = c("CAD","SuffHeartF"),type="marginal")

#Abssorbed  - After evidence
querygrain(grn1c.ev,nodes = c("CAD","SuffHeartF"),type="marginal")

#####################################################################################
#25 observations
####################################################################################
sim25<-simulate.grain(grn1c.ev,nsim=25,seed=100)
table(sim25$SuffHeartF)
table(sim25$CAD)

Heartfailurefreq<-table(sim25$SuffHeartF)
probHfail25<-Heartfailurefreq[1]/length(sim25$SuffHeartF)
probHfail25

CADfreq<-table(sim25$CAD)
probCAD25<-CADfreq[1]/length(sim25$CAD)
probCAD25

#####################################################################################
#500 observations
####################################################################################
sim500<-simulate.grain(grn1c.ev,nsim=500,seed=100)
table(sim500$SuffHeartF)
table(sim500$CAD)

Heartfailurefreq<-table(sim500$SuffHeartF)
probHfail500<-Heartfailurefreq[1]/length(sim500$SuffHeartF)
probHfail500

CADfreq<-table(sim500$CAD)
probCAD500<-CADfreq[1]/length(sim500$CAD)
probCAD500

################################The End##########################################