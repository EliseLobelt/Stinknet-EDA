##Code back up
#Libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(janitor)
library(tidyr)
library(vegan)

### ----data---- ###
data <- read_excel("/Users/elise/OneDrive/Documents/Desktop/Global Futures research/Cleaned_Data.xlsx")
data<- data %>% select(-"nitrB", -"phosB", -"ocarbB",-"ammonB")
head(data)
data[c("cnT","cnB")]<-sapply(data[c("cnT","cnB")],as.numeric)
data[data<0] <- 0
## round numbers in the columns clay-phosB
df <- data %>% mutate(across(where(is.numeric), round, 3))
## remove empty rows
df <- na.omit(df)


### Filter data based on Stinknet="TRUE" (sn_true) ###
df_true <- df[df$ONCPIL == "TRUE", ]
df_false <- df[df$ONCPIL == "FALSE", ]

### mosiac plot spread of stinknet presence and status of urban scenery ###
ggplot(df, aes(x = urban, fill = ONCPIL)) +
  geom_bar(position = "dodge") +
  labs(title = "Stinknet occurences in Urban/Non Urban areas",
       x = "Urban 0/1",
       y = "Count",
       fill = "Stinknet")

### average Amounts of nutrients in the soils where Sinknet is found ###
sn_nutrient <- df_true[,c("ammonT","nitrT","phosT")]
nutrient_averages <- colMeans(sn_nutrient)
print(nutrient_averages,3)

### Contigency table of broadclasses and stinknet ###
snbc_df = data.frame(df$broad_classes,df$ONCPIL)
snbc_df = table(df$broad_classes,df$ONCPIL) 
print(chisq.test(snbc_df))  
print(snur_df) #Where Stinknet: T/F
## chi squared test



#### ---initial data analysis--- ####

sndf <-df[,c("ammonT","nitrT", "phosT", "pH","conduct", 
             "bulkdens","ocarbT", "ONCPIL", "clay","sand")]
# create a dataset
stinknet <-sndf$ONCPIL 
condition <- rep(c("ammonT","nitrT","phosT") , 6)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(stinknet,condition,value)


#while stinknet is present
df_true <- df_true[,c("ammonT",
                      "nitrT",
                      "phosT",
                      "ocarbT")]
nutrient_averages <- colMeans(df_true)
print(nutrient_averages,3)

#while Stinknet is NOT present
df_false <- df_false[,c("ammonT",
                        "nitrT",
                        "phosT",
                        "ocarbT")]
nutrient_averages <- colMeans(df_false)
print(nutrient_averages,3)


#-------------------NMDS-----------------------------------------------

#create dataset:
sndf <-df[,c("ammonT", "nitrT", "phosT", "pH",
             "conduct", "bulkdens","ocarbT", "ONCPIL","clay","sand")]
#COLUMNS: such that the dataset includes only the columns with the needed metadata 
##(treatments, sites, etc.) and the actual species being analyzed.

#Then break into metadata and species data (AFTER remove rows/NAs):
SN<-sndf[,"ONCPIL"] #stinknet t/f
SOILSPECIE <-sndf[,c("clay","sand")] #different chemistry (soil)

#run NMDS:
Ord1<-metaMDS(SOILSPECIE, wascores=TRUE, k=3) #relationship between specie and site
Ord1
print(scores(Ord1, display="sites"))

#And record stress from results.
#wascores designates to use the weighted average for species scores.




### ---Vegan and NMDS--- ###

#create dataset:
soil <- read.table("/Users/elise/OneDrive/Documents/Desktop/Global Futures research/soil_df.txt")
#df <- df %>% select(-col1, -col2)
soil_df<- soil %>% select(-"nitrB", -"phosB", -"ocarbB",-"ammonB")
META <-df[,c("ONCPIL", "broad_classes")]
head(META)
head(soil_df)

###  Run NMDS:  ###
Ord1<-metaMDS(soil_df, wascores=TRUE, k=3) #relationship between specie and site
Ord1
print(scores(Ord1, display="sites"))


library(ggplot2)
plot(Ord1)
colfactorn<-factor(META$ONCPIL)
points(Ord1, display="sites",cex=1, pch=16, col=c(1,4,3,2)[colfactorn])
text(Ord1, display="species", cex=1)
ordiellipse(Ord1, colfactorn, col=c(1,4,3,2), kind="ehull", lwd=3)
legend(x="topleft", legend=levels(colfactorn), col=c(1,4,3,2), pch=16)

#Turn site scores into a dataset:
MDS1 = Ord1$points[,1]
MDS2 = Ord1$points[,2]
Ord1XY = data.frame(MDS1 = MDS1, MDS2 = MDS2)
print(Ord1XY)

###   ANOSIM:  ###
summary(anosim(Ord1XY, META$ONCPIL == "TRUE"))

###  or PERMANOVA: ###
print(adonis2(Ord1XY~META$ONCPIL, META))
print(adonis2(Ord1XY~META$ONCPIL*META$broad_classes, META))


#----------Distance matrix --------------------
#Create distance matrix:
DistRESULTS<-vegdist(soil_df)

#run PERMANOVA on it:
print(adonis2(DistRESULTS~META$ONCPIL, META))

#Followed by pairwise comparisons of treatment levels if significant 
##(have to do each treatment separately):
library(RVAideMemoire)
pairwise.perm.manova(DistRESULTS, META$ONCPIL) 
pairwise.perm.manova(DistRESULTS, META$broad_classes)


