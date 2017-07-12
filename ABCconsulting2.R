### Packages -----
require(XLConnect)
require(plyr)
require(dplyr)
require(tidyr)
require(corrplot)
require(ggplot2)
require(randomForest)
require(vegan)
require(ggrepel)
require(simputation)
require(HotDeckImputation)

### Data -----
wb = loadWorkbook('Book1.xlsx')
df = readWorksheet(wb, sheet = 'Random Sites 2008-2016 To R')

df <- plyr::rename(df,c('Station'='Station'
                        ,'Location'='Location'
                        ,'Year'='Year'
                        ,'StationType'='StationType'
                        ,'CSCI'='CSCI'
                        ,'MMI'='MMI'
                        ,'O.E'='OE'
                        ,'IBI'='IBI'
                        ,'D18'='D18'
                        ,'H20'='H2O'
                        ,'S2'='S2'
                        ,'BioticStructure'='% Biotic Structure'
                        ,'BufferLandscape'='% Buffer Landscape'
                        ,'Hydrology'='% Hydrology'
                        ,'OverallScore'='OverallScore'
                        ,'PhysicalStructure'='% Physical Structure'
                        ,'Eroded'='% Eroded'
                        ,'Stable'='% Stable'
                        ,'Vulnerable'='% Vulnerable'
                        ,'FastWater'='% Fast Water'
                        ,'SlowWater'='% Slow Water'
                        ,'ChannelAlteration'='Channel Alteration'
                        ,'EpifaunalSubstrate'='Epifaunal Substrate'
                        ,'SedimentDeposition'='Sediment Deposition'
                        ,'Alkalinity'='Alkalinity as CaCO3'
                        ,'DO'='DO'
                        ,'pH'='pH'
                        ,'Salinity'='Salinity'
                        ,'SpecificConductivity'='SpecificConductivity'
                        ,'Temperature'='Temperature'
                        ,'X.Slope'='Mean Slope'
                        ,'Discharge'='Discharge'
                        ,'WettedWidth'='Wetted Width'
                        ,'MicroalgaeThickness'='Microalgae Thickness'
                        ,'Macrophytes'='% Macrophytes'
                        ,'Macroalgae'='% Macroalgae'
                        ,'Cover'='% Cover Densiometer'
                        ,'CPOM'='% CPOM'
                        ,'sandfines'='% Sand and Fines'
                        ,'concrete.asphalt'='% Concrete and Asphalt'
                        ,'cobblegravel'='% Cobble and Gravel'
                        ,'Ammonia'='Ammonia as N'
                        ,'DOC'='Dissolved Organic Carbon'
                        ,'Hardness'='Hardness as CaCO3'
                        ,'NitrateNitriteN'='NitrateNitriteN'
                        ,'Nitrate'='Nitrate as N'
                        ,'Nitrite'='Nitrite as N'
                        ,'NitrogenTotal'='Total Nitrogen'
                        ,'TKN'='Total Kjeldahl Nitrogen'
                        ,'OrthoPhosphate'='Orthophosphate as P'
                        ,'Phosphorus'='Phosphorus as P'
                        ,'TSS'='TSS'
                        ,'TOC'='Total Organic Carbon'
                        ,'Chloride'='Chloride'
                        ,'Sulfate'='Sulfate'
                        ,'Arsenic'='Arsenic'
                        ,'Cadmium'='Cadmium'
                        ,'Chromium'='Chromium'
                        ,'Copper'='Copper'
                        ,'Iron'='Iron'
                        ,'Lead'='Lead'
                        ,'Mercury'='Mercury'
                        ,'Nickel'='Nickel'
                        ,'Selenium'='Selenium'
                        ,'Zinc'='Zinc'))

col_names <- names(dplyr::select(df,-c(Station,Location,Year,StationType)))
df[col_names] <- sapply(df[col_names],as.numeric)

df = dplyr::select(df,-c(MMI,OE,IBI,S2,D18,OverallScore,
                         DO,pH,Salinity,SpecificConductivity,Temperature,NitrateNitriteN)) %>%
  filter(Year!=2015)

df$categoryCSCI = ifelse(df$CSCI<=0.62,'Very Likely Altered',
                         ifelse(df$CSCI>0.62 & df$CSCI<=0.79,'Likely Altered',
                                ifelse(df$CSCI>0.79 & df$CSCI<=0.92,'Possibly Altered','Likely Intact')))

#data imputation
data_to_imp = dplyr::select(df,-c(Station,Location,Year,StationType,CSCI,H2O,categoryCSCI))
df_imputed = impute.NN_HD(data_to_imp,distance = 'eukl')

#checking imputation
df_mean = data.frame(
  mean_orig = apply(as.matrix(data_to_imp),2,function(x)mean(na.omit(x))),
  mean_imp = apply(as.matrix(df_imputed),2,mean)
)

#new data frame
df = cbind(
  Station = df$Station,
  Location = df$Location,
  Year = df$Year,
  StationType = df$StationType,
  CSCI = df$CSCI,
  H2O = df$H2O,
  df_imputed,
  categoryCSCI = df$categoryCSCI
)

#final data
df_csci = dplyr::select(df,-c(H2O)) %>% 
  drop_na()

df_h2o = dplyr::select(df,-c(CSCI,categoryCSCI)) %>% 
  drop_na()

### Exploratory Data Analysis -----

#considering correlations between variables
extracting_high_correlations <- function(data,threshold){
  correlations = cor(data,use = 'pairwise.complete.obs')
  
  cors_meeting_threshold = numeric()
  
  for(i in 1:dim(data)[2]){
    for(j in 1:dim(data)[2]){
      
      rowname = rownames(correlations)[i]
      colname = colnames(correlations)[j]
      correlation = correlations[i,j]
      
      if(abs(correlation) >= threshold & rowname < colname){
        cors_meeting_threshold = rbind(cors_meeting_threshold,c(rowname,colname,correlation))
      }else{
        NULL
      }
    }
  }
  
  cors_meeting_threshold = cors_meeting_threshold[order(cors_meeting_threshold[,3]),]
  list(correlations = correlations,cors_meeting_threshold = cors_meeting_threshold)
  
}
cor_func_output = extracting_high_correlations(dplyr::select(df_csci,-c(Station,Location,Year,StationType,categoryCSCI)),0.2)
View(cor_func_output$cors_meeting_threshold)

#plotting
dd = df_csci
for(var1 in names(dd)){
  for(var2 in names(dd)){
    
    if(var1 < var2){
      
      dftemp <- data.frame(v1 = dd[,var1], v2 = dd[,var2])
      
      jpeg(filename = paste0('~/Desktop/ConsultingProject2_ABCLabs_June2017/plots/',var1,'_',var2,'.jpeg'))
      plot(dftemp$v1, dftemp$v2, type = 'p')
      dev.off()
      
    }
    
  }
}

### Random Forest (Variable Importance) -----

rfvi <- function(data,response,percentage_train = 0.90,folds = 3,number_of_trees = 5000){
  
  #initializing results list
  results = list()
  
  #number of rows
  n = nrow(data)
  
  #mixing up the data
  scrambling = sample(x = 1:n,size = n,replace = FALSE)
  df_scramble = data[scrambling,]
  
  #training and validation sets
  break_point = ceiling(percentage_train*n) #due to limited data
  training = df_scramble[1:break_point,]
  validation = df_scramble[(break_point+1):n,]
  
  #cross-validation
  cv_param = rfcv(trainx = sqrt(training[,-which(names(training)%in%c(response))]),
                  trainy = training[,response],
                  cv.fold = folds) #due to limited data
  
  #optimum number of variables at each split
  results$optim_num_params = cv_param$n.var[which(cv_param$error.cv == min(cv_param$error.cv))]
  
  #training and testing model
  model = randomForest(x = sqrt(training[,-which(names(training)%in%c(response))]),
                       y = training[,response],
                       xtest = sqrt(validation[,-which(names(validation)%in%c(response))]),
                       ytest = validation[,response],
                       ntree = number_of_trees,
                       mtry = results$optim_num_params,
                       importance = TRUE)
  
  #variable importance
  df = data.frame(names = rownames(model$importance),
                  per_inc_mse = model$importance[,'%IncMSE']/model$importanceSD)
  results$dfvi = df[df$per_inc_mse >= 1,]
  
  #plotting variable importance
  results$plt = {
  ggplot(data = results$dfvi,aes(x = reorder(names,per_inc_mse),y = per_inc_mse))+
    geom_point(aes(size = per_inc_mse))+
    coord_flip()+
    labs(x = '',y = '% Increase MSE')+
    theme_bw()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.position = '')
  }
  
  return(results)
  
}

csciModel = rfvi(data = dplyr::select(df_csci,-c(Station,Location,Year,StationType,categoryCSCI)),response = 'CSCI')
csciModel$plt
csciModel$optim_num_params

h2oModel = rfvi(data = dplyr::select(df_h2o,-c(Station,Location,Year,StationType)),response = 'H2O')
h2oModel$plt
h2oModel$optim_num_params

### Non-metric Multidimensional Scaling -----

nmds <- function(data,distance = 'euclidean',k = 2){
  
  #initializing results
  results = list()
  
  d = dplyr::select(data,-c(Station,Location,Year,StationType,CSCI,H2O,categoryCSCI))
  
  #running the multidimensional scaling
  fit = metaMDS(comm = d,distance = distance,k = k)
  
  results$stress = fit$stress
  
  #data frames with results
  results$df_nmds = data.frame(Location = data$Location,
                       Category = data$categoryCSCI,
                       c1 = fit$points[,1],
                       c2 = fit$points[,2])
  results$df_nmds$Category = factor(results$df_nmds$Category,
                            levels=c('Very Likely Altered','Likely Altered','Possibly Altered','Likely Intact'),
                            ordered = T)
  results$df_nmds_species = as.data.frame(fit$species)
  
  #plot coloring
  cols = c('red','orange','yellow','blue')
  names(cols) = levels(results$df_nmds$Category)
  
  #plot mds solution
  results$plt = {
  ggplot()+ 
    geom_point(data = results$df_nmds,aes(x = c1,y = c2,shape = Location,colour = Category),size = 5)+
    geom_text_repel(data = results$df_nmds_species,aes(x = MDS1,y = MDS2,label = rownames(results$df_nmds_species)),size = 5)+
    scale_colour_manual(name = 'Category',values = cols)+
    labs(title = paste0('NMDS Ordination On PHAB Data (Stress=',round(results$stress,4),')'),x = 'MDS Coordinate 1',y = 'MDS Coordinate 2')+
    theme_bw()+
    theme(axis.text=element_text(size=14,color='black'),
          axis.title=element_text(size=14),
          legend.title=element_blank(),
          legend.text=element_text(size=14),
          legend.key=element_blank(),
          legend.position = 'top')
  }
  
  return(results)
  
}

nmdsModel = nmds(df)
nmdsModel$plt
nmdsModel$stress
