library(readxl)
Wine_example_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA6_Wine_example_mean.xlsx",col_names = FALSE) # mean
#Wine_example_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA6_Wine_example_std.xlsx",col_names = FALSE) # standard deviation


par(family = 'Kai')  
# mean
## Targeted PL
a1=as.numeric(Wine_example_mean[1,1:12])
a2=as.numeric(Wine_example_mean[2,1:12])
## Surrogate PL
#a1=as.numeric(Wine_example_mean[3,1:12])
#a2=as.numeric(Wine_example_mean[4,1:12])

# standard deviation
## Targeted PL
#a1=as.numeric(Wine_example_std[1,1:12])
#a2=as.numeric(Wine_example_std[2,1:12])
## Surrogate PL
#a1=as.numeric(Wine_example_std[3,1:12])
#a2=as.numeric(Wine_example_std[4,1:12])

data=as.matrix(rbind(a1,a2),2,12)
barplot(height = data,  
        names.arg = c('Full','AIC','BIC','SAIC','SBIC','MMA','Equal','K-data','RMA-op', 'RMA-2','RMA-5','RMA-n/2'), 
        family = 'Kai',  
        col = c('red', 'blue'), 
        border = '#ffffff', 
        # mean
        ## Targeted PL
        ylab = 'Relative Targeted Empirical Ranking Risk',  
        ylim = c(0, 2),
        ## Surrogate PL
#        ylab = 'Relative Surrogate Empirical Ranking Risk',  
#        ylim = c(0, 8),

        # standard deviation
        ## Targeted PL
#        ylab = 'Relative Standard Deviation of Targeted Empirical Ranking Risk',  
#        ylim = c(0, 4),
        ## Surrogate PL
#        ylab = 'Relative Standard Deviation of Surrogate Empirical Ranking Risk',  
#        ylim = c(0, 4),        
        horiz = FALSE,  
        
        legend.text = c('Red Wine', 'White Wine'),  
        beside = TRUE 
)

