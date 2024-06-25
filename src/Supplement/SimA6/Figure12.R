library(readxl)
Wine_example_compare_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Wine_example_compare_mean.xlsx",col_names = FALSE)
#Wine_example_compare_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Wine_example_compare_std.xlsx",col_names = FALSE)


par(family = 'Kai')  

# mean
## Targeted PL
a1=as.numeric(Wine_example_compare_mean[1,c(1,2)])
a2=as.numeric(Wine_example_compare_mean[2,c(1,2)])
## Surrogate PL
#a1=as.numeric(Wine_example_compare_mean[3,c(1,2)])
#a2=as.numeric(Wine_example_compare_mean[4,c(1,2)])

# standard deviation
## Targeted PL
#a1=as.numeric(Wine_example_compare_std[1,c(1,2)])
#a2=as.numeric(Wine_example_compare_std[2,c(1,2)])
## Surrogate PL
#a1=as.numeric(Wine_example_compare_std[3,c(1,2)])
#a2=as.numeric(Wine_example_compare_std[4,c(1,2)])

data=as.matrix(rbind(a1,a2),2,2)
barplot(height = data, 
        names.arg = c('Correlation','BIC'), 
        family = 'Kai', 
        col = c('red', 'blue'), 
        border = '#ffffff',  
        # mean
        ## Targeted PL
        ylab = 'Relative Targeted Empirical Ranking Risk',          
        ylim = c(0, 2), 
        ## Surrogate PL
#        ylab = 'Relative Surrogate Empirical Ranking Risk',  
#        ylim = c(0, 2.5),        
        
        # standard deviation
        ## Targeted PL
#        ylab = 'Relative Standard Deviation of Targeted Empirical Ranking Risk',  
#        ylim = c(0, 2),
        ## Surrogate PL
#        ylab = 'Relative Standard Deviation of Surrogate Empirical Ranking Risk',  
#        ylim = c(0, 2), 
        horiz = FALSE,  
        legend.text = c('Red Wine', 'White Wine'),  
        beside = TRUE  
)

