library(readxl)
Wine_example_compare_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Wine_example_compare_mean.xlsx",col_names = FALSE)
#Wine_example_compare_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Wine_example_compare_std.xlsx",col_names = FALSE)


par(family = 'Kai')  

# mean
## TPL
a1=as.numeric(Wine_example_compare_mean[1,c(1,2)])
a2=as.numeric(Wine_example_compare_mean[2,c(1,2)])
## SPL
#a1=as.numeric(Wine_example_compare_mean[3,c(1,2)])
#a2=as.numeric(Wine_example_compare_mean[4,c(1,2)])

# std
## TPL
#a1=as.numeric(Wine_example_compare_std[1,c(1,2)])
#a2=as.numeric(Wine_example_compare_std[2,c(1,2)])
## SPL
#a1=as.numeric(Wine_example_compare_std[3,c(1,2)])
#a2=as.numeric(Wine_example_compare_std[4,c(1,2)])

data=as.matrix(rbind(a1,a2),2,2)
barplot(height = data,  # 绘图数据（矩阵）
        names.arg = c('Correlation','BIC'),  # 柱子名称
        family = 'Kai',  # 中文字体
        col = c('red', 'blue'),  # 填充颜色
        border = '#ffffff',   # 轮廓颜色
        #        xlab = '频数',  # X轴名称
        ylab = 'Relative Standard Deviation of Surrogate Empirical Ranking Risk',  # Y轴名称
        #        main = '性别分布条形图',  # 主标题
        horiz = FALSE,  # 是否为水平放置
        ylim = c(0, 3), # Y轴取值范围
        legend.text = c('Red Wine', 'White Wine'),  # 图例文本
        beside = TRUE  # 是否平行排列
)

