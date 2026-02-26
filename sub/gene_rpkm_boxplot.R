#!/usr/bin/env Rscript

# 加载所需的包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(optparse)
library(ggprism)

# 解析命令行参数
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input data file [default %default]"),
  make_option(c("-g", "--group"), type="character", default=NULL,
              help="Group information file [default %default]"),
  make_option(c("-o", "--output"), type="character", default="output.pdf",
              help="Output PDF file [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=160,
              help="Output figure width in mm [default %default]"),
  # make_option(c("-x", "--xlabAngle"), type="logical", default=45,
  #               help="X lab set in angle [default %default]"),
  make_option(c("-h", "--height"), type="numeric", default=110,
              help="Output figure height in mm [default %default]")
)

# 创建解析器并设置 add_help_option=FALSE
opts_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opts = parse_args(opts_parser)

# 读取数据
data <- read.delim(opts$input, sep="\t", header=TRUE)

# 读取分组信息
group_data <- read.delim(opts$group, sep="\t", header=TRUE)

# 合并数据
data <- merge(data, group_data, by.x="SampleID", by.y="SampleID")

# 创建箱线图
p1 <- ggplot(data, aes(x=Group, y=RPKM, fill = Group)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA) +  # 添加箱线图层，不显示离群点
  scale_fill_brewer(palette="Set1") +  # 设置填充颜色
  geom_jitter(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.2), size=0.3) +  # 添加散点图层，并调整位置
  ylab('RPKM') +  # 设置 y 轴标签
  xlab('') +  # 清除 x 轴标签
  theme_prism() + # 使用 Prism 主题
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))

# 保存图形
ggsave(opts$output, p1, width = opts$width, height = opts$height, units = "mm")