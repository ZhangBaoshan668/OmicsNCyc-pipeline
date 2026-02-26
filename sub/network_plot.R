#!/usr/bin/env Rscript

# 加载所需的包
library(igraph)
library(dplyr)
library(Hmisc)
library(optparse)
library(purrr)
## 读入OTU/ASV表格，列为样本，行为物种
# 解析命令行参数
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input data file [default %default]"),
  make_option(c("-g", "--group"), type="character", default=NULL,
              help="Group information file [default %default]"),
  make_option(c("-o", "--output"), type="character", default="output.pdf",
              help="Output PDF file [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=16,
              help="Output figure width in mm [default %default]"),
  make_option(c("-h", "--height"), type="numeric", default=11,
              help="Output figure height in mm [default %default]")
)

# 创建解析器并设置 add_help_option=FALSE
opts_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opts = parse_args(opts_parser)

# 读取数据
otu_rare <- read.delim(opts$input, header = T,row.names = 1,stringsAsFactors = F)

# 读取分组信息
metadata <- read.table(opts$group, header=TRUE)

## 定义一些颜色
col_g <- "#C1C1C1"
cols <- c("#EA9527", "#5ECC6D", "#F16E1D" ,"#DEB99B" , "#5DAFD9", "#7ED1E4", "#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")

trt_id <- unique(metadata$Group)

# 3. 验证样本匹配情况
shared_samples <- intersect(metadata$SampleID, colnames(otu_rare))
if(length(shared_samples) == 0) stop("错误：元数据与OTU表样本完全不匹配")

# 4. 改进的分组方法（精确匹配）
split_otu <- lapply(trt_id, function(group){
  # 从元数据获取该组所有样本ID
  group_samples <- metadata$SampleID[metadata$Group == group]
  
  # 找出实际存在于OTU表中的样本
  valid_samples <- intersect(group_samples, colnames(otu_rare))
  
  if(length(valid_samples) == 0) return(NULL)
  
  # 提取数据并过滤全零OTU
  group_data <- otu_rare[, valid_samples, drop = FALSE]
  group_data[rowSums(group_data) > 0, , drop = FALSE]
}) %>% 
  setNames(trt_id) %>% 
  compact()  # 自动移除空分组

# 5. 检查结果
if(length(split_otu) == 0) stop("错误：未能创建任何分组")

# 使用实际组别命名结果列表
names(split_otu) <- trt_id

# 对分组的OTU数据列表(split_otu)中的每个组进行处理
g <- lapply(split_otu,function(x){
  occor<-WGCNA::corAndPvalue(t(x)/colSums(x),method = 'spearman')
  mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(t(x)/colSums(x))[2])
  ## R value
  occor.r<-occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0
  occor.r[is.na(occor.r)]=0
  g <-  graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  # 删除自相关
  g <- simplify(g)
  # 删除孤立节点
  g <- delete_vertices(g, which(degree(g)==0) )
  return(g)
})

# 定义计算Zi和Pi的函数
calculate_zi_pi <- function(graph) {
  mods <- V(graph)$modularity
  zi <- sapply(1:vcount(graph), function(i) {
    same_mod <- which(mods == mods[i])
    degrees <- degree(graph)[same_mod]
    if(length(degrees) < 2) return(0)
    (degree(graph, i) - mean(degrees)) / sd(degrees)
  })
  zi[is.na(zi)] <- 0
  
  pi <- sapply(1:vcount(graph), function(i) {
    neighbors_mod <- mods[neighbors(graph, i)]
    if(length(neighbors_mod) == 0) return(0)
    tab <- table(neighbors_mod)
    1 - sum((tab/sum(tab))^2)
  })
  
  return(data.frame(Zi = zi, Pi = pi))
}

# 遍历所有分组绘制网络图
for(i in 1:length(g)) {
  group_name <- names(g)[i]
  current_g <- g[[i]]
  
  # 跳过空网络
  if(vcount(current_g) == 0) {
    warning(paste("分组", group_name, "没有节点，跳过..."))
    next
  }
  
  cat("正在处理分组:", group_name, "\n")
  
  # 设置边的相关性和权重
  E(current_g)$correlation <- E(current_g)$weight
  E(current_g)$weight <- abs(E(current_g)$weight)
  
  # 计算模块性
  set.seed(007)
  fc <- cluster_fast_greedy(current_g)
  V(current_g)$modularity <- membership(fc)
  
  # 计算Zi和Pi值
  zp <- calculate_zi_pi(current_g)
  V(current_g)$Zi <- zp$Zi
  V(current_g)$Pi <- zp$Pi
  
  # 节点拓扑属性
  node_attrs <- data.frame(
    Group = group_name,
    Node = V(current_g)$name,
    Degree = degree(current_g),
    Betweenness = betweenness(current_g),
    Closeness = closeness(current_g),
    Eigenvector = eigen_centrality(current_g)$vector,
    Modularity = V(current_g)$modularity,
    Zi = V(current_g)$Zi,
    Pi = V(current_g)$Pi,
    Role = ifelse(V(current_g)$Zi > 2.5 & V(current_g)$Pi > 0.62, "Network hub",
                  ifelse(V(current_g)$Zi > 2.5, "Module hub",
                         ifelse(V(current_g)$Pi > 0.62, "Connector", "Peripheral")))
  )
  
  # 边属性
  edge_list <- as_edgelist(current_g)
  edge_attrs <- data.frame(
    Group = group_name,
    Source = edge_list[,1],
    Target = edge_list[,2],
    Weight = E(current_g)$weight,
    Correlation = E(current_g)$correlation,
    Direction = ifelse(E(current_g)$correlation > 0, "Positive", "Negative")
  )
  
  # 模块颜色设置
  modu_sort <- table(V(current_g)$modularity) %>% sort(decreasing = TRUE)
  top_num <- min(10, length(modu_sort))
  modu_name <- names(modu_sort[1:top_num])
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  
  # 顶点着色
  V(current_g)$color <- ifelse(V(current_g)$modularity %in% modu_name,
                               modu_cols[match(V(current_g)$modularity, modu_name)],
                               col_g)
  V(current_g)$frame.color <- V(current_g)$color
  
  # 边着色
  E(current_g)$color <- col_g
  for(mod in modu_name) {
    col_edge <- modu_cols[mod]
    otu_same_modu <- V(current_g)$name[V(current_g)$modularity == mod]
    same_modu_edges <- which(edge_list[,1] %in% otu_same_modu & edge_list[,2] %in% otu_same_modu)
    E(current_g)$color[same_modu_edges] <- col_edge
  }
  
  # 网络布局
  set.seed(42)
  current_layout <- layout_with_fr(current_g, niter=999, grid='nogrid')
  # current_layout <- layout_nicely(current_g)
# current_layout <-  layout_with_graphopt(current_g)
# current_layout <-  layout_with_kk(current_g)
# current_layout <-  layout_with_mds(current_g)
  
  # 绘制网络图
  pdf(file.path(opts$output,paste0("Network_", group_name, ".pdf")), width=opts$width, height=opts$height, paper="special")
  par(mar=c(0,0,3,0), font.main=4, cex.main=1)
  plot(current_g, 
       layout = current_layout,
       vertex.size = 2,
      #  vertex.label = V(current_g)$name,
       vertex.label = NA,
      #  vertex.label.cex = 0.7,
      #  vertex.label.color = "black",
       edge.width = 0.5 + 2*abs(E(current_g)$correlation),
       edge.color = E(current_g)$color,
       main = paste0(group_name, " Co-occurrence Network\n",
                     "Nodes: ", vcount(current_g), 
                     " | Edges: ", ecount(current_g),
                     "\nModularity: ", round(modularity(fc), 3),
                     " | Transitivity: ", round(transitivity(current_g), 3)))

  # # 添加图例
  legend_items <- c()
  legend_cols <- c()
  if(length(modu_name) > 0) {
    legend_items <- c(legend_items, paste("Module", modu_name))
    legend_cols <- c(legend_cols, modu_cols)
  }
  legend_items <- c(legend_items, c("Network hub", "Module hub", "Connector", "Peripheral"))
  legend_cols <- c(legend_cols, c("red", "blue", "green", "gray"))

  legend("topleft",
         legend = legend_items,
         col = legend_cols,
         pch = c(rep(15, length(modu_name)), rep(16, 4)),
         bty = "n",
         pt.cex = 1.5,
         cex = 0.8)
  
  dev.off()
  
  # 保存GraphML文件
  write_graph(current_g, file.path(opts$output,paste0("Network_", group_name, ".graphml")), format = "graphml")
  
  # 保存节点和边属性
  write.csv(node_attrs, file.path(opts$output,paste0("Node_attributes_", group_name, ".csv")), row.names = FALSE)
  write.csv(edge_attrs, file.path(opts$output,paste0("Edge_attributes_", group_name, ".csv")), row.names = FALSE)
}

# 合并所有分组的拓扑属性
all_nodes <- do.call(rbind, lapply(names(g), function(name) {
  if(vcount(g[[name]]) > 0) {
    read.csv(file.path(opts$output,paste0("Node_attributes_", name, ".csv")))
  }
}))

all_edges <- do.call(rbind, lapply(names(g), function(name) {
  if(vcount(g[[name]]) > 0) {
    read.csv(file.path(opts$output,paste0("Edge_attributes_", name, ".csv")))
  }
}))

write.csv(all_nodes, file.path(opts$output,"All_nodes_topology_attributes.csv"), row.names = FALSE)
write.csv(all_edges, file.path(opts$output,"All_edges_attributes.csv"), row.names = FALSE)

# 生成网络统计摘要
# 修改网络统计摘要的计算方式，确保处理负权重问题
network_stats <- data.frame(
  Group = names(g),
  Nodes = sapply(g, vcount),
  Edges = sapply(g, ecount),
  Average_degree = sapply(g, function(x) ifelse(vcount(x)>0, mean(degree(x)), NA)),
  Modularity = sapply(g, function(x) {
    if(vcount(x) == 0) return(NA)
    # 确保权重为非负
    E(x)$weight <- abs(E(x)$weight)
    modularity(cluster_fast_greedy(x))
  }),
  Transitivity = sapply(g, function(x) ifelse(vcount(x)>0, transitivity(x), NA)),
  Diameter = sapply(g, function(x) ifelse(vcount(x)>0, diameter(x, weights = NA), NA)), # 忽略权重计算直径
  Average_path_length = sapply(g, function(x) ifelse(vcount(x)>0, mean_distance(x, weights = NA), NA)) # 忽略权重计算平均路径长度
)


# 保存统计结果
write.csv(network_stats, file.path(opts$output, "Network_summary_statistics.csv"), row.names = FALSE)
