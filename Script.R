library(caret)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
library(ggplot2)
library(glmnet)
library(jsonlite)
library(parallel)
library(pROC)
library(readxl)
library(reshape2)
library(rpart)
library(rpart.plot)
library(writexl)
library(dplyr)

set.seed(281101)

d <- as.data.table(read_excel("processed_data_sampled.xlsx"))

iqr_d <- d[, .(iqr=IQR(tpm), mean_tpm=mean(tpm)), .(gene_id, gene_name)]
d <- d[gene_id %in% iqr_d[iqr>0.5 & mean_tpm > 5, ]$gene_id, ]

md <- unique(d[, -c("filename", "project_id", "gene_id", "gene_name", "gene_type", "tpm")])
setkey(md, "sample_type")

iqr_d <- acast(d, sample_id ~ gene_id, value.var = "tpm")
fit <- prcomp(iqr_d)
iqr_pca <- list(
  components=cbind(md[match(rownames(fit$x), md$sample_id), ], as.data.table(fit$x)),
  variance = data.frame(pc=1:length(fit$sdev), variance=(fit$sdev^2)/sum(fit$sdev^2)),
  loadings=cbind(parameter=rownames(fit$rotation), data.table(fit$rotation))
)

pdf("pca_plot.pdf", width=8, height=6)
for (cname in c("sample_type")) {
  p <- ggplot(iqr_pca$components, aes_string(x="PC1", y="PC2", col=cname)) + geom_point() + theme_bw() + ggtitle(sprintf("PCA plot colored by %s", cname)) + xlab(sprintf("PC1 (%0.3f%%)", iqr_pca$variance[1, ]$variance*100)) + ylab(sprintf("PC2 (%0.3f%%)", iqr_pca$variance[2, ]$variance*100))
  print(p)
}
dev.off()

iqr_d <- t(scale(iqr_d[, -1]))
pdf("heatmap.pdf", width=25, height=25)
ha <- HeatmapAnnotation(df=md[colnames(iqr_d), .(sample_type)])
p <- Heatmap(iqr_d, row_names_gp=gpar(fontsize=6), column_names_gp=gpar(fontsize=6), clustering_distance_columns="spearman", clustering_distance_rows="euclidean", column_dend_side="top", column_names_side="top", column_dend_height=unit(30, "mm"), row_dend_width=unit(30, "mm"), row_names_side="left", col=colorRamp2(c(min(iqr_d), 0, max(iqr_d)), c("blue", "white", "red")), show_row_names=F, top_annotation=ha)
print(p)
dev.off()

iqr_d <- dcast(d, sample_id + sample_type ~ gene_id, value.var="tpm")
rownames(iqr_d) <- iqr_d$sample_id
iqr_d$sample_id <- NULL
train_size <- floor(0.8*nrow(iqr_d))
training <- sample(rownames(iqr_d), train_size)
testing <- setdiff(rownames(iqr_d), training)
fit <- rpart(sample_type ~ ., data=iqr_d[training, ])
pdf("tree.pdf", width=8, height=6)
rpart.plot(fit, extra=2)
dev.off()

train_index<-sample(seq_len(nrow(iqr_d)), size = 0.8*nrow(iqr_d))
train_data<-iqr_d[train_index, ]
test_data<-iqr_d[-train_index, ]
tree_model<-rpart(sample_type ~., data = train_data, method = "class")
predictions<-predict(tree_model, test_data, type="class")
conf_matrix<-table(Predicted=predictions, Actual=test_data$sample_type)
tp<-conf_matrix[2,2]
tn<-conf_matrix[1,1]
fp<-conf_matrix[2,1]
fn<-conf_matrix[1,2]
results_df<-data.frame(True_Positive=tp, True_Negative=tn, False_Positive=fp, False_Negative=fn)
write_xlsx(results_df, "Decision_tree_results.xlsx")
