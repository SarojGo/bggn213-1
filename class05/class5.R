#' ---
#' title: "Bioinformatics Class 5"
#' author: "Barry"
#' output:
#'   html_document:
#'     code_folding: hide
#'     keep_md: yes
#' ---

# Bioinformatics Classs 5 
# Plots

x <- rnorm(1000, 0)

summary(x)

# lets see this data as a graph
boxplot(x)

# Good old histogram
hist(x)

# Section 1 from lab sheet
baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)

plot(baby, type="b", pch=15, cex=1, ylim=c(0,12), 
     lwd=3, main="Baby weight")

## Section 1B.
feat <- read.table("bggn213_05_rstats/feature_counts.txt", sep="\t",
                   header = TRUE)

par(mar=c(5, 11, 4, 2))  
barplot(feat$Count, names.arg=feat$Feature, horiz = TRUE, las=2)


# Section 2A
#read.table("bggn213_05_rstats/male_female_counts.txt", sep="\t", header=TRUE)

mfcount <- read.delim("bggn213_05_rstats/male_female_counts.txt")

barplot(mfcount$Count, col=rainbow(10))

mycols <- cm.colors(nrow(mfcount))
barplot(mfcount$Count, col=mycols)

# Section 2B
expr <- read.delim("bggn213_05_rstats/up_down_expression.txt")
plot(expr$Condition1, expr$Condition2, col=expr$State)

# How may genes went up down and around?
table(expr$State)

# Section 2C.
meth <- read.delim("bggn213_05_rstats/expression_methylation.txt")
plot(meth$promoter.meth, meth$gene.meth)

# make a gray to red color vector
expression.palette <- colorRampPalette(c("black","red"))(101) 
expression.low <- min(meth$expression) 
expression.high <- max(meth$expression) 

map.colors <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}

mycol=map.colors(meth$expression,c(expression.high, expression.low),expression.palette)

plot(meth$promoter.meth, meth$gene.meth, col=mycol, pch=15)

plot(meth$expression, meth$gene.meth, col=mycol, pch=15)

# double check color setting
library("colorspace")
chcols <- diverge_hcl(length(meth$expression))[rank(meth$expression)]
#chcols <- sequential_hcl(length(meth$expression))[rank(meth$expression)]

plot(meth$promoter.meth, meth$gene.meth, col=chcols, pch=15)
