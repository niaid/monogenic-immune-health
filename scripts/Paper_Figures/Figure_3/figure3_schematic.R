library(ggplot2)

FIG.OUT.PATH <- snakemake@output[[1]]

x <- rnorm(500, 0, 3)
y <- x + rnorm(500, 0, 1)

mat <- cbind(x,y)
pca <- prcomp(mat)
pcs <- pca$x

s <- svd(mat)
U <- s$u
V <- s$v
D <- diag(s$d) ##turn it into a matrix

DV <- sqrt(D) %*% t(V)

dat = data.frame(x=x, y = y)

pdf(FIG.OUT.PATH, height = 4, width =4)
ggplot(dat, aes(x = x, y = y)) + 
  geom_point(color = "steelblue") + 
  theme_classic() + 
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  geom_segment(aes(x=0, xend=DV[1,1], y=0, yend=DV[1,2]), size = 2, color = "darkblue",
               arrow = arrow(length = unit(0.5, "cm"))) + 
  geom_segment(aes(x=0, xend=DV[2,1], y=0, yend=DV[2,2]), size = 2, color = "blue",
               arrow = arrow(length = unit(0.5, "cm")))

ggplot(dat, aes(x = x, y = y)) + 
  geom_point(color = "red") + 
  theme_classic() + 
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  geom_segment(aes(x=0, xend=DV[1,1], y=0, yend=DV[1,2]), size = 2, color = "darkred",
               arrow = arrow(length = unit(0.5, "cm"))) + 
  geom_segment(aes(x=0, xend=DV[2,1], y=0, yend=DV[2,2]), size = 2, color = "salmon",
               arrow = arrow(length = unit(0.5, "cm")))
dev.off()

