source("code/CAI.R")

df = read.csv("output_data/PCA_output.csv")
df = merge(df, df2, by.y = "Name", by.x = "X")

par(mfrow = c(1,2))
summary(lm(data = df, PC1 ~ normalised))
plot(data = df, PC1 ~ normalised, 
     main = "R^2 = 0.1258",
     xlab = "nCAI")
abline(lm(data = df, PC1 ~ normalised), col = "red")

summary(lm(data = df, PC2 ~ normalised))
plot(data = df, PC2 ~ normalised,
     main = "R^2 = 0.6302",
     xlab = "nCAI")
abline(lm(data = df, PC2 ~ normalised), col = "red")
