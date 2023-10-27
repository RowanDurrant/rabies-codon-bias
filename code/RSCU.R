library("readxl")
df = as.data.frame(read_excel("RSCU.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]

hosts = c("Skunk", "Big brown bat", "Free-tailed bat")

for(i in 3:ncol(df)){
  print(colnames(df)[i])
  for(j in hosts){
    print(paste( j, mean(as.numeric(df[which(df$Host == j), i]))))
    
  }
  
}

for(k in hosts){
  print(paste(k, nrow(df[df$Host == k,])))
}

paste(df$`Accession no.`[df$Host == "Skunk"], collapse = ", ")
