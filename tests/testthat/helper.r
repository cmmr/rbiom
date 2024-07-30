
rare50 <- rarefy(hmp50)
# min50  <- as_rbiom(list(counts = hmp50$counts))

hmp5  <- hmp50[1:5]
rare5 <- rarefy(hmp5)
min5  <- as_rbiom(list(counts = hmp5$counts))

