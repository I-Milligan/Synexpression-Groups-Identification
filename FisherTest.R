library(tidyverse)

hits <- read_tsv("hits.txt", guess_max = 100000)
world <- read_tsv("GeneWorld.txt", guess_max = 100000)

world$hits <- world$ensembleID %in% hits$ensembleID

chroms <- unique(world$chromosome)

results <- list()

for(n0 in chroms){
  print(n0)
  world$ch <- world$chromosome==n0
  ans <- fisher.test(world$hits, world$ch)
  ans1 <- tibble(chromosome=n0,pvalue=ans$p.value,odds=ans$estimate)
  results[[n0]] <- ans1
}

all_results <- bind_rows(results)

qans <- p.adjust(all_results$pvalue,method="fdr")

all_results$qans <- qans
write.table(all_results, file = "nrf1_results", quote = FALSE, sep = "\t")
