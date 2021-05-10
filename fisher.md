## Fisher Test and FDR Verification

```R
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
```

**Results**
```
	  chromosome	              pvalue	        odds	        qans
1	    24	                    0.665309024   	1.326821209	  1
2	    1	                      0.18216189	    0	            1
3	    9	                      0.166067253   	1.944328675	  1
4	    17	                    0.142086337   	2.073582139	  1
5	    19	                    0.723682135   	0.472864582	  1
6	    8	                      0.733307957	    0.405067621	  1
7	    11	                    0.681985525	    1.216670178	  1
8	    21	                    0.722163431	    1.02745585	  1
9	    7	                      0.746722039	    1.148257961	  1
10	  6	                      0.300585339	    1.743573449	  1
11	  15	                    0.416254414	    0	            1
12	  12	                    1	              0.540682772	  1
13	  16	                    0.508753766	    1.29686217	  1
14	  5	                      0.37411421	    0.302210673	  1
15	  20	                    1	              0.783912415	  1
16	  10	                    0.707241933	    1.089263938	  1
17	  4	                      0.699876083   	1.122958761	  1
18	  25	                    0.408062313	    0	            1
19	  3	                      1	              0.724196497	  1
20	  18	                    0.681985525	    1.216670178	  1
 21	  22	                    0.677848404	    1.241454237	  1
22	  23	                    0.482458397	    1.392807437	  1
23	  2	                      0.739926493	    1.202008304	  1
24	  13	                    0.728904155	    0.424801564	  1
25	  14	                    0.009157349	    3.643292732	  0.375451308
26	  CHR_ALT_CTG9_2_10	      1	              0	            1
27	  CHR_ALT_CTG22_1_16	    1	              0	            1
28	  CHR_ALT_CTG22_1_13	    1	              0	            1
29	  CHR_ALT_CTG3_1_30	      1	              0	            1
30	  CHR_ALT_CTG22_1_6	      1	              0	            1
31	  CHR_ALT_CTG12_1_30	    1	              0	            1
32	  CHR_ALT_CTG17_1_16	    1	              0	            1
33	  CHR_ALT_CTG23_1_4	      1	              0	            1
34	  CHR_ALT_CTG6_1_12	      1	              0	            1
35	  CHR_ALT_CTG8_2_3	      1	              0	            1
36	  CHR_ALT_CTG3_1_31	      1	              0	            1
37	  CHR_ALT_CTG18_1_21	    1	              0	            1
38	  CHR_ALT_CTG11_1_14	    1	              0	            1
39	  KZ115982.1	            1	              0	            1
40	  CHR_ALT_CTG8_2_5	      1	              0	            1
```
