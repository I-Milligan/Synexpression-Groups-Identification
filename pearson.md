## Pearson Correlation

```python
### Create Pearson Correlation 
rawdataT = rawdata.T
pairwise = rawdataT.corr()
pairwise.columns = gene_mapping['gene'].to_list()
pairwise.index = gene_mapping['gene'].to_list()

# Specify Gene to look at
geneID = "etv2"
#Testing on ENSDARG00000052402
temp = pairwise[[geneID]]
temp = temp.sort_values(by=geneID,ascending=False)
print(temp.head(n=20))
```
