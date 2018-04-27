**Workflow to prepare files and data for loading into sleuth**

First, need to clean up trans_to_gene_key to make it into a three-column table with target_id, ens_gene. 

Get canary data from NCBI.

```bash
grep -v "Serinus canaria" trans_to_gene_key | cut -d ">" -f2,2 | perl -p -e 's/(^\S+)\s+.*gene:(\S+).*$/$1\t$2/' > t2g_key.clean
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 == "mRNA" || $1 == "ncRNA" || $1 == "misc_RNA") print $11, $16}' GCF_000534875.1_SCA1_feature_table.txt >> t2g_key.clean 
```

Verify t2g_key.clean is complete

```bash
branta:transcriptomes tim$ wc -l trans_to_gene_key 
  521808 trans_to_gene_key
branta:transcriptomes tim$ wc -l t2g_key.clean 
  521808 t2g_key.clean
```

Now need to make sleuth tables for each BioProject, which contain:

sample condition                          path
## 1 SRR493366  scramble ../results/SRR493366/kallisto
## 2 SRR493367  scramble ../results/SRR493367/kallisto
## 3 SRR493368  scramble ../results/SRR493368/kallisto
## 4 SRR493369   HOXA1KD ../results/SRR493369/kallisto
## 5 SRR493370   HOXA1KD ../results/SRR493370/kallisto
## 6 SRR493371   HOXA1KD ../results/SRR493371/kallisto

Make these by hand because I want to manually verify each experiment