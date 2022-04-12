#!/usr/bin/env Rscript
library(rtracklayer)
library(org.Hs.eg.db)
library(data.table)
args = commandArgs(trailing=TRUE)
fname = args[1]
oname = args[2]


gff = rtracklayer::readGFF(fname)
print(levels(gff$type))
gff = gff[gff$type == "gene",]
if (("gene_name" %in% colnames(gff)) & ("gene_id" %in% colnames(gff))) {
    df = as.data.frame(gff[c("seqid", "start", "end", "gene_name", "gene_id", "strand")])
} else if (("gene" %in% colnames(gff))) {
    df = as.data.frame(gff[c("seqid", "start", "end", "gene", "gene", "strand")])
} else {
    print("Error")
}
colnames(df) = c("seqid", "start", "end", "gene_name", "gene_id", "strand")
print(str(df))
df = df[!duplicated(df),]
is_human = "APOE" %in% df$gene_name
if (mean(df$gene_name == df$gene_id) == 1) {
    if (is_human) {
        gid = mapIds(org.Hs.eg.db, keys=df$gene_name, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
    } else {
        gid = mapIds(org.Mm.eg.db, keys=df$gene_name, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
    }
}

df$gene_name = gsub("[.][0-9]+$", "", df$gene_name)
df$gene_id = gsub("[.][0-9]+$", "", df$gene_id)
if (is_human) {
    gid = mapIds(org.Hs.eg.db, keys=df$gene_name, column="ALIAS", keytype="SYMBOL", multiVals="CharacterList")@listData
} else {
    gid = mapIds(org.Mm.eg.db, keys=df$gene_name, column="ALIAS", keytype="SYMBOL", multiVals="CharacterList")@listData
}
df1 = dplyr::bind_rows(lapply(gid, function(x) { data.frame(gene_name_2=x) }), .id="gene_name")
df2 = data.frame(gene_name=names(gid), gene_name_2=names(gid))
df12 = rbind(df1, df2)
df12 = df12[!is.na(df12$gene_name_2),]
df12 = df12[!duplicated(df12),]
df.final = merge(df, df12)[c("seqid", "start", "end", "gene_name_2", "gene_id", "strand")]
data.table::fwrite(df.final, file=oname, sep="\t", col.names=FALSE)
