library(docxtractr)

doc_tables <- read_docx(path = "~/Desktop/Supp_Table1.docx")
cpg_table <- docx_extract_tbl(docx = doc_tables, tbl_number = 1)
gene_cpg_list <- cpg_table$Closest.gene.name.and.ID

cpg_ids <- c()

for(gene_cpg in gene_cpg_list){
  split_gene_cpg <- unlist(strsplit(gene_cpg, " "))
  cpg_id <- split_gene_cpg[length(split_gene_cpg)]
  cpg_ids <- append(cpg_ids, cpg_id)
}

write (cpg_ids, file = "./data/DifferentiallyMethylatedProbes.txt")
