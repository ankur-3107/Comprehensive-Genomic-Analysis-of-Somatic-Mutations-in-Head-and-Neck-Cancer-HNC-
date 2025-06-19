
# Dataset Link -
# https://www.cbioportal.org/study/summary?id=hnc_mskcc_2016

# Head and Neck Cancer (HNC), HPV (Human Papillomavirus)

setwd("C:/Ankur_IIITD/cancer_genomics/HNC/hnc_mskcc_2016")
list.files()


### Clinical Summary and Survival Analysis ------------------------------------------------------------------------------------

library(survival)
library(survminer)
library(dplyr)

clin.patient <- read.delim("data_clinical_patient.txt", skip = 4)
clin.patient <- clin.patient %>%
  mutate(OS_MONTHS = as.numeric(OS_MONTHS),
         OS_STATUS = ifelse(grepl("DECEASED", OS_STATUS), 1, 0))

# Overall Survival (OS) Kaplan-Meier Plot
fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ HPV_STATUS, data = clin.patient)
ggsurvplot(fit, data = clin.patient, pval = TRUE,
           title = "Overall Survival by HPV Status")

# Patients with HPV_STATUS=Positive appear to have a generally better overall survival probability over time compared to those with HPV_STATUS=Negative.
# HPV_STATUS=Unknown and HPV_STATUS= groups has very limited data, making it hard to draw conclusions.
# p-value=0.057 is close to traditional significance level of of 0.05, suggests that there is a trend towards a statistically significant difference in overall survival between the groups.


### Mutation Landscape -----------------------------------------------------------------------------------------------------

library(maftools)
library(data.table)

# Calculate and plot VAF 
mutations <- fread("data_mutations.txt")
t_alt_count <- mutations$t_alt_count
t_ref_count <- mutations$t_ref_count
vaf <- t_alt_count / (t_ref_count + t_alt_count)
mutations <- cbind(mutations, vaf)
maf <- read.maf(maf = mutations)
plotVaf(maf = maf, vafCol = 'vaf')

# This plot highlights TP53 as the most frequently mutated gene in this dataset and suggests that CDKN2A and NOTCH1 mutations tend to occur at higher VAFs, potentially indicating they are more frequently clonal events compared to genes like NOTCH2 or NOTCH3.

# Data summary and oncoplot
plotmafSummary(maf)
oncoplot(maf = maf, top = 20, sortByAnnotation = TRUE)

# Missense Mutations are by far the most common type of variant, accounting for over 900 instances.
# Single Nucleotide Polymorphisms (SNPs) are overwhelmingly the most common type of variant, with over 1000 instances.
# C>T (Cytosine to Thymine) transitions are the most frequent single nucleotide variant class, with over 800 occurrences.
# Variants per sample plot shows a distribution where most samples have a relatively low number of variants. (Gives idea about total mutation burden i.e. TMB)
# The median number of missense mutations per sample is around 2, with some samples having significantly more (outliers).
# TP53 is the most frequently mutated gene, with over 60 mutations detected, affecting 46% of the samples. The different colored segments within the TP53 bar represent the different variant classifications found in that gene.
# The oncoplot strongly points to TP53, NOTCH1, KMT2D, CDKN2A, FAT1, and PIK3CA as key driver genes in our HNC dataset due to their high alteration frequencies.


# Lollipop plot to represent functional amino acid mutations
lollipopPlot(maf = maf, gene = "TP53", AACol = "MA:protein.change", showMutationRate = TRUE)

# Lollipop plot visualizes the distribution of somatic mutations within the TP53 gene. 
# TP53 protein has three main domain - 
#     1) P53_TAD (p53 Transactivation Domain): Located at the N-terminus (roughly amino acids 1-60). This domain is crucial for initiating transcription of target genes.
#     2) P53 (DNA Binding Domain): The large orange box, spanning approximately amino acids 100-300. This is the core domain responsible for binding to DNA and regulating gene expression. This domain is frequently targeted by mutations.
#     3) P53_tetramer (Oligomerization Domain): Located at the C-terminus (roughly amino acids 320-360). This domain is responsible for the formation of p53 tetramers, which is essential for its tumor suppressor activity.

# Finding driver genes
maf.drivers <- oncodrive(maf, AACol = "MA:protein.change", minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = maf.drivers, fdrCutOff = 0.5, useFraction = TRUE, labelSize = 0.5)

# Driver genes are the genes that are most likely to be oncogenic drivers by evaluating both the positional distribution (clustering) and the predicted functional consequence of their mutations across a cohort. 
# Genes appearing in the top-right (red coloured) - HRAS, CDKN2A, PIK3CA, PTEN, NFE2L2, FGFR3, PARP1 and SF3B1, are strong candidates for further investigation as key players in the cancer being studied.

# Adding and summarizing pfam domains
maf.pfam = pfamDomains(maf, AACol = 'MA:protein.change', top = 10)
maf.pfam$proteinSummary[, 1:7, with = FALSE]
maf.pfam$domainSummary[, 1:5, with = FALSE]

# Co-mutation and Mutual Exclusivity Analysis
somaticInteractions(maf, top = 20, pvalue = c(0.05, 0.01))

# Look for dark green cells with asterisks. These pairs of genes are frequently mutated together more often than expected by chance (eg. TP53-CDKN2A). They are might be involved in the same biological pathway, and mutations in both contribute to disease progression.
# Look for dark brown/orange cells with asterisks. These pairs of genes are rarely mutated together; if one is mutated, the other is often wild-type (unmutated).
# For co-occurring mutations, you might target both pathways simultaneously. While for mutually exclusive mutations, understanding which pathway is active (based on the mutated gene) can guide specific targeted therapies.

# Drug-gene interactions 
dgi = drugInteractions(maf, fontSize = 0.75)
tp53.dgi = drugInteractions(genes = "TP53", drugs = TRUE)
tp53.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

# Oncogenic signaling pathways
pws = pathways(maf, plotType = 'treemap')
plotPathways(maf, pathlist = pws)


### Pathway Enrichment Analysis of Mutated Gene -----------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load mutated genes
mut_genes <- unique(mutations$Hugo_Symbol)
entrez_ids <- bitr(mut_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG (Kyoto Encyclopedia of Genes and Genomes) Pathway Enrichment
kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa')
dotplot(kegg, showCategory = 10) + ggtitle("KEGG Pathway Enrichment") |
  barplot(kegg, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

# PI3K-Akt signaling pathway is a central and highly active pathway in your HNC data, indicating its critical role in the disease's development and progression.
# There is significant molecular overlap between your HNC samples and other common cancer types, as evidenced by the enrichment of multiple "cancer" pathways.
# EGFR pathway, particularly related to inhibitor resistance, is also highly implicated, suggesting potential therapeutic vulnerabilities or resistance mechanisms that need to be considered for HNC treatment.

# GO (Gene Ontology) Biological Process Enrichment
go <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
dotplot(go, showCategory = 10) + ggtitle("GO Biological Processes") |
  barplot(go, showCategory = 10) + ggtitle("GO Biological Process")

# HNC samples exhibit significant deregulation of protein phosphorylation and kinase activity, pointing to widespread alterations in cellular signaling networks.
# phosphatidylinositol 3-kinase/protein kinase B signal transduction (highly significant, large count) and regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction (significant, moderate count); independently confirm the high importance of the PI3K-Akt signaling pathway
# The strong enrichment of "gland development" and terms related to cell differentiation/fate suggests that developmental processes are aberrantly activated or disrupted in HNC samples, potentially contributing to tumorigenesis.



### CNA Landscape --------------------------------------------------------------------------

library(ggplot2)

cna <- fread("data_cna_hg19.seg")
# Convert 'chrom' to factor for proper ordering in plots
cna[, chrom := factor(chrom, levels = unique(chrom[order(as.numeric(chrom))]))]

# Calculate segment length
cna[, seg.length := loc.end - loc.start]

# Identify gains and losses based on seg.mean
gain_threshold <- 0.2
loss_threshold <- -0.2

cna[, CNA_type := "Neutral"]
cna[seg.mean > gain_threshold, CNA_type := "Gain"]
cna[seg.mean < loss_threshold, CNA_type := "Loss"]
cna[, CNA_type := factor(CNA_type, levels = c("Loss", "Neutral", "Gain"))]

ggplot(cna, aes(x = (loc.start + loc.end) / 2, y = seg.mean, color = CNA_type)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_segment(aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean), size = 1) +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("Loss" = "blue", "Neutral" = "grey", "Gain" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = gain_threshold, linetype = "dotted", color = "red") +
  geom_hline(yintercept = loss_threshold, linetype = "dotted", color = "blue") +
  labs(
    title = paste0("Genome-wide Copy Number Profile for ", unique(cna$ID)),
    x = "Genomic Position",
    y = "Segment Mean Log2 Ratio",
    color = "CNA Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Hide x-axis labels for genomic position
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.1, "cm"),
    strip.text.x = element_text(size = 8),
    legend.position = "bottom"
  )

# This plot clearly shows that sample P-0001467-T01-IM3 has a highly altered genome with numerous and significant copy number gains and losses across many chromosomes, which is characteristic of many types of cancer. 
# The specific regions of gain and loss can point to the activation of oncogenes or inactivation of tumor suppressor genes, which are critical in cancer development and progression.


ggplot(cna, aes(x = seg.length / 1e6, y = seg.mean, color = CNA_type)) + # Convert length to Mb
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Loss" = "blue", "Neutral" = "grey", "Gain" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Segment Length vs. Segment Mean Log2 Ratio",
    x = "Segment Length (Mb)",
    y = "Segment Mean Log2 Ratio",
    color = "CNA Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# This plot confirms the presence of a wide spectrum of copy number alterations in the sample, ranging from small, high-magnitude focal events (both gains and losses) to larger, lower-magnitude arm-level or chromosomal alterations.
# The abundance of both high-level amplifications and deep deletions, especially at shorter lengths, reinforces the idea of significant genomic instability and the likely presence of functionally important driver CNAs in this sample.



### TMB Comparison by HPV Status ------------------------------------------------------------------------

clin.sample <- read.delim("data_clinical_sample.txt", skip = 4)
tmb <- data.frame(SAMPLE_ID = clin.sample$SAMPLE_ID,
                  TMB_NONSYNONYMOUS = clin.sample$TMB_NONSYNONYMOUS,
                  HPV_STATUS = clin.patient$HPV_STATUS)

ggplot(tmb, aes(x = HPV_STATUS, y = TMB_NONSYNONYMOUS, fill = HPV_STATUS)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Tumor Mutational Burden by HPV Status", y = "TMB (log scale)")


### Circos Plot for Structural Variants (Gene Fusions) ----------------------------------------------------------------

library(circlize)
sv <- fread("data_sv.txt")
sv_data <- sv[!is.na(sv$Site1_Chromosome) & !is.na(sv$Site2_Chromosome), ]

chords <- data.frame(
  chr1 = paste0("chr", sv_data$Site1_Chromosome),
  pos1 = sv_data$Site1_Position,
  chr2 = paste0("chr", sv_data$Site2_Chromosome),
  pos2 = sv_data$Site2_Position
)

circos.clear()
circos.initializeWithIdeogram(species = "hg19")
circos.genomicLink(chords[,1:2], chords[,3:4], col = "red")

# This plot gives the information about inter- and intra- chromosomal translocations. Intrachromosomal translocations involve a segment moving to a different location on the same chromosome, while interchromosomal translocations involve a segment moving between two different chromosomes. 
# Inter-chromosomal translocations - chr1-chr12, chr12-chr15, chr11-chr22 and chr12-chr22.
# Intra-chromosomal translocations - chr3, chr4, chr6, chr7, chr9, chr12, chr17 and chrX.



### Mutational Signature Analysis -----------------------------------------------------------------

# Libraries
library(deconstructSigs)

# Input: mutation context (need VCF-like info)
# Assuming your data_mutations is properly formatted
sigs_input <- mut.to.sigs.input(
  mut.ref = mutations,
  sample.id = "Tumor_Sample_Barcode",
  chr = "Chromosome",
  pos = "Start_Position",
  ref = "Reference_Allele",
  alt = "Tumor_Seq_Allele2"
)

sigs_output <- whichSignatures(tumor.ref = sigs_input,
                               signatures.ref = signatures.cosmic,
                               sample.id = "P-0001467-T01-IM3",
                               contexts.needed = TRUE)

plotSignatures(sigs_output)

# This is the 96-triplet mutation profile.
# For specific tumor sample P-0001487-701-IM3; its observed mutational landscape is best explained by the activity of four distinct mutational processes: Signature 8 (most prominent), Signature 3, Signature 18, and Signature 7. 
# Presence of Signature 3 suggests potential homologous recombination deficiency in this tumor, which could have implications for treatment (e.g., sensitivity to PARP inhibitors).
# Presence of Signature 7 indicates prior exposure to UV radiation, suggesting this tumor might be a skin cancer (melanoma or squamous cell carcinoma) or has acquired this damage irrespective of its origin.
# Signature 8 and 18 contribute significantly, indicating other active mutational processes, whose precise biological mechanisms might require further investigation.
# "error = 0.266" value quantifies the overall discrepancy between the observed and reconstructed profiles. A smaller error indicates a better fit. (Quality of fit)


### Recurrent Fusions in Metastatic Disease --------------------------------------------------------------

sv.merge <- inner_join(sv, clin.sample %>% dplyr::select(SAMPLE_ID, SAMPLE_TYPE), by=c("Sample_Id" = "SAMPLE_ID"))

fusion_counts <- sv.merge %>%
  dplyr::mutate(fusion = paste(Site1_Hugo_Symbol, Site2_Hugo_Symbol, sep = "--")) %>%
  group_by(fusion, SAMPLE_TYPE) %>%
  tally()

ggplot(fusion_counts, aes(x=fusion, y=n, fill = SAMPLE_TYPE)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(title = "Recurrent Fusions by Sample Type")










