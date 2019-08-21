library(cgdsr)
library(tidyverse)
library(ggsci)


# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
# test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for thyroid cancer 
# 243 = panCancer Atlas
# 244 = TCGA provisional 
# thyroid = getCancerStudies(mycgds)[178,1]
# use study Papillary Thyroid Carcinoma (TCGA, Cell 2014) : thca_tcga_pub 
thyroid = getCancerStudies(mycgds)[180,1]
thyroid_case_list = getCaseLists(mycgds,thyroid)[1,1]

# Get available genetic profiles
thyroid_genetic_profile = getGeneticProfiles(mycgds,thyroid)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list
thyroid_mutations = getProfileData(mycgds,
                                  c('ANGPT2'),
                                  thyroid_genetic_profile,
                                  thyroid_case_list)

# thyroid_mutations = getMutationData(mycgds,
#                                     thyroid_case_list,
#                                     thyroid_genetic_profile,
#                                     c('KRAS','HRAS','NRAS'))

# Get clinical data for the case list
thyroid_clinical = getClinicalData(mycgds,thyroid_case_list)

# combine clinical and mutation data 
thyroid <- thyroid_mutations %>%
  mutate(ID = rownames(thyroid_mutations)) %>%
  left_join(thyroid_clinical %>%
              mutate(ID = rownames(thyroid_clinical)), by="ID") %>%
  mutate(BRAFV600E_RAS = ifelse(BRAFV600E_RAS=="","None",BRAFV600E_RAS)) %>%
  mutate(mut_type = case_when(
    grepl("BRAF",MUT_TUMORPORTAL_GENE_PROTEIN_CHANGE) == TRUE ~ "BRAF",
    grepl("HRAS",MUT_TUMORPORTAL_GENE_PROTEIN_CHANGE) == TRUE ~ "RAS",
    grepl("NRAS",MUT_TUMORPORTAL_GENE_PROTEIN_CHANGE) == TRUE ~ "RAS",
    grepl("KRAS",MUT_TUMORPORTAL_GENE_PROTEIN_CHANGE) == TRUE ~ "RAS",
    TRUE ~ "non-RAS/BRAF"
  ))

# plotting 

tcga_fig <- ggplot(thyroid, aes(mut_type,ANGPT2, fill=mut_type)) + 
  scale_fill_manual(values=rep("white",3)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0, dodge.width = 0.75),
              aes(fill = mut_type, col = mut_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="Mutation type") +
  scale_color_jama() 

t.test(ANGPT2 ~ mut_type, data = subset(thyroid, mut_type!="non-RAS/BRAF"))
t.test(ANGPT2 ~ mut_type, data = subset(thyroid, mut_type!="BRAF"))
t.test(ANGPT2 ~ mut_type, data = subset(thyroid, mut_type!="RAS"))
