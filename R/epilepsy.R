setwd("C:\\Users\\sguo\\Arrowhead Pharmaceuticals Inc\\Discovery Biology - Translationalgenetics\\Resource\\UkbUnivariateDiseaseResults-Olink-2024")

# Load the dplyr library
library(dplyr)
library(janitor)
library(ggplot2)
library(reshape2)
library(pheatmap)


data<-read.csv("UkbUnivariateDiseaseAssociations.csv")
data$pval[data$pval==0] <- 1.112537e-308
data$pval[data$pval<1.112537e-300] <- 1.112537e-300
input<-data %>% mutate(strength=-log(pval,10)*log2(hazard_ratio)) %>% select(phenotype,Assay,strength) %>% pivot_wider(names_from = phenotype,values_from =strength)

data1 <- input %>%
  clean_names()%>%
  mutate(
    cns = rowMeans(select(., epilepsy_recurrent_seizures_convulsions, other_peripheral_nerve_disorders), na.rm = TRUE),
    remain = rowMeans(select(., -c(epilepsy_recurrent_seizures_convulsions, other_peripheral_nerve_disorders, assay)), na.rm = TRUE),
    rate=cns/remain
  ) %>% select(assay,epilepsy_recurrent_seizures_convulsions,other_peripheral_nerve_disorders,cns,remain,rate) %>% 
  filter(cns>4) %>% 
  arrange(desc(rate)) 

data1

data2 <- input %>%
  clean_names()%>%
  mutate(
    cns = rowMeans(select(., epilepsy_recurrent_seizures_convulsions, other_peripheral_nerve_disorders), na.rm = TRUE),
    remain = rowMeans(select(., -c(epilepsy_recurrent_seizures_convulsions, other_peripheral_nerve_disorders, assay)), na.rm = TRUE),
    rate=cns/remain
  ) %>% select(assay,epilepsy_recurrent_seizures_convulsions,other_peripheral_nerve_disorders,cns,remain,rate) %>% 
  filter(cns>0) %>% 
  arrange(desc(epilepsy_recurrent_seizures_convulsions)) 

data2

write.xlsx(data2,file="UnivariateDiseaseAssociations_Christy_rank_by_epilepsy_strength.xlsx")


tmp <- t(input[match(unique(c(data1$assay[1:20],data2$assay[1:20])), input$Assay), ])
gene<-tmp[1,]
disease<-rownames(tmp)[-1]
tmp[1:5,1:5]
head(gene)
head(disease)
tmp <- tmp[-1, ]  # Remove the header row (if it's included in the matrix)
tmp <- apply(tmp, 2, function(x) as.numeric(trimws(x)))  # Convert to numeric and remove spaces
colnames(tmp)<-gene
rownames(tmp)<-disease

tmp[is.na(tmp)] <- 0

# Set rownames for the annotation to match the heatmap colnames
# Create a data frame for row annotations
row_annotations <- data.frame(
  Highlight = ifelse(
    rownames(tmp) %in% c("Epilepsy, recurrent seizures, convulsions", 
                         "Other peripheral nerve disorders"),
    "Highlighted",  # Label rows to be highlighted
    "Normal"         # Label other rows
  )
)

# Set row names for the annotation to match the heatmap row names
rownames(row_annotations) <- rownames(tmp)
# Define colors for the annotations
annotation_colors <- list(
  Highlight = c(Highlighted = "red", Normal = "grey")
)
# Generate the heatmap with row annotations
pheatmap(tmp, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         na_col = "white",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Clustered Heatmap of Assays vs Disease Conditions",
         display_numbers = FALSE,
         fontsize_row = 8,
         angle_col = 45,
         annotation_row = row_annotations,       # Add row annotations
         annotation_colors = annotation_colors)  # Define annotation colors

