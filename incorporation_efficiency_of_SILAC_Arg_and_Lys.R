#!/usr/bin/R

### R script to calculate the incorporation efficiency of isotope labelled amino acids in the proteins
### © Mani Mudaliar, 2015.
### This script is licensed under the BSD 3-Clause License
### Contact: Manikhandan.Mudaliar@glasgow.ac.uk
### Glasgow Polyomics, University of Glasgow, Glasgow G61 1QH UK
### If you are using this script, please cite https://github.com/Bioinfo1/assorted_R_scripts/blob/master/incorporation_efficiency_of_SILAC_Arg_and_Lys.R in your work (e.g. References section in your publication or thesis)
### There has been an increasing number of requests from postdocs and PhD students for analysing SILAC incorporation efficiency, and I hope this script will help them to perform the analysis themselves. Consequently, I have  kept the analysis script simple and tried to explain the background in clear terms.    

### Usage on a Linux/Mac terminal: ### R --no-save < incorporation_efficiency_of_SILAC_Arg_and_Lys.R > incorporation_efficiency_of_SILAC_Arg_and_Lys.log 2>&1
### (OR) interactively use on a Windows R Console

### Version 0.1 ### 14/06/2015 ###
######### Change Log ##########
### Version 0.1 ###

### References ###
### 1. A computational approach to correct arginine-to-proline conversion in quantitative proteomics, Nat Methods. 2009 Mar;6(3):184-5; doi:10.1038/nmeth0309-184
### 2. Gabriele Stohr and Andreas Tebbe, Book Chapter "8. Quantitative LC-MS of proteins" in "Protein and Peptide Analysis by LC-MS Experimental Strategies"; ISBN: 978-1-84973-182-9
### 3. Cox, J. and Mann, M. (2008) MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. Nat Biotechnol 26, 1367-72; doi:10.1038/nbt.1511 [MaxQuant (version 1.5.2.8)]

#########################################################################################################
####################################    General instructions      #######################################
#########################################################################################################
### 1. Because trypsin cleaves C-terminally to arginine (R) and lysine (K), isotope labelled arginine (Arg) and lysine (Lys) are used in the SILAC experiment as it guarantees that almost all the generated peptides in the experiment contain their isotopic counterparts.
### 2. Generally, heavy labelled arginine and lysine (R10[13C6, 15N4 arginine] and K8[13C6, 15N2 lysine]), medium labelled arginine and lysine (R6[13C6 arginine], K4[2H4 lysine]) and their natural light analogues (R0[12C, 14N, 1H  arginine], K0[12C, 14N, 1H  lysine]) are used in SILAC experiments.
### 3. Check the incorporation efficiency of the labelled amino acids into the proteins before exposing the labelled cells to treatment or other experimental conditions (before starting the SILAC experiment proper). 
### 4. The incorporation efficiency is calculated by comparing the heavy (or medium) to the light form of the isotopic label in the proteins harvested from the cells grown in isotope labelled medium. Cells grown in the normal (light) medium are not used. Do not mix the proteins from light labelled (normal) cells and the heavy (or medium) labelled cells as in the regular SILAC experiment. 
### 5. Harvest proteins from about one million cells grown in the isotope labelled medium. 
### 6. Use a sample of about one microgram digested proteins in a run on a 60-min gradient in a HPLC-MS (Dionex - Orbitrap Elite) system. 
### 7. To calculate the incorporation efficiency, analyse about one thousand (at least a few hundred) peptides from the LC-MS run. 
### 8. Check the conversion of Arginine (R) to Proline (P): The supplemented heavy Arg (R) might be metabolically converted into Pro (P) in certain cell types. This results in additional forms of isotopic peptides. Isotopic Pro may form up to 30–40% of Pro-containing peptides (doi:10.1038/nmeth0309-184), and hence it is very important to check the conversion of Arg to Pro. 
### 9. For best results, minimal label incorporation of 95% should be achieved. Remember the labelled amino acids commercially available guarantees only 98% purity! 
#########################################################################################################
######################   Labelling incorporation efficiency data analysis     ###########################
#########################################################################################################
### 1. Visually inspect the chromatogram and check whether the unlabelled light peptides almost disappeared. In principle, the SILAC isotopic labelled peptide pairs have the same physio-chemical properties. Generally, 13C and 15N stable isotope labelled arginine (R) and lysine (K) peptides co-elute in the reversed-phase HPLC. However, 2H (deuterium) based labelled amino acids peptides (e.g. Lys4) elute slightly earlier due to slightly more hydrophilic character of deuterated peptides compared to their non-deuterated counterparts.
### 2. Using an appropriate software (e.g. MaxQuant) identify the peptide sequences from the MS/MS data and calculate the peptide ratios for all analyzed isotopic peptide pairs.
### 3. Exclude missed cleavages to avoid peptides containing both arginine and lysine [MaxQuant -> Group-specific parameters -> General -> Max. missed cleavages = 0]
### 4. For calculating labelling efficiencies use the non-normalized ratios. Remember that you are analysing heavy/medium labelled samples only for incorporation tests, and the normalization acts on the assumption that most SILAC ratios are around 1 after combining light and heavy samples.
### 5. From the analysis output, extract the calculated SILAC ratios of the quantified peptides.
### 6. To get separate incorporation efficiencies for arginine and lysine, split the arginine and lysine containing peptides into two peptide populations.
### 7. Compute kernel density estimates for arginine and lysine peptides

################################################################################################################################################
####  Instructions to analyse MS/MS raw data using MaxQuant (version 1.5.2.8) for labelling incorporation efficiency ###########################
################################################################################################################################################
### 1. Load raw files: [MaxQuant -> Raw files -> Load -> select MS/MS raw file/s -> Open]
### 2. Set No Fractions: [MaxQuant -> Raw files = No fractions]
### 3. Set Type as Standard : [MaxQuant -> Group-specific parameters -> General -> Type = Standard]
### 4. Set SILAC multiplicity: [MaxQuant -> Group-specific parameters -> General -> Multiplicity = 2]
### 5. Set the correct SILAC labels used in the experiment: [MaxQuant -> Group-specific parameters -> General -> Labels -> Light labels -> No selection;  Heavy labels = Select appropriate labels e.g. Arg6 and Lys4]
### 6. Set Variable modifications : [MaxQuant -> Group-specific parameters -> General -> Variable modifications = Acetyl (Protein N-term) and Oxidation (M)]
### 7. Set enzyme as Trypsin: [Trypsin/P] 
### 8. Set Max. missed cleavages to 0: [MaxQuant -> Group-specific parameters -> General -> Max. missed cleavages = 0]
### 9. Set Main search peptide tolerance to 6: [MaxQuant -> Group-specific parameters -> Instrument -> Main search peptide tolerance = 6 ppm]
### 10. Set the correct fasta file: [MaxQuant -> Global parameters -> General -> Fasta files -> Add file -> Select the correct fasta file and ensure the Andromeda search engine has been configured to use this fasta file]
### 11. Set min. ratio count to 1: [MaxQuant -> Global parameters -> Protein quantification -> Min.ratio count = 1]
### 12. Find the MaxQuant output file 'peptides.txt' under 'combined\txt' folder. This file will be used in the analysis to calculate the incorporation efficiency.

#############################################################################################
####  Labelling incorporation efficiency analysis from MaxQuant (version 1.5.2.8) output ####
#############################################################################################

## 1. Define the path and other parameter specifications (Modify the paths and analysis name in this section to suit your analysis) 
print(paste("SILAC Labelling incorporation efficiency analysis started at: ", date()))  # Log date and time (Date and time when the job was started)
defWD <- "/Volumes/share1/mani/Projects/Proteomics/Silac/Analysis"; #### Specify the path to the working directory ### (In Windows, please use '\' instead of '/')
filePath <- "/Volumes/share1/mani/Projects/Proteomics/Silac/Analysis/combined/txt/peptides.txt"; #### Specify the path to the 'peptides.txt' file. The peptides.txt file contains the identified peptides, R count, K count, P count and Ratio H/L. 
analysisName <- "Silac_20150614"; #### Specify a name for your analysis 
#########################################################################################################

## 2. Read data
setwd(defWD);
getwd();     ### Check the working directory name for log 
peptide_table <- read.table(file = filePath, header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE, comment.char = "") ### Use peptides.txt from /combined/txt
### Inspect the input data
nrow(peptide_table);  ### check the number of rows present in peptides.txt file
str (peptide_table);  ### display the structure of the R object created from the peptides.txt file
names(peptide_table); ### print the column names in the peptides.txt file

## 3. Filter exclude contaminants, false (reverse) hits and peptides without Ratio H/L.
contaminants_removed <- peptide_table[peptide_table$Potential.contaminant != "+", ]; ### Filter exclude possible contaminants in the peptides.txt 
nrow(contaminants_removed);
reverse_match_removed <- contaminants_removed[contaminants_removed$Reverse != "+", ]; ### Filter exclude peptides identified from reverse match for FDR calculation in the peptides.txt file
nrow(reverse_match_removed);
quant_peptides <- reverse_match_removed[is.na(reverse_match_removed$Ratio.H.L) == FALSE, ]; ### Filter exclude peptides that do not have 'Ratio H/L' in the peptides.txt file
nrow(quant_peptides);

## 4. Subset peptides with arginine residues but without lysine residues
arginine_pep <- quant_peptides[quant_peptides$R.Count > 0, ]; ### Filter include peptides that have at least one arginine residue
nrow(arginine_pep);
arginine_pep <- arginine_pep[arginine_pep$K.Count == 0, ]; ### Filter exclude lysine containing peptides in the peptides with arginine
nrow(arginine_pep);

## 5. Subset peptides with lysine residues but without arginine residues
lysine_pep <- quant_peptides[quant_peptides$K.Count > 0, ]; ### Filter include peptides that have at least one lysine residue
nrow(lysine_pep);
lysine_pep <- lysine_pep[lysine_pep$R.Count == 0, ]; ### Filter exclude arginine containing peptides in the peptides with lysine
nrow(lysine_pep);

## 6. Compute kernel density estimates for arginine and lysine peptides
df <- density(1 - 1 / (quant_peptides$Ratio.H.L + 1));
eff <- df$x[which.max(df$y)];

df_arginine <- density(1 - 1 / (arginine_pep$Ratio.H.L + 1));
eff_arginine <- df_arginine$x[which.max(df_arginine$y)];

df_lysine <- density(1 - 1 / (lysine_pep$Ratio.H.L + 1));
eff_lysine <- df_lysine$x[which.max(df_lysine$y)];

all_pep <- paste("All peptides: ", round(eff * 100, 2), "%");
arg_pep <- paste("Isotopic Arginine : ", round(eff_arginine * 100, 2), "%")
lys_pep <- paste("Isotopic Lysine : ", round(eff_lysine * 100, 2), "%");

## 7. Plot incorporation efficiency
xmin <- min(c(min(df$x), min(df_arginine$x), min(df_arginine$x))); ### Calculate plot axis limits from density estimates
xmax <- max(c(max(df$x), max(df_arginine$x), max(df_arginine$x)));
ymax <- max(c(max(df$y), max(df_arginine$y), max(df_arginine$y)));

postscript(paste(analysisName, "_silac_incorporation_efficiency.eps", sep = ""), width = 10, height = 7, paper="a4", onefile = FALSE, family = "Helvetica", title = "", bg = "transparent", fg = "black", horizontal = FALSE, pointsize = 12, pagecentre = TRUE, print.it = FALSE, colormodel = "cmyk"); ### Save the plot in .eps format. Generally, all journals accept this file format, and this vector graphics is readable in Adobe Illustrator and Inkscape 

plot(density(1 - 1 / (quant_peptides$Ratio.H.L + 1)), xlim = c((xmin - 0.1), (xmax + 0.1)), col = "black", ylab = expression("Density  Function"), ylim = c(0, (ymax + 2)), main = "SILAC label incorporation", lwd = 2);
lines(density(1 - 1 / (arginine_pep$Ratio.H.L + 1)), lwd = 2, col = "red");
lines(density(1 - 1 / (lysine_pep$Ratio.H.L + 1)), lwd = 2, col = "green");
legend(x = "topleft", c(all_pep, arg_pep, lys_pep), col = c("black", "red", "green"), lwd = 2);
dev.off();
###########################################################################################################################
save(list = ls(all = TRUE), file= paste(analysisName, "_final.RData", sep = "")); ### Save R data for future use or for log purpose
print(paste("Analysis completed at: ", date()));
q()
