# The evaluation of tools used to predict the impact of missense mutations is hindered by two types of circularity

This is the code to reproduce all findings from the following paper:

Grimm, Dominik G., et al. "The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity." Human mutation 36.5 (2015): 513-523.
http://onlinelibrary.wiley.com/doi/10.1002/humu.22768/full

When using the code/data please cite out paper!



********************************
Scripts to reproduce all figures
********************************

Code by: Dominik Gerhard Grimm
Year: 2015
Group: Machine Learning and Computational Biology Research Group
Insitute: Max Planck Institute for Intelligent Systems and Max Planck Institute for Developmental Biology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

********************************
Required Python Packages
********************************

numpy (>= 1.8.1)

scipy (>= 0.10.0)

matplotlib (>= 1.3.1)

sklearn (>= 0.15.2)


********************************
Install Packages
********************************

Use the tool easy_install or pip to install the missing packages.


********************************
Run scripts
********************************

Running the scripts to reproduce all figures and table.

Go to the terminal/bash, move to the DataS1 folder and type in the following command:

>python start.py

This will start the reproduction of all results



********************************
Archive Folder
********************************

> Output: Contains all Main and Supplementary Figures and Tables after executing the start.py script.


> ToolScores: Data directory! Contains all variants and retrieved tool-scores, features, predicted labels and true labels for all investigated datasets


> Similiarities: Contains percentage of variants that can be found in similar proteins in HumVar/ExoVar for different protein similarity thresholds


> Scripts: Contains all nevessary Scripts to rerun all experiments and for plotting all figures. These scripts are executed when start.py is called


********************************
Header Explanation for Data files in ToolScores
********************************

The first row is the header!
Each row contains one variant and the tool-scores and predicted labels for different tools

Here is a description of the different coloumns:
------------------------------------------------

Column 1: True Label - The true label of this variant

Column 2: #RS-ID - if available the rs identifier for this variant 

Column 3: CHR - chromosome at which the variant is located

Column 4: Nuc-Pos - nucleotide position of the variant 

Column 5: REF-Nuc - the reference nucleotide

Column 6: ALT-Nuc - the alternative nucleotide

Column 7: MAF - minor allele frequence if available 

Column 8: Ensembl-Gene-ID - ensembl gene id for this variant

Column 9: Ensembl-Protein-ID - ensembl protein if for this variant 

Column 10: Ensembl-Transcript-ID - ensemble transcript id for this variant 

Column 11: UniProt-Accession - UniProt accession id 

Column 12: AA-Pos - amino acid position on the transcript 

Column 13: REF-AA - the reference amino acid 

Column 14: ALT-AA - the alternative amino acid for this variant 



Column 15: MutationTaster - the score retrived from the MutationTaster2 website for this variant 

Column 16: MutationTaster Predicted Label for this variant 

Column 17: MutationAssessor - the score retrived from the MutationAssessor website for this variant 

Column 18: MutationAssessor Predicted Label for this variant 

Column 19: PolyPhen2 - the score retrived from the PolyPhen2 website for this variant 

Column 20: PolyPhen2 Predicted Label for this variant 

Column 21: CADD - the score retrived from the CADD website for this variant 

Column 22: SIFT - the score retrived from the SIFT website for this variant 

Column 23: SIFT Predicted Label for this variant 

Column 24: LRT - the score retrived from the LRT website for this variant 

Column 25: LRT Predicted Label for this variant 

Column 26: FatHMM-U - the score retrived from the FatHMM-U website for this variant 

Column 27: FatHMM-U Predicted Label for this variant 

Column 28: FatHMM-W - the score retrived from the FatHMM-W website for this variant 

Column 29: FatHMM-W weighting feature ln(Wd) 

Column 30: FatHMM-W weighting feature ln(Wn) 

Column 31: FatHMM-W Predicted Label for this variant 

Column 32: GERP++ - the GERP++ score

Column 33: phyloP - the phyloP score

Column 34: Condel (PP2 + MutationAssessor + SIFT) - the Condel score 

Column 35: Condel Predicted Label for this variant 

Column 36: Condel+ (PP2 + MutationAssessor + SIFT + FatHMM-W) - the Condel+ score 

Column 37: Condel+ Predicted Label for this variant 

Column 38: Logit (PP2 + MutationAssessor + SIFT) - the Logit score 

Column 39: Logit Predicted Label for this variant 

Column 40: Logit+ (PP2 + MutationAssessor + SIFT + FatHMM-W) - the Logit+ score 

Column 41: Logit+ Predicted Label for this variant 

