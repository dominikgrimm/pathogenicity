'''
################################################################################################################
Manuscript Title: Evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity
################################################################################################################
Manuscripts Authors: Dominik G. Grimm, Chloe Agathe Azencott, Fabian Aicheler, Udo Gieraths, Daniel G. MacArthur, Kaitlin E. Samocha, David N. Cooper, Peter D. Stenson, Mark J. Daly, Jordan W. Smoller, Laramie E. Duncan, Karsten M. Borgwardt

Code by: Dominik Gerhard Grimm
Year: 2014/2015
Group: Machine Learning and Computational Biology Research Group
Insitute: Max Planck Institute for Intelligent Systems and Max Planck Institue for Developmental Biology

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
'''
import utils,os

def printTableS1(score_data=None):
	val = utils.Validation()
	f = open(os.path.abspath("Output/Supplementary/tableS1.csv"),'w')
	f.write(";;MT2;PP2;MASS;CADD;SIFT;LRT;FatHMM-U;FatHMM-W;Gerp++;phyloP\n")
	datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
	for i,dataset in enumerate(datasets):
		score_data.selectDataset(dataset)
		labels = score_data.getTrueLabels()
		string = dataset + ";TP;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['TP'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['TP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['TP']
		f.write(string + "\n")

		string = ";FP;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['FP'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['FP'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['FP']
		f.write(string + "\n")
		
		string = ";TN;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['TN'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['TN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['TN']
		f.write(string + "\n")
		
		string = ";FN;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['FN'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['FN'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['FN']
		f.write(string + "\n")
		
		string = ";AUC;"
		string += "%.2f"%val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('polyphen2'))['auc'] + ";"
		string += "%.2f"%val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('CADD'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('sift'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('lrt'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_u'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('gerp++'))['auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('phylop'))['auc'] 
		f.write(string + "\n")
		string = ";AUC-PR;"
		string += "%.2f"%val.getROCStats(labels,score_data.getScores('mutationtaster'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('polyphen2'))['pr_auc'] + ";"
		string += "%.2f"%val.getROCStats(labels,score_data.getScores('mutationassessor'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('CADD'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('sift'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('lrt'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_u'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_w'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('gerp++'))['pr_auc'] + ";"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('phylop'))['pr_auc']
		f.write(string + "\n")
		string = ";Accuracy;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['accuracy'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['accuracy'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['accuracy']
		f.write(string + "\n")
		string = ";F-Score;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['f1_score'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['f1_score'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['f1_score']
		f.write(string + "\n")
		string = ";MCC;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['mcc'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['mcc'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['mcc']
		f.write(string + "\n")
		
		string = ";Precision/Positive Predictive Value;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['precision'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['precision'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['precision']
		f.write(string + "\n")
		
		string = ";Recall/Sensitivity;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['recall'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['recall'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['recall']
		f.write(string + "\n")
		
		string = ";Specificity;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['specificity'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['specificity'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['specificity']
		f.write(string + "\n")
		
		string = ";Negative Predictive Value;"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationtaster'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('polyphen2'))['npv'] + ";"
		string += "%.2f"%val.getPredictionStats(labels,score_data.getPredictedLabels('mutationassessor'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('CADD'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('sift'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('lrt'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_u'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('gerp++'))['npv'] + ";"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('phylop'))['npv']
		f.write(string + "\n")
	f.close()
