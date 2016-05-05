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

def printTableS2(score_data=None):
	val = utils.Validation()
	mv_predictor = utils.ProteinMajorityVote()
	lr_predictor = utils.LogisticRegression()
	f = open(os.path.abspath("Output/Supplementary/tableS2.csv"),'w')
	f.write(";;FatHMM-W;Logistic Regression over the features ln(Wn) & ln(Wd);Protein Majority Vote (MV)\n")
	datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
	for i,dataset in enumerate(datasets):
		score_data.selectDataset(dataset)
		labels = score_data.getTrueLabels()
		
		print "\tTraining Logistic Regression on Features ln(Wn) and ln(Wd)"
		lr = lr_predictor.run(true_labels=labels,features=score_data.getFatHMMFeatures(),folds=10)
		print "\tPerforming a Protein Majority Vote for dataset: " + dataset
		mv =  mv_predictor.getMV4Dataset(true_labels=labels,proteins=score_data.getUniprotIDs(),folds=10)
		
		string = dataset + ";AUC;"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc'] + ";"
		string += "%.2f (%.2f)"% (lr['auc'],lr['auc_std']) + ";" 
		string += "%.2f (%.2f)"% (mv['auc'],mv['auc_std'])
		f.write(string + "\n")
		
		string = ";AUC-PR;"
		string += "%.2f"% val.getROCStats(labels,score_data.getScores('fathmm_w'))['pr_auc'] + ";"
		string += "%.2f (%.2f)"% (lr['auc_pr'],lr['auc_pr_std']) + ";"
		string += "%.2f (%.2f)"% (mv['auc_pr'],mv['auc_pr_std'])
		f.write(string + "\n")
		
		string = ";Accuracy;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['accuracy'] + ";"
		string += "%.2f (%.2f)"% (lr['accuracy'],lr['accuracy_std']) + ";"
		string += "%.2f (%.2f)"% (mv['accuracy'],mv['accuracy_std'])
		f.write(string + "\n")
		
		string = ";F1-Score;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['f1_score'] + ";"
		string += "%.2f (%.2f)"% (lr['f1_score'],lr['f1_score_std']) + ";"
		string += "%.2f (%.2f)"% (mv['f1_score'],mv['f1_score_std'])
		f.write(string + "\n")
		
		string = ";MCC;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['mcc'] + ";"
		string += "%.2f (%.2f)"% (lr['mcc'],lr['mcc_std']) + ";"
		string += "%.2f (%.2f)"% (mv['mcc'],mv['mcc_std'])
		f.write(string + "\n")
		
		string = ";Precision;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['precision'] + ";"
		string += "%.2f (%.2f)"% (lr['precision'],lr['precision_std']) + ";"
		string += "%.2f (%.2f)"% (mv['precision'],mv['precision_std'])
		f.write(string + "\n")
		
		string = ";Recall;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['recall'] + ";"
		string += "%.2f (%.2f)"% (lr['recall'],lr['recall_std']) + ";"
		string += "%.2f (%.2f)"% (mv['recall'],mv['recall_std'])
		f.write(string + "\n")
		
		string = ";Negative Predictive Value;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['npv'] + ";"
		string += "%.2f (%.2f)"% (lr['npv'],lr['npv_std']) + ";"
		string += "%.2f (%.2f)"% (mv['npv'],mv['npv_std'])
		f.write(string + "\n")
		
		string = ";Specificity;"
		string += "%.2f"% val.getPredictionStats(labels,score_data.getPredictedLabels('fathmm_w'))['specificity'] + ";"
		string += "%.2f (%.2f)"% (lr['specificity'],lr['specificity_std']) + ";"
		string += "%.2f (%.2f)"% (mv['specificity'],mv['specificity_std'])
		f.write(string + "\n")
	f.close()
