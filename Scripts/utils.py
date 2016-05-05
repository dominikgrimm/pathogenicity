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
import numpy as np
import scipy.stats as stats
import os
from sklearn.metrics import roc_curve,auc,precision_recall_curve
from sklearn import cross_validation
from sklearn import linear_model

import sklearn.metrics as metrics

_humvar_ = os.path.abspath("ToolScores/humvar_tool_scores.csv")
_exovar_ = os.path.abspath("ToolScores/exovar_tool_scores.csv")
_varibench_selected_ = os.path.abspath("ToolScores/varibench_selected_tool_scores.csv")
_predictSNP_selected_ = os.path.abspath("ToolScores/predictSNP_selected_tool_scores.csv")
_swissvar_selected_ = os.path.abspath("ToolScores/swissvar_selected_tool_scores.csv")

dataset_map = {'humvar':_humvar_,
           'exovar':_exovar_,
           'varibench_selected':_varibench_selected_,
           'swissvar_selected':_swissvar_selected_,
           'predictSNP_selected':_predictSNP_selected_}

similarity_data = {'varibench_selected':os.path.abspath("Similarities/varibench_selected_similarities.csv"),
                   'predictSNPSelected':os.path.abspath("Similarities/predictSNPSelected_similarities.csv"),
                   'swissvar_selected':os.path.abspath("Similarities/swissvar_selected_similarities.csv"),
                   }

class DataScores():

    def __init__(self):
        self.__dataset_names = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
        self.__data_map = {}
        self.__dataMatrix = []
        self.__file_index = [14,16,18,20,21,23,25,27,28,29,31,32,33,35,37,39]
        self.__tool_names = np.array(['mutationtaster','mutationassessor','polyphen2',
                     'CADD','sift','lrt','fathmm_u','fathmm_w','wd','wn','gerp++','phylop',
                     'condel_old','condel_new','logit_condel_old','logit_condel_new'])
        self.__real_labels = []
        self.__real_labels_map = {}
        self.__predicted_labels = []
        self.__predicted_labels_map = {}
        self.__uniprot_map = {}
        self.__uniprot_ids = []
        self.__category_map = {}
        self.__category_variant_map = {}
        self.__selectedDataset = None
        self.__imputed = False
        self.__n_total_proteins = 0
        self.__n_mixed_proteins = 0
        self.__n_del_proteins = 0
        self.__n_neutral_proteins = 0
        self.__n_total_variants = 0
        self.__n_mixed_variants = 0
        self.__n_del_variants = 0
        self.__n_neutral_variants = 0

    '''
    Private member methods
    '''
    def __compute_ratio(self,ratio=None,protein=None,a1=None,a2=None):
        ind_ratio_del = np.where(ratio>a1)[0]
        r_mixed = ratio[ind_ratio_del]
        protein_tmp = protein[ind_ratio_del]
        ind_ratio_neutral = np.where(r_mixed<a2)[0]
        pst = np.unique(protein_tmp[ind_ratio_neutral]).shape[0]
        mixed = np.unique(protein_tmp[ind_ratio_neutral])
        return [pst,r_mixed[ind_ratio_neutral].shape[0],mixed]

    '''
    Public member methods
    '''
    def load_data(self):
        for dataset in self.__dataset_names:
            self.__filename = dataset_map[dataset]
            f = open(self.__filename)
            self.__dataMatrix = []
            self.__real_labels = []
            self.__predicted_labels = []
            self.__uniprot_ids = []
            for i,line in enumerate(f):
                if i==0:
                    continue
                else:
                    sv = line.strip().split(",")
                    self.__real_labels.append(int(sv[0]))
                    data_line = np.ones(len(self.__file_index))
                    #self.__uniprot_ids.append(sv[10])
                    self.__uniprot_ids.append(sv[8])
                    predicted_line = np.ones(len(self.__file_index))*np.nan
                    for j,index in enumerate(self.__file_index):
                        data_line[j] = float(sv[index])
                        if index==20 or index==31 or index==32 or index==28 or index==29:
                            continue
                        else:
                            if index==27:
                                predicted_line[j] = float(sv[index+3])
                            else:
                                predicted_line[j] = float(sv[index+1])
                    self.__dataMatrix.append(data_line)
                    self.__predicted_labels.append(predicted_line)
            self.__uniprot_map[dataset] = np.array(self.__uniprot_ids)
            self.__predicted_labels = np.array(self.__predicted_labels)
            self.__dataMatrix = np.array(self.__dataMatrix)
            self.__real_labels = np.array(self.__real_labels)
            self.__data_map[dataset] = self.__dataMatrix
            self.__predicted_labels_map[dataset] = self.__predicted_labels
            self.__real_labels_map[dataset] = self.__real_labels
            f.close()

    def selectDataset(self,dataset=None):
        self.__dataMatrix = self.__data_map[dataset]
        self.__predicted_labels = self.__predicted_labels_map[dataset]
        self.__real_labels = self.__real_labels_map[dataset]
        self.__uniprot_ids = self.__uniprot_map[dataset]

    def getUniprotIDs(self):
        return self.__uniprot_ids

    def getToolNames(self):
        return np.array(['mutationtaster','mutationassessor','polyphen2',
                         'CADD','sift','lrt','fathmm_u','fathmm_w','gerp++','phylop',
                 'condel_old','condel_new','logit_condel_old','logit_condel_new'])
    
    def getRealNames(self):
        return np.array(['MutationTaster2','MutationAssessor','PolyPhen-2',
                     'CADD','SIFT','LRT','FATHMM-U','FATHMM-W','GERP++','PhyloP',
                     'Condel','Condel+','Logit','Logit+'])

    def getFatHMMFeatures(self):
        return self.__dataMatrix[:,[8,9]]

    def getScores(self,tool_name=None):
        ind = np.where(tool_name==self.__tool_names)[0]
        if self.__imputed==False:
            if tool_name=="fathmm_w" or tool_name=="fathmm_u" or tool_name=="sift" or tool_name=="lrt":
                return self.__dataMatrix[:,ind]*(-1)
            else:
                return self.__dataMatrix[:,ind]
        else:
            return self.__dataMatrix[:,ind]

    def getPredictedLabels(self,tool_name=None):
        ind = np.where(tool_name==self.__tool_names)[0]
        labels = self.__predicted_labels[:,ind]
        #return labels[:,np.newaxis]
        return labels.ravel()

    def getTrueLabels(self):
        return self.__real_labels

    def getDatasets(self):
        return self.__dataset_names
    
    def getSimilarityData(self,dataset=None):
        if dataset in similarity_data:
            f = open(similarity_data[dataset])
            data_matrix = []
            for i,line in enumerate(f):
                if i==0:
                    continue
                sv = line.strip().split(" ")
                row = []
                row.append(sv[5])
                row.append(sv[7])
                data_matrix.append(row)
            f.close()
            data_matrix = np.array(data_matrix)
            return data_matrix
        else:
            return None

    def loadCategories(self):
        [proteins_unique,unique_ind,inverse_ind] = np.unique(self.__uniprot_ids,return_inverse=True,return_index=True)
        fmatrix = np.zeros((self.__real_labels.shape[0],2))
        for i,protein in enumerate(proteins_unique):
            n = np.where(protein==self.__uniprot_ids)[0]
            deleterious = (self.__real_labels[n]==1).sum()
            neutral = (self.__real_labels[n]==-1).sum()
            fmatrix[n,:] = np.array([deleterious,neutral])
        ind_del = np.where(fmatrix[:,0]>0)[0]
        del_mat = fmatrix[ind_del]
        ind_mixed = np.where(del_mat[:,1]>0)[0]
        mixed_mat = del_mat[ind_mixed]
        ind_only_del = np.where(del_mat[:,1]==0)[0]
        #only_del = del_mat[ind_only_del]
        ind_neutral = np.where(fmatrix[:,0]==0)[0]
        neutral_mat = fmatrix[ind_neutral]
        ind_only_neutral = np.where(neutral_mat[:,1]>0)[0]
        #only_neutral = neutral_mat[ind_only_neutral]
        #Store statistics
        self.__n_total_proteins = float(proteins_unique.shape[0])
        self.__n_mixed_proteins = np.unique(self.__uniprot_ids[ind_del][ind_mixed]).shape[0]
        self.__n_del_proteins = np.unique(self.__uniprot_ids[ind_del][ind_only_del]).shape[0]
        self.__n_neutral_proteins = np.unique(self.__uniprot_ids[ind_neutral][ind_only_neutral]).shape[0]
        self.__n_total_variants = self.__uniprot_ids.shape[0]
        self.__n_mixed_variants = (self.__uniprot_ids[ind_del][ind_mixed]).shape[0]
        self.__n_del_variants = (self.__uniprot_ids[ind_del][ind_only_del]).shape[0]
        self.__n_neutral_variants = (self.__uniprot_ids[ind_neutral][ind_only_neutral]).shape[0]
        #compute ratios
        mixed_ratio = mixed_mat[:,1]/mixed_mat.sum(axis=1)
        mixed_proteins = self.__uniprot_ids[ind_del][ind_mixed]
        unique_mixed_proteins = np.unique(mixed_proteins)
        [cat2,variants2,proteins2] = self.__compute_ratio(mixed_ratio,mixed_proteins,0.1,0.9)
        [cat3,variants3,proteins3] = self.__compute_ratio(mixed_ratio,mixed_proteins,0.2,0.8)
        [cat4,variants4,proteins4] = self.__compute_ratio(mixed_ratio,mixed_proteins,0.3,0.7)
        [cat5,variants5,proteins5] = self.__compute_ratio(mixed_ratio,mixed_proteins,0.4,0.6)
        
        #retrive single variant proteins
        #1. Deleterious single variants
        t1 = np.where(fmatrix[:,0]==1)[0]
        f1 = fmatrix[t1]
        t2 = np.where(f1[:,1]==0)[0]
        single_variant_proteins = np.unique(self.__uniprot_ids[t1][t2])
        #2. Neutral single variants
        t1 = np.where(fmatrix[:,0]==0)[0]
        f1 = fmatrix[t1]
        t2 = np.where(f1[:,1]==1)[0]
        single_variant_proteins = np.concatenate([single_variant_proteins,np.unique(self.__uniprot_ids[t1][t2])])

        #retrive pure proteins 
        tmp = np.union1d(unique_mixed_proteins,single_variant_proteins)
        pure_proteins = np.unique(np.setdiff1d(proteins_unique,tmp))
        pure_proteins = np.unique(np.union1d(pure_proteins,single_variant_proteins))
        #store data 
        #self.__category_map['pure'] = pure_proteins
        #self.__category_map['single_variant'] = single_variant_proteins
        #self.__category_map['cat1'] = mixed_proteins
        #self.__category_map['cat2'] = proteins2
        #self.__category_map['cat3'] = proteins3
        #self.__category_map['cat4'] = proteins4
        #self.__category_map['cat5'] = proteins5
        self.__category_map['pure'] = pure_proteins
        self.__category_map['single_variant'] = single_variant_proteins
        self.__category_map['cat1'] = unique_mixed_proteins
        self.__category_map['cat2'] = proteins2
        self.__category_map['cat3'] = proteins3
        self.__category_map['cat4'] = proteins4
        self.__category_map['cat5'] = proteins5
        self.__category_variant_map['cat1'] = mixed_proteins.shape[0]
        self.__category_variant_map['cat2'] = variants2
        self.__category_variant_map['cat3'] = variants3
        self.__category_variant_map['cat4'] = variants4
        self.__category_variant_map['cat5'] = variants5

    def getVariantFractions(self):
        return [self.__n_total_variants,self.__n_mixed_variants,self.__n_del_variants,self.__n_neutral_variants]

    def getProteinFractions(self):
        return [self.__n_total_proteins,self.__n_mixed_proteins,self.__n_del_proteins,self.__n_neutral_proteins]

    def getProteinsCategorie(self,category=None):
        return self.__category_map[category]
    
    def getVariantsCategorie(self,category=None):
        return self.__category_variant_map[category]

    def getData4Categorie(self,category=None,tool_name=None):
        unique_acc = self.__category_map[category]
        ind = (np.reshape(self.__uniprot_ids,(self.__uniprot_ids.shape[0],1))==unique_acc).nonzero()[0]
        scores = self.getScores(tool_name)
        return [self.__real_labels[ind],scores[ind]]

class Validation():
    
    def getROCStats(self,true_labels=None,ranking=None):
        validation = {}
        ind = np.where(~np.isnan(ranking))[0]
        not_enough_data = False
        if ind.shape[0]<ranking.shape[0]*0.05:
            not_enough_data = True
        if ind.shape[0]==0:
            not_enough_data = True
        if not_enough_data:
            validation['auc'] = np.nan
            validation['fpr'] = None
            validation['tpr'] = None
            validation['thresholds'] = None
            validation['pr_auc'] = None
            validation['precission_ranking'] = None
            validation['recall_ranking'] = None
            validation['pr_thresholds'] = None
            return validation
        
        ranking = ranking[ind]
        true_labels = true_labels[ind]

        fpr,tpr,thresholds = roc_curve(true_labels,ranking)
        roc_auc = auc(fpr,tpr)
        validation['auc'] = roc_auc
        validation['fpr'] = fpr
        validation['tpr'] = tpr
        validation['thresholds'] = thresholds
        try:
            precission_ranking,recall_ranking,pr_thresholds = precision_recall_curve(true_labels,ranking)
            pr_auc = auc(recall_ranking,precission_ranking)
            validation['pr_auc'] = pr_auc
            validation['precission_ranking'] = precission_ranking
            validation['recall_ranking'] = recall_ranking
            validation['pr_thresholds'] = pr_thresholds
        except:
            validation['pr_auc'] = None
            validation['precission_ranking'] = None
            validation['recall_ranking'] = None
            validation['pr_thresholds'] = None
        return validation
    
    def getTTestPval(self,true_labels=None,ranking=None):
        ind = np.where(~np.isnan(ranking))[0]
        ranking = ranking[ind]
        true_labels = true_labels[ind]
        cases = np.where(true_labels==1)[0]
        controls = np.where(true_labels==-1)[0]
        if controls.shape[0]==0 or cases.shape[0]==0:
            validation = {}
            validation['cases'] = 0
            validation['controls'] = 0
            validation['cases_mean'] = 0
            validation['controls_mean'] = 0
            validation['cases_variance'] = 0
            validation['controls_variance'] = 0
            validation['df'] = 0
            validation['tstat'] = 0
            validation['pv'] = 1
            return validation
        else:
            ts,pv = stats.ttest_ind(ranking[cases].flatten(),ranking[controls].flatten())
            validation = {}
            validation['cases'] = cases.flatten().shape[0]
            validation['controls'] = controls.flatten().shape[0]
            validation['cases_mean'] = ranking[cases].flatten().mean()
            validation['controls_mean'] = ranking[controls].flatten().mean()
            validation['cases_variance'] = ranking[cases].flatten().var(ddof=1)
            validation['controls_variance'] = ranking[controls].flatten().var(ddof=1)
            validation['df'] = cases.flatten().shape[0]+controls.flatten().shape[0]-2.0
            validation['tstat'] = ts
            validation['pv'] = pv
            return validation

    def getPredictionStats(self,true_labels=None,predictions=None):
        ind = np.where(~np.isnan(predictions))[0]
        validation = {}
        if ind.shape[0] == 0:
            validation['TP'] = np.nan
            validation['FP'] = np.nan
            validation['TN'] = np.nan
            validation['FN'] = np.nan
            validation['accuracy'] = np.nan
            validation['f1_score'] = np.nan
            validation['mcc'] = np.nan
            validation['precision'] = np.nan
            validation['recall'] = np.nan
            validation['npv'] = np.nan
            validation['specificity'] = np.nan
        else:
            true_labels = true_labels[ind]
            predictions = predictions[ind]
            validation = {}
            validation['TP'] = 0
            validation['FP'] = 0
            validation['TN'] = 0
            validation['FN'] = 0
            for i,pred in enumerate(predictions):
                if pred>0 and true_labels[i]>0:
                    validation['TP'] += 1
                elif pred>0 and true_labels[i]<0:
                    validation['FP'] += 1
                elif pred<0 and true_labels[i]<0:
                    validation['TN'] += 1
                elif pred<0 and true_labels[i]>0:
                    validation['FN'] += 1
            validation['TP'] = int(validation['TP'])
            validation['FP'] = int(validation['FP'])
            validation['TN'] = int(validation['TN'])
            validation['FN'] = int(validation['FN'])
            validation['accuracy'] = metrics.accuracy_score(true_labels,predictions)
            validation['f1_score'] = metrics.f1_score(true_labels,predictions)
            validation['mcc'] = metrics.matthews_corrcoef(true_labels,predictions)
            validation['precision'] = metrics.precision_score(true_labels,predictions)
            validation['recall'] = metrics.recall_score(true_labels,predictions)
            validation['npv'] = float(validation['TN'])/float(validation['TN'] + validation['FN'])
            validation['specificity'] = float(validation['TN'])/float(validation['FP'] + validation['TN'])
        return validation

class ProteinMajorityVote():

    def getMV4Dataset(self,true_labels=None,proteins=None,folds=10):
        skf = cross_validation.StratifiedKFold(true_labels,n_folds=10,shuffle=True)
        auc_l = []
        auc_pr_l = []
        accuracy_l = []
        f1_score_l = []
        mcc_l = []
        precision_l = []
        recall_l = []
        npv_l = []
        specificity_l = []
        for train_index, test_index in skf:
            protein_train, protein_test = proteins[train_index], proteins[test_index]
            y_train, y_test = true_labels[train_index], true_labels[test_index]
            
            #get unique proteins 
            [proteins_unique,unique_ind,inverse_ind] = np.unique(protein_train,return_inverse=True,return_index=True)

            mv_map = {}
            for i,protein in enumerate(proteins_unique):
                ind = np.where(protein==protein_train)[0]
                deleterious = (y_train[ind]==1).sum()
                neutral = (y_train[ind]==-1).sum()
                mv_map[protein] = float(deleterious)/float(deleterious+neutral)
            #create test_scores
            [proteins_unique,unique_ind,inverse_ind] = np.unique(protein_test,return_inverse=True,return_index=True)
            test_scores = np.zeros((y_test.shape[0],1))
            test_predictions = np.zeros((y_test.shape[0],1))
            for i,protein in enumerate(proteins_unique):
                ind = np.where(protein==protein_test)
                if protein in mv_map:
                    test_scores[ind] = mv_map[protein]
                    if mv_map[protein]>0.5:
                        test_predictions[ind] = 1
                    else:
                        test_predictions[ind] = -1
                else:
                    test_scores[ind] = 0.5
                    if np.random.binomial(1,0.5)==1:
                        test_predictions[ind] = 1
                    else:
                        test_predictions[ind] = -1
            fpr,tpr,thresholds = roc_curve(y_test,test_scores)
            auc_l.append(auc(fpr,tpr))
            precission_ranking,recall_ranking,pr_thresholds = precision_recall_curve(y_test,test_scores)
            auc_pr_l.append(auc(recall_ranking,precission_ranking))
            val = Validation()
            validation_test = val.getPredictionStats(y_test,test_predictions)
            accuracy_l.append(validation_test['accuracy'])
            f1_score_l.append(validation_test['f1_score'])
            mcc_l.append(validation_test['mcc'])
            precision_l.append(validation_test['precision'])
            recall_l.append(validation_test['recall'])
            npv_l.append(validation_test['npv'])
            specificity_l.append(validation_test['specificity'])

        auc_l = np.array(auc_l)
        auc_pr_l = np.array(auc_pr_l)
        accuracy_l = np.array(accuracy_l)
        f1_score_l = np.array(f1_score_l)
        mcc_l = np.array(mcc_l)
        precision_l = np.array(precision_l)
        recall_l = np.array(recall_l)
        npv_l = np.array(npv_l)
        specificity_l = np.array(specificity_l)

        validation = {}
        validation['auc'] = auc_l.mean()
        validation['auc_pr'] = auc_pr_l.mean()
        validation['auc_std'] = auc_l.std()
        validation['auc_pr_std'] = auc_pr_l.std()
        validation['accuracy'] = accuracy_l.mean()
        validation['accuracy_std'] = accuracy_l.std()
        validation['f1_score'] = f1_score_l.mean()
        validation['f1_score_std'] = f1_score_l.std()
        validation['mcc'] = mcc_l.mean()
        validation['mcc_std'] = mcc_l.std()
        validation['precision'] = precision_l.mean()
        validation['precision_std'] = precision_l.std()
        validation['recall'] = recall_l.mean()
        validation['recall_std'] = recall_l.std()
        validation['npv'] = npv_l.mean()
        validation['npv_std'] = npv_l.std()
        validation['specificity'] = specificity_l.mean()
        validation['specificity_std'] = specificity_l.std()
        return validation
class LogisticRegression():

    def run(self,true_labels=None,features=None,folds=10):
        c_values = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2]
        #remove missing values
        ind = np.where(~np.isnan(features[:,0]))[0]
        true_labels = true_labels[ind]
        features = features[ind,:]
        skf = cross_validation.StratifiedKFold(true_labels,n_folds=10,shuffle=True)
        auc_l = []
        auc_pr_l = []
        accuracy_l = []
        f1_score_l = []
        mcc_l = []
        precision_l = []
        recall_l = []
        npv_l = []
        specificity_l = []
        for train_index, test_index in skf:
            train_y, test_y = true_labels[train_index], true_labels[test_index]
            train_features, test_features = features[train_index], features[test_index]
            
            #create subtrainingsets for internal parameter estimation
            best_auc = 0
            best_c = 0
            skf2 = cross_validation.StratifiedKFold(train_y,n_folds=10,shuffle=True)
            for j in xrange(len(c_values)):
                subtrain_y, subtest_y = train_y[list(skf2)[0][0]], train_y[list(skf2)[0][1]]
                subtrain_features, subtest_features = train_features[list(skf2)[0][0]], train_features[list(skf2)[0][1]]
                logit = linear_model.LogisticRegression(C = c_values[j], penalty="l1")
                logit.fit(subtrain_features,subtrain_y)
                z = logit.decision_function(subtest_features)
                fpr,tpr,thresholds = roc_curve(subtest_y,z)
                tmp_auc = auc(fpr,tpr)
                if tmp_auc > best_auc:
                    best_auc = tmp_auc
                    best_c = c_values[j]

            logit = linear_model.LogisticRegression(C = best_c, penalty="l1")
            logit.fit(train_features,train_y)
            z = logit.decision_function(test_features)
            fpr,tpr,thresholds = roc_curve(test_y,z)
            auc_l.append(auc(fpr,tpr))
            precission_ranking,recall_ranking,pr_thresholds = precision_recall_curve(test_y,z)
            auc_pr_l.append(auc(recall_ranking,precission_ranking))
            val = Validation()
            validation_test = val.getPredictionStats(test_y,logit.predict(test_features))
            accuracy_l.append(validation_test['accuracy'])
            f1_score_l.append(validation_test['f1_score'])
            mcc_l.append(validation_test['mcc'])
            precision_l.append(validation_test['precision'])
            recall_l.append(validation_test['recall'])
            npv_l.append(validation_test['npv'])
            specificity_l.append(validation_test['specificity'])
            
        auc_l = np.array(auc_l)
        auc_pr_l = np.array(auc_pr_l)
        accuracy_l = np.array(accuracy_l)
        f1_score_l = np.array(f1_score_l)
        mcc_l = np.array(mcc_l)
        precision_l = np.array(precision_l)
        recall_l = np.array(recall_l)
        npv_l = np.array(npv_l)
        specificity_l = np.array(specificity_l)

        validation = {}
        validation['auc'] = auc_l.mean()
        validation['auc_pr'] = auc_pr_l.mean()
        validation['auc_std'] = auc_l.std()
        validation['auc_pr_std'] = auc_pr_l.std()
        validation['accuracy'] = accuracy_l.mean()
        validation['accuracy_std'] = accuracy_l.std()
        validation['f1_score'] = f1_score_l.mean()
        validation['f1_score_std'] = f1_score_l.std()
        validation['mcc'] = mcc_l.mean()
        validation['mcc_std'] = mcc_l.std()
        validation['precision'] = precision_l.mean()
        validation['precision_std'] = precision_l.std()
        validation['recall'] = recall_l.mean()
        validation['recall_std'] = recall_l.std()
        validation['npv'] = npv_l.mean()
        validation['npv_std'] = npv_l.std()
        validation['specificity'] = specificity_l.mean()
        validation['specificity_std'] = specificity_l.std()
        return validation
