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
import pylab as pl
import numpy as np
from matplotlib import rc
import utils
import os

def plotFigure(score_data=None,release=True):
    pl.ion()
    
    font_size = 13
    
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    
    color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D','#FFF79A','#998675','#A4A4A4','#5AFF00','#29A3A3']
    
    hatch = pl.Rectangle((0,0),1,1,fill=None,hatch="///")
    circle = pl.Rectangle((0,0),1,1,fill=None,hatch="ooo")
    
    #load data
    val = utils.Validation()
    mv_predictor = utils.ProteinMajorityVote()
    lr_predictor = utils.LogisticRegression()
    
    biased_map = {'mutationtaster':['humvar','exovar','varibench',
                    'varibench_selected','predictSNP_selected','swissvar_selected'],
              'mutationassessor':['humvar','exovar','varibench'],
              'polyphen2':['humvar','exovar','varibench'],
              'sift':[''],
              'fathmm_u':[''],
              'fathmm_w':['humvar','exovar','varibench'],
              'gerp++':[''],
              'phylop':[''],
              'CADD':['']
    }
    
    datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
    n_datasets = len(datasets)
    mt_biased = np.zeros(n_datasets)
    mt = np.zeros(n_datasets)
    mass_biased = np.zeros(n_datasets)
    mass = np.zeros(n_datasets)
    pp2_biased = np.zeros(n_datasets)
    pp2 = np.zeros(n_datasets)
    fathmmw = np.zeros(n_datasets)
    fathmmw_biased = np.zeros(n_datasets)
    fathmmw_type2_biased = np.zeros(n_datasets)
    sift = np.zeros(n_datasets)
    lrt = np.zeros(n_datasets)
    fathmmu = np.zeros(n_datasets)
    gerp = np.zeros(n_datasets)
    phylop = np.zeros(n_datasets)
    cadd = np.zeros(n_datasets)
    mv = np.zeros(n_datasets)
    mv_biased = np.zeros(n_datasets)
    features = np.zeros(n_datasets)
    features_biased = np.zeros(n_datasets)
    features_type2_biased = np.zeros(n_datasets)

    for i,dataset in enumerate(datasets):
        score_data.selectDataset(dataset)
        labels = score_data.getTrueLabels()
        print "\tPerforming a Logistic Regression over the weighting features of FatHMM-W for dataset: " + dataset
        lr_value = lr_predictor.run(true_labels=labels,features=score_data.getFatHMMFeatures(),folds=10)['auc']
        print "\tComputing AUC values for dataset: " + dataset
        if dataset in biased_map['mutationtaster']:
            mt_biased[i] = val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc']
        else:
            mt[i] = val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc']
        if dataset in biased_map['mutationassessor']:
            mass_biased[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
        else:
            mass[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
        if dataset in biased_map['polyphen2']:
            pp2_biased[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
        else:
            pp2[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
        if dataset in biased_map['fathmm_w']:
            fathmmw_biased[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
            features_biased[i] = lr_value
        else:
            if dataset == "varibench_selected" or dataset=="predictSNP_selected" or dataset=="swissvar_selected":
                fathmmw_type2_biased[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
                features_type2_biased[i] = lr_value
            else:
                fathmmw[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
                features[i] = lr_value
        fathmmu[i] = val.getROCStats(labels,score_data.getScores('fathmm_u'))['auc']
        sift[i] = val.getROCStats(labels,score_data.getScores('sift'))['auc']
        lrt[i] = val.getROCStats(labels,score_data.getScores('lrt'))['auc']
        gerp[i] = val.getROCStats(labels,score_data.getScores('gerp++'))['auc']
        phylop[i] = val.getROCStats(labels,score_data.getScores('phylop'))['auc']
        cadd[i] = val.getROCStats(labels,score_data.getScores('CADD'))['auc']
        print "\tPerforming a Protein Majority Vote for dataset: " + dataset
        mv_biased[i] = mv_predictor.getMV4Dataset(true_labels=labels,proteins=score_data.getUniprotIDs(),folds=10)['auc']
    
    pl.figure(figsize=(15,5))
    fig1 = pl.subplot(111)

    width = 0.05

    x=np.arange(n_datasets)

    tool_names = np.array(['FatHMM-W','MutationTaster-2','PolyPhen-2','MutationAssessor','CADD','SIFT','LRT','FatHMM-U','Gerp++','phyloP','Features ln(Wn), ln(Wd)','Protein Majority Vote','Potentially Type 1 Biased','Potentially Type 2 Biased'])

    spines_to_remove = ['top','right','bottom']
    ax = fig1.get_axes()
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.grid(True)

    t0 = fig1.bar(x-width/2.0-5*width-0.1,fathmmw,width=width,color=color_t[5])
    t1 = fig1.bar(x-width/2.0-4*width-0.08,mt,width=width,color=color_t[0])
    t2 = fig1.bar(x-width/2.0-3*width-0.06,pp2,width=width,color=color_t[1])
    t3 = fig1.bar(x-width/2.0-2*width-0.04,mass,width=width,color=color_t[8])
    t4 = fig1.bar(x-width/2.0-width-0.02,cadd,width=width,color=color_t[9])
    t5 = fig1.bar(x-width/2.0,sift,width=width,color=color_t[2])
    t6 = fig1.bar(x+width/2.0+0.02,lrt,width=width,color=color_t[3])
    t7 = fig1.bar(x+width/2.0+width+0.04,fathmmu,width=width,color=color_t[4])
    t8 = fig1.bar(x+width/2.0+2*width+0.06,gerp,width=width,color=color_t[6])
    t9 = fig1.bar(x+width/2.0+3*width+0.08,phylop,width=width,color=color_t[7])
    t10 = fig1.bar(x+width/2.0+4*width+0.1,features,width=width,color=color_t[11])
    t11 = fig1.bar(x+width/2.0+5*width+0.12,mv,width=width,color=color_t[10])

    #fig1.text(-0.05,1.02,"b",fontsize=15,fontweight="bold",va="top",transform=fig1.transAxes)

    light_grey = np.array([float(248)/float(255)]*3)
    light_grey = "#FFFFFF"
    almost_black = '#262626'
    legend = ax.legend([t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,hatch,circle],tool_names,frameon=True, scatterpoints=1,
                prop={'size':11},ncol=7,loc="upper center",fancybox=True,bbox_to_anchor=(0.5, 1.02))
    legend.get_frame().set_alpha(0.5)
    rect = legend.get_frame()
    rect.set_facecolor(light_grey)
    rect.set_linewidth(0.0)
    # Change the legend label colors to almost black, too
    texts = legend.texts
    for t in texts:
        t.set_color(almost_black)

    fig1.bar(x-width/2.0-5*width-0.1,fathmmw_biased,width=width,color=color_t[5],hatch="/o/o/")
    fig1.bar(x-width/2.0-5*width-0.1,fathmmw_type2_biased,width=width,color=color_t[5],hatch="ooo")
    fig1.bar(x-width/2.0-4*width-0.08,mt_biased,width=width,color=color_t[0],hatch="///")
    fig1.bar(x-width/2.0-3*width-0.06,pp2_biased,width=width,color=color_t[1],hatch="///")
    fig1.bar(x-width/2.0-2*width-0.04,mass_biased,width=width,color=color_t[8],hatch="///")
    fig1.bar(x+width/2.0+4*width+0.1,features_biased,width=width,color=color_t[11],hatch="/o/o/")
    fig1.bar(x+width/2.0+4*width+0.1,features_type2_biased,width=width,color=color_t[11],hatch="ooo")
    fig1.bar(x+width/2.0+5*width+0.12,mv_biased,width=width,color=color_t[10],hatch="ooo")

    pl.xticks(x,['HumVar','ExoVar','VariBenchSelected','predictSNPSelected','SwissVarSelected'],fontsize=font_size)
    fig1.set_ylabel("AUC")
    fig1.set_ylim(0.5,1.06)
    fig1.set_xlim(-0.5,n_datasets-0.5)
    pl.yticks([0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
    pl.subplots_adjust(left=0.06,bottom=0.07,right=0.99,top=0.99,wspace=0.05)

    if release:
        pl.savefig(os.path.abspath('Output/Figures/Figure1.tiff'),dpi=300)
    else:
        pl.savefig(os.path.abspath('Output/Figures/Figure1.pdf'))
        pl.savefig(os.path.abspath('Output/Figures/Figure1.tiff'),dpi=300)
        pl.savefig(os.path.abspath('Output/Figures/Figure1.jpg'))
    pl.close()
