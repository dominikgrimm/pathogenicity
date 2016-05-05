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
import utils,os
import matplotlib as mpl

def plotFigure(score_data=None,release=True):
    pl.ion()

    font_size = 13
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['font.size'] = font_size
    mpl.rcParams['patch.edgecolor'] = 'black'
    color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D',
           '#FFF79A','#998675','#A4A4A4','#5AFF00','#29A3A3','#F53DD6','#F2800D','#3399FF']
    
    hatch = pl.Rectangle((0,0),1,1,fill=None,hatch="///")
    circle = pl.Rectangle((0,0),1,1,fill=None,hatch="ooo")
    
    #load data
    val = utils.Validation()
    biased_map = {'mutationtaster':['humvar','exovar','varibench',
                    'varibench_selected','swissvar_selected'],
              'mutationassessor':['humvar','exovar','varibench'],
              'polyphen2':['humvar','exovar','varibench'],
              'logit_condel_old':['humvar','exovar','varibench'],
              'logit_condel_new':['humvar','exovar','varibench'],
              'condel_old':['humvar','exovar','varibench'],
              'condel_new':['humvar','exovar','varibench'],
              'sift':[''],
              'fathmm_u':[''],
              'fathmm_w':['humvar','exovar','varibench'],
              'gerp++':[''],
              'phylop':[''],
              'CADD':['']
    }
    
    datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
    n = len(datasets)
    mass_biased = np.zeros(n)
    mass = np.zeros(n)
    pp2 = np.zeros(n)
    pp2_biased = np.zeros(n)
    sift = np.zeros(n)
    logit_biased = np.zeros(n)
    logit = np.zeros(n)
    condel_biased = np.zeros(n)
    condel = np.zeros(n)
    
    fathmmw = np.zeros(n)
    fathmmw_biased = np.zeros(n)
    fathmmw_type2_biased = np.zeros(n)
    logit_p_biased = np.zeros(n)
    logit_p_type2_biased = np.zeros(n)
    logit_p = np.zeros(n)
    condel_p_biased = np.zeros(n)
    condel_p_type2_biased = np.zeros(n)
    condel_p = np.zeros(n)

    for i,dataset in enumerate(datasets):
        score_data.selectDataset(dataset)
        labels = score_data.getTrueLabels()
        print "\tComputing AUC values for dataset: " + dataset
        if dataset in biased_map['mutationassessor']:
            mass_biased[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
        else:
            mass[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
        if dataset in biased_map['polyphen2']:
            pp2_biased[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
        else:
            pp2[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
        if dataset in biased_map['logit_condel_old']:
            logit_biased[i] = val.getROCStats(labels,score_data.getScores('logit_condel_old'))['auc']
        else:
            logit[i] = val.getROCStats(labels,score_data.getScores('logit_condel_old'))['auc']
        if dataset in biased_map['condel_old']:
            condel_biased[i] = val.getROCStats(labels,score_data.getScores('condel_old'))['auc']
        else:
            condel[i] = val.getROCStats(labels,score_data.getScores('condel_old'))['auc']
        if dataset in biased_map['fathmm_w']:
            fathmmw_biased[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
        else:
            if dataset == "varibench_selected" or dataset=="predictSNP_selected" or dataset=="swissvar_selected":
                fathmmw_type2_biased[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
            else:
                fathmmw[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
        if dataset in biased_map['condel_new']:
            condel_p_biased[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
        else:
            if dataset == "varibench_selected" or dataset=="predictSNP_selected" or dataset=="swissvar_selected":
                condel_p_type2_biased[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
            else:
                condel_p[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
        if dataset in biased_map['logit_condel_new']:
            logit_p_biased[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
        else:
            if dataset == "varibench_selected" or dataset=="predictSNP_selected" or dataset=="swissvar_selected":
                logit_p_type2_biased[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
            else:
                logit_p[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
        sift[i] = val.getROCStats(labels,score_data.getScores('sift'))['auc']
    
    pl.figure(figsize=(10,5))
    fig1 = pl.subplot(111)

    width = 0.07

    x=np.arange(n)

    tool_names = np.array(['FatHMM-W','PolyPhen-2','MutationAssessor','SIFT','Condel','Logit','Condel+','Logit+','Type 1 Biased','Type 2 Biased'])

    spines_to_remove = ['top','right','bottom']
    ax = fig1.get_axes()
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.grid(True)

    t0 = fig1.bar(x-width/2.0-3*width-0.06,fathmmw,width=width,color=color_t[5])
    t1 = fig1.bar(x-width/2.0-2*width-0.04,pp2,width=width,color=color_t[1])
    t2 = fig1.bar(x-width/2.0-width-0.02,mass,width=width,color=color_t[8])
    t3 = fig1.bar(x-width/2.0,sift,width=width,color=color_t[2])
    t4 = fig1.bar(x+width/2.0+0.02,condel,width=width,color=color_t[11])
    t5 = fig1.bar(x+width/2.0+width+0.04,logit,width=width,color=color_t[12])
    t6 = fig1.bar(x+width/2.0+2.0*width+0.06,condel_p,width=width,color=color_t[13])
    t7 = fig1.bar(x+width/2.0+3.0*width+0.08,logit_p,width=width,color=color_t[14])

    light_grey = np.array([float(248)/float(255)]*3)
    light_grey = "#FFFFFF"
    almost_black = '#262626'
    legend = ax.legend([t0,t1,t2,t3,t4,t5,t6,t7,hatch,circle],tool_names,frameon=True, scatterpoints=1,
            prop={'size':font_size},ncol=5,loc="upper center",fancybox=True,bbox_to_anchor=(0.5,1.0))
    rect = legend.get_frame()
    rect.set_facecolor(light_grey)
    rect.set_linewidth(0.0)
    texts = legend.texts
    for t in texts:
        t.set_color(almost_black)

    fig1.bar(x-width/2.0-3*width-0.06,fathmmw_biased,width=width,color=color_t[5],hatch="/o/o/")
    fig1.bar(x-width/2.0-3*width-0.06,fathmmw_type2_biased,width=width,color=color_t[5],hatch="ooo")
    fig1.bar(x-width/2.0-2*width-0.04,pp2_biased,width=width,color=color_t[1],hatch="///")
    fig1.bar(x-width/2.0-width-0.02,mass_biased,width=width,color=color_t[8],hatch="///")
    fig1.bar(x+width/2.0+0.02,condel_biased,width=width,color=color_t[11],hatch="///")
    fig1.bar(x+width/2.0+width+0.04,logit_biased,width=width,color=color_t[12],hatch="///")
    fig1.bar(x+width/2.0+2.0*width+0.06,condel_p_biased,width=width,color=color_t[13],hatch="/o/o/")
    fig1.bar(x+width/2.0+3.0*width+0.08,logit_p_type2_biased,width=width,color=color_t[14],hatch="ooo")
    fig1.bar(x+width/2.0+2.0*width+0.06,condel_p_type2_biased,width=width,color=color_t[13],hatch="ooo")
    fig1.bar(x+width/2.0+3.0*width+0.08,logit_p_biased,width=width,color=color_t[14],hatch="/o/o/")

    pl.xticks(x,['HumVar','ExoVar','VaribenchSelected','predictSNPSelected','SwissVarSelected'],fontsize=font_size)
    pl.ylabel("AUC")
    fig1.set_ylim(0.5,1.03)
    fig1.set_xlim(-0.5,n-0.5)
    fig1.text(-0.05,1.02,"a",fontsize=15,fontweight="bold",va="top",transform=fig1.transAxes)
    pl.yticks([0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
    pl.subplots_adjust(left=0.07,bottom=0.07,right=0.99,top=0.97,wspace=0.05)
    
    if release:
        pl.savefig(os.path.abspath('Output/Figures/Figure4.tiff'),dpi=300)
    else:
        pl.savefig(os.path.abspath('Output/Figures/Figure4.pdf'))
        pl.savefig(os.path.abspath('Output/Figures/Figure4.tiff'),dpi=300)
        pl.savefig(os.path.abspath('Output/Figures/Figure4.jpg'))
    pl.close()
