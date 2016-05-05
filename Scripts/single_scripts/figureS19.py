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
import numpy as sp
import matplotlib as mpl
import os,utils
from matplotlib.font_manager import FontProperties

def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    ax = axes or pl.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()


def plotFigure(score_data=None,release=True):
    pl.ion()
    font_size = 10

    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 0.8
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'black'
    
    color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D',
           '#FFF79A','#998675','#A4A4A4','#5AFF00','#29A3A3','#F53DD6','#F2800D','#3399FF']
    
    val = utils.Validation()
    categories = ['all','pure','cat1','cat2','cat3','cat4','cat5']

    
    pl.figure(figsize=(10,7))
    pl.subplot(311)
    condel = sp.zeros(7)
    condel_p = sp.zeros(7)
    logit = sp.zeros(7)
    logit_p = sp.zeros(7)
    score_data.selectDataset('varibench_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            condel[i] = val.getROCStats(labels,score_data.getScores('condel_old'))['auc']
            condel_p[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
            logit[i] = val.getROCStats(labels,score_data.getScores('logit_condel_old'))['auc']
            logit_p[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_old')
            condel[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_new')
            condel_p[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_old')
            logit[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_new')
            logit_p[i] = val.getROCStats(labels,scores)['auc']
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(condel,'o-',color=color_t[11])
    pl.plot(condel_p,'h--',color=color_t[13])
    pl.plot(logit,'x-',color=color_t[12])
    pl.plot(logit_p,'<--',color=color_t[14])
    pl.ylim(0.5,1)
    pl.xlim(-.1,6.1)
    pl.yticks([0.5,0.6,0.7,0.8,0.9,1.0])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(7),['All','Pure',']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    leg = pl.legend(['Condel','Condel+','Logit','Logit+'],loc='upper center',fancybox=True, ncol=5,prop={'size':10},numpoints=1)
    leg.get_frame().set_alpha(0.2)
    leg.get_frame().set_edgecolor("none")
    pl.title("VariBenchSelected")
    remove_border()
    
    pl.subplot(312)
    condel = sp.zeros(7)
    condel_p = sp.zeros(7)
    logit = sp.zeros(7)
    logit_p = sp.zeros(7)
    score_data.selectDataset('predictSNP_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            condel[i] = val.getROCStats(labels,score_data.getScores('condel_old'))['auc']
            condel_p[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
            logit[i] = val.getROCStats(labels,score_data.getScores('logit_condel_old'))['auc']
            logit_p[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_old')
            condel[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_new')
            condel_p[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_old')
            logit[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_new')
            logit_p[i] = val.getROCStats(labels,scores)['auc']
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(condel,'o-',color=color_t[11])
    pl.plot(condel_p,'h--',color=color_t[13])
    pl.plot(logit,'x-',color=color_t[12])
    pl.plot(logit_p,'<--',color=color_t[14])
    pl.ylim(0.5,1)
    pl.xlim(-.1,6.1)
    pl.yticks([0.5,0.6,0.7,0.8,0.9,1.0])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(7),['All','Pure',']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    pl.title("predictSNPSelected")
    remove_border()
    
    pl.subplot(313)
    condel = sp.zeros(7)
    condel_p = sp.zeros(7)
    logit = sp.zeros(7)
    logit_p = sp.zeros(7)
    score_data.selectDataset('swissvar_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            condel[i] = val.getROCStats(labels,score_data.getScores('condel_old'))['auc']
            condel_p[i] = val.getROCStats(labels,score_data.getScores('condel_new'))['auc']
            logit[i] = val.getROCStats(labels,score_data.getScores('logit_condel_old'))['auc']
            logit_p[i] = val.getROCStats(labels,score_data.getScores('logit_condel_new'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_old')
            condel[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='condel_new')
            condel_p[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_old')
            logit[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='logit_condel_new')
            logit_p[i] = val.getROCStats(labels,scores)['auc']
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(condel,'o-',color=color_t[11])
    pl.plot(condel_p,'h--',color=color_t[13])
    pl.plot(logit,'x-',color=color_t[12])
    pl.plot(logit_p,'<--',color=color_t[14])
    pl.ylim(0.5,1)
    pl.xlim(-.1,6.1)
    pl.yticks([0.5,0.6,0.7,0.8,0.9,1.0])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(7),['All','Pure',']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    pl.title("SwissVarSelected")
    remove_border()
    
    pl.subplots_adjust(left=0.06,bottom=0.11,right=0.98,top=0.91,wspace=0.03,hspace=0.55)
    
    if release:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS19.pdf'))
    else:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS19.pdf'))
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS19.tiff'),dpi=300)
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS19.jpg'))
    pl.close()
