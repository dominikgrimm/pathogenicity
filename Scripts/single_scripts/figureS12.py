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
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'black'

    color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D','#FFF79A','#998675','#A4A4A4']
    
    val = utils.Validation()
    categories = ['cat1','cat2','cat3','cat4','cat5']

    n = len(categories)
    mt = sp.zeros(n)
    pp2 = sp.zeros(n)
    mass = sp.zeros(n)
    cadd = sp.zeros(n)
    sift = sp.zeros(n)
    lrt = sp.zeros(n)
    fathmmu = sp.zeros(n)
    fathmmw = sp.zeros(n)
    gerp = sp.zeros(n)
    phylop = sp.zeros(n)
    
    score_data.selectDataset('varibench_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            mt[i] = val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc']
            pp2[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
            mass[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
            cadd[i] = val.getROCStats(labels,score_data.getScores('CADD'))['auc']
            sift[i] = val.getROCStats(labels,score_data.getScores('sift'))['auc']
            fathmmu[i] = val.getROCStats(labels,score_data.getScores('fathmm_u'))['auc']
            fathmmw[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
            gerp[i] = val.getROCStats(labels,score_data.getScores('gerp++'))['auc']
            phylop[i] = val.getROCStats(labels,score_data.getScores('phylop'))['auc']
            lrt[i] = val.getROCStats(labels,score_data.getScores('lrt'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationtaster')
            mt[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='polyphen2')
            pp2[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationassessor')
            mass[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='CADD')
            cadd[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='sift')
            sift[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_u')
            fathmmu[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_w')
            fathmmw[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='gerp++')
            gerp[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='phylop')
            phylop[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='lrt')
            lrt[i] = val.getROCStats(labels,scores)['auc']
    
    pl.figure(figsize=(12,8))
    fig = pl.subplot(311)
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(fathmmw,'o-',color=color_t[0])
    pl.plot(fathmmu,'h-',color=color_t[1])
    pl.plot(mt,'x-',color=color_t[2])
    pl.plot(mass,'<-',color=color_t[8])
    pl.plot(pp2,'>-',color=color_t[3])
    pl.plot(cadd,'D-',color=color_t[9])
    pl.plot(sift,'1-',color=color_t[4])
    pl.plot(lrt,'2-',color=color_t[5])
    pl.plot(gerp,'3-',color=color_t[6])
    pl.plot(phylop,'4-',color=color_t[7])
    pl.ylim(0.53,0.7)
    pl.xlim(-.1,4.1)
    pl.yticks([0.53,0.6,0.7])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(5),[']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    leg = pl.legend(['FatHMM-W','FatHMM-U','MT2','MASS','PP2','CADD','SIFT','LRT','GERP++','pyhloP'],loc='upper right',fancybox=True, ncol=5,prop={'size':10},numpoints=1)
    leg.get_frame().set_alpha(0.2)
    leg.get_frame().set_edgecolor("none")
    fig.text(-0.05,1.02,"a",fontsize=15,fontweight="bold",va="top",transform=fig.transAxes)
    pl.title("VariBenchSelected")
    remove_border()
    
    mt = sp.zeros(n)
    pp2 = sp.zeros(n)
    mass = sp.zeros(n)
    cadd = sp.zeros(n)
    sift = sp.zeros(n)
    lrt = sp.zeros(n)
    fathmmu = sp.zeros(n)
    fathmmw = sp.zeros(n)
    gerp = sp.zeros(n)
    phylop = sp.zeros(n)
    score_data.selectDataset('predictSNP_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            mt[i] = val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc']
            pp2[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
            mass[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
            cadd[i] = val.getROCStats(labels,score_data.getScores('CADD'))['auc']
            sift[i] = val.getROCStats(labels,score_data.getScores('sift'))['auc']
            fathmmu[i] = val.getROCStats(labels,score_data.getScores('fathmm_u'))['auc']
            fathmmw[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
            gerp[i] = val.getROCStats(labels,score_data.getScores('gerp++'))['auc']
            phylop[i] = val.getROCStats(labels,score_data.getScores('phylop'))['auc']
            lrt[i] = val.getROCStats(labels,score_data.getScores('lrt'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationtaster')
            mt[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='polyphen2')
            pp2[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationassessor')
            mass[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='CADD')
            cadd[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='sift')
            sift[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_u')
            fathmmu[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_w')
            fathmmw[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='gerp++')
            gerp[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='phylop')
            phylop[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='lrt')
            lrt[i] = val.getROCStats(labels,scores)['auc']
    fig = pl.subplot(312)
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(fathmmw,'o-',color=color_t[0])
    pl.plot(fathmmu,'h-',color=color_t[1])
    pl.plot(mt,'x-',color=color_t[2])
    pl.plot(mass,'<-',color=color_t[8])
    pl.plot(pp2,'>-',color=color_t[3])
    pl.plot(cadd,'D-',color=color_t[9])
    pl.plot(sift,'1-',color=color_t[4])
    pl.plot(lrt,'2-',color=color_t[5])
    pl.plot(gerp,'3-',color=color_t[6])
    pl.plot(phylop,'4-',color=color_t[7])
    pl.ylim(0.55,0.8)
    pl.xlim(-.1,4.1)
    pl.yticks([0.55,0.6,0.7,0.8])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(5),[']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    fig.text(-0.05,1.02,"b",fontsize=15,fontweight="bold",va="top",transform=fig.transAxes)
    pl.title("predictSNPSelected")
    remove_border()
    
    mt = sp.zeros(n)
    pp2 = sp.zeros(n)
    mass = sp.zeros(n)
    cadd = sp.zeros(n)
    sift = sp.zeros(n)
    lrt = sp.zeros(n)
    fathmmu = sp.zeros(n)
    fathmmw = sp.zeros(n)
    gerp = sp.zeros(n)
    phylop = sp.zeros(n)
    score_data.selectDataset('swissvar_selected')
    score_data.loadCategories()
    for i,cat in enumerate(categories):
        if cat=="all":
            labels = score_data.getTrueLabels()
            mt[i] = val.getROCStats(labels,score_data.getScores('mutationtaster'))['auc']
            pp2[i] = val.getROCStats(labels,score_data.getScores('polyphen2'))['auc']
            mass[i] = val.getROCStats(labels,score_data.getScores('mutationassessor'))['auc']
            cadd[i] = val.getROCStats(labels,score_data.getScores('CADD'))['auc']
            sift[i] = val.getROCStats(labels,score_data.getScores('sift'))['auc']
            fathmmu[i] = val.getROCStats(labels,score_data.getScores('fathmm_u'))['auc']
            fathmmw[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
            gerp[i] = val.getROCStats(labels,score_data.getScores('gerp++'))['auc']
            phylop[i] = val.getROCStats(labels,score_data.getScores('phylop'))['auc']
            lrt[i] = val.getROCStats(labels,score_data.getScores('lrt'))['auc']
        else:
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationtaster')
            mt[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='polyphen2')
            pp2[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='mutationassessor')
            mass[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='CADD')
            cadd[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='sift')
            sift[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_u')
            fathmmu[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_w')
            fathmmw[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='gerp++')
            gerp[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='phylop')
            phylop[i] = val.getROCStats(labels,scores)['auc']
            [labels,scores] = score_data.getData4Categorie(category=cat,tool_name='lrt')
            lrt[i] = val.getROCStats(labels,scores)['auc']
    fig = pl.subplot(313)
    font = FontProperties()
    font.set_weight('bold')
    pl.plot(fathmmw,'o-',color=color_t[0])
    pl.plot(fathmmu,'h-',color=color_t[1])
    pl.plot(mt,'x-',color=color_t[2])
    pl.plot(mass,'<-',color=color_t[8])
    pl.plot(pp2,'>-',color=color_t[3])
    pl.plot(cadd,'D-',color=color_t[9])
    pl.plot(sift,'1-',color=color_t[4])
    pl.plot(lrt,'2-',color=color_t[5])
    pl.plot(gerp,'3-',color=color_t[6])
    pl.plot(phylop,'4-',color=color_t[7])
    pl.ylim(0.55,0.73)
    pl.xlim(-.1,4.1)
    pl.yticks([0.55,0.6,0.7,0.75])
    pl.grid(axis='y')
    pl.ylabel("AUC")
    pl.xticks(sp.arange(5),[']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]',],fontsize=font_size,rotation=90)
    remove_border()
    fig.text(-0.05,1.02,"c",fontsize=15,fontweight="bold",va="top",transform=fig.transAxes)
    pl.title("SwissVarSelected")

    pl.subplots_adjust(left=0.06,bottom=0.11,right=0.98,top=0.91,wspace=0.03,hspace=0.5)
    
    if release:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS12.pdf'))
    else:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS12.pdf'))
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS12.tiff'),dpi=300)
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS12.jpg'))
    pl.close()
