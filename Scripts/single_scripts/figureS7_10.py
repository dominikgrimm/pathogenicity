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
import matplotlib as mpl
import numpy as sp
import os

def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    ax = axes or pl.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    
    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    
    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()

def plot(score_data=None,did=None,release=True):
    pl.ion()
    
    font_size = 11
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['font.size'] = font_size
    mpl.rcParams['patch.edgecolor'] = 'white'

    labels = ['Pathogenic Only Proteins','Neutral Only Proteins','Mixed Proteins']
    [total,n_mixed,n_del,n_neutral] = score_data.getProteinFractions()
    #print n_del, n_neutral, n_mixed, total
    
    total_protein = total

    
    fracs = sp.array([float(n_del)/float(total),float(n_neutral)/float(total),float(n_mixed)/float(total)])
    pl.figure(figsize=(14,5))
    fig = pl.subplot(221,aspect=True)
    patches, texts, autotexts = fig.pie(fracs,labels=labels,autopct="%1.1f%%",shadow=True,colors = ['#DBDBDB','#888888','#AAAAAA'],explode=(0,0,0.1),radius=1)
    for i,text in enumerate(texts):
        text.set_fontsize(9)
        autotexts[i].set_fontsize(9)
    fig.text(-0.7,1.05,"a",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
    
    fig = pl.subplot(223,aspect=True)
    labels = ['Variants in Pathogenic Only Protein','Variants in Neutral Only Proteins','Variants in Mixed Proteins']
    
    [total,n_mixed,n_del,n_neutral] = score_data.getVariantFractions()
    #print n_del, n_neutral, n_mixed, total
    total_variant = total
    
    fracs = sp.array([float(n_del)/float(total),float(n_neutral)/float(total),float(n_mixed)/float(total)])
    
    patches, texts, autotexts = fig.pie(fracs,labels=labels,autopct="%1.1f%%",shadow=True,colors = ['#DBDBDB','#888888','#AAAAAA'],explode=(0,0,0.1),radius=1)
    for i,text in enumerate(texts):
        text.set_fontsize(9)
        autotexts[i].set_fontsize(9)
    fig.text(-0.7,1.05,"c",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
    
    cat1 = score_data.getProteinsCategorie(category='cat1')
    cat2 = score_data.getProteinsCategorie(category='cat2')
    cat3 = score_data.getProteinsCategorie(category='cat3')
    cat4 = score_data.getProteinsCategorie(category='cat4')
    cat5 = score_data.getProteinsCategorie(category='cat5')
    data = [sp.unique(cat1).shape[0],sp.unique(cat2).shape[0],sp.unique(cat3).shape[0],sp.unique(cat4).shape[0],sp.unique(cat5).shape[0]]
    
    fig = pl.subplot(222)
    pl.bar(sp.arange(5)-0.4,data,width=0.8,color='#AAAAAA',alpha=0.6)
    pl.grid(axis='y', color='white', linestyle='-', lw=1)
    #remove_border(left=False)
    remove_border()
    ratio = sp.array(data)/float(total_protein)
    for x, y, r in zip(sp.arange(5), data, ratio):
        pl.annotate("%1.1f%%" % (r*100.0), (x, y + 1), ha='center')
    pl.ylabel("#Mixed Proteins")
    pl.xticks(sp.arange(5),[']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]'])
    fig.text(-0.2,1.05,"b",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
    
    cat1 = score_data.getVariantsCategorie(category='cat1')
    cat2 = score_data.getVariantsCategorie(category='cat2')
    cat3 = score_data.getVariantsCategorie(category='cat3')
    cat4 = score_data.getVariantsCategorie(category='cat4')
    cat5 = score_data.getVariantsCategorie(category='cat5')
    data = [cat1,cat2,cat3,cat4,cat5]
    fig = pl.subplot(224)
    pl.bar(sp.arange(5)-0.4,data,width=0.8,color="#AAAAAA",alpha=0.6)
    pl.grid(axis='y', color='white', linestyle='-', lw=1)
    remove_border()
    ratio = sp.array(data)/float(total_variant)
    for x, y, r in zip(sp.arange(5), data, ratio):
        pl.annotate("%1.1f%%" % (r*100.0), (x, y + 1), ha='center')
    pl.ylabel("#Variants in Mixed Proteins")
    pl.xticks(sp.arange(5),[']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]'])
    fig.text(-0.2,1.05,"d",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
    
    pl.subplots_adjust(left=0.08,bottom=0.06,right=0.99,top=0.9,wspace=0.45)
    
    if release:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS' + str(did) + ".pdf"))
    else:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS' + str(did) + ".pdf"))
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS' + str(did) + ".tiff"),dpi=300)
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS' + str(did) + ".jpg"))
    pl.close()


def plotFigure(score_data=None,release=True):
    dataset_map = {'humvar':7,
           'exovar':8,
           'swissvar_selected':10,
           'predictSNP_selected':9,}
    datasets = ['humvar','exovar','predictSNP_selected','swissvar_selected']
    for i,dataset in enumerate(datasets):
        print "\tComputing Categories and plotting Supplementary Fig. S" + str(dataset_map[dataset])
        score_data.selectDataset(dataset)
        score_data.loadCategories()
        plot(score_data,did=dataset_map[dataset],release=release)

