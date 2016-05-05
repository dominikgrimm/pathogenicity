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
import os

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
    
    data = score_data.getSimilarityData("varibench_selected")
    data1 = score_data.getSimilarityData("predictSNPSelected")
    data2 = score_data.getSimilarityData("swissvar_selected")
    
    data = sp.array(data,dtype="float")*100.0
    data1 = sp.array(data1,dtype="float")*100.0
    data2 = sp.array(data2,dtype="float")*100.0

    x = sp.arange(0,8,1)
    x_labels = sp.arange(0,8,1)
    labels = sp.arange(30,101,10)[::-1]

    width = 0.3
    pl.figure(figsize=(10,6))
    fig = pl.subplot(311)
    pl.bar(x-0.35,data[:,0],color=color_t[0],width=width)
    pl.bar(x+0.05,data[:,1],color=color_t[3],width=width)
    #pl.legend(['Deleterious Variants','Neutral Variants'],fancybox=True,loc='lower right',ncol=2,fontsize=10)
    pl.xticks(x_labels,labels)
    pl.xlim(-1,8)
    pl.xlabel("Protein Similarity in %")
    pl.ylabel("Variants in Similar Proteins %",fontsize=8)
    pl.grid()
    remove_border()   
    pl.title("VariBenchSelected")
    
    ax = fig.get_axes()
    legend = ax.legend(['Deleterious Variants in Similar Proteins','Neutral Variants in Similar Proteins'],frameon=True, 
                        scatterpoints=1,prop={'size':9},ncol=2,
                        loc="upper center",fancybox=True,bbox_to_anchor=(0.5, 1.6))
    legend.get_frame().set_alpha(0.5)
    rect = legend.get_frame()
    rect.set_linewidth(0.0)

    pl.subplot(312)
    pl.bar(x-0.35,data1[:,0],color=color_t[0],width=width)
    pl.bar(x+0.05,data1[:,1],color=color_t[3],width=width)
    pl.xticks(x_labels,labels)
    pl.xlim(-1,8)
    pl.xlabel("Protein Similarity in %")
    pl.ylabel("Variants in Similar Proteins %",fontsize=8)
    pl.grid()
    remove_border()   
    pl.title("predictSNPSelected")

    pl.subplot(313)
    pl.bar(x-0.35,data2[:,0],color=color_t[0],width=width)
    pl.bar(x+0.05,data2[:,1],color=color_t[3],width=width)
    pl.xticks(x_labels,labels)
    pl.xlim(-1,8)
    pl.xlabel("Protein Similarity in %")
    pl.ylabel("Variants in Similar Proteins %",fontsize=8)
    pl.grid()
    pl.subplots_adjust(left=0.06,bottom=0.11,right=0.98,top=0.91,wspace=0.03,hspace=0.73)
    remove_border()   
    pl.title("SwissVarSelected")

    if release:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS13.pdf'))
    else:
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS13.pdf'))
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS13.tiff'),dpi=300)
        pl.savefig(os.path.abspath('Output/Supplementary/FigureS13.jpg'))
    pl.close()
