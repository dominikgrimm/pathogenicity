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
	
color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D','#FFF79A','#998675','#A4A4A4']
font_size = 11

hatch = pl.Rectangle((0,0),1,1,fill=None,hatch="///")
circle = pl.Rectangle((0,0),1,1,fill=None,hatch="ooo")

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

def plotBar(data=None,color_id=None,figure_id=None,name=None,flag=False):
	ax = pl.subplot(figure_id)
	width = 0.8
	x=sp.arange(7)
	if not (name=="VaribenchSelected"):
		pl.bar(x-0.4,data,width=width,color=color_t[color_id],hatch="/o/o/")
	else:
		pl.bar(x-0.4,data,width=width,color=color_t[color_id],hatch="ooo")

	tmp = data.copy()
	tmp[1::] = 0
	pl.xticks(x,['All','Pure',']0.0,1.0[','[0.1,0.9]','[0.2,0.8]','[0.3,0.7]','[0.4,0.6]'],fontsize=font_size,rotation=90)
	ln = sp.log10(len(name))
	pl.text(3.5-ln,0.95,name)
	if flag:
		remove_border(left=False)
		pl.yticks([0.5,0.6,0.7,0.8,0.9,1.0])
		pl.grid(axis='y')
		pl.tick_params(axis='y',which="both",labelleft='off',left='off')
	else:
		pl.ylabel("AUC")
		remove_border()
		pl.yticks([0.5,0.6,0.7,0.8,0.9,1.0])
		pl.grid(axis='y')
	pl.ylim(0.5,1)
	pl.xlim(-0.5,7.5)
	return ax


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
	mpl.rcParams['font.size'] = font_size
	mpl.rcParams['patch.edgecolor'] = 'black'
	
	val = utils.Validation()
	categories = ['all','pure','cat1','cat2','cat3','cat4','cat5']

	humvar = sp.zeros(7)
	exovar = sp.zeros(7)
	varibench_selected = sp.zeros(7)
	predictSNP_selected = sp.zeros(7)
	swissvar_selected = sp.zeros(7)
	
	datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
	for dataset in datasets:
		score_data.selectDataset(dataset)
		score_data.loadCategories()
		for i,cat in enumerate(categories):
			if cat=="all":
				labels = score_data.getTrueLabels()
				if dataset=='humvar':
					humvar[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
				elif dataset=='exovar':
					exovar[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
				elif dataset=='swissvar_selected':
					swissvar_selected[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
				elif dataset=='predictSNP_selected':
					predictSNP_selected[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
				else:
					varibench_selected[i] = val.getROCStats(labels,score_data.getScores('fathmm_w'))['auc']
			else:
				[labels,scores] = score_data.getData4Categorie(category=cat,tool_name='fathmm_w')
				if dataset=='humvar':
					humvar[i] = val.getROCStats(labels,scores)['auc']
				elif dataset=='exovar':
					exovar[i] = val.getROCStats(labels,scores)['auc']
				elif dataset=='predictSNP_selected':
					predictSNP_selected[i] = val.getROCStats(labels,scores)['auc']
				elif dataset=='swissvar_selected':
					swissvar_selected[i] = val.getROCStats(labels,scores)['auc']
				else:
					varibench_selected[i] = val.getROCStats(labels,scores)['auc']
	pl.figure(figsize=(10,7))
	font = FontProperties()
	font.set_weight('bold')
	
	plotBar(humvar,0,321,"HumVar",flag=False)
	ax = plotBar(exovar,1,322,"ExoVar",flag=True)
	plotBar(varibench_selected,3,323,"VariBenchSelected",flag=False)
	plotBar(predictSNP_selected,4,324,"predictSNPSelected",flag=True)
	plotBar(swissvar_selected,5,325,"SwissVarSelected",flag=False)

	rect = pl.Rectangle((0,0),1,1,fill=None)
	
	leg = ax.legend([rect,hatch,circle],['FatHMM-W','Type 1 Biased','Type 2 Biased'],loc='upper center', bbox_to_anchor=(0.0, 1.15),fancybox=True, ncol=5,prop={'size':10},numpoints=1)
	leg.get_frame().set_alpha(0.2)
	leg.get_frame().set_edgecolor("none")
	
	pl.subplots_adjust(left=0.065,bottom=0.11,right=0.99,top=0.94,wspace=0.03,hspace=0.5)
	
	if release:
		pl.savefig(os.path.abspath('Output/Supplementary/FigureS11.pdf'))
	else:
		pl.savefig(os.path.abspath('Output/Supplementary/FigureS11.pdf'))
		pl.savefig(os.path.abspath('Output/Supplementary/FigureS11.tiff'),dpi=300)
		pl.savefig(os.path.abspath('Output/Supplementary/FigureS11.jpg'))
	pl.close()
