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
import utils,os
import matplotlib as mpl

def plot_roc(results=None,name=None,legend=None,color=None,style='-'):
	auc = results['auc']
	tpr = results['tpr']
	fpr = results['fpr']
	legend.append(name + ": %.2f (AUC)"%auc)
	pl.plot(fpr,tpr,color,linewidth=2.5,linestyle=style)
	return legend

def plot_roc_pr(results=None,name=None,legend=None,color=None,style='-'):
	auc = results['pr_auc']
	tpr = results['recall_ranking']
	fpr = results['precission_ranking']
	legend.append(name + ": %.2f (AUC-PR)"%auc)
	pl.plot(tpr,fpr,color,linewidth=2.5,linestyle=style)
	return legend

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
	mpl.rcParams['patch.edgecolor'] = 'black'
	
	color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D','#FFF79A','#998675','#A4A4A4','#5AFF00']
	
	#load data
	val = utils.Validation()

	datasets = ['humvar','exovar','varibench_selected','predictSNP_selected','swissvar_selected']
	fig_names = ['FigureS2','FigureS3','FigureS4','FigureS5','FigureS6']
	for i,dataset in enumerate(datasets):
		score_data.selectDataset(dataset)
		labels = score_data.getTrueLabels()
		print "\tCreating ROC and ROC-PR curves for dataset:  " + dataset + " (" + fig_names[i] + ")"
		legend = []
		pl.figure(figsize=(12,6))

		fig = pl.subplot(121)
		pl.grid(True)
		fig.set_xlim([0,1])
		fig.set_ylim([0,1])
	
		spines_to_remove = ['top','right']
		ax = fig.get_axes()
		for spine in spines_to_remove:
			ax.spines[spine].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
	
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('fathmm_w')),
				  name="FatHMM-W",
			  	  legend=legend,
			  	  color=color_t[5])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('mutationtaster')),
					name="MutationTaster2",
					legend=legend,
					color=color_t[0])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('polyphen2')),
					name="PolyPhen2",
					legend=legend,
					color=color_t[1])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('mutationassessor')),
					name="MutationAssessor",
					legend=legend,
					color=color_t[8])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('CADD')),
					name="CADD",
					legend=legend,
					color=color_t[9])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('sift')),
					name="Sift",
					legend=legend,
					color=color_t[2])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('lrt')),
					name="LRT",
					legend=legend,
					color=color_t[3])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('fathmm_u')),
					name="FatHMM-U",
					legend=legend,
					color=color_t[4])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('gerp++')),
					name="Gerp++",
					legend=legend,
					color=color_t[6])
		legend = plot_roc(results=val.getROCStats(labels,score_data.getScores('phylop')),
					name="PhyloP",
					legend=legend,
					color=color_t[7])
		
		
		leg = ax.legend(legend,'lower right',numpoints=1,prop={'size':12},fancybox=True)
		leg.get_frame().set_alpha(0.5)
		fig.set_xlabel("False Positive Rate or (1-Specificity)")
		fig.set_ylabel("True Positive Rate or (Sensitivity)")
		fig.plot([0,1],[0,1],'--',color='#ACACAC')
		fig.text(-0.1,1.05,"a",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
		
		fig = pl.subplot(122)
		legend = []
		pl.grid(True)
		fig.set_xlim([0,1])
		fig.set_ylim([0,1])
		
		spines_to_remove = ['top','right']
		ax = fig.get_axes()
		for spine in spines_to_remove:
			ax.spines[spine].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('fathmm_w')),
				  name="FatHMM-W",
			  	  legend=legend,
			  	  color=color_t[5])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('mutationtaster')),
					name="MutationTaster-2",
					legend=legend,
					color=color_t[0])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('polyphen2')),
					name="PolyPhen-2",
					legend=legend,
					color=color_t[1])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('mutationassessor')),
					name="MutationAssessor",
					legend=legend,
					color=color_t[8])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('CADD')),
					name="CADD",
					legend=legend,
					color=color_t[9])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('sift')),
					name="SIFT",
					legend=legend,
					color=color_t[2])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('lrt')),
					name="LRT",
					legend=legend,
					color=color_t[3])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('fathmm_u')),
					name="FatHMM-U",
					legend=legend,
					color=color_t[4])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('gerp++')),
					name="Gerp++",
					legend=legend,
					color=color_t[6])
		legend = plot_roc_pr(results=val.getROCStats(labels,score_data.getScores('phylop')),
					name="phyloP",
					legend=legend,
					color=color_t[7])
		
		
		leg = ax.legend(legend,'lower right',numpoints=1,prop={'size':12},fancybox=True)
		leg.get_frame().set_alpha(0.5)
		fig.set_xlabel("Recall")
		fig.set_ylabel("Precision")
		fig.text(-0.1,1.05,"b",fontsize=14,fontweight='bold',va='top',transform=fig.transAxes)
		
		pl.subplots_adjust(left=0.07,bottom=0.08,right=0.99,top=0.93,wspace=0.15,hspace=0.15)
		
		if release:
			pl.savefig(os.path.abspath('Output/Supplementary/' + fig_names[i] + '.pdf'))
		else:
			pl.savefig(os.path.abspath('Output/Supplementary/' + fig_names[i] + '.pdf'))
			pl.savefig(os.path.abspath('Output/Supplementary/' + fig_names[i] + '.tiff'),dpi=300)
			pl.savefig(os.path.abspath('Output/Supplementary/' + fig_names[i] + '.jpg'))
		pl.close()
