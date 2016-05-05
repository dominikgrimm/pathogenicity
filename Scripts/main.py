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


#import all necessary packages
import utils
import single_scripts.figure1 as figure1
import single_scripts.figure2 as figure2
import single_scripts.figure3 as figure3
import single_scripts.figure4 as figure4
import single_scripts.figureS2_6 as figureS2_6
import single_scripts.figureS7_10 as figureS7_10
import single_scripts.figureS11 as figureS11
import single_scripts.figureS12 as figureS12
import single_scripts.figureS13 as figureS13
import single_scripts.figureS14_18 as figureS14_18
import single_scripts.figureS19 as figureS19
import single_scripts.tableS1 as tableS1
import single_scripts.tableS2 as tableS2

def generateFiguresAndTables(release=False):
    print
    #load data
    print "*************************************************************************************************************"
    print "*Evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity*"
    print "*************************************************************************************************************"
    print
    print "Script to reproduce all figures from this manuscript"
    print
    print "Loading Data..."
    score_data = utils.DataScores()
    score_data.load_data()

    print "Creating Figure1..."
    figure1.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Figures S2-S6..."
    figureS2_6.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Table S1..."
    tableS1.printTableS1(score_data)
    
    print "Creating Supplementary Table S2..."
    tableS2.printTableS2(score_data)
    
    print "Creating Figure2..."
    figure2.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Figures S7-S10..."
    figureS7_10.plotFigure(score_data,release=release)
    
    print "Creating Figure3..."
    figure3.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Figure S11..."
    figureS11.plotFigure(score_data,release=release)
    
    print "Creating FigureS12..."
    figureS12.plotFigure(score_data,release=release)
    
    print "Creating FigureS13..."
    figureS13.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Figures S14-S18..."
    figureS14_18.plotFigure(score_data,release=release)
    
    print "Creating Figure4..."
    figure4.plotFigure(score_data,release=release)
    
    print "Creating Supplementary Figure S19..."
    figureS19.plotFigure(score_data,release=release)
    
    print "Done"
    print 
