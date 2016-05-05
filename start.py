'''
################################################################################################################
Manuscript Title: The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity
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

#Check if all packages are installed
import imp 

modules = ['scipy','numpy','pylab','matplotlib','sklearn']
missing = False
for module in modules:
	try:
		imp.find_module(module)
	except ImportError:
		missing = True
		print
		print "[ERROR] Package " + module + " is missing. Please check dependencies.txt and install packages (e.g. pip)"
		print
if missing==True:
	quit()

#set path
import sys,os
sys.path.append(os.path.join(os.getcwd(),"Scripts"))

import Scripts.main as main

if __name__ in "__main__":
	release = True
	main.generateFiguresAndTables(release=release)
