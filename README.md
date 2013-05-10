Lake Analyzer v3.3.1
====

Lake Analyzer is a MATLAB code package for calculating physical indices of water column stability from temperature and (if applicable) wind measurements. 
Citing Lake Analyzer
-------
The package has a user manual and a publication: Read JS, DP Hamilton, ID Jones, K Muraoka, LA Winslow, R Kroiss, CH Wu, E Gaiser. 2011. Derivation of lake mixing and stratification indices from high-resolution lake buoy data. Environmental Modelling and Software. 26: 1325-1336.

Contributing to Lake Analyzer
--------
**Lake Analyzer is in active development**, and we are hoping to add flux estimates to the code, including evaporation and the gas transfer velocity. Users of the code base are encouraged to fork this code, modify it, and send a pull request to integrate changes made to the algorithms. At this point, the MATLAB Lake Analyzer code and web portal (see below) will continue to focus on physical derivatives. See http://lakeanalyzer.gleon.org/ for the web portal for these tools (MATLAB not required)

In order to contribute to this code, we recommend the following workflow: 

1) *fork* this repository to your own personal github account

2) *clone* the github repository to your computer:

	$git clone https://github.com/username/Lake-Analyzer.git

3) modify code or add new functionality, save the code

4) add the repository master to a remote master called *upstream*

	$cd Lake-Analyzer

	$git remote add upstream https://github.com/GLEON/Lake-Analyzer.git

5) before pushing your changes to your repository, pull in the current version of the GLEON master:

	$git fetch upstream

6) merge these differences with your own *master* version:

	$git merge upstream/master

7) push your changes to your github repository, in addition to changes made by pulling in the GLEON master:

	$git push

8) submit a pull request to GLEON master using your account at github.com