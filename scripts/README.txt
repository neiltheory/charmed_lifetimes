A guide to using the scripts in this folder. All pathways will need to be reconsidered with regards to input/output files...

=======PREPARING A FILE FOR USE WITH TMVA========

1) Make a ‘reduced’ file from the full data file turbo_2015_data.root. This will contain the events within a certain TAU and M range as defined in mkDataSetFromTTree.C.

2) Using sWeights.C make a fit to the reduced data and output a file with the calculated sWeights

3) Using the python script addAngularVars_neil.py, create angular variables for each entry in the full turbo_2015_data.root 

4) Using mergetrees.py merge the newly created angular variables with the full turbo_2015_data.root file using the command: python mergetrees.py <output file name> <input file name 1> <input TTreename 1> <input file name 2> <input TTree name 2>

5) Using TTree_stripper.C, 'strip' the merged file from the process in the previous step so that it only contain the entries within the same ranges that were used to create the sWeights

6) Again using mergetrees.py, merge the ‘stripped’ full data with the sWeights to produce a file that has all the data for all the variables (including the new add angular variables) and all the sWeights. Each entry in this file will now lie in the ranges that were used to make the fit and obtain the sWeights. The output file will now be suitable for use with the TMVA package.


=======TMVA========

7) run the TMVA package using the BDTG functionality with any variables (in this case we have included the angular variables which we went to all the effort to make in step 3). It will output a file (TMVA.root) and an .xml file (the .xml file can be used to assign BDTG weights for all candidates in the original data file, even the ones outside the ranges used in the training in the TMVA package.)

8) using the TMVA GUI one can (using button 5a) find the optimal cut value of the BDTG weight by inputting the background and signal yeilds as calculated in the original fit (these values were outputted from the script sWeights.C, that was run in step 2).

=======REMOVE THE BACKGROUND======
9) using the TMVA.root file (which contains the original data and the newly calculated BDTG weights) as an input, use the lifetimebinning.C script to remove some background components of the data by making the appropriate cut (as calculated in step 8) on the BDTG variable, then divide the remaning candidates into equal sized bins in the TAU variable and, for each bin, perform the same background+signal fit to ascertain the quantity of signal candidates and remaining background candidates in each bin. Finaly the lifetimebinning.C script will populate a new histogram with the signal yeild number for each bin and save it.

======CALCULATE LIFETIME AND ERROR=====
10) Use the output from the lifetimebinning.c script as the input to the ExpFitTAU_SigOnly and create a blinded fit of the exponential decay of the TAU variable. The blinding constant and it's error are ouput to screen. 1/(blinding constant) gives the (blinded) lifetime and the (blinding error)/((blinding constant)^2) is the error in the lifetime.

=====CALCULATE LIFETIME AND ERROR FOR OTHER RESONANCES====
(i.e. Calculate lifetime and error for other particles in the data that may lie outside the mass range of the data originally used to train the TMVA package to produce the BDTG weights.)

9) If another decay resonance is known to happen, say, one tenth of the time that the originally considered decay happens, then it is possible to calculate BDTG weights for all candidates in the data using the BDTG weights (.xml) output file. Use addmva.py to calculate the weights for all the candidates, even ones not included in the TMVA training range. An example of how to use it is given in testmvamerge.sh. This will produce an ntuple containing the BDTG weight for each candidate, which can then be merged with the original full data file using the mergetrees.py script (as was done in steps 4 & 6).

10) Find the optimal cut value of the BDTG weight by inputting, say, 10% of the signal yeild in the TMVA GUI (button 5a) when running the GUI with the TMVA.root file produced earlier, this will give a crude but workable optimal cut value for the expeted reduced sigal.

11) using the lifetimebinning_Xi.c script, define the range the new (mass) range of interest and produce an output histogram of the signal candidates by the same process as that used in step 9.

12) use ExpFitTAU_SigOnly to make a blinded fit to the new TAU variable in the new particles mass range.
