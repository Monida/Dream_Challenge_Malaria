README_Subch1.txt 

The purpose of this README file is to describe
a) the files included in this folder
b) the lay out of the data and define each variable

########################
# a) DATASETS INCLUDED #
########################

The files listed in this folder include the following

SubCh1_TrainingData.csv
SubCh1_TestData.csv

#################################
# b) LAYOUT/DESCRIPTION OF DATA #
#################################

SubCh1_TrainingData.csv and SubCh1_TestData.csv are the training and test data sets for the DREAM1 subchallenge 1 datasets, respectively. The training data set will be used to build a model to predict the IC50 of the test data sets. The training data set consists of 30 parasite isolates (or lines) whilst the test set consists of 25 parasite isolates. 

Sample_Names = Unique identifier for each sample

DHA_IC50 = Dihydroartemisinin IC50 (nM). IC50 of DHA, the metabolically active form of artemisinin. Located at the very last column of the dataset.

Isolate = Parasite sample name. Each isolate is a parasite sample isolated from a single patient. 

Timepoint = Estimated time after invasion which transcription sample was collected. Two options are available here 6hr or 24hr. Malaria has a 48hr life cycle in the blood stage and has very different development/transcription/biology at these two timepoints.

Treatment = Whether the transcription sample was either perturbed with 5nM DHA or perturbed with DMSO (our control, listed as UT in dataset).

BioRep = Biological Replicate. Indicates which biological replicate. Options available are Brep1 or Brep2.

The transcription data was collected using a custom Agilent microarray. The transcription data set consists of MAL genes followed by PF3D7 genes with transcription values associated with each gene where available. The MAL genes are 92 non-coding RNAs while the PF3D7 genes are protein coding genes. An in-depth description about the microarray and its contents can be found here https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187595 . We recommend using plasmodb.org/plasmo for descriptions about individual gene function. 
