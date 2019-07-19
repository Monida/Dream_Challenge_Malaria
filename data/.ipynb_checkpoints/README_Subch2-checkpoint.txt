README_Subch2.txt 

The purpose of this README file is to describe
a) the files included in this folder
b) the lay out of the data and define each variable

########################
# a) DATASETS INCLUDED #
########################

Sub-Challenge 2 Datasets

SubCh2_TrainingData.csv 
SubCh2_TestData.csv 

#################################
# b) LAYOUT/DESCRIPTION OF DATA #
#################################

SubCh2_TrainingData.csv and SubCh1_TestData.csv are the training and test data sets for the sub-challenge 2 datasets, respectively. The training data set will be used to build a model to predict the categorical clearance rates of the test data sets. The training data set consists of 1043 in vivo parasite isolates.  This training dataset comes from a published source (http://science.sciencemag.org/content/347/6220/431.long) from the Bozdech lab at Nanyang University in Singapore. Each sample has a unique identifier that starts with GSM that is used to match the transcription data with the metadata. The transcription values for this dataset are all log2 normalized against a co-hybed reference sample in a custom array.
The test data consists of in vitro transcription samples from 32 Southeast Asian parasite isolates. The goal of this sub-challenge is to predict the categorical clearance rate of these 32 isolates. This transcription dataset was collected using a custom Agilent microarray by the Ferdig lab. The transcription data set consists of MAL genes followed by PF3D7 genes with transcription values associated with each gene where available. The MAL genes are 92 non-coding RNAs while the PF3D7 genes are protein coding genes. An in-depth description about the microarray and its contents can be found here https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187595 . We recommend using plasmodb.org/plasmo for descriptions about individual gene function.


ClearanceRate = Categorical classification of DHA resistance measured in patients at the clinic. Slow stands for slow clear lines, or resistant lines. Fast stands for fast clear lines, or sensitive lines. Located in the last column for the dataset. Variable that needs prediction for Sub-challenge 2.

Sample_Name = Unique identifier for each sample.

Country = Country of origin. Where the malaria parasite was isolated from the patient. All samples from the test data were collected from the Thailand-Myanmar border.

Asexual_Stage = The stage of the parasite based on its transcriptional profile. Transcriptional profiles compared against the 3D7 sample from Bozdech et al (2003) PLoS Biol.
Asexual.stage..hpi. = The stage of the parasite based on its transcriptional profile, provided by the authors. For more information, please view the Mok et al (2015) Science paper.

Kmeans.Grp = Transcriptional based clustering of training data. Grp A corresponds to early ring stage parasites, Grp B corresponds to late stage parasites, Grp C consists of a mixture between early stage parasites and gametocytes. For more information, please visit the Mok et al (2015) Science paper.

Isolate = Parasite sample name. Each isolate is a parasite sample isolated from a single patient. 

Timepoint = Estimated time after invasion which transcription sample was collected. Two options are available here 6hr or 24hr. Malaria has a 48hr life cycle in the blood stage and has very different development/transcription/biology at these two timepoints.

Treatment = Whether the transcription sample was either perturbed with 5nM DHA or perturbed with DMSO (our control, listed as UT in dataset).

BioRep = Biological Replicate. Indicates which biological replicate. Options available are Brep1 or Brep2.


