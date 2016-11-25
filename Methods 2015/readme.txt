Four datasets : Enzyme.mat, IC.mat, GPCR.mat and NR.mat
Main Function : PredictingInteractionForS2orS3.m

Each mat file contains the following variables:

DTI --> drug-target interaction matrix
d_s --> chemical-structure-based similarity matrix
t_s --> sequence-alignment-based similarity matrix
d_ATC--> ATC-based similarity matrix
t_s_Class--> FC-based similarity matrix
DrugName--> all names of drugs corresponding to the rows of DTI
TargetName--> all names of targets corresponding to the columns of DTI

NOTE: Users should run the main function in MATLAB . 