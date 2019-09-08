#######################
# IMPORT LIBRARIES    #
#######################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from matplotlib import colors
import numpy as np
import operator
from sklearn.kernel_ridge import KernelRidge

#######################
# IMPORT DATASETS     #
#######################

GO=pd.read_csv('../data/GO_Processes_In_Resistance.csv')
PF_genes=pd.read_csv('../data/PFalciparum_Gene_GOterms.csv')
train=pd.read_csv('../data/SubCh1_TrainingData.csv')
test=pd.read_csv('../data/flat_test.csv',index_col=0)
IC50=pd.read_csv('../data/flat_IC50.csv',index_col=0)
genes_column_labels_df=pd.read_csv('../data/genes_column_labels.csv')
isolates_test=[line.rstrip('\n') for line in open('../data/isolates_test.txt')]


#############################
# 1. FIND IMPORTANT GENES   #
#############################

# Identify GO processes involved in artemisinin resistance

Gene_to_Go_Process_df=pd.DataFrame({'Gene_ID':PF_genes['Gene ID'],\
                                               'Curated_Go_Process_IDs':PF_genes['Curated GO Process IDs']})

Gene_to_Go_Process_df.dropna(inplace=True)
Gene_to_Go_Process_df['Curated_Go_Process_IDs']=\
    Gene_to_Go_Process_df['Curated_Go_Process_IDs']\
    .str.split(';', n=-1)

GO_Process_to_Gene_list=[]

for idx, row in Gene_to_Go_Process_df.iterrows():
    for GO_Process in row['Curated_Go_Process_IDs']:
        GO_Process_to_Gene_list.append([GO_Process,row['Gene_ID']])

GO_Process_to_Gene_df=pd.DataFrame(GO_Process_to_Gene_list, 
                                   columns=['Curated_Go_Process_IDs','Gene_ID'])

GO_Process_to_Gene_df.sort_values('Curated_Go_Process_IDs', 
                                 inplace=True)
GO_Process_to_Gene_df.reset_index(drop=True,inplace=True)

GO_Process_to_Gene_df=GO_Process_to_Gene_df\
.groupby('Curated_Go_Process_IDs')['Gene_ID']\
.apply(list).to_frame()

GO_Process_to_Gene_df.reset_index(inplace=True)

Gene_Count=[len(gene_list) for gene_list in list(GO_Process_to_Gene_df["Gene_ID"])]
GO_Process_to_Gene_df['Gene_Count']=Gene_Count

GO=pd.merge(GO,GO_Process_to_Gene_df,left_on='GO Process ID',
             right_on='Curated_Go_Process_IDs')
GO.drop('Curated_Go_Process_IDs', axis=1, inplace=True)

# Explore the distribution of Gene_Coun

plt.hist(GO['Gene_Count'].values,bins=50)
plt.title('Gene_Count distribution: All GO Processes')

plt.hist(GO[GO['IsInvolvedInArtemisinResistance']==0]['Gene_Count'].values,bins=50)
plt.title('Gene_Count distribution: GO Processes not involved in Artemisin Resistance')

plt.hist(GO[GO['IsInvolvedInArtemisinResistance']==1]['Gene_Count'].values,bins=50)
plt.title('Gene_Count distribution: GO Processes involved in Artemisin Resistance')


# Filter the GO dataframe

GO=GO[GO['IsInvolvedInArtemisinResistance']==1].reset_index(drop=True)
GO=GO[GO['Gene_Count']<=25].reset_index(drop=True)

# Create a list of unique genes involved in artemisinin resistance

result = set()
for row in range(len(GO)):
    result = result | set(GO["Gene_ID"][row])
result = list(result)
print("Total of genes including after filtering: {}".format(len(result)))


##########################################
# 2. FILTER AND FLATTEN TRAIN DATAFRAME  #
##########################################

# Create list of important genes
important_genes = set(list(train.columns)) & set(result)
important_genes = list(important_genes)

# Filter train dataframe to keep only important genes
columns_to_keep = list(train.columns)[0:5] + \
[list(train.columns)[-1]] + important_genes
filtered_train = train[columns_to_keep]
filtered_train.head(10)

# Aggregate and calculate mean across bioreplicates
aggregated_means_train = filtered_train.groupby(["Isolate","Timepoint","Treatment"]).mean()
aggregated_means_train.head(10)

# Generate a gene dictionary with the following format
# key: geneX
# values: 'geneX_24HR_DHA','geneX_24HR_UT','geneX_6HR_DHA','geneX_6HR_UT'

gene_dict={}
isolate_comb_list=['_24HR_DHA','_24HR_UT','_6HR_DHA','_6HR_UT']
for gene in important_genes:
    gene_dict[gene]=[gene+comb for comb in isolate_comb_list]

# Create a list of isolates
isolates_list=aggregated_means_train.index.get_level_values(0)
isolates_list=list(dict.fromkeys(isolates_list))

# Flatten train data frame
flat_train=pd.DataFrame()
genes=[]
for gene in important_genes:
    temp_list=[]
    for isolate in isolates_list:
        temp_list.append(aggregated_means_train.loc[isolate][gene].tolist())
    flat_train=pd.concat([flat_train,pd.DataFrame(data=temp_list,columns=gene_dict[gene])],axis=1)

flat_train.index=isolates_list
flat_train.head(10)

####################
# TRAIN KRR MODELS #
####################

# Select validation isolates
val_isolates=['isolate_07','isolate_11','isolate_16','isolate_18','isolate_27']

# Helper functions

def getGenesColumnLabels(genes):
#Takes a list of genes and returns the genes labels
    genes_list=[]
    for gene in genes:
        genes_list=genes_list+list(genes_column_labels_df[gene][0:4])
    return genes_list

def findTrainIsolates(val_isolates,df):
#Takes a list of isolates chosen to be the validation set and the train df.
#Returns the isolates not in val_isolates
    train_isolates=list(set(df.index)-set(val_isolates))
    train_isolates.sort()
    return train_isolates

def filterDfByIsolates(isolates,df):
#Takes a list of isolates and a dataframe
#Returns the rows of the dataframe whose index are in isolates
    filtered_df=df.loc[isolates]
    return filtered_df

def filterDfByGenes(genes,df):
#Takes a list of genes e.g.['PF3D7_1360200','PF3D7_1212800'] and a dataframe
#Returns the columns of the dataframe whose labels are in genes
    genes_to_filter=getGenesColumnLabels(genes)
    filtered_df=df[genes_to_filter]
    return filtered_df

def trainModelKRR(genes,val_isolates,dfx,dfy):
# Takes a list of genes to include in the model, a list of validation isolates, a X dataframe and a y dataframe
# Returns the score of the model. The score is R^2. 
    train_isolates=findTrainIsolates(val_isolates,dfx)
    
    X=filterDfByIsolates(train_isolates,dfx)
    X=filterDfByGenes(genes,X)
    
    y=filterDfByIsolates(train_isolates,dfy)
    
    model=KernelRidge(alpha=1)
    model=model.fit(X.values,y.values)
    
    X_val=filterDfByIsolates(val_isolates,dfx)
    X_val=filterDfByGenes(genes,X_val)
    
    y_val=filterDfByIsolates(val_isolates,dfy).values
    
    score=model.score(X_val,y_val)
    
    return score

# Train model multiple times to find the genes that generate the best model's scores

selected_genes=[]
selected_scores=[0]
N=100
for i in range(N):
    best_gene_so_far=''
    best_score_so_far=selected_scores[-1]
    
    for gene in important_genes:
        
        if gene not in selected_genes:
            genes_to_use=selected_genes+[gene]
            
            model_score=trainModelKRR(genes_to_use,val_isolates,flat_train,IC50)
            
            if model_score > best_score_so_far:
                best_score_so_far=model_score
                best_gene_so_far=gene

    
    if best_score_so_far>selected_scores[-1]:
        selected_genes.append(best_gene_so_far)
        selected_scores.append(best_score_so_far)

# Plot gene vs r2
# Note each y point is an accumulative score i.e. score calculated with the current + the previous found gene

x=np.arange(len(selected_scores))
y=selected_scores
plt.plot(x,y)
plt.xticks(x,['']+selected_genes,rotation='vertical')
plt.show

# Keep only the necesary genes
selected_genes=selected_genes[0:4]

# Train model with selected genes

def tryModel(genes,val_isolates,dfx,dfy):
    model_info={}
    
    train_isolates=findTrainIsolates(val_isolates,dfx)
    
    X=filterDfByIsolates(train_isolates,dfx)
    X=filterDfByGenes(genes,X)
    
    y=filterDfByIsolates(train_isolates,dfy)
    
    model=KernelRidge(alpha=1)
    model=model.fit(X.values,y.values)
    
    X_val=filterDfByIsolates(val_isolates,dfx)
    X_val=filterDfByGenes(genes,X_val)
    
    y_val=filterDfByIsolates(val_isolates,dfy).values
    
    y_pred=model.predict(X_val)
    
    score=model.score(X_val,y_val)
    
    model_info['y']=y.values
    model_info['y_val']=y_val
    model_info['X']=X.values
    model_info['X_val']=X_val
    model_info['y_pred']=y_pred
    model_info['score']=score
    return model_info


try_model=tryModel(selected_genes,val_isolates,flat_train,IC50)
y=[try_model['y_val'][i][0] for i in range(len(try_model['y_val']))]
y_pred=try_model['y_pred']

plt.scatter(y,y_pred)


# Train final model with all selected genes and all isolates

X=flat_train[getGenesColumnLabels(selected_genes)]
X_test=test[getGenesColumnLabels(selected_genes)]
y=IC50

model=KernelRidge(alpha=1)
model=model.fit(X.values,y.values)

#######################
# MAKE PREDICTIONS    #
#######################

IC50_pred=model.predict(X_test)

#######################
# WRITE RESULTS       #
#######################

temp_list=[]
for value in IC50_pred:
    temp_list.append(value[0])
IC50_pred=temp_list
submission=pd.DataFrame({'Isolate':isolates_test,'Predicted_IC50':IC50_pred})

submission.to_csv('../submissions/lylat_SubCh1_model5_150819.csv')
