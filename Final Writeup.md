# **Subchallenge 1**

For SubChallenge 1 we used Kernel Ridge Regression to predict the amount of the IC50 values for the isolates. Before, that the important genes, meaning the genes that are biologically related to artemisinin resistant, were identified and used as features.

In summary these are the steps that were taken:

1. Find important genes
2. Filter and Flatten the train data
3. Find genes that generete the best model
4. Train final model with all the isolates and the selected genes.
   
## Find important genes
After reviewing the literature suggested for the SubChallenge1, we manually listed the biological processes that have been found to be resistant to artemisinin. Later, we match these processes to the columns **Go Process ID** and **Go Process Desc** of the **PFalciparum_Gene_GOterms** (from now on PF_genes data frame) and created the **GO** data frame to which we manually added a column specifying with 1 that the process is "important" for artemisinin resistance and with 0 that it is "not important".

From the **PF_genes** dataframe, we extracted the genes that are related to each GO Process that appear in the **GO** dataframe and added them as a new column Gene_ID. We also added the column Gene_Count, which is the Gene_ID count.

After that, we filtered the GO processes **GO** dataframe to remove all rows that are not important for artemisinin resistance (those labeled with 0), as well as the rows that have a Gene_count <=25 (i.e. the GO with more than 25 genes were considered ubiquituos and not specific to drug resistance).

Finally, the all the genes that appeared in all the filtered GO processes were considered as the "important" genes and used for next steps. In total we found 489 important genes.

## Filter and flatten the train data

The **train** dataframe was filtered by leaving only the columns of genes that made it to the important genes.

The next step was to "flatten" the **train** data frame. To do that we first grouped the **train** dataframe by Isolate, Timepoint and Treatment and calculated the mean of the different bioreplicate to form and aggregated dataframe:

```python
aggregated_means_train = train.groupby(["Isolate","Timepoint","Treatment"]).mean()
aggregated_means_train.head(10)
```
Next we concatenated and transposed each group of rows that represent each isolate-timepoint-treatment-gene condition combination, and transformed each combination into a single column feature. Therefore instead of having as the rows: 

isolate1-24Hr-DHA 

isolate1-24Hr-UT 

isolate1-6Hr-DHA 

isolate1-6Hr-UT

...

isolate3-24Hr-DHA 

isolate30-24Hr-UT 

isolate30-6Hr-DHA 

isolate30-6Hr-UT

and as columns: 

Gene1, Gene2, Gene3..., Gene489


we "flattened" the data frame and only had one row per isolate:


isolate1 

...

isolate30

and as columns, we had each timepoint-treatment-gene condition combination:

gene1-24Hr-DHA, gene1-24Hr-UT, gene1-6Hr-DHA, gene1-6Hr-UT, ..., gene489-24Hr-DHA, gene489-24Hr-UT, gene489-6Hr-DHA, gene489-6Hr-UT. 

## Find genes that generate the best model

After reading the paper from [Jinyu Chen and Louxin Zhang](https://www.biorxiv.org/content/biorxiv/early/2019/07/11/697896.full.pdf), which surveyed different models for drug response prediction, we decided to use the Kernel ridge regression approach [KRR](https://scikit-learn.org/stable/modules/kernel_ridge.html), implementing using the Scikit Learn library.

To find the genes that improved the model's score the model was trained multiple times, following these steps:

1. The first time, the model was trained with only one gen at a time and then the gene that yielded the best score was selected.

2. The second time, the model was trained with a pair of genes at a time. The pair of genes consisted on the gene found in the previous round plus another gene from the list that was left.

3. This loop was ran 100 times, every time the model was trained the set of selected genes consisted on the previous genes found to generate the best score plus a new gene.

After training the model multiple times, four genes were found to be necesary to reach a good score of 0.99 r2: PF3D7_1323000, PF3D7_1362200, PF3D7_1329100, PF3D7_0818900. Therefore, these were the selected genes to train the final model. 

Each time the model was trained, it was validated agains a group of 5 isolates manually selected. Given that the range of the IC50 of the 30 isolates was wide and there were too few sample points, we decide to select an isolate that would represente every group of 6 isolates sorted in ascending order; to make sure the validation set was a good representation of the whole isolates set. 

These were the selected isolates:
* isolate_07
* isolate_27
* isolate_18
* isolate_11
* isolate_16

## Train final model with all the isolates and the selected genes

Finally we trained the KRR model with all the isolates and the selected genes. 


## **Subchallenge 2**
For this challenge, we decided to use a tree ensemble approach to predict a parasite by their clearance rate. Namely, we used the implementation of gradient boosted trees from the [XGBoost python library](https://xgboost.readthedocs.io/en/latest/tutorials/model.html).

## Data preprocessing
We removed 9 samples from the training dataset which did not have a clearance rate (rows 74, 235, 497, 501, 629, 836, 849, 909, and 956 in the original table). Also, we decided to remove the "country", "asexual stage", and "kmeans group" columns so that we could train a model solely on gene expression levels. For training and validation of our models, we used random splits of 90% of the data for training and 10% for validation.

## Identification and extraction of relevant features
When we first trained XGBoost models with our preprocessed table, we noticed that the model performance (Area Under the ROC curve) on the validation set varied between random splits from ~0.71 up to ~0.92. We decided to look at the feature importance scores that XGBoost assigned to the genes in the best performing model (AUC = 0.92 in validation set). We picked the top 100 most important genes according to XGBoost and used them for training subsequent models that we used in our final submission. The top 100 genes we selected can be found in the file "XGBoost_top_100_genes.txt" (syn20684655).

## Final model: Voting system with XGBoost trees
For our final model, we decided to train ten XGBoost models using only the top 100 genes we identified above. To select the ten XGBoost models, we ran a loop that would train an XGBoost model from scratch using a random parition of the training data. The loop would check if the model performance on the validation set was acceptable. An acceptable model would be one whose best validation AUC during training was >= 0.92). Below we outline the specific parameters we used for all XGBoost models we trained and which we fine tuned manually.

```python
# XGBoost parameters
num_round = 200
param = {
    'booster': 'gbtree',
    'max_depth': 5,
    'eta': 0.1,
    'lambda': 1,
    'alpha': 1,
    'colsample_bytree': 0.5,
    'objective': 'binary:logistic',
    'eval_metric': 'auc'
    }
```
We were able to find 10 acceptable XGBoost models after 316 iterations of our loop.

Below we outline the procedure we used to make final predictions for each isolate on the test dataset:

1. For a given isolate in the test set, we extract the expression values of the top 100 genes across all of its variants (that is, Biological replicates, treatments and timepoints). 
2. Using one XGBoost model at a time, we compute the inference score (probability of having "SLOW" clearance rate and thus artemisinin resistant) for each of the isolate variants. If the maximum inference score across all variants is > 0.5, then the isolate is labeled as "SLOW", or "FAST" otherwise.
3. Once we have labeled the same isolate by the ten XGBoost models, we assign the most frequent (or voted) label to the isolate. The final probability assigned to the isolate is the maximum of the scores given by the XGBoost models that voted on the selected label.


All the associated code and model binaries can be found in our Project Files section (folders syn20684656,syn20684638, respectively).