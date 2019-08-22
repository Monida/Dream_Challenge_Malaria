# **Subchallenge 2**
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
2. Using one XGBoost model at a time, we compute the inference score (probability of having "SLOW" clearance rate and thus artemisin resistant) for each of the isolate variants. If the maximum inference score across all variants is > 0.5, then the isolate is labeled as "SLOW", or "FAST" otherwise.
3. Once we have labeled the same isolate by the ten XGBoost models, we assign the most frequent (or voted) label to the isolate. The final probability assigned to the isolate is the maximum of the scores given by the XGBoost models that voted on the selected label.

All the associated code and model binaries can be found in our Project Files section (folders syn20684656,syn20684638, respectively).