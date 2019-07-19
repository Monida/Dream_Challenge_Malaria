import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

# Load data
train = pd.read_csv("../data/SubCh2_TrainingData.csv")

# Filter out rows with Clearance rate = 'nan'
rows_to_keep = [row for row in range(len(train)) if str(train['ClearanceRate'][row]) != 'nan']
train = train.iloc[rows_to_keep,:].reset_index(drop=True)
gene_names = list(train.columns)[4:4956]

# Encode Clearance Rates as labels
y_raw = list(train['ClearanceRate'])
labelE = LabelEncoder().fit(y_raw)
y = labelE.transform(y_raw)

# Build feature matrix X
X = train.iloc[:,4:4956].values

# Split data into training and validation
X_train, X_val, y_train, y_val = train_test_split(X,y,test_size=0.1,random_state=8)
print("Training size: {}, Validation size:{}".format(len(y_train),len(y_val)))

# Build DMatrixes for XGBoost
dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=gene_names)
dval = xgb.DMatrix(X_val, label=y_val, feature_names=gene_names)

# Declare XGBoost parameters
num_round = 200
param = {
    'booster': 'gbtree',
    'max_depth': 5,
    'eta': 0.1,
    'lambda': 1,
    'alpha': 1,
    'colsample_bytree': 0.5,
    'objective': 'binary:logistic',
    'eval_metric': 'auc',
    'verbose': 2
    }

# Prepare XGBoost run
evallist = [(dtrain, 'training'), (dval, 'validation')]

# Train model with early stopping (stop if validation score doesn't improve after 10 rounds)
bst = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=30)

print("Finished")