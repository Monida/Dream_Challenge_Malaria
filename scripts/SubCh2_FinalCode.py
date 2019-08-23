#######################
# IMPORT LIBRARIES    #
#######################
import xgboost as xgb
import pandas as pd
import numpy as np
import pickle
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

#######################
# PREPROCESS DATA     #
#######################
train = pd.read_csv("../data/SubCh2_TrainingData.csv")

# Filter out rows with clearance rate = 'nan'
rows_to_keep = [row for row in range(len(train)) if str(train['ClearanceRate'][row]) != 'nan']
train = train.iloc[rows_to_keep,:].reset_index(drop=True)
gene_names = list(train.columns)[4:4956]

# Encode Clearance Rates as labels
y_raw = list(train['ClearanceRate'])
labelE = LabelEncoder().fit(y_raw)
y = labelE.transform(y_raw)

# Build feature matrix X
X = train.iloc[:,4:4956].values

#######################
# FEATURE EXTRACTION  #
#######################

# Split data into training and validation
X_train, X_val, y_train, y_val = train_test_split(X,y,test_size=0.1)
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

feature_gains = bst.get_score(importance_type='gain')
feature_gains_df = pd.DataFrame()
feature_gains_df['Gene'] = list(feature_gains.keys())
feature_gains_df['XGBoostGain'] = list(feature_gains.values())

important_genes = list(feature_gains_df.sort_values(by="XGBoostGain", ascending=False).head(100)['Gene'])

# Save list of important genes
f = open('xgboost_SubCh2_July21st_Top100_genes.txt', 'w')
f.write("\n".join(important_genes)+"\n")
f.close()

#######################
# TRAIN XGBOOST MODELS#
#######################

num_models = 10
success = 0
thr = 0.92
new_X = train[important_genes].values
iterations = 0
while success < num_models:
    X_train, X_val, y_train, y_val = train_test_split(new_X,y,test_size=0.1)

    # Build DMatrixes for XGBoost
    dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=important_genes)
    dval = xgb.DMatrix(X_val, label=y_val, feature_names=important_genes)

    # Declare XGBoost parameters
    num_round = 50
    param = {
        'booster': 'gbtree',
        'max_depth': 5,
        'eta': 0.1,
        'lambda': 1,
        'alpha': 1,
        'colsample_bytree': 0.5,
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'verbose': 0
        }

    # Prepare XGBoost run
    evallist = [(dtrain, 'training'), (dval, 'validation')]

    # Train model with early stopping (stop if validation score doesn't improve after 10 rounds)
    new_bst = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=30, verbose_eval=False)
    
    iterations += 1
    
    if new_bst.best_score >= thr:
        print("Found a good model after {} iterations!".format(iterations))
        # Save model
        with open('xgboost_SubCh2_FINAL_MODEL_August15th_{}.pkl'.format(success + 1), 'wb') as xgboost_model_file:
            pickle.dump(new_bst, xgboost_model_file)
        success += 1

#######################
# MAKE PREDICTIONS    #
#######################

test = pd.read_csv('../data/SubCh2_TestData.csv')
test_isolates = list(set(test['Isolate']))
test_isolates.sort()

# Make list of 10 XGBoost models
xgb_models = []
    
for i in range(10):
    f = open('xgboost_SubCh2_FINAL_MODEL_August15th_{}.pkl'.format(i+1), 'rb')
    new_bst = pickle.load(f)
    f.close()
    xgb_models.append(new_bst)

# results = {'model-#': {'isolate-#': {'Clearance': label, 'Probability': p}}}
results = {}
for i in range(len(xgb_models)):
    new_bst = xgb_models[i]
    model = 'model_{}'.format(i+1)
    results[model] = {}
    for isolate in test_isolates:
        tmp = test[test['Isolate']==isolate].copy().reset_index(drop=True)
        tmp_X = tmp[important_genes].values
        dtmp = xgb.DMatrix(tmp_X, feature_names=important_genes)
        tmp_y = new_bst.predict(dtmp)
        predicted = max(tmp_y)
        results[model][isolate] = {}
        if predicted < 0.5:
            results[model][isolate]['Clearance'] = 'FAST'
            results[model][isolate]['Probability'] = 1 - predicted
        else:
            results[model][isolate]['Clearance'] = 'SLOW'
            results[model][isolate]['Probability'] = predicted

results_df = pd.DataFrame.from_dict(results)
isolates = list(results_df.index)
models = list(results_df.columns)

def getAnswer(isolate):
    slow_count = 0
    slow_probs = []
    fast_count = 0
    fast_probs = []
    for model in models:
        clearance = results_df[model][isolate]['Clearance']
        if clearance == 'SLOW':
            slow_count += 1
            slow_probs.append(results_df[model][isolate]['Probability'])
        else:
            fast_count += 1
            fast_probs.append(results_df[model][isolate]['Probability'])
    # Vote and decide
    if slow_count > 5:
        return 'SLOW',np.max(slow_probs), slow_count
    elif slow_count == 5:
        if np.mean(slow_probs) > np.mean(fast_probs):
            return 'SLOW', np.max(slow_probs), slow_count
        else:
            return 'FAST', np.max(fast_probs), fast_count
    else:
        return 'FAST', np.max(fast_probs), fast_count

#######################
# WRITE RESULTS       #
#######################

f = open('../submissions/lylat_SubCh2_model2_150819.txt','w')
f.write("Isolate\tPredicted_Categorical_Clearance\tProbability\n")
for isolate in isolates:
    label, prob, count = getAnswer(isolate)
    f.write("\t".join([isolate, label, str(prob)]))
    f.write("\n")
f.close()
