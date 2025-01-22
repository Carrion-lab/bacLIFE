# evaluation of a model fit using mutual information input features
from pandas import read_csv
import argparse
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import chi2
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
from sklearn.neighbors import KNeighborsClassifier

# load the dataset
def load_dataset(filename):
	# load the dataset as a pandas DataFrame
	data = read_csv(filename, sep = '\t', header = 0,index_col=0)
	col = 'Sample'
	data = data.loc[:, data.columns != col]
	# retrieve numpy array
	dataset = data.values
	# split into input (X) and output (y) variables
	X = dataset[:, :-1]
	y = dataset[:,-1]
	# format all fields as string
	#X = X.astype(str)
	return data, X, y

 
# prepare target
def prepare_targets(y_train, y_test):
	le = LabelEncoder()
	le.fit(y_train)
	y_train_enc = le.transform(y_train)
	y_test_enc = le.transform(y_test)
	return y_train_enc, y_test_enc



# feature selection
def select_features(X_train, y_train, X_test):
	fs = SelectKBest(score_func=chi2, k=5000) #score_func=mutual_info_classif
	fs.fit(X_train, y_train)
	X_train_fs = fs.transform(X_train)
	X_test_fs = fs.transform(X_test)
	
	return X_train_fs, X_test_fs, fs


# Get columns to keep and create new dataframe with those only
def write_feature_names(fs, data, filename):
  cols = fs.get_support(indices=True)
  features_df_new = data.iloc[:,cols]
  features_df_new = features_df_new.columns.values.tolist()
  f = open(filename, 'w')
  for i in range(len(features_df_new)):
    f.write(features_df_new[i])
    f.write('\n')
  f.close()
 
#Function extract first element of every sublist 
def Extract(lst,position):
    return list(list(zip(*lst))[position])



#Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--lifestyle', dest = 'lifestyle', help = 'lifestyle', required = True)
parser.add_argument('-u', '--unknown', dest = 'input_unknown', help = 'input_unknown', required = True)
parser.add_argument('-i', '--inputtable', dest = 'input_table', help = 'input_table', required = True)
args = parser.parse_args()

print('\n')
print('\n')
print('-----Running random forest for class: %s-----' % args.lifestyle )

# load the dataset
data, X, y = load_dataset(args.input_table)
# split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1)


# prepare output data
y_train_enc, y_test_enc = y_train, y_test
# feature selection
X_train_fs, X_test_fs, fs = select_features(X_train, y_train_enc, X_test)
feature_file = str('classifier/features_%s' % (args.lifestyle))
write_feature_names(fs, data,feature_file)
print('Features selected for modelling saved as: %s' % feature_file )


##Join train and test
X_cv = np.vstack([X_train_fs, X_test_fs])
y_cv = np.concatenate([y_train_enc, y_test_enc])



##Model
clf = RandomForestClassifier(
    n_estimators=50,
    class_weight='balanced'
)


##DummyModel
#clf = KNeighborsClassifier(n_neighbors = 1, metric = 'jaccard')


###cross validation2
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

cv = StratifiedKFold(n_splits=5)
fig, ax = plt.subplots()
for i, (train, test) in enumerate(cv.split(X_cv, y_cv)):
    clf.fit(X_cv[train], y_cv[train])
    viz = RocCurveDisplay.from_estimator(
        clf,
        X_cv[test],
        y_cv[test],
        name="ROC fold {}".format(i),
        alpha=0.3,
        lw=1,
        ax=ax,
    )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)

ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", alpha=0.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(
    mean_fpr,
    tprs_lower,
    tprs_upper,
    color="grey",
    alpha=0.2,
    #label=r"$\pm$ 1 std. dev.",
)

ax.set(
    xlim=[-0.05, 1.05],
    ylim=[-0.05, 1.05],
    title="Receiver operating characteristic example",
)
ax.legend(loc="lower right")
plot_file = str('classifier/ROC_%s') % (args.lifestyle)
plt.savefig(plot_file)
print('ROC plot saved as: %s' % plot_file )





###Classify Unknow lifestyle bacteria
model_all_data = clf.fit(X_cv, y_cv)
data = read_csv(args.input_unknown, sep = '\t', header = 0,index_col=0)
# retrieve numpy array
col = "Sample"
X_new = data.loc[:, data.columns != col]
X_new = X_new.values
X_new_fs = fs.transform(X_new)
ynew = model_all_data.predict_proba(X_new_fs)

data['0'] = Extract(ynew, 0)
data['1'] = Extract(ynew, 1)
data = data[['Sample', '0', '1']]
file = str('classifier/predicted_%s.csv' % (args.lifestyle))
data.to_csv(file)
print('Predictions on unknown genomes saved as: %s' % file )




