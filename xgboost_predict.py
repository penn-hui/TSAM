from __future__ import division
import numpy as np
import scipy.io as sio  
from scipy import stats
import xgboost as xgb
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import math
import csv
import pickle
from sklearn.metrics import roc_auc_score
import sys

def countGC(s, start,end):
#'compute the GC counts and GC contents of a sequences from (start,end)'
    GCcounts=len(s[start:end].replace('A', '').replace('T', ''))
    GCcontent=GCcounts/len(s)
    return GCcounts, GCcontent

def ncf(s,order):
#nucleotide composition features and position nucletide binary features
#s--sequence
#order--1,2,3: A,AA,AAA
     t_list = ["A","G","C","T"]
     L=len(s)
     if order==1:
         nc=t_list
     elif order==2:
         nc=[m+n for m in ["A","G","C","T"] for n in ["A","G","C","T"]]
     elif order==3:
         nc=[m+n+k for m in ["A","G","C","T"] for n in ["A","G","C","T"] for k in ["A","G","C","T"]]
     nc_f=np.zeros((1,4**order))
     
     pos_fea=np.zeros((1,1))
     for i in range(0,L-order+1):
         pos=np.zeros((1,4**order))
         for j in range(0,len(nc)):
             if s[i:i+order]==nc[j]:
                 nc_f[0][j]=nc_f[0][j]+1
                 pos[0][j]=1
                 pos_fea=np.hstack((pos_fea,pos))
         
     n=len(pos_fea[0])
     pos_fea=pos_fea[0][1:n]
     return nc_f, pos_fea
 
def evaluate_performance(predict, real, pos_lab, neg_lab, Type):
    if Type=='cls': # classification
       TN=0;
       TP=0;
       FP=0;
       FN=0;
       for i in range(len(predict)):
           #TN p=r=n_l
           if predict[i]==real[i]==neg_lab:
               TN=TN+1
           elif predict[i]==real[i]==pos_lab:
               TP=TP+1
           elif predict[i]>real[i]:
               FP=FP+1
           elif predict[i]<real[i]:
               FN=FN+1
       Specificity=TN/(TN+FP);
       Recall=TP/(TP+FN);
       Precision=TP/(TP+FP);
       Accuracy=(TP+TN)/(TP+FP+TN+FN);
       F1=2*Recall*Precision/(Recall+Precision);
       MCC=(TP*TN-FP*FN)/math.sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP))
         
       return Specificity, Recall, Precision, Accuracy, F1, MCC
    elif Type=='reg':
        spr = stats.spearmanr(real, predict)
        return spr

def TM_cal(s):
    # computing tm features, s should be 30
   TM_region=np.zeros((1,4))
   s_20=s[4:24]
   s_30=s
   tm_sub1=mt.Tm_NN(Seq(s_20[2:7]))
   tm_sub2=mt.Tm_NN(Seq(s_20[7:15]))
   tm_sub3=mt.Tm_NN(Seq(s_20[15:20]))
   tm_sub4=mt.Tm_NN(Seq(s_30[0:30]))
   
   TM_region[0][0]=tm_sub1
   TM_region[0][1]=tm_sub2
   TM_region[0][2]=tm_sub3
   TM_region[0][3]=tm_sub4

   return TM_region

def analysisImprotances(Improtance,featureNum):
    L=len(Improtance)
    keyMatrix=np.zeros((featureNum,L))
    impMatrix=np.zeros((featureNum,L))
    for i in range(0,L):
        impro=Improtance[i]
        key=impro.keys()
        for j in range(0,len(key)):
            fea_ind=int(key[j][1:len(key[j])])
            fea_imp=impro[key[j]]
            keyMatrix[j,i]=fea_ind
            impMatrix[j,i]=fea_imp
            
    return keyMatrix, impMatrix

#def xgboostPredict(seqs,features,featype):
    


if __name__ == '__main__':
    featype = int(sys.argv[2])
    pretype = int(sys.argv[1])
#    featype = int('1')
#    pretype = int('1')
    matfn='test_seqs.mat'
    data=sio.loadmat(matfn)
    seqs=data['test_seqs']
    cut_fea=data['cut_fea']

    ######extract features########
    Seqs=seqs.tolist()
    L=len(Seqs)
    seqs_20=[]
    seqs_30=[]
    for i in range(0,L):
        seq=Seqs[i][1][0]
        seqs_20.append(seq[4:24])
        seqs_30.append(seq)
    
    TM_region=TM_cal(seqs_30[0])
    for i in range(1,L):
        tm_region=TM_cal(seqs_30[i])
        TM_region=np.vstack((TM_region,tm_region))
        
    nc1, pos_fea1=ncf(seqs_30[0],1)
                     
    for i in range(1,L):
        nc_f,pos_fea=ncf(seqs_30[i],1)
        nc1=np.vstack((nc1,nc_f))
        pos_fea1=np.vstack((pos_fea1,pos_fea))    
    
    nc2, pos_fea2=ncf(seqs_30[0],2)                 
    for i in range(1,L):
        nc_f, pos_fea=ncf(seqs_30[i],2)
        nc2=np.vstack((nc2,nc_f))
        pos_fea2=np.vstack((pos_fea2,pos_fea))
        
    nc3, a=ncf(seqs_30[0],3)                 
    for i in range(1,L):
        nc_f, a=ncf(seqs_30[i],3)
        nc3=np.vstack((nc3,nc_f))
        
    GC_fea=np.zeros((L,2))
    for i in range(0,L):
        seq=seqs_30[i]
        GCcounts, GCcontent=countGC(seq, 0,30)
        GC_fea[i][0]=GCcounts
        GC_fea[i][1]=GCcontent
    
    if featype==1:
        Feature_all=np.hstack((TM_region,nc1,nc2,nc3,GC_fea,pos_fea1,pos_fea2,cut_fea))

    elif featype==2:
        Feature_all=np.hstack((TM_region,nc1,nc2,nc3,GC_fea,pos_fea1,pos_fea2))
    elif featype==3:
        feature_all=np.hstack((TM_region,nc1,nc2,nc3,GC_fea,pos_fea1,pos_fea2,cut_fea))
        Feature_all=feature_all[:,0:675]
    fea_apply=Feature_all
    ###### save features #########
    file_path='test_fea'+str(pretype)+'_'+str(featype)+'.csv'
    with open(file_path,'wb') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for line in fea_apply:
            writer.writerow(line)
            
    modelfile='finalized_model_'+str(pretype)+'_'+str(featype)+'.model'
    #loaded_model = pickle.load(open(modelfile, 'rb'))
    bst = xgb.Booster({'nthread':8}) #init model
    bst.load_model(modelfile) # load data
    x_test=fea_apply
    y_test=np.zeros((len(x_test[:,0]),1))
    data_test = xgb.DMatrix(x_test, y_test)
    y_pred = bst.predict(data_test)
    Y_pred=np.zeros((len(y_pred),1))
    
    for i in range(0,len(y_pred)):
        Y_pred[i][0]=y_pred[i]
        
    file_path='xgboost_pre_score'+str(pretype)+'_'+str(featype)+'.csv'    

    with open(file_path,'wb') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for line in Y_pred:
            writer.writerow(line)        
