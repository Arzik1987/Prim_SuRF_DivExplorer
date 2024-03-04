import numpy as np
import pandas as pd


class PRIM_nominal_pat:
    def __init__(self, mass_min = 0.1, threshold = 1):
        self.mass_min = mass_min
        self.threshold = threshold

    def fit(self, X, y):
        self.X_ = X.copy()
        self.y_ = y.copy()
        self._get_initial_restrictions(X)
        self.N_ = len(y)
        
        hgh = -1   
        self.boxes_ = []
        self.quals_ = []
        self.mass_ = []
        i = 1
        
        while hgh != -np.inf and i < 100 and np.sum(self.y_)/len(self.y_) < self.threshold:
            self.boxes_.append(self.box_.copy())
            self.quals_.append(np.sum(self.y_)/len(self.y_))
            self.mass_.append(len(self.y_)/self.N_ )
            i = i + 1
            hgh = self._peel_one()
            
        self.X_ = None
        self.y_ = None
        self.box_ = None
        self.N_ = None
        
        return self
    
    def _peel_one(self):
        hgh = -np.inf
        cn, valind = -1, -1
        
        for i in range(0, self.X_.shape[1]):
            for j in range(0, len(self.box_[i])):
                retain = self.X_.iloc[:,i] != self.box_[i][j]
                if np.count_nonzero(retain)/self.N_ > self.mass_min:
                    tar = self._target_fun(retain)
                    if tar > hgh:
                        hgh = tar
                        inds = retain
                        cn = i
                        valind = j
        
        if hgh != -np.inf:
            self.X_ = self.X_[inds]
            self.y_ = self.y_[inds]
            self.box_[cn] = np.delete(self.box_[cn], valind)
        
        return hgh
    
    def _get_initial_restrictions(self, X):
        self.box_ = []
        for i in X.columns:
            self.box_ += [pd.unique(X[i])]
    
    def _target_fun(self, retain):
        # return np.sum(self.y_[retain])/np.count_nonzero(retain) - np.sum(self.y_[~retain])/np.count_nonzero(~retain)
        bsize = np.count_nonzero(~retain)
        Bsize = len(retain)
        Bby = np.sum(self.y_[retain])/np.count_nonzero(retain)
        By = np.sum(self.y_)/len(retain)
        return (Bsize - bsize)*(Bby - By)/bsize
    
    def get_res(self):
        return self.boxes_, self.quals_, self.mass_


class PRIM_nominal_fpr:
    def __init__(self, mass_min = 0.1, threshold = 1):
        self.mass_min = mass_min
        self.threshold = threshold

    def fit(self, X, y, ypred):
        self.X_ = X.copy()
        self.y_ = y.copy()
        self.ypred_ = ypred.copy()
        self._get_initial_restrictions(X)
        self.N_ = len(y)
        
        hgh = -1   
        self.boxes_ = []
        self.quals_ = []
        self.mass_ = []
        i = 1
        
        while hgh != -np.inf and i < 100:
            self.boxes_.append(self.box_.copy())
            tn = self.y_ == 0
            self.quals_.append(np.sum(self.ypred_[tn])/len(self.ypred_[tn]))
            self.mass_.append(len(self.y_)/self.N_ )
            i = i + 1
            hgh = self._peel_one()
            
        self.X_ = None
        self.y_ = None
        self.box_ = None
        self.N_ = None
        
        return self
    
    def _peel_one(self):
        hgh = -np.inf
        cn, valind = -1, -1
        
        for i in range(0, self.X_.shape[1]):
            for j in range(0, len(self.box_[i])):
                retain = self.X_.iloc[:,i] != self.box_[i][j]
                if np.count_nonzero(retain)/self.N_ > self.mass_min:
                    tar = self._target_fun(retain)
                    if tar > hgh:
                        hgh = tar
                        inds = retain
                        cn = i
                        valind = j
        
        if hgh != -np.inf:
            self.X_ = self.X_[inds]
            self.y_ = self.y_[inds]
            self.ypred_ = self.ypred_[inds]
            self.box_[cn] = np.delete(self.box_[cn], valind)
        
        return hgh
    
    def _get_initial_restrictions(self, X):
        self.box_ = []
        for i in X.columns:
            self.box_ += [pd.unique(X[i])]
    
    def _target_fun(self, retain):
        bsize = np.count_nonzero(~retain)
        Bsize = len(retain)
        Bbtn = self.y_[retain] == 0
        Btn = self.y_ == 0
        Bby = sum(self.ypred_[retain][Bbtn])/np.count_nonzero(Bbtn)
        By = np.sum(self.ypred_[Btn])/np.count_nonzero(Btn)
        return (Bsize - bsize)*(Bby - By)/bsize
    
    def get_res(self):
        return self.boxes_, self.quals_, self.mass_


class PRIM_nominal_fnr:
    def __init__(self, mass_min = 0.1, threshold = 1):
        self.mass_min = mass_min
        self.threshold = threshold

    def fit(self, X, y, ypred):
        self.X_ = X.copy()
        self.y_ = y.copy()
        self.ypred_ = ypred.copy()
        self._get_initial_restrictions(X)
        self.N_ = len(y)
        
        hgh = -1   
        self.boxes_ = []
        self.quals_ = []
        self.mass_ = []
        i = 1
        
        while hgh != -np.inf and i < 100:
            self.boxes_.append(self.box_.copy())
            tp = self.y_ == 1
            self.quals_.append(np.sum(self.ypred_[tp] == 0)/len(self.ypred_[tp]))
            self.mass_.append(len(self.y_)/self.N_ )
            i = i + 1
            hgh = self._peel_one()
            
        self.X_ = None
        self.y_ = None
        self.box_ = None
        self.N_ = None
        
        return self
    
    def _peel_one(self):
        hgh = -np.inf
        cn, valind = -1, -1
        
        for i in range(0, self.X_.shape[1]):
            for j in range(0, len(self.box_[i])):
                retain = self.X_.iloc[:,i] != self.box_[i][j]
                if np.count_nonzero(retain)/self.N_ > self.mass_min:
                    tar = self._target_fun(retain)
                    if tar > hgh:
                        hgh = tar
                        inds = retain
                        cn = i
                        valind = j
        
        if hgh != -np.inf:
            self.X_ = self.X_[inds]
            self.y_ = self.y_[inds]
            self.ypred_ = self.ypred_[inds]
            self.box_[cn] = np.delete(self.box_[cn], valind)
        
        return hgh
    
    def _get_initial_restrictions(self, X):
        self.box_ = []
        for i in X.columns:
            self.box_ += [pd.unique(X[i])]
    
    def _target_fun(self, retain):
        bsize = np.count_nonzero(~retain)
        Bsize = len(retain)
        Bbtp = self.y_[retain] == 1
        Btp = self.y_ == 1
        Bby = sum(self.ypred_[retain][Bbtp] == 0)/np.count_nonzero(Bbtp)
        By = np.sum(self.ypred_[Btp] == 0)/np.count_nonzero(Btp)
        return (Bsize - bsize)*(Bby - By)/bsize
    
    def get_res(self):
        return self.boxes_, self.quals_, self.mass_
