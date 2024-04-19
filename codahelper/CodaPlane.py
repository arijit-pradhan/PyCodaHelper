import numpy as np

class CodaPlane:
    def __init__(self, origin, normal):
        self.__origin = origin
        self.__normal = normal
        
        # check statements
        if (not (isinstance(self.__origin, np.ndarray))) or (not(isinstance(self.__normal, np.ndarray))):
            raise Exception("Both normal and origin have to be numpy arrays ")
        self.__dim = origin.shape[0]
        if(self.__origin.ndim != self.__origin.ndim):
            raise Exception("origin and normal different dimensional nd-array")
        if(self.__origin.ndim > 1) or (self.__normal.ndim > 1):
            raise Exception("Both origin and normal should be 1 or 0 dimensional nd-array")
        if(self.__normal.shape[0] != self.__dim):
            raise Exception("Mismatch in dimension of normal and origin")
        
        self.__origin.shape = (self.__dim, 1)
        self.__normal.shape = (self.__dim, 1)
        
    def EvalPlaneEquation(self, pt):
        # checks
        if(not(isinstance(pt, np.ndarray))):
            raise Exception("pt has to be a numpy arrays")
        if(pt.shape[0] != self.__dim):
            raise Exception("Pt has to be {}-dimensional".format(self.__dim))
        pt.shape = (self.__dim, 1)
        
        d = self.__normal.T @ self.__origin
        return (self.__normal.T @ pt) - d
    
    def ReflectPoint(self, pt):
        # checks
        if(not(isinstance(pt, np.ndarray))):
            raise Exception("pt has to be a numpy arrays")
        if(pt.shape[0] != self.__dim):
            raise Exception("Pt has to be {}-dimensional".format(self.__dim))
        pt.shape = (self.__dim, 1)
        
        xr = np.zeros((self.__dim,1))
        lam = ( (self.__normal.T @ self.__origin) - (self.__normal.T @ pt) ) / (self.__normal.T @ self.__normal);
        xr = pt + 2*lam*self.__normal
        return xr