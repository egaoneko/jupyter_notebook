import numpy as np

def solve (a, b):
    return np.linalg.lstsq(a, b) if isinstance(a, int) else np.linalg.lstsq(a, b)

def is_superfluous(L, i):
    zero_like = 1e-14
    A = np.transpose(np.array(L[:i] + L[i+1:]))
    b = np.array(L[i])
    u = np.transpose(np.array([solve(A,b)[0]]))

    residual = np.squeeze(np.asarray((np.transpose(np.asmatrix(b)) - np.asmatrix(A)*np.asmatrix(u))))
    return np.dot(residual, residual) < zero_like

def is_independent(L):
    for idx, val in enumerate(L):
        if(sc.is_superfluous(L,idx)):
            return False
    return True
