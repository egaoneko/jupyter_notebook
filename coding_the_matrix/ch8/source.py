import numpy as np
import random
from numbers import Number

class GF2(object):
    MODULO = 2

    def __init__(self, val):
        self.val = int(val) % self.MODULO

    def __add__(self, val):
        return self.__class__((self.val + int(val)) % self.MODULO)
    def __sub__(self, val):
        return self.__class__((self.val - int(val)) % self.MODULO)
    def __mul__(self, val):
        return self.__class__((self.val * int(val)) % self.MODULO)
    def __div__(self, other):
        if other == 0: raise ZeroDivisionError
        return self
    __truediv__ = __div__
    def __rdiv__(self,other): return other
    __rtruediv__ = __rdiv__
    __radd__ = __add__
    __rsub__ = __add__
    __rmul__ = __mul__
    def __int__(self):
        return self.val
    def __repr__(self):
        return "%s(%d)" % (self.__class__.__name__, self.val)
    def __float__(self):
        return float(self.val)
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.val == other.val
        return False
    def __ne__(self, other):
        return not self.__eq__(other)
    
class GF2Utils(object):
    
    @staticmethod
    def getVector(vec):
        vector = np.array(vec)
        Gv = np.vectorize(GF2)
        return Gv(vector)
    
    @staticmethod
    def getMatrix(mat):
        matrix = np.matrix(mat)
        Gv = np.vectorize(GF2)
        return Gv(matrix)
    
    @staticmethod
    def randGF2():
        return GF2(random.randint(0,1))
    

def solve (a, b):
    return np.linalg.lstsq(a, b) if isinstance(a, int) else np.linalg.lstsq(a, b)

def is_superfluous(L, i):
    if len(L) == 1:
        return False
    
    is_gf2 = isinstance(L[0][0], GF2)
    
    if isinstance(L, np.ndarray):
        L = L.tolist()
    
    if is_gf2:
        return _is_superfluous_for_gf2(L, i)
    else:
        return _is_superfluous(L, i)

def _is_superfluous(L, i):
    zero_like = 1e-14
    A = np.transpose(np.array(L[:i] + L[i+1:]))
    b = np.array(L[i])
    u = np.transpose(np.array([solve(A,b)[0]]))
    residual = np.squeeze(np.asarray((np.transpose(np.asmatrix(b)) - np.asmatrix(A)*np.asmatrix(u))))
    return np.dot(residual, residual) < zero_like

def _is_superfluous_for_gf2(L, i):
    A = np.transpose(GF2Utils.getVector(L[:i] + L[i+1:]))
    b = GF2Utils.getVector(L[i])
    u = np.transpose(GF2Utils.getVector([solve(A,b)[0]]))
    residual = np.squeeze(np.asarray((np.transpose(np.asmatrix(b)) - np.asmatrix(A)*np.asmatrix(u))))
    au = np.squeeze(np.asarray(np.transpose(np.asmatrix(A)*np.asmatrix(u)))).tolist()
    b = np.squeeze(np.asarray(np.asmatrix(b))).tolist()
    return au == b

def is_independent(L):
    for idx, val in enumerate(L):
        if(is_superfluous(L,idx)):
            return False
    return True

def subset_basis_shrink(T):
    S = T[:]

    for v in T:
        size = len(S)-1
        if(is_superfluous(S,size)):
            S = S[:size]
    return S

def subset_basis_glow(T):
    S = []

    for v in T:
        S.append(v)
        if(not is_independent(S)):
            S.pop()
    return S

def subset_basis(T):
    return subset_basis_glow(T)

def superset_basis(T, L):
    super_basis_list = T[:]

    # step2: add the vec from L by grow style
    for i in range(len(L)):
        temp = super_basis_list+[L[i]]
        if is_independent(temp):
            super_basis_list.append(L[i])
    return super_basis_list

#def exchange(S, A, z):
#    w = []
#
#    S_set = convertListInList2TupleInSet(S)
#    A_set = convertListInList2TupleInSet(A)
#    S_int_A = list(A_set)
#    S_sub_A = list(S_set - A_set)
#    S_uni_Z = S_int_A + [tuple(z)] + S_sub_A
#    prev_idx = len(S_uni_Z) - len(S_sub_A)
#    for i, v in enumerate(S_sub_A):
#        idx = prev_idx + i
#        if sc.is_superfluous(S_uni_Z, idx):
#            w.append(S_uni_Z[idx])
#    return w

def exchange(S, A, z):
    n_R = []
    n_S = S[:]
    n_S.append(z)

    for i in range(len(n_S)):
        if is_superfluous(n_S, i) and n_S[i] not in A and n_S[i] != z:
            n_R.append(n_S[i])
    return n_R

def rank(L):
    return len(subset_basis(L))

#def is_independent(L):
#    return len(L) == rank(L)

def is_invertible(M):
    m = np.matrix(M)
    t_m = np.transpose(m)

    r_rank = rank(m.tolist())
    # print(r_rank)

    c_rank = rank(t_m.tolist())
    # print(c_rank)

    return (r_rank == c_rank) and is_independent(t_m.tolist())

def find_matrix_inverse(A):
    if not is_invertible(A):
        return None
    A_m = np.matrix(A)
    # print(A)
    I = np.identity(rank(A_m.tolist()))
    # print(b)
    vec = solve(A_m, I)[0]
    # print(vec)
    return vec.tolist()

def generate_echelon_form(M):
    m = np.array(M, dtype="float64")
    row, col = m.shape
    
    i = np.identity(row)
    
    new_rowlist = []
    new_m_rowlist = []
    rows_left = set(range(row))
    
    for c in range(col):
        rows_with_nonzero = [r for r in rows_left if m[r][c] != 0]
        
        if rows_with_nonzero != []:
            pivot = rows_with_nonzero[0]
            
            new_rowlist.append(m[pivot])
            rows_left.remove(pivot)
            new_m_rowlist.append(i[pivot])
            
            for r in rows_with_nonzero[1:]:
                multiplier = m[r][c] / m[pivot][c]
                m[r] -= multiplier * m[pivot]
                i[r] -= multiplier * i[pivot]
                
    for r in rows_left:
        new_rowlist.append(m[r])
        
    for r in rows_left:
        new_m_rowlist.append(i[r])
    
    return np.matrix(new_rowlist), np.matrix(new_m_rowlist)
