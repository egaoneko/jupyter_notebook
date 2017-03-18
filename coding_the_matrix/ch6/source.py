import numpy as np

def solve (a, b):
    return np.linalg.lstsq(a, b) if isinstance(a, int) else np.linalg.lstsq(a, b)

def is_superfluous(L, i):
    zero_like = 1e-14

    if len(L) == 1:
        return False
    A = np.transpose(np.array(L[:i] + L[i+1:]))
    b = np.array(L[i])
    u = np.transpose(np.array([solve(A,b)[0]]))

    residual = np.squeeze(np.asarray((np.transpose(np.asmatrix(b)) - np.asmatrix(A)*np.asmatrix(u))))
    return np.dot(residual, residual) < zero_like

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

def exchange(S, A, z):
    w = []

    S_set = convertListInList2TupleInSet(S)
    A_set = convertListInList2TupleInSet(A)
    S_int_A = list(A_set)
    S_sub_A = list(S_set - A_set)
    S_uni_Z = S_int_A + [tuple(z)] + S_sub_A
    prev_idx = len(S_uni_Z) - len(S_sub_A)
    for i, v in enumerate(S_sub_A):
        idx = prev_idx + i
        if sc.is_superfluous(S_uni_Z, idx):
            w.append(S_uni_Z[idx])
    return w

def exchange_another(S, A, z):
    n_R = []
    n_S = S[:]
    n_S.append(z)

    for i in range(len(n_S)):
        if sc.is_superfluous(n_S, i) and n_S[i] not in A and n_S[i] != z:
            n_R.append(n_S[i])
    return n_R
