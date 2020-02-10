import matplotlib.pyplot as plt
from numbers import Number

plt.style.use('seaborn-whitegrid')

def plot(L, scale=4, dot_size = 5):
    xAxis = []
    yAxis = []
    
    for pt in L:
        if isinstance(pt, Number):
            x,y = pt.real, pt.imag
        else:
            if isinstance(pt, tuple) or isinstance(pt, list):
                x,y = pt
            else:
                raise ValueError
        xAxis.append(x)
        yAxis.append(y)
    plt.xlim(-scale, scale)
    plt.ylim(-scale, scale)
    plt.scatter(xAxis, yAxis, s = dot_size)
    plt.show()