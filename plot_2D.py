import os, sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter
import numpy as np

def plot_2D():
    fig = plt.figure()
    fig2, ax = plt.subplots()
    # plt.tick_params(axis='bsoth', which='major')


    path_serial_output = "/Users/liudmilakaragyaur/Documents/Eurohack2020_local"

    filename_serial_output = os.path.join(path_serial_output, 'output_2D_new.csv')
    # data_serial_output = np.genfromtxt(filename_serial_output, delimiter=",", names=["n_dofs", "L2_err", "n_iter"]) 

    data_serial_output = np.genfromtxt(filename_serial_output, delimiter=",", names=["n_dofs", "L2_err"]) 

    filename_serial_output_old = os.path.join(path_serial_output, 'output_old_2D.csv')
    data_serial_output_old = np.genfromtxt(filename_serial_output_old, delimiter=",", names=["n_dofs", "L2_err", "n_iter"]) 


    y = quadr(data_serial_output_old["n_dofs"])

    print
    print(y)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.grid(True)
    plt.locator_params(nbins=6)
    plt.loglog(data_serial_output["n_dofs"], data_serial_output["L2_err"],  color = "green" ) # GR "-*",
    plt.loglog(data_serial_output_old["n_dofs"], data_serial_output_old["L2_err"],  color = "blue" ) # GR "-*",
    plt.loglog(data_serial_output_old["n_dofs"], y, color = "grey" ) # GR "-*",

    
    # fig2.suptitle("Error in L2-norm", fontsize=15)


    plt.xlabel(r"Problem size", fontsize=10)
    plt.ylabel(r"BiCGStab iter", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    # plt.ylim((0, 350))
    # plt.yticks(np.arange(10e-3, 10e-5, 10e-5))
    
    outname= "execution_time"

    # plt.legend([r'SPMV - MacBook Pro 2.2 GHz Intel Core i7', r'MatrixFree - OpenMP', r'MatrixFree - GPU'], fontsize = 15)
    plt.legend([r'GPU code', r'CPU code', r'Linear'], fontsize = 15)

    plt.savefig('%s_2D.pdf'%outname, dpi=1000)
    plt.show()



def quadr(x):
    res = []
    for n in range(np.shape(x)[0]):
        res.append(1/(x[n]))
        # res.append(np.sqrt(x[n]))
    return np.array(res)
    

plot_2D()

