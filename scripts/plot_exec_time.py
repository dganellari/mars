import os, sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter
import numpy as np

def plot_2D():
    fig = plt.figure()
    fig2, ax = plt.subplots()
    # plt.tick_params(axis='bsoth', which='major')


    path_serial_mycode = "/Users/liudmilakaragyaur/Documents/Eurohack2020_local"
    path_serial_newcode = "/Users/liudmilakaragyaur/Documents/Eurohack2020_local"
    path_GPU_newcode = "/Users/liudmilakaragyaur/Documents/Eurohack2020_local"

    
    filename_serial_mycode = os.path.join(path_serial_mycode, 'timing_spmv_CPU.csv')
    data_serial_mycode = np.genfromtxt(filename_serial_mycode, delimiter=",", names=["n_dofs", "seconds"]) 


    filename_serial_newcode = os.path.join(path_serial_newcode, 'timing_mf_OMP.csv')
    data_serial_newcode = np.genfromtxt(filename_serial_newcode, delimiter=",", names=["n_dofs", "seconds"]) 

    filename_GPU_newcode = os.path.join(path_GPU_newcode, 'timing_mf_GPU_release.csv')
    data_GPU_newcode = np.genfromtxt(filename_GPU_newcode, delimiter=",", names=["n_dofs", "seconds"]) 


    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.grid(True)
    plt.locator_params(nbins=6)
    # plt.plot(data_serial_mycode["n_dofs"], data_serial_mycode["seconds"],  color = "green" ) # GR "-*",
    plt.loglog(data_serial_newcode["n_dofs"], data_serial_newcode["seconds"],  color = "red" ) # RB "--",
    plt.loglog(data_GPU_newcode["n_dofs"], data_GPU_newcode["seconds"],color = "blue" ) # unif

    
    fig2.suptitle("Comparison  of execution time of BiCGStab", fontsize=15)


    plt.xlabel(r"Problem size", fontsize=10)
    plt.ylabel(r"Execution time (in seconds)", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylim((0, 350))
    # plt.yticks(np.arange(10e-3, 10e-5, 10e-5))
    
    outname= "execution_time"

    # plt.legend([r'SPMV - MacBook Pro 2.2 GHz Intel Core i7', r'MatrixFree - OpenMP', r'MatrixFree - GPU'], fontsize = 15)
    plt.legend([r'MatrixFree - OpenMP', r'MatrixFree - GPU'], fontsize = 15)

    plt.savefig('%s_3D.pdf'%outname, dpi=1000)
    plt.show()




plot_2D()




