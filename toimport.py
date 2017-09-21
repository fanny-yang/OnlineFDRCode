import numpy as np
import os


proc_list = ['Bonferroni', 'LORD++', 'LORD', 'mem-LORD++', 'mem-LORD++ (abstinent)']


def saveres(direc, filename, mat, ext = 'dat', verbose = True):
    filename = "%s.%s" % (filename, ext)
    if not os.path.exists(direc):
        os.makedirs(direc)
    savepath = os.path.join(direc, filename)
    np.savetxt(savepath, mat, fmt='%.3e', delimiter ='\t')
    if verbose:
        print("Saving results to %s" % savepath)

def compute_memFDR(totrej_vec, falrej_vec, penw_vec, mempar):
    
    NUMHYP = len(totrej_vec)
    memory_vec = np.zeros(NUMHYP)
    memFDR_vec = np.zeros(NUMHYP)
    for j in range(NUMHYP):
        memory_vec[0:j] = [mempar**(j-i) for i in range(j)]
        memFDR_num = np.sum(falrej_vec*np.array(memory_vec)*penw_vec)
        memFDR_denom = np.sum(totrej_vec*np.array(memory_vec)*penw_vec)
        if memFDR_denom > 0:
            memFDR_vec[j] = np.true_divide(memFDR_num, max(1, memFDR_denom))
        else:
            memFDR_vec[j] = 0
    return memFDR_vec

def compute_memTDR(totrej_vec, Hypo, penw_vec, mempar):
    
    NUMHYP = len(totrej_vec)
    memory_vec = np.zeros(NUMHYP)
    memBDR_vec = np.zeros(NUMHYP)
    correj_vec = np.array(Hypo)*np.array(totrej_vec)
    for j in range(NUMHYP):
        memory_vec[0:j] = [mempar**(j-i) for i in range(j)]
        memBDR_num = np.sum(correj_vec*np.array(memory_vec)*penw_vec)
        memBDR_denom = np.sum(Hypo*np.array(memory_vec)*penw_vec)
        if memBDR_denom > 0:
            memBDR_vec[j] = np.true_divide(memBDR_num, memBDR_denom)
        else:
            memBDR_vec[j] = 0
    return memBDR_vec
    
def str2list(string, type = 'int'):
    str_arr =  string.split(',')
    if type == 'int':
        str_list = [int(char) for char in str_arr]
    elif type == 'float':
        str_list = [float(char) for char in str_arr]
    return str_list
