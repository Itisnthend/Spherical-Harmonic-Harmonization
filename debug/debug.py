# -*- coding: utf-8 -*-
"""
8/8/2020, 11:34:27 AM
@author: qianwcontact@gmail.com
"""

import os
import shutil
import sys

import numpy as np
import pandas as pd
import scipy.io as scio
import nibabel as nib


def save_as_nii(data, objDir):

    #print(temp.shape)
    originalDir = '/ifs/loni/faculty/shi/spectrum/qwang/HarmExperiment/HCPDevelopment/HCD0259243_V1_MR/FA.nii.gz'
    original = nib.load(originalDir)
    affine = original.affine
    img = nib.Nifti1Image(data, affine)
    nib.save(img, objDir)




if __name__ == "__main__":
    
    rootDir = "/ifs/loni/faculty/shi/spectrum/qwang"
    #rootDir = '/Users/itisnthend/INI/qwang'
    objDir = os.path.join(rootDir, 'HarmExperiment/HCPDevelopment/mapFunction.nii.gz')
    #subjects_list = ["002_S_0413"]

    alphaDir = os.path.join(rootDir, 'HarmExperiment/HCPDevelopment/mapFunction.mat')
    mapFunction = scio.loadmat(alphaDir)['mapFunction']

    
    #save_as_nii(mapFunction, objDir)
