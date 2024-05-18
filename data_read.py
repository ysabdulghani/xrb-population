import os
import sys
import glob


def readDatesFile(fileName):
    with open(fileName) as date_f:
        dates = date_f.read().splitlines()
    return dates

def findProductsFolders(instrument,dates):
    if instrument == 'maxi':
        if len(dates) == 4:
            if dates[0] != dates[1]:
                softProductsFolder = glob.glob('products_'+dates[0]+'_to_'+dates[1])[0]
            else:
                softProductsFolder = glob.glob('products_'+dates[0])[0]
            if dates[2] != dates[3]:
                transProductsFolder = glob.glob('products_'+dates[2]+'_to_'+dates[3])[0]
            else:
                transProductsFolder = glob.glob('products_'+dates[2])[0]
        else:
            if dates[0] != dates[1]:
                softProductsFolder = glob.glob('products_'+dates[0]+'_to_'+dates[1])[0]
            else:
                softProductsFolder = glob.glob('products_'+dates[0])[0]
            transProductsFolder = None
    elif instrument == 'xrt' or instrument == 'rxte_pca':
        if len(dates) == 2:
            softProductsFolder = glob.glob('products_'+dates[0])[0]
            transProductsFolder = glob.glob('products_'+dates[1])[0]
        elif len(dates) == 1:
            softProductsFolder = glob.glob('products_'+dates[0])[0]
            transProductsFolder = None
    else:
        raise ValueError('Cannot find products folder!')

    return (softProductsFolder,transProductsFolder)

def findSavedModel():
    rootDir = os.getcwd()
    for dirName, subdirList, fileList in os.walk(rootDir):
        for fname in fileList:
            if 'PLonlyFlux' in fname and '.xcm' in fname and '2023' not in fname and 'fake' in fname:
                return fname
    return None





