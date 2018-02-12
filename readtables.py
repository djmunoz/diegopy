import numpy as np
import pandas as pd
from string import split


def readapj(filename):
    """
    Read an ApJ-format table

    """

    f = open(filename)
    header_cols=[]
    col_locs=[]
    data=[]
    for kk,line in enumerate(f.readlines()):
        if (line[0] == '#'):
            if not(' ' in line):
                continue
            if (kk <= 1):
                #extent = len(split(split(line,'- ')[1],' '))
                extent = line.index(split(split(line,'-')[1])[0])+1
                while (line[extent] != ' '):
                    extent+=1
            print extent
            header_cols.append(split(line[extent:])[2])
            #col_locs.append("".join(split(line[1:])[:2]))
            col_locs.append("".join(line[1:extent]))
        else:
            ncols= len(header_cols)
            data_line=[]
            byte_list=[]
            for col in range(ncols):
                byte_list.append([int(split(col_locs[col],'-')[0]),int(split(col_locs[col],'-')[1])+1])
                data_line.append(line[int(split(col_locs[col],'-')[0])-1:int(split(col_locs[col],'-')[1])+1])   
            data.append(np.asarray(data_line))
    #print col_locs
    #print byte_list
    data = np.asarray(data)
    df = pd.DataFrame(data[:,:ncols],    # values
                      columns=header_cols[:ncols])

    return df


def readtex(filename):
    """
    Read a table in tex format
    
    """
    f = open(filename)
    header_cols=[]
    col_locs=[]
    data=[]
    
    for line in f.readlines():
        return None
