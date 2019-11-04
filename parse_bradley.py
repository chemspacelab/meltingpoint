
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

import misc
from chemhelp import cheminfo

@misc.memory.cache
def main():


    xsl = pd.read_excel("data/BradleyMeltingPointDataset.xlsx")
    col_names = [col for col in xsl]

    data = xsl.loc[xsl["donotuse"] != "x"]

    data = data.drop(columns="donotuse")
    data = data.drop(columns="donotusebecause")
    data = data.drop(columns="name")
    data = data.drop(columns="csid")
    data = data.drop(columns="link")
    data = data.drop(columns="source")

    # print(data)

    data.to_csv("data/melting_bradley.csv", sep=",", index=False)


    return data


def clean_data(df):

    smiles = df.iloc[1]

    data = {}

    for index, row in df.iterrows():

        smi = row.smiles
        value = row.mpC + 273.15

        molobj, status = cheminfo.smiles_to_molobj(smi)

        if molobj is None:
            print("error:", smi)
            continue

        smi = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        if smi not in data:
            data[smi] = []

        data[smi].append(value)

    misc.save_obj("data/melting_bradley_clean", data)

    return



if __name__ == "__main__":
    pd = main()
    data = clean_data(pd)

