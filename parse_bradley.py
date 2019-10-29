
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def main():


    xsl = pd.read_excel("data/BradleyMeltingPointDataset.xlsx")
    col_names = [col for col in xsl]
    print(col_names)

    data = xsl.loc[xsl["donotuse"] != "x"]

    data = data.drop(columns="donotuse")
    data = data.drop(columns="donotusebecause")
    data = data.drop(columns="name")
    data = data.drop(columns="csid")
    data = data.drop(columns="link")
    data = data.drop(columns="source")

    print(data)

    data.to_csv("data/melting_bradley.csv", sep=",", index=False)

    return


if __name__ == "__main__":
    main()

