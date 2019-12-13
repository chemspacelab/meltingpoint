
import tqdm

import numpy as np
import os
import re
from pprint import pprint
import json
import urllib.request


def get_json(url):
    """
    download json from url
    ignore byte that cannot be converted to utf-8
    """

    # print("wget", url)

    response = urllib.request.urlopen(url)
    data = response.read()
    values = json.loads(data.decode("utf-8", "ignore"))

    return values


def save_json(name, obj, indent=4):
    with open(name + ".json", 'w') as f:
        json.dump(obj, f, indent=indent)


def load_json(name):

    if not name.split(".")[-1] == "json":
        name = name + ".json"

    with open(name, 'r') as f:
        content = f.read()
        content = json.loads(content)
        return content

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def clean_unit(unit):

    unit = unit.replace('deg', '').replace('°','').replace('º','')
    unit = unit.strip()

    return unit


def parse_range(datain):
    """
    if range is provided, take the mean
    """

    n_minus = datain.count("-")
    data = datain.split("-")
    data = [x.strip() for x in data]
    n_split = len(data)
    data = list(filter(None, data))
    data = [float(x) for x in data]

    if len(data) == 1:
        return None

    if len(data) != 2:
        print("error in parse_range", datain)
        quit()

    mean = np.mean(data)

    if n_split == 4:
        mean *= -1

    return mean

def parse_value(data):
    """

    dict?

    """

    data = data.lower()

    if "less than" in data:
        return None

    if "greater than" in data:
        return None

    if "<" in data:
        return None

    if ">" in data:
        return None

    data = data.replace("degress Celcius", "deg c")
    data = data.replace(" to ", " - ")

    data = re.sub(r'\([^)]*\)', '', data)

    pattern_value = re.compile('(?:\-?\d+\.?\d*)(?:\s?\-\s?\-?\d+\.?\d*)?')
    values = pattern_value.findall(data)

    # Find unit
    unit_patterns = [
        'deg [ckf]',
        '°\s?[ckf]',
    ]
    units = []

    for unit in unit_patterns:
        pattern_unit = re.compile(unit)
        units += pattern_unit.findall(data)

    units = [clean_unit(unit) for unit in units]

    n_units = len(units)
    n_values = len(values)

    if n_units != n_values:
        return None
        print(data, values, units)

    return values, units


def parse_line(data):

    name = data["Name"]

    if "RADIOACTIVE" in name:
        return None

    try:
        names = data["Link"]["ToCompound"]["ByCID"]
    except KeyError:
        names = []

    values = data["Data"]

    collect_values = []

    for value in values:

        value = value["Value"]

        if "Number" in value and "Unit" in value:
            hm = str(value["Number"]) + " " + value["Unit"]
            value["StringWithMarkup"] = [{"String": hm}]

        if "StringWithMarkup" in value:
            value = value["StringWithMarkup"][0]["String"]
            value = value.lower()

            ignore_keywords = ["decomp", "anhydrous form"]

            if any(x in value for x in ignore_keywords):
                continue

            value = parse_value(value)
            if value is None:
                return None

            for val, unt in zip(*value):
                collect_values.append((val, unt))

        elif "Number" in value and "Unit" in value:

            number = value["Number"]
            unit = value["Unit"].lower().replace("dec", "deg")
            unit = clean_unit(unit)

            for val in number:
                collect_values.append((val, unit))

        else:
            return None

    if len(values) == 0:
        return None

    shortcut = {
        "name": name,
        "names": names,
        "values": collect_values
    }

    return shortcut


def download_pubchem_molecule(name, scr=None):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{name}/JSON"
    data = get_json(url)

    if scr is not None:
        save_json(scr + name, data)

    return data


def get_pubchem_molecule(name, scr="_pubchem_/"):

    name = str(name)

    filename = scr + name

    if os.path.isfile(filename + ".json"):
        data = load_json(filename)

    else:
        data = download_pubchem_molecule(name, scr=scr)

    try:
        data = data["Record"]["Section"]
        headings = [x["TOCHeading"] for x in data]
        idx = headings.index("Names and Identifiers")
        section = data[idx]["Section"]
        headings = [x["TOCHeading"] for x in section]
        idx = headings.index("Computed Descriptors")
        section = section[idx]["Section"]
        descs = [x["Description"] for x in section]
        idxs = [i for i,x in enumerate(descs) if "SMILES" in x]
        idx = idxs[0]
        smiles = section[idx]["Information"][0]['Value']['StringWithMarkup'][0]["String"]
    except:
        smiles = None

    return smiles


def yield_data_rows(filenames):

    for filename in filenames:

        print(filename)

        data = load_json(filename)
        data = data["Annotations"]["Annotation"]

        for row in data:
            yield row

def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', action='store', help='', metavar="FILE", nargs="+")
    parser.add_argument('--scratch', action='store', help='', metavar="FILE")
    args = parser.parse_args()

    data = yield_data_rows(args.json)

    manual_list = []
    automatic_list = []

    for inx in data:

        x = parse_line(inx)

        if x is None:
            manual_list.append(inx)
            continue

        # if len(x["names"]) == 0:
        #     manual_list.append(inx)
        #     continue

        idxs = x["names"]
        chemname = x["name"]
        data = x["values"]

        smiles = None

        for idx in idxs:
            smiles = get_pubchem_molecule(idx)
            if smiles is not None:
                break

        data = list(filter(None, data))

        if len(data) == 0:
            manual_list.append(inx)
            continue

        # smiles, values, units
        # TODO Convert
        kelvins = []
        for val, un in data:

            if isfloat(val):
                newval = float(val)
            else:
                newval = parse_range(val)

            if un == "k":
                pass
            elif un == "c":
                newval += 273.15
            elif un == "f":
                newval += 459.67
                newval *= 5.0/9.0
            else:
                print("error: unknown unit", un)
                manual_list.append(inx)
                continue

            automatic_list.append((smiles, chemname.lower(), newval))

    # TODO Merge on keys, remove null

    name2smiles = {}
    smiles2value = {}

    null_list = []

    for datapoint in automatic_list:

        smiles, name, value = datapoint

        if smiles is None:
            null_list.append(datapoint)
            continue

        if smiles not in smiles2value:
            smiles2value[smiles] = []

        smiles2value[smiles].append(value)

        if name not in name2smiles:
            name2smiles[name] = smiles

    failed_smiles = []

    for datapoint in null_list:

        smiles, name, value = datapoint

        if name in name2smiles:
            smiles = name2smiles[name]
        else:
            failed_smiles.append(datapoint)
            continue

        smiles2value[smiles].append(value)

    #

    if args.scratch is not None:
        save_json(args.scratch + "/" + "success", automatic_list)
        save_json(args.scratch + "/" + "failed", manual_list)
        save_json(args.scratch + "/" + "values", smiles2value)
        save_json(args.scratch + "/" + "failed_smiles", failed_smiles)

    print()
    print("errors:", len(manual_list))
    print("success:", len(automatic_list))
    print()
    print("datapoints:", len(smiles2value.keys()))

    return


if __name__ == '__main__':
    main()
