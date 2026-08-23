import sys
from pathlib import Path
import pandas as pd
from rdkit import Chem
from biosynfoni import Biosynfoni

# use python validate.py <path-to-coconut-sdf> to set file path to the data set
if len(sys.argv) != 2:
    print("Usage: python compare.py <path-to-coconut-sdf>")
    sys.exit(1)

data_path = sys.argv[1]

cases = {
    "no_filter": {
        "intrasub_overlap": False,
        "intersub_overlap": False
    },
    "inter_filter": {
        "intrasub_overlap": False,
        "intersub_overlap": True
    },
    "intra_filter": {
        "intrasub_overlap": True,
        "intersub_overlap": False
    }
}

for case_name, settings in cases.items():

    data = []

    supplier = Chem.SDMolSupplier(
        data_path,
        sanitize=False,
        removeHs=False
    )

    limit = 738824 #if the limit is greater than the Data it stops at the End

    for i, mol in enumerate(supplier):
        if i >= limit:
            break

        if mol is None:
            continue

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        Chem.RemoveStereochemistry(mol)

        bsf = Biosynfoni(
            mol,
            intrasub_overlap=settings["intrasub_overlap"],
            intersub_overlap=settings["intersub_overlap"]
        )

        countFingerprint = bsf.fingerprint

        smiles = Chem.MolToSmiles(
            mol,
            canonical=False,
            isomericSmiles=False
        )

        row = [smiles] + list(countFingerprint)
        data.append(row)

    fp_size = len(data[0]) - 1
    columns = ["smiles"] + [f"count_{i}" for i in range(fp_size)]

    df_out = pd.DataFrame(data, columns=columns)

    df_out.to_csv(
        f"dataPython_{case_name}.csv",
        index=False
    )

    print(f"Created dataPython_{case_name}.csv")


