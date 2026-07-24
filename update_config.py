"""Update config.py with Lung and Breast cancer definitions."""
import re

with open("config.py") as f:
    c = f.read()

# Add Lung + Breast cancer types
old = """"Cervical": [
        "Cervical Cancer", "Cervical Squamous Cell Carcinoma",
        "Cervical Adenocarcinoma",
    ],
}"""
new = """"Cervical": [
        "Cervical Cancer", "Cervical Squamous Cell Carcinoma",
        "Cervical Adenocarcinoma",
    ],
    "Lung": [
        "Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma",
        "Non-Small Cell Lung Cancer", "Small Cell Lung Cancer",
        "Lung Cancer",
    ],
    "Breast": [
        "Breast Invasive Carcinoma", "Breast Invasive Ductal Carcinoma",
        "Breast Invasive Lobular Carcinoma", "Breast Cancer",
    ],
}"""
c = c.replace(old, new)

# Add driver genes
old2 = """"Cervical":    ["PIK3CA", "EP300", "FBXW7", "STK11", "ERBB2",
                     "MAPK1", "PTEN", "KRAS"],
}"""
new2 = """"Cervical":    ["PIK3CA", "EP300", "FBXW7", "STK11", "ERBB2",
                     "MAPK1", "PTEN", "KRAS"],
    "Lung":        ["TP53", "KRAS", "EGFR", "STK11", "KEAP1",
                     "NF1", "BRAF", "PIK3CA", "ALK", "MET"],
    "Breast":      ["TP53", "PIK3CA", "PTEN", "BRCA1", "ERBB2",
                     "GATA3", "CDH1", "RB1", "NF1", "MAP3K1"],
}"""
c = c.replace(old2, new2)

with open("config.py", "w") as f:
    f.write(c)
print("Config updated: Lung + Breast added")
