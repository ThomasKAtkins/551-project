# code by Masha Larina

import pandas as pd

abolish_binding = [455, 456, 459, 474, 475, 486, 490, 493, 499]
enhance_binding = [439, 452, 470, 484, 498, 501]

values = [0 for _ in range(1274)]
labels = ["NA" for _ in range(1274)]

for i in abolish_binding:
    values[i] = 1
    labels[i] = "-"

for i in enhance_binding:
    values[i] = 1
    labels[i] = "+"

out = pd.DataFrame({"seq_position": list(range(1, 1275)), "value" : values, "affinity_label" : labels})
out.to_csv("data\\binding_affinity\\binding_info.csv", index=False)
