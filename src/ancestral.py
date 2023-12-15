# ML

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# fraction of ancestral at non-mutation positions locations
peptide_df = pd.read_csv("data/NetMHCpan/peptides.csv")

df_vals = peptide_df[peptide_df["ancestral"]==True]

xvals = list(df_vals["pos"])
yvals = list(df_vals["freq"])

idx = np.argsort(yvals)

sorted_yvals = np.array(yvals)[idx]
sorted_xvals = np.array(xvals)[idx]


plt.scatter(xvals, sorted_yvals, s=1, color="blue")

# order x values by increasing y values
abolish_binding = [455, 456, 459, 474, 475, 486, 490, 493, 499]
enhance_binding = [439, 452, 470, 484, 498, 501]

abolish_binding_2 = [list(sorted_xvals).index(x) for x in abolish_binding]
enhance_binding_2 = [list(sorted_xvals).index(x) for x in enhance_binding]

plt.scatter(abolish_binding_2, [sorted_yvals[x] for x in abolish_binding_2], color="red", label="abolish binding", s=20)
plt.scatter(enhance_binding_2, [sorted_yvals[x] for x in enhance_binding_2], color="green", label="enhance binding", s=20)
plt.legend()
plt.title("Ancestral frequency at mutation-sensitive RBD sites")
plt.ylabel("Frequency of ancestral peptide")
plt.xlabel("Sequence position, ordered by increasing ancestral frequency")

plt.show()

# # fraction of ancestral at non-mutation positions locations
# peptide_df = pd.read_csv("data/NetMHCpan/peptides.csv")

# df_vals = peptide_df[peptide_df["ancestral"]==True]

# xvals = list(df_vals["pos"])
# yvals = list(df_vals["freq"])


# plt.scatter(df_vals["pos"], df_vals["freq"], s=1, color="blue")

# # order x values by increasing y values
# abolish_binding = [455, 456, 459, 474, 475, 486, 490, 493, 499]
# enhance_binding = [439, 452, 470, 484, 498, 501]

# plt.scatter(abolish_binding, [list(df_vals["freq"])[x] for x in abolish_binding], color="red", label="abolish binding", s=20)
# plt.scatter(enhance_binding, [list(df_vals["freq"])[x] for x in enhance_binding], color="green", label="enhance binding", s=20)
# plt.legend()
# plt.title("Ancestral frequency at mutation-sensitive RBD sites")
# plt.ylabel("Frequency of ancestral peptide")
# plt.xlabel("Sequence position")

# plt.show()






