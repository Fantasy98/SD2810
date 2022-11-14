#%%
from scipy.io import loadmat
Data = loadmat("dlmqip.mat")
print(Data.keys())

# %%
qdlm = Data["Qdlm"]
print(qdlm.shape)
# %%
