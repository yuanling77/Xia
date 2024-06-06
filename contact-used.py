import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB_small
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import pylab
import seaborn as sns
import math

#You can adjust the following parameters based on your file name and the number of frames in the trajectory. 
ori_pdb='output.pdb'   #Initial PDB structure name
trr='P132H-run.trr'             #Name of simulated trajectory file
d=4.0                               #Defined standard distance for atomic contact (nm)
f=8000                             #The frame number of trajectory
r=611                               #Number of residues, including ligands
vframe=55000                  #Cutoff for the number of contacts

u = mda.Universe(ori_pdb, trr)
Mpro = u.select_atoms ('resid 1-610')
nirma = u.select_atoms('resid 611')
atom=Mpro.atoms[1000]
atom.residue
atom.residue.resid
atom.name
n_Mpro=len(Mpro)
n_nirma=len(nirma)

print('Mpro has {} atoms and nirma has {} atoms'.format(n_Mpro, n_nirma))
def atom_contact_matrix (dist, frames, res_num):
  amino=np.zeros(dtype=float,shape=(1,res_num))
  dists = []
  for ts in u.trajectory:
    if ts.frame < frames:
      dist_arr = distances.distance_array(Mpro.positions, nirma.positions)
      for i in range(n_nirma):
        for j in range(n_Mpro):
          if dist_arr[j][i]<=dist:
            amino[0,Mpro.atoms[j].residue.resid]=1
  return amino          
        #dists.append(dist_arr)

amino_con_ID=atom_contact_matrix (d, f, r)
amino_contact=[]
for i in range(len(amino_con_ID[0])):
  if amino_con_ID[0,i]==1:
    amino_contact.append(i)
    


atoms_contact = Mpro.select_atoms('resid ' + ' '.join(str(resid) for resid in amino_contact))    
dist_lig = distances.distance_array(nirma.positions, atoms_contact.positions)
atom_res_value=np.zeros(dtype=float,shape=(len(dist_lig),len(amino_contact)))
atom_res_sum =pd.DataFrame(atom_res_value, columns=[resid for resid in amino_contact])

for ts in u.trajectory:
  if ts.frame < f:
    #atoms_contact=Mpro.select_atoms(amino_contact)
    atoms_contact=Mpro.select_atoms('resid' + ' '+' '.join(str(resid) for resid in amino_contact))
    dist_lig = distances.distance_array(nirma.positions, atoms_contact.positions)
    atom_res_value=np.zeros(dtype=float,shape=(len(dist_lig),len(amino_contact)))
    atom_res_contact =pd.DataFrame(atom_res_value, columns=[resid for resid in amino_contact])
    for i in range(len(dist_lig)):
      for j in range(len(dist_lig[0])):
        if dist_lig[i][j] <= d:
          atom_res_contact.loc[[i],[atoms_contact[j].resid]]+=1
    atom_res_sum=atom_res_sum+atom_res_contact 
    
print('The max number of contact is:)',atom_res_sum)  
  
numb=[14, 7, 4, 37, 10, 39, 40, 13, 41, 15, 42, 43, 16, 44, 45, 5, 38, 17, 3, 2, 36, 1, 6, 18, 46, 8, 19, 47, 48, 20, 49, 21, 50, 27, 29, 62, 63, 61, 28, 59, 58, 60, 22, 9, 23, 51, 24, 25, 53, 54, 52, 26, 56, 55, 57, 32, 66, 64, 65, 11, 67, 30, 12, 31, 34, 33, 35]
atom_res_sum_new=pd.DataFrame(atom_res_value, columns=[resid for resid in amino_contact])
for i in range(len(numb)):
  atom_res_sum_new.loc[i]=atom_res_sum.loc[numb[i]-1]
  
  
  
atom_res_sum_new.to_excel('atom_res_con.xlsx', index=True)

zp=pd.read_excel('atom_res_con.xlsx',index_col=0)
sns.heatmap(zp, annot=False, cmap='rocket_r', xticklabels=True, yticklabels=True, vmin=0, vmax=vframe)
plt.savefig('heatmap.png', dpi=600)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Heatmap')
plt.xlabel('Columns')
plt.ylabel('Rows')
plt.show()    
