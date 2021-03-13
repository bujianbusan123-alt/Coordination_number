from ase.io import read
from ase.data import covalent_radii as CR
import sys
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

def get_coordination_numbers(atoms, covalent_percent=1.25):
    """Returns an array of coordination numbers and an array of existing bonds determined by
    distance and covalent radii.  By default a bond is defined as 120% of the combined radii
    or less. This can be changed by setting 'covalent_percent' to a float representing a 
    factor to multiple by (default = 1.2).
    If 'exclude' is set to an array,  these atomic numbers with be unable to form bonds.
    This only excludes them from being counted from other atoms,  the coordination
    numbers for these atoms will still be calculated,  but will be unable to form
    bonds to other excluded atomic numbers.
    """
    fig = plt.figure()
    cm_all = open('cm_all.txt', 'w')
    for atom in atoms:
        atoms = read(atom)
        cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, cm11, cm12 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        cm = {'cm1':cm1, 'cm2': cm2, 'cm3': cm3, 'cm4':cm4,
              'cm5': cm5, 'cm6': cm6, 'cm7': cm7, 'cm8':cm8,
              'cm9': cm9, 'cm10': cm10, 'cm11': cm11, 'cm12': cm12}
    #for 

        # Get all the distances
        distances = np.divide(atoms.get_all_distances(mic=True), covalent_percent)
    
        # Atomic Numbers
        numbers = atoms.numbers
        print(len(numbers))
        #numbers = len(atoms)
        # Coordination Numbers for each atom
        cn = []
        cr = np.take(CR, numbers)
        # Array of indices of bonded atoms.  len(bonded[x]) == cn[x]
        bonded = []
        indices = list(range(len(atoms)))
        for i in indices:
            bondedi = []
            for ii in indices:
                # Skip if measuring the same atom
                if i == ii:
                    continue
                if (cr[i] + cr[ii]) >= distances[i,ii]:
                    bondedi.append(ii)
            # Add this atoms bonds to the bonded list
            bonded.append(bondedi)
        for i in bonded:
            cn.append(len(i))
        print(len(cn))
        for j in cn:
            if j == 1:
                cm['cm1'] += 1
            elif j == 2:
                cm['cm2'] += 1
            elif j == 3:
                cm['cm3'] += 1
            elif j == 4:
                cm['cm4'] += 1
            elif j == 5:
                cm['cm5'] += 1
            elif j == 6:
                cm['cm6'] += 1
            elif j == 7:
                cm['cm7'] += 1
            elif j == 8:
                cm['cm8'] += 1
            elif j == 9:
                cm['cm9'] += 1
            elif j == 10:
                cm['cm10'] += 1
            elif j == 11:
                cm['cm11'] += 1
            elif j == 12:
                cm['cm12'] += 1
        for indices, values in cm.items():
            cm[indices] = float('%5.2f' % (values/len(numbers)* 100)) #设置小数点位数
        print(cm) #切记返回有时候不能采用
        plt.plot([x for x in range(1, 13)], [values for indices, values in cm.items()], label=atom)
        cm_all.write(atom + '\n')
        cm_all.write(str(cm) + '\n')
        print('\n')
        
    plt.ylabel('N/Ntol(%)', fontdict={'size':15}, fontweight='medium')
    plt.xlabel('Coordination Numbers', fontdict={'size':15}, fontweight='medium')
    plt.legend(prop={'size':10}) 
    plt.tight_layout()
    plt.savefig('cn.png', dpi=600)
    plt.show()
    cm_all.close()
   # print(bonded)
if __name__ == '__main__':
    atoms = ['np.xyz', 'long-1nm.xyz',  'long-15nm.xyz', 'long-2nm.xyz']
    get_coordination_numbers(atoms, covalent_percent=1.25)
