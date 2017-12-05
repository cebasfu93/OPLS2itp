import numpy as np
from optparse import OptionParser

masses={'H':'1.00800', 'N':'14.01000', 'C':'12.01000', 'O':'16.00000', 'S':'32.06000'}
name_ndx={}
#Opciones del script
parser=OptionParser()
parser.add_option("-i","--input", action="store", type="string", dest="InputFile", default='test.opls', help="Path of the .opls input file.")

(options, args) = parser.parse_args()

inp_name=options.InputFile
resname=inp_name[:3].upper()
inp_file=np.genfromtxt(inp_name, delimiter='\n', dtype='string')

itp=open(resname+".itp", "w")

############Atom types###########
itp.write("[ atomtypes ]\n;name   bond_type     mass     charge   ptype    sigma         epsilon\n")
for i  in range(len(inp_file)):
    if "OPLSAA"in inp_file[i]:
        ini=i+4
    if "Stretch"in inp_file[i]:
        fin=i-1
N_atoms=fin-ini
types=[]
for i in range(N_atoms):
    atom_types=inp_file[ini+i].split()
    types.append(atom_types[3])
unique_types=list(set(types))
for j in range(len(unique_types)):
    for i in range(N_atoms):
        atom_types=inp_file[ini+i].split()
        if atom_types[3]==unique_types[j]:
            itp.write(" "+atom_types[3].ljust(9) + atom_types[3].ljust(12) + "0.00000  0.00000   A\t" + atom_types[5].ljust(14) + atom_types[6].ljust(12) + "\n")
            break

############Molecule type###########
itp.write("\n[ moleculetype ] \n; molname         nrexcl \n "+ resname.ljust(18) + "3\n\n")

############Atoms###########
itp.write("[ atoms ] \n; id  at_type     res_nr  res_name  at_name  cg_nr  charge    mass \n")
for i in range(N_atoms):
    atom_types=inp_file[ini+i].split()
    name_ndx[atom_types[0]]=str(i+1)
    itp.write(str(i+1).rjust(6) + atom_types[3].rjust(5)+ "1".rjust(6)+resname.rjust(6)+atom_types[0].rjust(6) + str(i+1).rjust(5)+atom_types[4].rjust(13)+ masses[atom_types[3][0]].rjust(13)+"\n")

############Stretch###########
itp.write("\n[ bonds ] \n; i     j       funct   length  force_constant\n")
ini=fin+2
for i  in range(len(inp_file)):
    if "Bending"in inp_file[i]:
        fin=i
        break
N_bonds=fin-ini
for i in range(N_bonds):
    bonds=inp_file[ini+i].split()
    itp.write(name_ndx[bonds[0]].rjust(6) + name_ndx[bonds[1]].rjust(7)+"1".rjust(4)+ bonds[3].rjust(14)+bonds[2].rjust(14)+"\n")

############Bending###########
itp.write("\n[ angles ] \n; i     j       k       funct   angle   force_constant\n")
ini=fin+1
for i in range(len(inp_file)):
    if "proper Torsion" in inp_file[i]:
        fin=i
        break
print inp_file[ini]
print inp_file[fin]
N_angles=fin-ini
for i in range(N_angles):
    angles=inp_file[ini+i].split()
    itp.write(name_ndx[angles[0]].rjust(6)+ name_ndx[angles[1]].rjust(7)+ name_ndx[angles[2]].rjust(7)+ "1".rjust(7)+angles[4].rjust(14)+angles[3].rjust(14)+"\n")

############Proper-Dihedrals###########
itp.write("\n[ dihedrals ] ; propers \n; i      j       k       l       func   phase   kd    pn\n")
ini=fin+1
for i in range(len(inp_file)):
    if "improper Torsion" in inp_file[i]:
        fin=i
        break
N_dih=fin-ini
for i in range(N_dih):
    dihedrals=inp_file[ini+i].split()
    itp.write(name_ndx[dihedrals[0]].rjust(6)+ name_ndx[dihedrals[1]].rjust(7)+ name_ndx[dihedrals[2]].rjust(7)+name_ndx[dihedrals[3]].rjust(7)+"9".rjust(7)+"0.0".rjust(9)+"0.0".rjust(10)+"0".rjust(4)+"\n")

############Improper-Dihedrals###########
itp.write("\n[ dihedrals ] ; impropers \n; i      j       k       l       func   phase   kd    pn\n")
ini=fin+1
fin=len(inp_file)
N_imdih=fin-ini
for i in range(N_imdih):
    im_dihedrals=inp_file[ini+i].split()
    itp.write(name_ndx[im_dihedrals[0]].rjust(6)+ name_ndx[im_dihedrals[1]].rjust(7)+ name_ndx[im_dihedrals[2]].rjust(7)+name_ndx[im_dihedrals[3]].rjust(7)+"4".rjust(7)+"0.0".rjust(9)+"0.0".rjust(10)+"0".rjust(4)+"\n")

itp.close()
