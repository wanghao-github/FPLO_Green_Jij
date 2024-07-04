import numpy as np
from math import pi
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3d


def gen_kpath(file_path):
    full_path = file_path + 'syml'
    with open(full_path,'r') as f:
        lines = f.readlines()
        
    num_high_sym_points = int(lines[0].strip())
    num_kpts_on_path_values = list(map(int, lines[1].strip().split()))
    num_kpts_on_path = np.array(num_kpts_on_path_values)

    high_sym_point_dict = {}
    k_symbol = []
    
    for i in range(num_high_sym_points):
        key = str(lines[2+i].strip()[0])
        value = np.array(list(map(float,(lines[2+i].strip()[1:].split()))))
        high_sym_point_dict[key] = value
        k_symbol.append(key)
    num_high_sym_path = num_high_sym_points - 1

    kpath = {}
    for i in range(num_high_sym_path):
        kpath[i] = np.zeros((num_kpts_on_path[i],3))
    
    for i in range(num_high_sym_path):
        # print(i)
        for j in range(num_kpts_on_path[i]):
            # print(j)
            for k in range(3):
                # print(k)
                if i != range(num_high_sym_path):
                    # print(kpath[i][j][k])
                    kpath[i][j][k] = j * ((high_sym_point_dict[k_symbol[i+1]][k]-high_sym_point_dict[k_symbol[i]][k])/num_kpts_on_path[i]) + high_sym_point_dict[k_symbol[i]][k]
    kpoint_array = np.vstack([kpath[key] for key in sorted(kpath.keys())])
    
    return kpoint_array


def read_wanbandtb(file_path):
    full_path = file_path + '+wanbandtb'
    
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        kpoints_line = lines[1::2]
        energy_line = lines[2::2]
        
        kpts_num = len(kpoints_line)
        bands_num = len(energy_line[0].split()) - 1
        
        kpoints = np.zeros((kpts_num, 3))
        bands = np.zeros((kpts_num, bands_num))
        kdist = np.zeros(kpts_num)
        
        
        for i, line in enumerate(kpoints_line):
            values = line.split()
            for j in range(3):
                kpoints[i][j] = float(values[j+1])
        
        for k, line in enumerate(energy_line):
            values = line.split()
            kdist[k] = values[0]
            for n in range(bands_num):
                bands[k][n] = float(values[n+1])
        
        return kpoints, bands ,kdist

    except FileNotFoundError:
        print(f"File not found: {full_path}")
        return None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None



    
def read_hamdata(file_path):
    full_path = file_path + '+hamdata'
    with open(full_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    nwan_value = None
    lattice_vectors = []
    tij_hij_blocks = []
    wancenters = []
    block_lengths = []
    p_wfpairs_list = []
    for i, line in enumerate(lines):
        if 'nwan:' in line:
            if i + 1 < len(lines):
                nwan_value = int(lines[i + 1].strip())

        if 'lattice_vectors:' in line:
            for j in range(1, 4):
                if i + j < len(lines):
                    vector = list(map(float, lines[i + j].strip().split()))
                    lattice_vectors.append(vector)
                    
        if 'wancenters:' in line:
            for j in range(1, int(nwan_value) + 1):
                if i + j < len(lines):
                    center = list(map(float, lines[i + j].strip().split()))
                    wancenters.append(center)
                    
        # if 'Tij, Hij:' in line:
        if 'Tij, Hij, Sij:' in line:
            block = []
            # p_wfpairs = np.empty((3,nwan_value, nwan_value,))
            try:
                index1, index2 = map(int, lines[i + 1].strip().split())

                for j in range(2, len(lines) - i):
                    if 'end Tij, Hij:' in lines[i + j]:
                        break
                    try:
                        block.append(list(map(float, lines[i + j].strip().split())))
                    except ValueError as e:
                        # print(f"Error parsing data at line {i + j}: {lines[i + j].strip()}")
                        continue
                    
                if block:
                    tij_hij_blocks.append((index1, index2, block))
                    block_lengths.append((index1, index2,len(block)))
                    
                    for irs in range(len(block)):
                        p_wfpairs_list.append((block[irs][:3],index1, index2, block[irs][3],block[irs][4],irs)) 
                    
                    
            except ValueError as e:
                # print(f"Error parsing indices at line {i + 1}: {lines[i + 1].strip()}")
                continue
          
    lattice_vectors = np.array(lattice_vectors)
    lattice_moduli = np.linalg.norm(lattice_vectors, axis=1)
    p_wfpairs = np.array(p_wfpairs_list, dtype=object)
    return nwan_value, lattice_vectors, lattice_moduli, tij_hij_blocks, wancenters,block_lengths,p_wfpairs

def nrpts_mapping(supercell_size,lattice_vectors,wancenters):
    lattice_vectors

def is_hermitian(matrix):
    return np.allclose(matrix, np.conj(matrix.T), atol=1e-4,)

def process_tij_hij_blocks(tij_hij_blocks, lattice_vectors, lattice_moduli, nwan_value, wancenters,p_wfpairs):
    # unique_r_vectors = set()
    unique_r_vectors = []
    inv_lattice_vectors = np.linalg.inv(lattice_vectors)
    r_vector_count = {}
    
    # for p_wfpairs
    
    for (index1, index2, block) in tij_hij_blocks:
        for row in block:
            r = np.array(row[:3])
            si = np.array(wancenters[index1-1])#+np.array([0,13.217662,0])
            sj = np.array(wancenters[index2-1])#+np.array([7.6312,0,0])
            R = r - (sj - si)
            # R = r
            # T=R
            T = np.dot(R, inv_lattice_vectors)
            # T = np.dot(lattice_vectors,R)#
            # T_int = tuple(T.astype(int))
            # T_int = T
            T_rounded = np.rint(T)
            # 将四舍五入后的浮点数转换为整数
            T_int = T_rounded.astype(int)
            # unique_r_vectors.append(T_int)
            unique_r_vectors.append((T_int[0],T_int[1],T_int[2],index1,index2))
            # if T_int in r_vector_count:
            #     r_vector_count[T_int] += 1
            # else:
            #     r_vector_count[T_int] = 1
    unique_r_vectors = np.array(unique_r_vectors)

    data_dict = {}
    
    for r_vector in unique_r_vectors:
        r1, r2, r3 = r_vector[:3]
        for i in range(1, nwan_value + 1):
            for j in range(1, nwan_value + 1):
                data_dict[(r1, r2, r3, i, j)] = (0.0, 0.0)
    
    for (index1, index2, block) in tij_hij_blocks:
        for row in block:
            r = np.array(row[:3])
            si = np.array(wancenters[index1-1])#+np.array([0,13.217662,0])
            sj = np.array(wancenters[index2-1])#+np.array([7.6312,0,0])
            R = r - (sj - si)
            T = np.dot(R, inv_lattice_vectors)
            T_rounded = np.rint(T)
            T_int = T_rounded.astype(int)
            real = row[3]
            img = row[4]

            key = (T_int[0], T_int[1], T_int[2], index1, index2)
          
            # data_dict[key] = (real, img)
            if key in data_dict:
                existing_real, existing_img = data_dict[key]
                data_dict[key] = (existing_real + real, existing_img + img)
            else:
                data_dict[key] = (real, img)
    processed_data = [(*k, *v) for k, v in data_dict.items()]
    
    r_vector_count_list = sorted(r_vector_count.items())
    
    return processed_data, unique_r_vectors,r_vector_count_list

def save_to_wannier90_centres_xyz(file_path,data,nwan):
    full_path = file_path + 'wannier90_centres.xyz'
    header = "Generated by Hao Wang's script\n"
    with open(full_path, 'w', encoding='utf-8') as f:
        f.write(f"{nwan}\n")
        f.write(header)
        for line in data:
            f.write(
                f"X          {line[0]:18.8f} {line[1]:18.8f} {line[2]:18.8f}\n"
                )
        f.write(f"Bi         2.367314417586326     -4.100308848749840     29.667930445390219 \n")
        f.write(f"Bi        -2.367314417586326      4.100308848749840    -29.667930445390219 \n")

def save_to_wannier_hr(file_path, data,nwan,nrpts):
    full_path = file_path + 'wannier90_hr.dat'
    header = "Generated by Hao Wang's script\n"
    with open(full_path, 'w', encoding='utf-8') as f:
        f.write(header)
        f.write(f"{nwan}\n")
        f.write(f"{nrpts}\n")

        # Write nrpts 1s, 15 per line
        ones_line = '    1' * 15
        for _ in range(nrpts // 15):
            f.write(f"{ones_line}\n")
        if nrpts % 15 > 0:
            f.write(f"{'    1' * (nrpts % 15)}\n")
            
        # for line in data:
        #     f.write(
        #         f"{line[0]:5f} {line[1]:5f} {line[2]:5f} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
        #     )
        # for line in data:
        #     f.write(
        #         f"{int(line[0]):5d} {int(line[1]):5d} {int(line[2]):5d} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
        #     )
        for line in data:
            f.write(
                f"{line[0]:5.0f} {line[1]:5.0f} {line[2]:5.0f} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
            )
            
            
            
    # full_path1 = file_path + 'wannier90_hr1.dat'
    # header = "Generated by Hao Wang's script\n"
    # with open(full_path1, 'w', encoding='utf-8') as f1:
    #     f1.write(header)
    #     f1.write(f"{nwan}\n")
    #     f1.write(f"{nrpts}\n")

    #     # Write nrpts 1s, 15 per line
    #     ones_line = '    1' * 15
    #     for _ in range(nrpts // 15):
    #         f1.write(f"{ones_line}\n")
    #     if nrpts % 15 > 0:
    #         f1.write(f"{'    1' * (nrpts % 15)}\n")
            
    #     # for line in data:
    #     #     f.write(
    #     #         f"{line[0]:5d} {line[1]:5d} {line[2]:5d} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
    #     #     )
    #     for line in data:
    #         f1.write(
    #             f"{int(line[0]):5d} {int(line[1]):5d} {int(line[2]):5d} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
    #         )
        # for line in data:
        #     f1.write(
        #         f"{line[0]:12.0f} {line[1]:12.0f} {line[2]:12.0f} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
        #     )

# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\mos2\+hamdata'#klist6
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\mos2\bulk\+hamdata'#klist7
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\bi\+hamdata' #klist5
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\Fe\automode\+hamdata'#klist9
file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\CrI3\z\\'#klist3
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\Fe\fplo18\\'
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\Fe\conv\+hamdata'   #klist10
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\Fe\\'#klist8
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\BHZ\+hamdata'
# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\SrMn2Bi2\\'#klist11
nwan, lattice_vectors, lattice_moduli, tij_hij_blocks, wancenters ,block_lengths, p_wfpairs= read_hamdata(file_path)

if nwan:
    print(f'nwan: {nwan}')
else:
    print('nwan: not found')

if lattice_vectors.size != 0:
    print('lattice_vectors:')
    print(lattice_vectors)
else:
    print('lattice_vectors: not found')

processed_data, unique_r_vectors,r_vector_count_list = process_tij_hij_blocks(tij_hij_blocks, lattice_vectors, lattice_moduli, nwan, wancenters,p_wfpairs)

# unique_r_vectors_set = set(unique_r_vectors)
# nrpts = len(unique_r_vectors_set)
# # print(nrpts)

print(r_vector_count_list)

unique_r_vectors_set = set(map(tuple, unique_r_vectors[:, :3]))
nrpts = len(unique_r_vectors_set)

# unique_r_vectors_set = set(map(tuple, unique_r_vectors))
# nrpts = len(unique_r_vectors_set)

sorted_data = sorted(processed_data, key=lambda x: (x[0], x[1], x[2]))

save_to_wannier_hr(file_path, sorted_data,nwan,nrpts)
save_to_wannier90_centres_xyz(file_path,wancenters,nwan)

# print(f"Number of unique R vectors: {len(unique_r_vectors)}")
# print("Unique R vectors:")

# for R in unique_r_vectors:
#     print(R)
# print(nrpts)





def read_wannier90_hr(file_path):
    full_path = file_path + 'wannier90_hr.dat'
    # with open(r"C:\Users\wangh\Documents\HaoWang\Hao_code\from_Hr_plot_band_python_Tprime_MoS2\WF1_hr.dat","r") as f:
    with open(full_path,'r') as f:
        lines=f.readlines()
        num_wann=int(lines[1]); nrpts=int(lines[2])
        lines_nrpts=int(np.ceil(nrpts/15.0))
        ndegen_list = []
        for i in range(lines_nrpts):
            for j in range(len(lines[3+i].split())):
                ndegen_list.append(int(lines[3+i].split()[j]))
        ndegen = np.array(ndegen_list)
        # HmnR = []
        HmnR_np = np.zeros((num_wann**2*nrpts,7))
        for i in range(num_wann**2*nrpts): 
            for j in range(7):
                HmnR_np[i][j] = (float(lines[3+lines_nrpts+i].split()[j]))
        # HmnR_np = np.array(HmnR)
        # HmnR_np = HmnR_np.reshape((num_wann**2*nrpts,7))
        # print(np.shape(HmnR_np))
        HmnR_np_iR = np.zeros((num_wann,num_wann,nrpts),dtype=complex)
        irvec = np.zeros((nrpts,3))
        for ir in range(nrpts):
            for n in range(num_wann):
                for m in range(num_wann):
                    HmnR_np_iR[m,n,ir] = complex(HmnR_np[ir*num_wann**2+n*num_wann+m][5],HmnR_np[ir*num_wann**2+n*num_wann+m][6])
            irvec[ir][0]=HmnR_np[ir*num_wann**2][0]
            irvec[ir][1]=HmnR_np[ir*num_wann**2][1]
            irvec[ir][2]=HmnR_np[ir*num_wann**2][2]

    return num_wann, ndegen, irvec, HmnR_np_iR
lat=lattice_vectors      
# lat=np.array([[5.415957,0,0],[0,5.415957,0],[0,0,5.415957]])       
        
  
kpath1, band_wanbandtb ,kdist= read_wanbandtb(file_path)
num_wann, ndegen, irvec, HmnR_np_iR = read_wannier90_hr(file_path)
# kpath1=np.loadtxt('klist11.txt')

# band_structure2 = np.loadtxt('band_data11.txt')
# band_structure2.reshape((len(kpath1),-1))
# band_structure2 = band_structure2[:,1:]

rep = np.linalg.inv(lat)
# new_k = np.dot(kpath1,rep)
def get_bands_from_ham(num_wann, kpath,irvec,ndegen,lat,lattice_moduli,HmnR_np_iR):
    band_structure = np.zeros((len(kpath),num_wann))
    for i in range(len(kpath)):
        Ham_bulk = np.zeros((num_wann,num_wann),dtype=complex)
        for iR in range(nrpts):
            ia = irvec[iR][0]
            ib = irvec[iR][1]
            ic = irvec[iR][2]
            
            R = ia * lat[0,:] +ib * lat[1,:] + ic * lat[2,:]
            
            kdotR = np.dot(R,kpath[i])

            factor = np.exp(1j*kdotR*2*np.pi/lattice_moduli[0])
            Ham_bulk[:,:] = Ham_bulk[:,:] +(HmnR_np_iR[:,:,iR] * factor/ndegen[iR])
            
        matrix_to_check = Ham_bulk[:,:]
        if is_hermitian(matrix_to_check):
            pass
            # print("The matrix is Hermitian.")
            # print("r is abc ", ia,ib,ic)
        else:
            pass
            # print("The matrix is not Hermitian.")
        eigen_value, eigen_vector = np.linalg.eig(Ham_bulk)
        sorted_eig = np.sort(np.real(eigen_value))
        band_structure[i][:] = sorted_eig
    return band_structure
x_1d=np.arange(len(kpath1))
plt.figure(figsize=(10,12))
# plt.xlim(0,440)


band_structure = get_bands_from_ham(num_wann, kpath1,irvec,ndegen,lat,lattice_moduli,HmnR_np_iR)
plt.ylim(-10,3)
for i in range(num_wann):
    plt.plot(kdist,band_wanbandtb[:,i], '-', color='red', linewidth=4)
    plt.plot(kdist,band_structure[:,i],'o', color='blue', markersize=2)
# plt.savefig(r'C:\Users\wangh\Documents\HaoWang\Hao_code\from_Hr_plot_band_python_Tprime_MoS2\from_Hr_T_prime_MoS2_2.png')
plt.savefig(file_path + 'biband.png')
plt.show()
print(len(r_vector_count_list))
# print(block_lengths)
# # print(len(block_lengths))
# print(p_wfpairs)

# def save_to_p_wfpairs(data):
#     with open('pwfpairs', 'w', encoding='utf-8') as f:
#         for line in data:
#             f.write(
#                 f"{line[0]} {line[1]:5d} {line[2]:5d} {int(line[3]):8d}\n"
#             )
# save_to_p_wfpairs(p_wfpairs)

scale = lattice_moduli[0]

print(p_wfpairs)
p_wfpairs1 = p_wfpairs[:, [0, 1, 2, 5]]

# 获取2, 3, 4, 5列
p_hopsym = p_wfpairs[:, [1, 2, 3, 4, 5]]

print(p_wfpairs1[0,:])
# def reassemble_Ham(nwan,kapth,p_hopsym,block_lengths,p_wfpairs1,scale):
#     p_H = np.zeros((nwan, nwan), dtype=complex)
#     ci = 1j
#     for q in range(len(kapth)):
#         for iwan in range(nwan):
#             for jwan in range(nwan):
#                 p_nwfpairs_value = next((value for iw, jw, value in block_lengths if iw == iwan and jw == jwan), 0)
#                 for irs in range(p_nwfpairs_value):
#                     condition1 = (p_wfpairs1[:, 1] == iwan) & (p_wfpairs1[:, 2] == jwan) & (p_wfpairs1[:, 3] == irs)
#                     r = p_wfpairs1[condition1]
#                     condition2 = (p_hopsym[:, 0] == iwan) & (p_hopsym[:, 1] == jwan) & (p_hopsym[:, 4] == irs)
#                     cH = p_hopsym[condition2]
#                     qr = np.dot(kapth[q],r ) * scale
#                     p_H[iwan, jwan] += cH * np.exp(ci * qr)
# def reassemble_Ham(nwan, kapth, p_hopsym, block_lengths, p_wfpairs1, scale):
#     p_H = np.zeros((nwan, nwan), dtype=complex)
#     ci = 1j
    
#     for q in range(len(kapth)):
#         for iwan in range(nwan):
#             for jwan in range(nwan):
#                 p_nwfpairs_value = next((value for iw, jw, value in block_lengths if iw == iwan and jw == jwan), 0)
#                 for irs in range(p_nwfpairs_value):
#                     condition1 = (p_wfpairs1[:, 1] == iwan) & (p_wfpairs1[:, 2] == jwan) & (p_wfpairs1[:, 3] == irs)
#                     r = p_wfpairs1[condition1]
                    
#                     if r.size > 0:
#                         r_value = r[0, 0]  # 获取匹配行的第一个值
                        
#                         condition2 = (p_hopsym[:, 0] == iwan) & (p_hopsym[:, 1] == jwan) & (p_hopsym[:, 4] == irs)
#                         cH = p_hopsym[condition2]
                        
#                         if cH.size > 0:
#                             cH_real = cH[0, 2]  # 获取匹配行的第三列作为实部
#                             cH_imag = cH[0, 3]  # 获取匹配行的第四列作为虚部
#                             cH_complex = cH_real + ci * cH_imag  # 构建复数
#                             qr = np.dot(kapth[q], r_value) * scale
#                             p_H[iwan, jwan] += cH_complex * np.exp(ci * qr)
    
#     return p_H

# # 调用函数
# p_H = reassemble_Ham(nwan, kpath1, p_hopsym, block_lengths, p_wfpairs1, scale)
# print(p_H)

print(block_lengths)
p_nwfpairs = np.zeros((5, 5), dtype=int)
for x, y, value in block_lengths:
    p_nwfpairs[x-1, y-1] = value
    
print(p_wfpairs)
print(block_lengths)