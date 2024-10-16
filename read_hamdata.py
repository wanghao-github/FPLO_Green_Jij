# hamiltonian_processor.py

import numpy as np
import matplotlib.pyplot as plt

class HamiltonianProcessor:
    def __init__(self, file_path):
        
        self.file_path = file_path
        self.Ham_bulk = None
        self.nwan, self.lattice_vectors, self.centering, self.lattice_moduli, self.tij_hij_blocks, self.wancenters, self.block_lengths, self.p_wfpairs = self.read_hamdata()
        self.nrpts, self.processed_data, self.unique_r_vectors, self.r_vector_count_list = self.process_tij_hij_blocks()
        # self.MI, self.Mx, self.My, self.Mz = self.pauli_block_all()
        
        
    def gen_kpath(self):
        full_path = self.file_path + 'syml'
        with open(full_path, 'r') as f:
            lines = f.readlines()
        num_high_sym_points = int(lines[0].strip())
        num_kpts_on_path_values = list(map(int, lines[1].strip().split()))
        num_kpts_on_path = np.array(num_kpts_on_path_values)
        high_sym_point_dict = {}
        k_symbol = []

        for i in range(num_high_sym_points):
            key = str(lines[2 + i].strip()[0])
            value = np.array(list(map(float, (lines[2 + i].strip()[1:].split()))))
            high_sym_point_dict[key] = value
            k_symbol.append(key)
        num_high_sym_path = num_high_sym_points - 1
        kpath = {}
        for i in range(num_high_sym_path):
            kpath[i] = np.zeros((num_kpts_on_path[i], 3))
        for i in range(num_high_sym_path):
            for j in range(num_kpts_on_path[i]):
                for k in range(3):
                    if i != range(num_high_sym_path):
                        kpath[i][j][k] = j * ((high_sym_point_dict[k_symbol[i + 1]][k] - high_sym_point_dict[k_symbol[i]][k]) / num_kpts_on_path[i]) + high_sym_point_dict[k_symbol[i]][k]
        kpoint_array = np.vstack([kpath[key] for key in sorted(kpath.keys())])
        return kpoint_array

    # def monkhorst_pack(size,gamma_cnenter=False):
    #     kpts = np.indices(size).transpose((1,2,3,0)).reshape((-1,3))
    #     asize = np.array(size)
    #     shift = 


    def read_wanbandtb(self):
        full_path = self.file_path + '+wanbandtb'
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
                    kpoints[i][j] = float(values[j + 1])
            for k, line in enumerate(energy_line):
                values = line.split()
                kdist[k] = values[0]
                for n in range(bands_num):
                    bands[k][n] = float(values[n + 1])
            return kpoints, bands, kdist

        except FileNotFoundError:
            print(f"File not found: {full_path}")
            return None, None
        except Exception as e:
            print(f"An error occurred: {e}")
            return None, None

    def read_hamdata(self):
        full_path = self.file_path + '+hamdata'
        with open(full_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        nwan_value = None
        lattice_vectors = []
        tij_hij_blocks = []
        wancenters = []
        centering = []
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

            if 'centering:' in line:
                for j in range(1, 4):
                    if i + j < len(lines):
                        center = list(map(float, lines[i + j].strip().split()))
                        centering.append(center)
            
            if 'wancenters:' in line:
                for j in range(1, int(nwan_value) + 1):
                    if i + j < len(lines):
                        center = list(map(float, lines[i + j].strip().split()))
                        wancenters.append(center)

            if 'Tij, Hij, Sij:' in line:
                block = []
                try:
                    index1, index2 = map(int, lines[i + 1].strip().split())

                    for j in range(2, len(lines) - i):
                        if 'end Tij, Hij:' in lines[i + j]:
                            break
                        try:
                            block.append(list(map(float, lines[i + j].strip().split())))
                        except ValueError as e:
                            continue
                    if block:
                        tij_hij_blocks.append((index1, index2, block))
                        block_lengths.append((index1, index2, len(block)))
                        for irs in range(len(block)):
                            p_wfpairs_list.append((block[irs][:3], index1, index2, block[irs][3], block[irs][4], irs))
                except ValueError as e:
                    continue

        lattice_vectors = np.array(lattice_vectors)
        
        centering =np.array(centering)
        conventional_matrix = np.dot(np.linalg.inv(centering),lattice_vectors)
        lattice_moduli = np.linalg.norm(conventional_matrix, axis=1)
        p_wfpairs = np.array(p_wfpairs_list, dtype=object)
        return nwan_value, lattice_vectors, centering, lattice_moduli, tij_hij_blocks, wancenters, block_lengths, p_wfpairs

    def process_tij_hij_blocks(self):
        unique_r_vectors = []
        inv_lattice_vectors = np.linalg.inv(self.lattice_vectors)
        r_vector_count = {}

        for (index1, index2, block) in self.tij_hij_blocks:
            for row in block:
                r = np.array(row[:3])
                si = np.array(self.wancenters[index1 - 1])
                sj = np.array(self.wancenters[index2 - 1])
                R = r - (sj - si)
                T = np.dot(R, inv_lattice_vectors)
                T_rounded = np.rint(T)
                T_int = T_rounded.astype(int)
                unique_r_vectors.append((T_int[0], T_int[1], T_int[2], index1, index2))
        unique_r_vectors = np.array(unique_r_vectors)

        data_dict = {}

        for r_vector in unique_r_vectors:
            r1, r2, r3 = r_vector[:3]
            for i in range(1, self.nwan + 1):
                for j in range(1, self.nwan + 1):
                    data_dict[(r1, r2, r3, i, j)] = (0.0, 0.0)

        for (index1, index2, block) in self.tij_hij_blocks:
            for row in block:
                r = np.array(row[:3])
                si = np.array(self.wancenters[index1 - 1])
                sj = np.array(self.wancenters[index2 - 1])
                R = r - (sj - si)
                T = np.dot(R, inv_lattice_vectors)
                T_rounded = np.rint(T)
                T_int = T_rounded.astype(int)
                real = row[3]
                img = row[4]
                key = (T_int[0], T_int[1], T_int[2], index1, index2)
                if key in data_dict:
                    existing_real, existing_img = data_dict[key]
                    data_dict[key] = (existing_real + real, existing_img + img)
                else:
                    data_dict[key] = (real, img)
        processed_data = [(*k, *v) for k, v in data_dict.items()]
        r_vector_count_list = sorted(r_vector_count.items())
        unique_r_vectors_set = set(map(tuple, unique_r_vectors[:, :3]))
        nrpts = len(unique_r_vectors_set)
        return nrpts, processed_data, unique_r_vectors, r_vector_count_list

    def save_to_wannier90_centres_xyz(self):
        full_path = self.file_path + 'wannier90_centres.xyz'
        header = "Generated by Hao Wang's script\n"
        with open(full_path, 'w', encoding='utf-8') as f:
            f.write(f"{self.nwan}\n")
            f.write(header)
            for line in self.wancenters:
                f.write(
                    f"X          {line[0]:18.8f} {line[1]:18.8f} {line[2]:18.8f}\n"
                )
            f.write(f"Bi         2.367314417586326     -4.100308848749840     29.667930445390219 \n")
            f.write(f"Bi        -2.367314417586326      4.100308848749840    -29.667930445390219 \n")

    def save_to_wannier_hr(self):
        full_path = self.file_path + 'wannier90_hr.dat'
        header = "Generated by Hao Wang's script\n"
        with open(full_path, 'w', encoding='utf-8') as f:
            f.write(header)
            f.write(f"{self.nwan}\n")
            f.write(f"{self.nrpts}\n")

            ones_line = '    1' * 15
            for _ in range(self.nrpts // 15):
                f.write(f"{ones_line}\n")
            if self.nrpts % 15 > 0:
                f.write(f"{'    1' * (self.nrpts % 15)}\n")
            for line in self.processed_data:
                f.write(
                    f"{line[0]:5.0f} {line[1]:5.0f} {line[2]:5.0f} {int(line[3]):5d} {int(line[4]):5d} {line[5]:18.8f} {line[6]:18.8f}\n"
                )

    def read_wannier90_hr(self):
        full_path = self.file_path + 'wannier90_hr.dat'
        with open(full_path, 'r') as f:
            lines = f.readlines()
            num_wann = int(lines[1])
            nrpts = int(lines[2])
            lines_nrpts = int(np.ceil(nrpts / 15.0))
            ndegen_list = []
            for i in range(lines_nrpts):
                for j in range(len(lines[3 + i].split())):
                    ndegen_list.append(int(lines[3 + i].split()[j]))
            ndegen = np.array(ndegen_list)
            HmnR_np = np.zeros((num_wann ** 2 * nrpts, 7))
            for i in range(num_wann ** 2 * nrpts):
                for j in range(7):
                    HmnR_np[i][j] = (float(lines[3 + lines_nrpts + i].split()[j]))
            HmnR_np_iR = np.zeros((num_wann, num_wann, nrpts), dtype=complex)
            irvec = np.zeros((nrpts, 3))
            for ir in range(nrpts):
                for n in range(num_wann):
                    for m in range(num_wann):
                        HmnR_np_iR[m, n, ir] = complex(HmnR_np[ir * num_wann ** 2 + n * num_wann + m][5], HmnR_np[ir * num_wann ** 2 + n * num_wann + m][6])
                irvec[ir][0] = HmnR_np[ir * num_wann ** 2][0]
                irvec[ir][1] = HmnR_np[ir * num_wann ** 2][1]
                irvec[ir][2] = HmnR_np[ir * num_wann ** 2][2]

        return num_wann, ndegen, irvec, HmnR_np_iR

    def is_hermitian(self, matrix):
        return np.allclose(matrix, np.conj(matrix.T), atol=1e-4, )

    def assembly_ham(self,klist, irvec, ndegen, lat, lattice_moduli, HmnR_np_iR):
        for i in range(len(klist)):
            Ham_bulk = np.zeros((self.nwan, self.nwan), dtype=complex)
            for iR in range(self.nrpts):
                ia = irvec[iR][0]
                ib = irvec[iR][1]
                ic = irvec[iR][2]
                R = ia * lat[0, :] + ib * lat[1, :] + ic * lat[2, :]
                kdotR = np.dot(R, klist[i])
                factor = np.exp(1j * kdotR * 2 * np.pi / lattice_moduli[0])
                Ham_bulk[:, :] = Ham_bulk[:, :] + (HmnR_np_iR[:, :, iR] * factor / ndegen[iR])
        self.Ham_bulk = Ham_bulk[:, :]
        print('Ham_bulk is OK')
        return Ham_bulk

    def get_bands_from_ham(self, kpath, irvec, ndegen, lat, lattice_moduli, HmnR_np_iR):
        band_structure = np.zeros((len(kpath), self.nwan))
        for i in range(len(kpath)):
            Ham_bulk = np.zeros((self.nwan, self.nwan), dtype=complex)
            for iR in range(self.nrpts):
                ia = irvec[iR][0]
                ib = irvec[iR][1]
                ic = irvec[iR][2]
                R = ia * lat[0, :] + ib * lat[1, :] + ic * lat[2, :]
                kdotR = np.dot(R, kpath[i])
                factor = np.exp(1j * kdotR * 2 * np.pi / lattice_moduli[0])
                Ham_bulk[:, :] = Ham_bulk[:, :] + (HmnR_np_iR[:, :, iR] * factor / ndegen[iR])
            matrix_to_check = Ham_bulk[:, :]
            if self.is_hermitian(matrix_to_check):
                pass
            else:
                pass
            eigen_value, eigen_vector = np.linalg.eig(Ham_bulk)
            sorted_eig = np.sort(np.real(eigen_value))
            band_structure[i][:] = sorted_eig
        return band_structure

    def plot_bands(self, kdist, band_wanbandtb, band_structure):
        plt.figure(figsize=(10, 12))
        plt.ylim(-12, 10)
        for i in range(self.nwan):
            plt.plot(kdist, band_wanbandtb[:, i], '-', color='red', linewidth=4)
            plt.plot(kdist, band_structure[:, i], 'o', color='blue', markersize=2)
        plt.savefig(self.file_path + 'biband.png')
        plt.show()
        
        
        
    def pauli_block_all(self):
        if self.Ham_bulk is None:
            raise ValueError("Ham_bulk is not initialized. Call assembly_ham first.")
        self.MI = (self.Ham_bulk[::2, ::2] + self.Ham_bulk[1::2, 1::2]) / 2.0
        self.Mx = (self.Ham_bulk[::2, 1::2] + self.Ham_bulk[1::2, ::2]) / 2.0
        # Note that this is not element wise product with sigma_y but dot product
        self.My = (self.Ham_bulk[::2, 1::2] - self.Ham_bulk[1::2, ::2]) * 0.5j
        self.Mz = (self.Ham_bulk[::2, ::2] - self.Ham_bulk[1::2, 1::2]) / 2.0
        return self.MI, self.Mx, self.My, self.Mz
    
    
