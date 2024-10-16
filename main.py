# main.py

from read_hamdata import HamiltonianProcessor

# file_path = r'C:\Users\wangh\OneDrive\Desktop\Codes\fplohamdatatest\Fe\fplo22\Fe_4sp_3d_cutoff_25\\'
file_path = r'C:\Users\wangh\OneDrive\Desktop\Works\Orbital Hall\fcc Pt\FPLO\\'
processor = HamiltonianProcessor(file_path)

kpath, band_wanbandtb, kdist = processor.read_wanbandtb()

processor.save_to_wannier_hr()

num_wann, ndegen, irvec, HmnR_np_iR = processor.read_wannier90_hr()

band_structure = processor.get_bands_from_ham(kpath, irvec, ndegen, processor.lattice_vectors, processor.lattice_moduli, HmnR_np_iR)

processor.plot_bands(kdist, band_wanbandtb, band_structure)
# processor.save_to_wannier90_centres_xyz()
# processor.save_to_wannier_hr()

# processor.assembly_ham(kpath, irvec, ndegen, processor.lattice_vectors, processor.lattice_moduli, HmnR_np_iR)
# processor.pauli_block_all()