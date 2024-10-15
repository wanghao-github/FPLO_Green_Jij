import numpy as np
import matplotlib.pyplot as plt
import sys

class GreenFunction:
    
    def fermi(e, mu, width=0.01):
        x = (e - mu) / width
        return np.where(x < np.log(sys.float_info.max), 1 / (1.0 + np.exp(x)), 0.0)

    def get_density_matrix(self):
        rho = np.zeros((self.nbasis, self.nbasis), dtype=complex)
        if self.is_orthogonal:
            for ik, _ in enumerate(self.kpts):
                rho += (
                    (self.get_evecs(ik) * fermi(self.evals[ik], self.efermi))
                    @ self.get_evecs(ik).T.conj()
                    * self.kweights[ik]
                )
        else:
            for ik, _ in enumerate(self.kpts):
                rho += (
                    (self.get_evecs(ik) * fermi(self.evals[ik], self.efermi))
                    @ self.get_evecs(ik).T.conj()
                    @ self.get_Sk(ik)
                    * self.kweights[ik]
                )
        return rho
    