import math

import numpy as np

class GaussianTypeOrbital(object):
    def __init__(self, center, alpha, normalization = 0):
        # add normalization
        self.center = center
        self.alpha = alpha

        if normalization == 0:
            self.normalization = (2 * self.alpha / np.pi) ** 0.75
        else: self.normalization = normalization

    def __mul__(self, other):
        if type(other) is GaussianTypeOrbital:
            return product_of_gaussian_type_orbitals(self, other)
        elif type(other) is SumOfGaussianTypeOrbitals:
            return product_of_sum_gaussian_type_orbitals(SumOfGaussianTypeOrbitals([self]), other)
        else:
            return GaussianTypeOrbital(self.center, self.alpha, self.normalization * other)

    def __add__(self, other):
        if type(other) is GaussianTypeOrbital:
            return SumOfGaussianTypeOrbitals([self, other])


    def overlap_integral(self, orbital):
        product = self * orbital
        return product.normalization * ((np.pi / product.alpha) ** 1.5)

    def ke_integral(self, orbital):
        overlap = self.overlap_integral(orbital)
        intermediate = self.alpha * orbital.alpha / (self.alpha + orbital.alpha)
        return intermediate * (3 - 2 * intermediate * np.linalg.norm(self.center - orbital.center) ** 2) * overlap

    def nuc_integral(self, orbital, nucleus_center):
        def F0(x):
            if x == 0: return 1
            else: return (0.5 * (np.pi / x) ** 0.5) * math.erf(x ** 0.5)
        product_center = ((self.center * self.alpha) + (orbital.center * orbital.alpha)) / (self.alpha + orbital.alpha)
        intermediate = self.alpha + orbital.alpha
        overlap = self.overlap_integral(orbital)
        return 2 * (intermediate / np.pi) ** 0.5 * F0(intermediate * np.linalg.norm(product_center - nucleus_center) ** 2) * overlap



class SumOfGaussianTypeOrbitals(object):
    def __init__(self, gaussians):
        self.gaussians = gaussians

    def center(self, center):
        for gaussian in self.gaussians:
            gaussian.center = center

    # maybe wrap these so its not repeated?
    def overlap_integral(self, orbital):
        overlap = np.float64(0)
        for gaussian1 in self.gaussians:
            for gaussian2 in orbital.gaussians:
                overlap += gaussian1.overlap_integral(gaussian2)

        return overlap

    def ke_integral(self, orbital):
        ke = np.float64(0)
        for gaussian1 in self.gaussians:
            for gaussian2 in orbital.gaussians:
                ke += gaussian1.ke_integral(gaussian2)

        return ke

    def nuc_integral(self, orbital, nucleus_center):
        result = np.float64(0)
        for gaussian1 in self.gaussians:
            for gaussian2 in orbital.gaussians:
                result += gaussian1.nuc_integral(gaussian2, nucleus_center)

        return result

# GaussianTypeOrbital * GaussianTypeOrbital -> GaussianTypeOrbital
def product_of_gaussian_type_orbitals(orbital1, orbital2):
    product_alpha = orbital1.alpha + orbital2.alpha
    product_center = ((orbital1.center * orbital1.alpha) + (orbital2.center * orbital2.alpha)) / product_alpha
    product_normalization = orbital1.normalization * orbital2.normalization * np.exp(
        -orbital1.alpha * orbital2.alpha / product_alpha * (np.linalg.norm(orbital1.center - orbital2.center) ** 2))
    return GaussianTypeOrbital(product_center, product_alpha, product_normalization)

# SumOfGaussianTypeOrbitals + SumOfGaussianTypeOrbitals -> SumOfGaussianTypeOrbitals
def sum_of_gaussian_type_orbitals(orbital1, orbital2):
    return SumOfGaussianTypeOrbitals(orbital1.gaussians + orbital2.gaussians)

# SumOfGaussianTypeOrbitals * SumOfGaussianTypeOrbitals -> SumOfGaussianTypeOrbitals
def product_of_sum_gaussian_type_orbitals(orbital1, orbital2):
    sum_gaussians = []
    for gaussian1 in orbital1.gaussians:
        for gaussian2 in orbital2.gaussians:
            sum_gaussians += gaussian1 * gaussian2







sum = SumOfGaussianTypeOrbitals([
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 0.168856)*0.444635,
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 0.623913)*0.535328,
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 3.42525)*0.154329])
sum2 = SumOfGaussianTypeOrbitals([
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 0.168856)*0.444635,
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 0.623913)*0.535328,
    GaussianTypeOrbital(center = np.array((0, 0, 0)), alpha = 3.42525)*0.154329])

sum2.center(np.array([0,0,1.4]))

print (sum.overlap_integral(sum2))
print (sum.ke_integral(sum))
print (sum2.nuc_integral(sum2, np.array([0,0,1.4])))