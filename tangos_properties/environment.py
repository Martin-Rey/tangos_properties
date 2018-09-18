from tangos.properties.pynbody import PynbodyHaloProperties
from tangos.properties.pynbody.centring import centred_calculation


class DensityAtRadius(PynbodyHaloProperties):

    @staticmethod
    def _number_radii():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_radius(cls):
        return cls._number_radii()

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        import numpy as np
        import pynbody.analysis.cosmology

        mass = particle_data.d['mass'].sum().in_units("1 Msol")
        volume = 4 * np.pi / 3 * (self.get_radius() * existing_properties['max_radius']) ** 3
        rho_m = pynbody.analysis.cosmology.rho_M(particle_data)
        return mass/volume, mass/volume/rho_m

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(self.get_radius() * existing_properties['max_radius'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["shrink_center", "max_radius"]


class DensityAt4Radii(DensityAtRadius):
    names = "rho_4rmax", "overdensity_4rmax"

    @staticmethod
    def _number_radii():
        return 4.0


class DensityAt8Radii(DensityAtRadius):
    names = "rho_8rmax", "overdensity_8rmax"

    @staticmethod
    def _number_radii():
        return 8.0


class DensityAt10Radii(DensityAtRadius):
    names = "rho_10rmax", "overdensity_10rmax"

    @staticmethod
    def _number_radii():
        return 10.0


class DensityAtScale(PynbodyHaloProperties):

    @staticmethod
    def _number_Mpc():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_scale(cls):
        return cls._number_Mpc() * 1000

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        import numpy as np
        import pynbody.analysis.cosmology

        mass = particle_data.d['mass'].sum().in_units("1 Msol")
        volume = 4 * np.pi / 3 * (self.get_scale()) ** 3
        rho_m = pynbody.analysis.cosmology.rho_M(particle_data)
        return mass/volume, mass/volume/rho_m

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(self.get_scale(),
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["shrink_center", "max_radius"]


class DensityAt15Mpc(DensityAtScale):
    names = "rho_15Mpc", "overdensity_15Mpc"

    @staticmethod
    def _number_Mpc():
        return 15.0


class DensityAt10Mpc(DensityAtScale):
    names = "rho_10Mpc", "overdensity_10Mpc"

    @staticmethod
    def _number_Mpc():
        return 10.0


class DensityAt4Mpc(DensityAtScale):
    names = "rho_4Mpc", "overdensity_4Mpc"

    @staticmethod
    def _number_Mpc():
        return 4.0


class DensityAt8Mpc(DensityAtScale):
    names = "rho_8Mpc", "overdensity_8Mpc"

    @staticmethod
    def _number_Mpc():
        return 8.0


class DensityAt1Mpc(DensityAtScale):
    names = "rho_1Mpc", "overdensity_1Mpc"

    @staticmethod
    def _number_Mpc():
        return 1.0


class DensityAt2Mpc(DensityAtScale):
    names = "rho_2Mpc", "overdensity_2Mpc"

    @staticmethod
    def _number_Mpc():
        return 2.0
