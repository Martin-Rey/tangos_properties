from tangos.properties.pynbody import PynbodyHaloProperties
from tangos.properties.pynbody.centring import centred_calculation


class DensityAtRadius(PynbodyHaloProperties):
    """
    Modular class to calculate the density in a sphere at a given radius.
    Examples of simple class children are given below.
    """

    @staticmethod
    def _get_radius_def():
        """
        Implement your own method depending on how radius is defined.
        A string corresponding to a halo property can be supplied, e.g. "max_radius" or "nfw_rs" or else.
        A single number is assumed to correspond as a conversion factor between your radius unit and
        the unit in which radii are internally stored
        """
        raise NotImplementedError("This is meant to be an abstract class")

    @staticmethod
    def _number_radii():
        """
        Implement your own method to define the radius of the final sphere as
        _number_radii times _get_radius_def
        """
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_radius(cls, halo):
        if cls._get_radius_def is str:
            return cls._number_radii() * halo[cls._get_radius_def()]
        elif cls._get_radius_def is float:
            return cls._number_radii() * cls._get_radius_def()
        else:
            raise KeyError("Radius definition not supported")

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        import numpy as np
        import pynbody.analysis.cosmology

        mass = particle_data.d['mass'].sum().in_units("1 Msol")
        volume = 4 * np.pi / 3 * (self.get_radius(existing_properties)) ** 3
        rho_m = pynbody.analysis.cosmology.rho_M(particle_data)
        return mass/volume, mass/volume/rho_m

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(self.get_radius(existing_properties), existing_properties['shrink_center'])

    def requires_property(self):
        if self._get_radius_def is str:
            return ["shrink_center", self._get_radius_def()]
        else:
            return ["shrink_center"]


class DensityAt4Rmax(DensityAtRadius):
    names = "rho_4rmax", "overdensity_4rmax"

    @staticmethod
    def _get_radius_def():
        return "max_radius"

    @staticmethod
    def _number_radii():
        return 4.0


class DensityAt10Mpc(DensityAtRadius):
    names = "rho_10Mpc", "overdensity_10Mpc"

    @staticmethod
    def _get_radius_def():
        return 1000.0

    @staticmethod
    def _number_radii():
        return 10.0


class DensityAt1NFWRs(DensityAtRadius):
    names = "rho_1nfw_rs", "overdensity_1nfw_rs"

    @staticmethod
    def _get_radius_def():
        return "nfw_rs"

    @staticmethod
    def _number_radii():
        return 1.0
