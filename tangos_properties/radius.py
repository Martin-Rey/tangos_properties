from tangos.properties.pynbody import PynbodyHaloProperties
from tangos.properties.pynbody.centring import centred_calculation


class Radius(PynbodyHaloProperties):

    @staticmethod
    def _get_overdensity_contrast():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_contrast(cls):
        return cls._get_overdensity_contrast()

    @staticmethod
    def _get_reference_definition():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_rhodef(cls):
        return cls._get_reference_definition()

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        import pynbody.analysis as analysis
        return analysis.halo.virial_radius(particle_data, overden=self.get_contrast(), rho_def=self.get_rhodef())

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['max_radius'] * 3,
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["shrink_center", "max_radius"]


class Radius200m(Radius):
    names = "r200m"

    @staticmethod
    def _get_overdensity_contrast():
        return 200

    @staticmethod
    def _get_reference_definition():
        return 'matter'


class Radius200c(Radius):
    names = "r200c"

    @staticmethod
    def _get_overdensity_contrast():
        return 200

    @staticmethod
    def _get_reference_definition():
        return 'critical'


class RadiusVirial(Radius):
    names = "rvirial"

    @staticmethod
    def _get_overdensity_contrast():
        return 178

    @staticmethod
    def _get_reference_definition():
        return 'matter'


class Radius500m(Radius):
    names = "r500m"

    @staticmethod
    def _get_overdensity_contrast():
        return 500

    @staticmethod
    def _get_reference_definition():
        return 'matter'


class Radius500c(Radius):
    names = "r500c"

    @staticmethod
    def _get_overdensity_contrast():
        return 500

    @staticmethod
    def _get_reference_definition():
        return 'critical'