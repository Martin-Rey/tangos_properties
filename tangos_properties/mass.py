from tangos.properties.pynbody import PynbodyHaloProperties
from tangos.properties.pynbody.centring import centred_calculation


class HaloFinderMass(PynbodyHaloProperties):
    names = "halo_finder_mass"

    def calculate(self, halo, existing_properties):
        return halo['mass'].sum().in_units('1 Msol')


class Mass(PynbodyHaloProperties):

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        return particle_data['mass'].sum().in_units("1 Msol")


class Mass200m(Mass):
    names = "M200m"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r200m'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r200m"]


class Mass200c(Mass):
    names = "M200c"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r200c'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r200c"]


class VirialMass(Mass):
    names = "Mvirial"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['rvirial'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["rvirial"]


class Mass500m(Mass):
    names = "M500m"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r500m'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r500m"]


class Mass500c(Mass):
    names = "M500c"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r500c'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r500c"]
