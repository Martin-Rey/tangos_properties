from tangos.properties.pynbody import PynbodyHaloProperties
from tangos.properties.pynbody.centring import centred_calculation


class Spin(PynbodyHaloProperties):

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        import pynbody.analysis.angmom as angmom

        with angmom.sideon(particle_data, return_transform=True,
                           cen_size=self._simulation.get("approx_resolution_kpc", 0.1)*10.,
                           disk_size=existing_properties['max_radius']):
            s = angmom.spin_parameter(particle_data.d)

        return s


class SpinB01Virial(Spin):
    names = "spinB01_virial"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['rvirial'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["rvirial", "max_radius", "shrink_center"]


class SpinB01R200m(Spin):
    names = "spinB01_r200m"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r200m'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r200m", "max_radius", "shrink_center"]


class SpinB01R200c(Spin):
    names = "spinB01_r200c"

    def region_specification(self, existing_properties):
        import pynbody
        return pynbody.filt.Sphere(existing_properties['r200c'],
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["r200c", "max_radius", "shrink_center"]
