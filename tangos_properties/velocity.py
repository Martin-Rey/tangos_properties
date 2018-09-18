from tangos.properties import LivePropertyCalculation


class VirialVelocity(LivePropertyCalculation):
    names = "Vvirial"

    def calculate(self, _, existing_properties):
        import numpy as np
        import pynbody.units
        import pynbody.array
        M = pynbody.array.SimArray(existing_properties['Mvirial'], "Msol")
        R = pynbody.array.SimArray(existing_properties['rvirial'], "kpc")

        return (np.sqrt(pynbody.units.G * M / R)).in_units("km s**-1")

    def requires_property(self):
        return ['rvirial', 'Mvirial']


class Velocity200c(LivePropertyCalculation):
    names = "V200c"

    def calculate(self, _, existing_properties):
        import numpy as np
        import pynbody.units
        import pynbody.array
        M = pynbody.array.SimArray(existing_properties['M200c'], "Msol")
        R = pynbody.array.SimArray(existing_properties['r200c'], "kpc")

        return (np.sqrt(pynbody.units.G * M / R)).in_units("km s**-1")

    def requires_property(self):
        return ['r200c', 'M200c']


# Super simple velocity profile DM only
# Would really ought to be using the pynbody rotation curves tool but somehow it requires 'eps' array

class VelocityProfile(LivePropertyCalculation):
    names = "dm_vcirc_profile"

    def calculate(self, _, existing_properties):
        import numpy as np
        import pynbody.units
        import pynbody.array

        M = pynbody.array.SimArray(existing_properties['dm_mass_profile'], "Msol")
        R = pynbody.array.SimArray(existing_properties['rbins_profile'], "kpc")
        v_profile = (np.sqrt(pynbody.units.G * M / R)).in_units("km s**-1")
        return v_profile

    def requires_property(self):
        return ['dm_mass_profile', 'rbins_profile']