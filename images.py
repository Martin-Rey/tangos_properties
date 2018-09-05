from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.centring import centred_calculation


class DarkMatterImages(PynbodyPropertyCalculation):

    @classmethod
    def plot_extent(cls):
        return cls._number_radii()

    @classmethod
    def plot_xlabel(cls):
        return "x/max_radius"

    @classmethod
    def plot_ylabel(cls):
        return "y/max_radius"

    def plot_clabel(self):
         return r"M$_{\odot}$ kpc$^{-2}$"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    @staticmethod
    def _number_radii():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def _plot_size(cls, existing_properties):
        import pynbody.array
        return pynbody.array.SimArray(cls._number_radii() * existing_properties['max_radius'], "kpc")

    def region_specification(self, db_data):
        import pynbody.filt
        side = self._plot_size(db_data)
        xcenter = db_data['shrink_center'][0]
        ycenter = db_data['shrink_center'][1]
        zcenter = db_data['shrink_center'][2]

        x1 = xcenter - side/2
        y1 = ycenter - side/2
        z1 = zcenter - side/2
        x2 = xcenter + side/2
        y2 = ycenter + side/2
        z2 = zcenter + side/2
        return pynbody.filt.Cuboid(x1, y1=y1, z1=z1, x2=x2, y2=y2, z2=z2)

    def _render_projected(self, f, size):
        import pynbody.plot
        im = pynbody.plot.sph.image(f, qty='rho', width=size, units="Msol kpc^-2", noplot=True)
        return im

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        im_z = self._render_projected(particle_data.dm, self._plot_size(existing_properties))
        with particle_data.rotate_y(90):
            im_x = self._render_projected(particle_data.dm, self._plot_size(existing_properties))
        return im_z, im_x


class DMCloseImages(DarkMatterImages):
    names = "dm_closeup_image_z", "dm_closeup_image_x"

    @staticmethod
    def _number_radii():
        return 4.0


class DMEnvironmentImages(DarkMatterImages):
    names = "dm_environment_z", "dm_environment_x"

    @staticmethod
    def _number_radii():
        return 15.0
