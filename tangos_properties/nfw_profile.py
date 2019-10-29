from tangos.properties.pynbody.spherical_region import SphericalRegionPropertyCalculation
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties import LivePropertyCalculation
import numpy as np


class HaloDensityProfile(SphericalRegionPropertyCalculation):

    names = "dm_density_profile", "dm_mass_profile", "rbins_profile", "number_counts_profile"

    def plot_x0(self):
        return self.plot_xdelta()/2

    def plot_xdelta(self):
        return self._simulation.get("approx_resolution_kpc", 0.1)

    def plot_xlabel(self):
        return "r/kpc"

    def plot_ylabel(self):
        return r"$\rho/M_{\odot}\,kpc^{-3}$", r"$M/M_{\odot}$"

    def _get_profile(self, halo, maxrad):
        import pynbody

        delta = self.plot_xdelta()
        nbins = int(maxrad / delta)
        maxrad = delta * (nbins + 1)

        pro = pynbody.analysis.profile.Profile(halo, type='log', ndim=3,
                                               min=self.plot_x0(), max=maxrad, nbins=nbins)

        rho_a = pro['density']
        mass_a = pro['mass_enc']
        rbins_a = pro['rbins']
        n_a = pro['n']

        rho_a = rho_a.view(np.ndarray)
        mass_a = mass_a.view(np.ndarray)
        rbins_a = rbins_a.view(np.ndarray)
        n_a = n_a.view(np.ndarray)

        return rho_a, mass_a, rbins_a, n_a

    @centred_calculation
    def calculate(self, data, existing_properties):

        dm_a, dm_b, dm_r, dm_n = self._get_profile(data.dm, existing_properties["max_radius"])

        return dm_a, dm_b, dm_r, dm_n


class NFWParameters(HaloDensityProfile):
    """ In this class, concentration/ overdensity should be live properties.
        Theoeretical profile should be called from pynbody.
        Knowledge of dm_density_profile should be avoided at all cost.
        Fit could be restrained in inner region and outer regions for physical reasons.
    """

    names = "nfw_rhos", "nfw_rs", "nfw_goodness_fit", "nfw_reduced_goodness_fit", "nfw_cov", "nfw_rms"

    def requires_property(self):
        return ['dm_density_profile', "rbins_profile", "number_counts_profile", "max_radius", "shrink_center"]

    @centred_calculation
    def calculate(self, halo, existing_properties):
        from pynbody.analysis.theoretical_profiles import NFWprofile

        halo_rmasked, halo_profile_masked, err = self._prepare_profile(halo, existing_properties)

        rhos_s_guess = np.mean(halo_profile_masked)
        rs_guess = np.median(halo_rmasked)

        # NFW fitting with leastsquare
        parameters, cov = NFWprofile.fit(halo_rmasked, halo_profile_masked, profile_err=err,
                                         use_analytical_jac=True, guess=[rhos_s_guess, rs_guess])

        rhos_s = parameters[0]
        rs = parameters[1]
        theory = NFWprofile.profile_functional_static(halo_rmasked, rhos_s, rs)

        chi_square = (((halo_profile_masked - theory) / err) ** 2).sum()
        reduced_chi_square = chi_square / (len(halo_rmasked) - len(parameters))

        rms = (1 / len(halo_rmasked)) * (((np.log10(halo_profile_masked) - np.log10(theory)) ** 2).sum())

        return rhos_s, rs, chi_square, reduced_chi_square, cov, rms

    def _prepare_profile(self, halo, existing_properties):
        """ Makes sure that arrays are non zeros and masked """

        rbins = existing_properties['rbins_profile']
        rho = existing_properties['dm_density_profile']
        bin_counts = existing_properties['number_counts_profile']
        err = rho / np.sqrt(bin_counts)

        r_convergence = self._get_power_radius(halo, rbins, rho, bin_counts)

        r_relaxation = self._get_r_relaxation(existing_properties)

        r, profile, err = self._extract_non_zero(rbins, rho, err)

        return self._physical_mask(r, profile, err, r_convergence, r_relaxation)

    @staticmethod
    def _get_power_radius(halo,  r, profile, bin_counts):
        """ Determine the convergence radius from Eq 20 with 0.6 t_0 in  Power et al 2003."""
        import pynbody.analysis.cosmology

        number = 0
        i = 0

        while number < 0.6 or np.isnan(number):
            rconv = r[i]
            enclosed_number_particles = bin_counts[:i].sum()
            density_contrast = np.mean(profile[:i]) / pynbody.analysis.cosmology.rho_crit(halo, unit="Msol kpc**-3")
            number = np.sqrt(200) / 8 * enclosed_number_particles / np.log(enclosed_number_particles) / \
                     np.sqrt(density_contrast)
            i += 1

        return rconv

    @staticmethod
    def _get_r_relaxation(existing_properties):
        return 0.6 * existing_properties['max_radius']

    @staticmethod
    def _extract_non_zero(r, profile, err):
        halo_r_nonzero = r[np.where(profile > 0)]
        halo_profile_nonzero = profile[np.where(profile > 0)]
        err_nonzero = err[np.where(profile > 0)]
        return halo_r_nonzero, halo_profile_nonzero, err_nonzero

    @staticmethod
    def _physical_mask(r, profile, err, rmin, rmax):

        r0_mask = np.amin(r)
        rmax_mask = np.amax(r)

        if rmin > r0_mask:
            r0_mask = rmin

        if rmax < rmax_mask:
            rmax_mask = rmax

        rmasked = r[np.where(r > r0_mask)[0][0]:np.where(r < rmax_mask)[0][-1]]
        profile_masked = profile[np.where(r > r0_mask)[0][0]:np.where(r < rmax_mask)[0][-1]]
        err_masked = err[np.where(r > r0_mask)[0][0]:np.where(r < rmax_mask)[0][-1]]

        return rmasked, profile_masked, err_masked


class NFWConcentration(LivePropertyCalculation):
    names = "nfw_concentration"

    def calculate(self, _, existing_properties):
        return existing_properties['rvirial'] / existing_properties['nfw_rs']

    def requires_property(self):
        return ['rvirial', 'nfw_rs']


class NFWConcentrationPlus1Sigma(LivePropertyCalculation):
    names = "nfw_concentration_plus1s"

    def calculate(self, _, existing_properties):
        sigma = np.sqrt(existing_properties['nfw_cov'][1][1])
        return existing_properties['rvirial'] / (existing_properties['nfw_rs'] - sigma)

    def requires_property(self):
        return ['rvirial', 'nfw_rs']


class NFWConcentrationMinusSigma(LivePropertyCalculation):
    names = "nfw_concentration_minus1s"

    def calculate(self, _, existing_properties):
        sigma = np.sqrt(existing_properties['nfw_cov'][1][1])
        return existing_properties['rvirial'] / (existing_properties['nfw_rs'] + sigma)

    def requires_property(self):
        return ['rvirial', 'nfw_rs', 'nfw_cov']


class NFWVelocityConcentration(LivePropertyCalculation):
    names = "nfw_concentration_vel"

    def no_proxies(self):
        return True

    def calculate(self, _, existing_properties):
        import numpy as np
        import pynbody.analysis.theoretical_profiles
        import scipy.optimize as so

        vmax = np.amax(existing_properties['dm_vcirc_profile'])
        vvirial = existing_properties.calculate("V200c()")
        ratio = (vmax / vvirial) ** 2

        if ratio < 1:
            # Pathological halo that can't be described by NFW
            return None

        xmax = 2.163    # Klypin 2016 for NFW

        def f(c):
            return c / pynbody.analysis.theoretical_profiles.NFWprofile._helper_function(c) \
                   / xmax * pynbody.analysis.theoretical_profiles.NFWprofile._helper_function(xmax)

        def root(c):
            return f(c) - ratio

        return so.newton(root, 10)

    def requires_property(self):
        return ['r200c', 'dm_vcirc_profile', "M200c"]


class NFWHalfMassConcentration(LivePropertyCalculation):
    names = "nfw_concentration_halfmass"

    def no_proxies(self):
        return True

    def calculate(self, _, existing_properties):
        import numpy as np
        import scipy.optimize as so

        r200c = existing_properties['r200c']
        r_half_mass = existing_properties['half_mass_radius']
        ratio = r_half_mass / r200c

        def f(c):   # Lokas and Mamon 2000, eq 28
            return 0.6082 - 0.1843 * np.log10(c) - 0.1011 * (np.log10(c) ** 2) + 0.03918 * (np.log10(c) ** 3)

        def root(c):
            return f(c) - ratio

        return so.newton(root, 1)

    def requires_property(self):
        return ['r200c', 'half_mass_radius', "M200c"]