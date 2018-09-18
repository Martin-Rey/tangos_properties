from tangos.properties import LivePropertyCalculation


class Mergers(LivePropertyCalculation):
    requires_particle_data = False
    names = "mergers_z, mergers_r, mergers_halos"

    def calculate(self, _, existing_properties):
        return self._get_mergers_of_major_progenitor(existing_properties)

    @staticmethod
    def _get_mergers_of_major_progenitor(input_halo):
        """Given a halo, return the redshifts and ratios of all mergers on the major progenitor branch

        :parameter input_halo - the halo to consider the history of
        :type input_halo tangos.core.Halo

        :returns redshift, ratio, halo - arrays of the redshifts, ratios (1:X) and halo DB objects for the mergers
        """
        import numpy as np
        redshift = []
        ratio = []
        halo = []
        while input_halo is not None:
            mergers = tangos.relation_finding.MultiHopMostRecentMergerStrategy(input_halo).all()
            if len(mergers) > 0:
                for m in mergers[1:]:
                    redshift.append(mergers[0].timestep.next.redshift)
                    halo.append((mergers[0], m))
                    ratio.append(float(mergers[0].NDM) / m.NDM)
                input_halo = mergers[0]
            else:
                input_halo = None

        return np.array(redshift), np.array(ratio), halo