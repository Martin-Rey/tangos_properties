from tangos.properties import LivePropertyCalculation
from tangos.live_calculation import NoResultsError
from tangos.properties import PropertyCalculation
from tangos import relation_finding


class Mergers(PropertyCalculation):
    requires_particle_data = False
    names = "mergers_z", "mergers_ratios"

    @classmethod
    def no_proxies(self):
        return True

    def calculate(self, _, existing_properties):
        return self._get_mergers_of_major_progenitor(existing_properties)

    @staticmethod
    def _get_mergers_of_major_progenitor(input_halo):
        import numpy as np
        redshift = []
        ratio = []
        while input_halo is not None:
            mergers = relation_finding.MultiHopMostRecentMergerStrategy(input_halo).all()
            if len(mergers) > 0:
                for m in mergers[1:]:
                    redshift.append(mergers[0].timestep.next.redshift)
                    ratio.append(float(mergers[0].NDM) / m.NDM)
                input_halo = mergers[0]
            else:
                input_halo = None

        return np.array(redshift), np.array(ratio)


class NumberMajorMergers(LivePropertyCalculation):
    requires_particle_data = False
    names = "number_major_mergers"

    def calculate(self, _, halo):
        import numpy as np
        ratios = halo['mergers_ratios']
        major_mergers = np.where(ratios <= self._get_major_merger_ratio_def())[0]
        return len(major_mergers)


    @staticmethod
    def _get_major_merger_ratio_def():
        return 4

    def requires_property(self):
        return ["mergers_ratios"]


class RedshiftLastMajorMergers(LivePropertyCalculation):
    requires_particle_data = False
    names = "z_last_major_mergers"

    def calculate(self, _, halo):
        import numpy as np
        ratios = halo['mergers_ratios']

        if len(ratios) == 0:
            return 0.0
        else:
            major_mergers = np.where(ratios <= self._get_major_merger_ratio_def())
            return len(major_mergers)


    @staticmethod
    def _get_major_merger_ratio_def():
        return 4

    def requires_property(self):
        return ["mergers_z", "mergers_ratios"]