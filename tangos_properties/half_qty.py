from tangos.properties import PropertyCalculation
from tangos.live_calculation import NoResultsError


""" This property calculates the time/redshift at which a quantity is equal to a percentage of its final value.
    For example, the half mass redshift could be calculated by overriding this class with 
    _get_qty_name() as finder_mass and _get_target() as 0.5. 
    This property could be refactored as a live calculation for more modular use, similar to at(). """


class TimeAtWhichQuantityIsXXPerCent(PropertyCalculation):
    requires_particle_data = False

    @staticmethod
    def _get_qty_name():
        raise NotImplementedError("This is meant to be an abstract class")

    @staticmethod
    def _get_target_value():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_qty(cls):
        return cls._get_qty_name()

    @classmethod
    def get_target(cls):
        return cls._get_target_value()

    @classmethod
    def no_proxies(self):
        return True

    @staticmethod
    def _find_linear_intersect(y1, y0, x1, x0, target):
        slope = (y1 - y0) / (x1 - x0)
        return x0 + (target - y0) / slope

    def calculate(self, _, existing_properties):
        # TODO All these catch blocks are ugly but this is the only to make sure that the property will spit out intelligble errors
        parent = existing_properties
        qty_in_parent = parent[self.get_qty()]
        overall_target = self.get_target() * qty_in_parent

        try:
            child = parent.calculate("earlier(1)")
        except NoResultsError:
            raise NoResultsError("Could not find a descendant to %r" % parent)

        try:
            qty_in_child = child[self.get_qty()]
        except KeyError:
            raise KeyError(self.get_qty() + " does not exist in %r " % child)

        while qty_in_child > overall_target:
            parent = child
            qty_in_parent = qty_in_child

            try:
                child = parent.calculate("earlier(1)")
            except NoResultsError:
                raise NoResultsError("Could not find a descendant to %r" % child)

            try:
                qty_in_child = child[self.get_qty()]
            except KeyError:
                raise KeyError(self.get_qty() + " does not exist in %r " % child)

        return self._find_linear_intersect(qty_in_parent, qty_in_child, parent.calculate("z()"),
                                           child.calculate("z()"), overall_target), \
               self._find_linear_intersect(qty_in_parent, qty_in_child, parent.calculate("a()"),
                                           child.calculate("a()"), overall_target)

    def requires_property(self):
        return [self._get_qty_name()]


class HalfVirialMass(TimeAtWhichQuantityIsXXPerCent):
    names = "z_halfMvir", "a_halfMvir"

    @staticmethod
    def _get_qty_name():
        return "Mvirial"

    @staticmethod
    def _get_target_value():
        return 0.5


class HalfMass(TimeAtWhichQuantityIsXXPerCent):
    names = "z_halfMass", "a_halfMass"

    @staticmethod
    def _get_qty_name():
        return "halo_finder_mass"

    @staticmethod
    def _get_target_value():
        return 0.5


class HalfM200c(TimeAtWhichQuantityIsXXPerCent):
    names = "z_halfM200c", "a_halfM200c"

    @staticmethod
    def _get_qty_name():
        return "M200c"

    @staticmethod
    def _get_target_value():
        return 0.5


class FourPerCentM200c(TimeAtWhichQuantityIsXXPerCent):
    names = "z_fourpercentM200c", "a_fourpercentM200c"

    @staticmethod
    def _get_qty_name():
        return "M200c"

    @staticmethod
    def _get_target_value():
        return 0.04


class TenPerCentM200c(TimeAtWhichQuantityIsXXPerCent):
    names = "z_tenpercentM200c", "a_tenpercentM200c"

    @staticmethod
    def _get_qty_name():
        return "M200c"

    @staticmethod
    def _get_target_value():
        return 0.1


class TwentyPerCentM200c(TimeAtWhichQuantityIsXXPerCent):
    names = "z_twentypercentM200c", "a_twentypercentM200c"

    @staticmethod
    def _get_qty_name():
        return "M200c"

    @staticmethod
    def _get_target_value():
        return 0.2
