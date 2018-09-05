from tangos.properties import PropertyCalculation
from tangos.live_calculation import NoResultsError


class HalfValue(PropertyCalculation):
    requires_particle_data = False

    @staticmethod
    def _get_qty_name():
        raise NotImplementedError("This is meant to be an abstract class")

    @classmethod
    def get_qty(cls):
        return cls._get_qty_name()

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
        overall_target = 0.5 * qty_in_parent

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


class HalfVirialMass(HalfValue):
    names = "z_halfMvir", "a_halfMvir"

    @staticmethod
    def _get_qty_name():
        return "Mvirial"


class HalfMass(HalfValue):
    names = "z_halfMass", "a_halfMass"

    @staticmethod
    def _get_qty_name():
        return "halo_finder_mass"


class HalfM200c(HalfValue):
    names = "z_halfM200c", "a_halfM200c"

    @staticmethod
    def _get_qty_name():
        return "M200c"
