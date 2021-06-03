try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:
    bdist_wheel = None


def universal_wheel(dist, attr, value):
    value = getattr(dist, 'universal_wheel', None)
    if value is None:
        dist.universal_wheel = True

    base_bdist_wheel = dist.cmdclass.get('bdist_wheel', bdist_wheel)

    if base_bdist_wheel is None:
        return

    class UniversalBdistWheel(base_bdist_wheel):
        def get_tag(self):
            rv = base_bdist_wheel.get_tag(self)
            if not dist.universal_wheel:
                return rv
            return ('py3', 'none',) + rv[2:]

    dist.cmdclass['bdist_wheel'] = UniversalBdistWheel
