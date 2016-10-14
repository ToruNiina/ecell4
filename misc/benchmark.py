import time

from ecell4 import *

radius = 0.005
D = 1

def singlerun(f, L, num, max_steps, min_duration=0.0):
    """
    Parameters
    ----------
    f : Factory
    L : Real
        A size of the World
    num : Real
        The number of molecules
    max_steps : Integer
        The maximum number of steps
    min_duration : Real, optional
        The minimum duration
    """
    m = NetworkModel()
    m.add_species_attribute(Species("A", str(radius), str(D)))
    w = f.create_world(ones() * L)
    w.bind_to(m)
    w.add_molecules(Species("A"), num)
    sim = f.create_simulator(w)
    sim.initialize()
    tstart = time.time()
    i, telapsed = 0, 0.0
    while i < max_steps or telapsed < min_duration:
        sim.step()
        telapsed = time.time() - tstart
        i += 1
    return telapsed / sim.t()

def run(num_trials, *args):
    retval = []
    for _ in range(num_trials):
        retval.append(singlerun(*args))
    return retval

def matrix_sizes(L, N, r):
    N = int(min(L / (2 * r), max(3, cbrt(N))))
    return Integer3(N, N, N)

def partitioned_factory_maker(ftype, *args):
    def create_factory(L, num):
        rng = GSLRandomNumberGenerator()
        rng.seed(0)
        return ftype(matrix_sizes(L, num, radius), rng, *args)
    return create_factory

def non_partitioned_factory_maker(ftype, *args):
    def create_factory(L, num):
        rng = GSLRandomNumberGenerator()
        rng.seed(0)
        return ftype(*(args + (rng, )))
    return create_factory

def savedata(filename, num, data):
    with open(filename, "a") as fout:
        line = "{}\t{}".format(num, "\t".join([str(t) for t in data]))
        fout.write(line)
        fout.write("\n")
        print("{} => {}".format(filename, line))

def plotdata(ax, filename, label=None, c="k", marker="s", line=(0, 1.0)):
    import numpy
    data = numpy.loadtxt(filename)
    # data = numpy.log10(data)
    data = numpy.array([(row[0], numpy.mean(row[1: ]), numpy.std(row[1: ])) for row in data]).T

    ax.errorbar(data[0], data[1], data[2], fmt='o', color=c, marker=marker, mew=0, label=label)
    # ax.plot(data[0], data[0] + data[1][0] - data[0][0], '--', color=c)
    if line is not None:
        x = numpy.logspace(0.5, 6.5, 5)
        # x = data[0]
        ax.plot(x, numpy.power(x, line[1]) * (data[1][line[0]] / numpy.power(data[0][line[0]], line[1])), '--', color=c)


if __name__ == "__main__":
    import numpy
    import os
    import os.path

    def profile(filename, ns, create_factory, fixed_volume, one_particle_per_step, max_steps, min_duration):
        if os.path.isfile(filename):
            os.remove(filename)

        numarray = numpy.logspace(*ns).astype(int)
        for num in numarray:
            if fixed_volume:
                L = cbrt(40.0)  # 3.42
            else:
                L = cbrt(num / 60.0 * 1.0)  # 100nM

            if one_particle_per_step:
                max_steps_ = num * max_steps
            else:
                max_steps_ = max_steps

            savedata(filename, num, run(5, create_factory(L, num), L, num, max_steps_, min_duration))

    def profileall(solvers, ftypes=None):
        if ftypes is None:
            ftypes = solvers.keys()

        for ftype in ftypes:
            for fixed_volume in (True, False):
                if (ftype, fixed_volume) in (("Spatiocyte", False), ):
                    ns = (1.0, 5.0, 9)
                else:
                    ns = (1.0, 6.0, 11)
                create_factory, one_particle_per_step, c, marker = solvers[ftype]
                filename = "{}-{}.tsv".format(ftype.replace(".", "-"), "volume" if fixed_volume else "conc")
                profile(filename, ns, create_factory, fixed_volume, one_particle_per_step, max_steps, min_duration)

    def plotall(outputfilename, solvers, ftypes=None):
        if ftypes is None:
            ftypes = sorted(tuple(solvers.keys()))

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        plt.subplots_adjust(left = 0.1, right = 0.6)
        ax.set_xscale("log")
        ax.set_yscale("log")
        for ftype in ftypes:
            for fixed_volume in (True, False):
                create_factory, one_particle_per_step, c, marker = solvers[ftype]
                filename = "{}-{}.tsv".format(ftype.replace(".", "-"), "volume" if fixed_volume else "conc")
                label = "{} ({})".format(ftype, "volume" if fixed_volume else "conc")
                if (ftype, fixed_volume) == ("eGFRD", True):
                    plotdata(ax, filename, label, c, "^" if fixed_volume else "v", (0, 5.0 / 3.0))
                elif ftype == "Spatiocyte":
                    plotdata(ax, filename, label, c, "^" if fixed_volume else "v", (5, 1.0))
                else:
                    if fixed_volume:
                        plotdata(ax, filename, label, c, "^", None)
                        # plotdata(ax, filename, label, c, "^")
                    else:
                        plotdata(ax, filename, label, c, "v")
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] for h in handles]  # remove the errorbars
        ax.set_xlim(1.0, 6.0)
        ax.set_xlim(10.0 ** 0.5, 10.0 ** 6.5)
        ax.set_xlabel("N [# particles]")
        ax.set_ylabel("time [sec]")
        ax.grid()
        ax.legend(handles, labels, loc='upper left', numpoints=1, shadow=True, bbox_to_anchor=(1.0, 1.0))
        plt.savefig(outputfilename)
        plt.show()

    max_steps = 10
    min_duration = 1.0

    solvers = {
        "Mesoscopic": (non_partitioned_factory_maker(meso.MesoscopicFactory, 0.1), True, "b", "o"),
        "Mesoscopic relaxed": (non_partitioned_factory_maker(meso.MesoscopicFactory, 0.3), True, "navy", "o"),
        "BD": (partitioned_factory_maker(bd.BDFactory, 1e-5), False, "k", "^"),
        "BD relaxed": (partitioned_factory_maker(bd.BDFactory, 1e-3), False, "gray", "^"),
        "BD eGFRD": (partitioned_factory_maker(egfrd.BDFactory, 1e-5), False, "silver", "v"),
        "eGFRD": (partitioned_factory_maker(egfrd.EGFRDFactory), True, "r", "v"),
        "Spatiocyte": (non_partitioned_factory_maker(spatiocyte.SpatiocyteFactory, radius), False, "g", "o"),
        }

    # profileall(solvers)
    plotall("benchmark.png", solvers)
