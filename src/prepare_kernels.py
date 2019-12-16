

def dump_distances_and_kernels(scr):

    # TODO Properties should be read by scr!!

    # properties
    print("Saving properties")
    with open(scr + 'properties.csv', 'r') as f:
        properties = f.readlines()
        properties = [x.split()[0] for x in properties]
        properties = [float(x) for x in properties]
        properties = np.array(properties)

    print(properties.shape)

    misc.save_npy(scr + "properties", properties)

    # Prepare distances
    representation_names = ["cm", "bob", "slatm"] # + ["avgslatm"]
    for name in representation_names:
        print("Distance", name)
        representations = misc.load_npy(scr + "repr." + name)
        print(representations.shape)
        dist = generate_l2_distances(representations)
        misc.save_npy(scr + "dist." + name, dist)

        dist = None
        del dist

    # Prepare fchl kernels
    if False:
        print("Generating fchl18 kernel")
        start = time.time()
        reps = misc.load_npy(scr + "repr." + "fchl18")
        print("shape:", reps.shape)
        sigmas, kernels = get_fchl18_kernels(reps, return_sigmas=True)
        end = time.time()
        print("time:", end-start)
        misc.save_npy(scr + "fchl18." + "sigmas", sigmas)
        misc.save_npy(scr + "kernels." + "fchl18", kernels)

        reps = None
        del reps
        kernels = None
        del kernels

    if False:
        print("Generating fchl19 kernel")
        reps = misc.load_npy(scr + "repr." + "fchl19")
        print("shape:", reps.shape)
        atoms = misc.load_obj(scr + "atoms")
        start = time.time()
        sigmas, kernels = get_fchl19_kernels(reps, atoms, return_sigmas=True)
        end = time.time()
        print("time:", end-start)
        misc.save_npy(scr + "fchl19." + "sigmas", sigmas)
        misc.save_npy(scr + "kernels." + "fchl19", kernels)

    if True:
        print("Generating fingerprint kernel")
        representations_fp = misc.load_obj(scr + "repr.fp")
        kernel = get_fp_kernel(representations_fp)
        misc.save_npy(scr + "kernel.fp", kernel)

    return



def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0)
    parser.add_argument('--get-kernels', action='store_true', help='')

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    if args.get_kernels:
        dump_distances_and_kernels(args.scratch)

    return


