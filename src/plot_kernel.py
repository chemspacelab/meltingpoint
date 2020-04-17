
import numpy as np
import misc
import views


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('--kernel', action='store', help='', metavar="file")
    parser.add_argument('--dist', action='store', help='', metavar="file")

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    properties = misc.load_npy(args.scratch + "properties")

    if args.kernel:
        filename = ".".join(args.kernel.split(".")[:-1])
        kernel = misc.load_npy(filename)

        views.pca_with_properties(kernel, properties, filename)

    if args.dist:
        filename = ".".join(args.dist.split(".")[:-1])
        kernel = misc.load_npy(filename)

        k_sigma = 5
        k_lambda = 10**-6
        kernel /= (2*k_sigma**2)
        kernel = np.exp(-kernel)
        diag_kernel = kernel[np.diag_indices_from(kernel)]
        # kernel[np.diag_indices_from(kernel)] = diag_kernel + k_lambda

        check = np.max(kernel)
        print(check)

        views.pca_with_properties(kernel, properties, filename)

    return

if __name__ == "__main__":
    main()
