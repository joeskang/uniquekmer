"""
    main driver program
"""
import logging
import os
import time
from collections import deque, OrderedDict
import arg_reader as ar
from ambiguous_split import split_sequence, rb_tree_ambs
from exception_manager import get_time, InsufficientArguments, print_memory_usage
from file_manager import read_byte_numpy, reads, msgunpack_dict, append_file_name, json_it, unjson_it, chr_splits, write_array_to_byte, read_byte_array
from kasai import kasai, inverse1
from true_address import true_address_with_sort, true_address_no_sort
from rb_tree import RedBlackTree


def driver():
    """
        part 0: read genome and s_array files
        part 1: compute lcp array with kasai algorithm
        part 2: find groups of unambiguous and ambiguous chars in genome (unam/ambs)
        part 3: with provided length and distance, use unam/ambs to find valid start addresses of k-mers

    :return:

    """
    try:

        # read arguments
        args = ar.sa_lcp_read_args()
        if not args.outfile:
            raise InsufficientArguments

        # set logging for quiet/verbose
        VERBOSE = not args.verbose
        log = logging.info
        logging.basicConfig(level=logging.WARNING if VERBOSE else logging.INFO, format="%(message)s")

        iprint = logging.warning

        # part 0
        # keep s_array
        genome, past, s_array, start = _part_0(print=log, args=args)
        # part 1
        past, lcp, inv_suff = _part_1(genome=genome, past=past, s_array=s_array, args=args, print=log)
        # part 2
        au, past = _part_2(past=past, args=args, print=log)
        # part 3
        _part_3(lcp=lcp, au=au, args=args, past=past, print=log, inv_suff=inv_suff)

        # print total time that passed
        iprint("time elapsed: %s", str(time.time() - start))
        iprint("Program end.\n")

    except Exception as e:
        raise


def _part_0(args=None, print=print):
    try:
        start = time.time()
        # _____________________________________________
        print('\n_____________________________________')
        print('PART 0: READ ARGS AND GENOME/SA FILES')
        print('_____________________________________\n')
        # _____________________________________________

        past = start
        print('reading SA...\n')

        # read suffix array from bytes to ints
        # reading with numpy then converting to 1-D array much slower than array.array
        # however, array cannot read files larger than ~3GB

        s_array = read_byte_numpy(filename=args.SA)
        print('SA read.\n')
        past = get_time(past, print=print)
        print('reading genome...\n')

        # read with Reads instead
        # ! genome has ambs
        #genome = reads(filename=args.genome)
        chrs, genome = chr_splits(filename=args.genome)


        json_it(data=chrs,filename=append_file_name(args.outfile+"json_chrs"))

        print('genome read.\n')
        past = get_time(past, print=print)

        # TODO: change below line as necessary
        # args.LCPfile = '../data/lcp_pickle'

        return genome, past, s_array, start

    except Exception as e:
        raise


def _part_1(genome, past, s_array, args=None, print=print):
    try:
        # check if args.LCPfile exists
        # if it does, read the pickle file instead of calculating new lcp
        inv_suff = []

        # ___________________________________________
        print('\n_____________________________________')
        print('PART 1: COMPUTE LCP ARRAY')
        print('_____________________________________\n')
        # ____________________________________________

        # if user has specified a LCP pickle file that already exists
        if args.lcpfile and os.path.isfile(path=args.lcpfile):

            print("uniques file exists: ")
            #print(args.lcpfile, '\n')

            # TODO: change this as necessary
            #   hopefully start_uniques will be pickled/jsom/msgpacked in the future
            #lcp = unpickle_dict(filename=args.lcpfile)
            lcp = unjson_it(args.lcpfile)
            # find out what format the lcp was pickled
            if type(lcp) == dict or type(lcp) == OrderedDict:
                key = lcp.keys()[0]
                value = lcp[key]

                if key > value and (100 >= value >= 20):
                    print("old lcp pickle was in format sa:lcp")
                    lcp = deque(lcp.values())

                elif key < value and (100 >= key >= 20):
                    print("old lcp pickle was in format lcp:sa")
                    lcp = deque(lcp.keys())

                else:
                    print("not sure what's going on here for sa_lcp dict")
                    raise KeyboardInterrupt

                s_array = deque(s_array)

            elif type(lcp) == list:
                print('LCP file read as list')

            print("uniques unpacked\n")
            past = get_time(past, print=print)
            print("Computing Unique Start Lengths")

            # combine sa and lcp to form a dict with keys: sa, values: unique_starts
            # TODO: creating OrderedDict consumes too much memory

            filename = append_file_name('json_lcp')
            if args.outfile:
                filename = args.outfile

            if args.inverse:
                inv_suff = unjson_it(args.inverse)
            else:
                inv_suff = inverse1(s_array=s_array)

        else:
            if args.inverse:
                inv_suff = unjson_it(args.inverse)
            inv_suff, lcp = kasai(genome=genome,inv_suff=inv_suff, s_array=s_array, print=print)
            past = time.time()

            # convert suffix array (list) to suffix array (deque) for increased efficiency
            print('Completed.')

            # json it
            filename = append_file_name('json_lcp')
            if args.outfile:
                filename = append_file_name(args.outfile + 'json_lcp')

            print('json\'ing lcp array to %s', filename)
            json_it(data=lcp, filename=filename)

            print('LCP json\'ed!')
            past = get_time(past, print=print)

        return past, lcp, inv_suff

    except Exception as e:
        raise


def _part_2(past, args=None, print=print):
    try:
        # _________________________________
        print('\n_____________________________________')
        print('PART 2: FIND UNAM AND AMBS')
        print('_____________________________________\n')
        # _________________________________
        au = None

        if not args.ambs:
            ambs, unam = split_sequence(filename=args.genome)
            au = rb_tree_ambs(ambs, unam)

            print('args and unam successfully split\n')
            past = get_time(past, print=print)

        else:
            print("Ambs unam file exists: ", args.ambs if args.ambs else args.unam)

            a_u_dict = msgunpack_dict(args.ambs if args.ambs else args.unam)

            print('AMBS and UNAM Unpacked\n')
            past = get_time(past, print=print)

        assert au
        return au, past

    except Exception as e:
        raise


def _part_3(lcp, au:RedBlackTree, past, inv_suff, args=None, print=print):
    try:
        # ____________________________________
        print('\n_____________________________________')
        print('PART 3: VALIDATE STARTING ADDRESSES')
        print('_____________________________________\n')
        # ____________________________________

        true_addresses = []
        # tops = []
        unique_starts = []

        for tup in true_address_with_sort(
            lcp=lcp,
            au=au,
            top=args.length,
            bot=args.low,
            distance=args.distance,
            inv_suff=inv_suff
        ):
            true_addresses.append(tup[0])
            # tops.append(tup[1])
            unique_starts.append(tup[1])

        print('valid addresses calculated\n')
        past = get_time(past, print=print)
        # d_sa = within_distance(in_dict=d_sa, top=args.length, distance=args.distance)
        print('addresses within %s%s', str(args.distance), ' calculated')
        past = get_time(past, print=print)

        filename = append_file_name(filename=args.outfile + 'true_addresses')
        # MSGPACK DOES NOT PRESERVE ORDER
        print('saving true addresses as byte file')
        write_array_to_byte(filename=filename,byte_arr=true_addresses)
        print('saving tops as byte file')
        # filename = append_file_name(filename=args.outfile + 'tops')
        # write_array_to_byte(filename=filename, byte_arr=tops)
        print('saving unique starts as byte file')
        filename = append_file_name(filename=args.outfile + 'unique_starts')
        write_array_to_byte(filename=filename, byte_arr=unique_starts)
        # json_it(d_sa, filename)

        print('default dict msgpack\'ed\n')
        get_time(past, print=print)

        # delete lcp file if final file successfully written
        # don't delete an input file the user has specified
        # if not args.lcpfile:
        #    filename = append_file_name(args.lcpfile + 'json_lcp')
        #    os.remove(filename)

        # write_dictionary(in_dict=d_sa, filename='../default_d_sa_json')

        # print('wrote dictionary to json file\n')

        return

    except InsufficientArguments as e:
        print("Insufficient number of arguments passed!")

    except MemoryError:
        print_memory_usage()
        raise


if __name__ == "__main__":
    print("running program")
    driver()
