#!/usr/bin/env python

import os
import tempfile
import signal


def cleanup():
    if main.fntmp != None and os.path.exists(main.fntmp):
        # print main.fntmp
        os.unlink(main.fntmp)
        main.fntmp = None

    if main.fntmp2 != None and os.path.exists(main.fntmp2):
        # print main.fntmp2
        os.unlink(main.fntmp2)
        main.fntmp2 = None

    if main.fnmwks != None:
        for f in main.fnmwks:
            # print f
            os.unlink(f)
        main.fnmwks = None


def signal_handler(sig, frame):
    print 'Got abort signal...'
    cleanup()
    sys.exit(0)


def main():
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        main.fntmp = fntmp = tempfile.mktemp()
        os.system('./02_par_merge.sh %s' % fntmp)

        # fake merged mwks to get the output from 03_par_collect_PS_firing.sh
        main.fnmwks = fnmwks = []
        # actual commands to run for each mwk file
        cmds = {}
        for l0 in open(fntmp, 'rt').readlines():
            l = l0.strip()
            f = l.split()[2]
            fnmwks.append(f)
            cmds[f] = l + ' && '

            if not os.path.exists(f):
                open(f, 'wt').close()


        main.fntmp2 = fntmp2 = tempfile.mktemp()
        os.system('./03_par_collect_PS_firing.sh %s' % fntmp2)
        for l0 in open(fntmp2, 'rt').readlines():
            l = l0.strip()
            f = l.split()[1]

            if f not in cmds: continue
            cmds[f] += l

        # print the completed job list
        for f in fnmwks:
            print cmds[f]
    finally:
        cleanup()

main.fntmp = None
main.fntmp2 = None
main.fnmwks = None


if __name__ == '__main__':
    # assumes no arguments
    main()
