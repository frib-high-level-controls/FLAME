
import logging
_log = logging.getLogger(__name__)

import os, shutil
datadir = os.path.join(os.path.dirname(__file__), 'data')

from collections import defaultdict, OrderedDict

from . import GLPSParser

def default_mangler(elem):
    return r'$(PREF)%(name)s:'%elem

class Generator(object):
    """Helper to generate supporting files for FLAME VIOC
    """

    def __init__(self, inp, mangler=default_mangler, indir=datadir):
        self._M, self._I = mangler, indir
        self._input = inp
        P = GLPSParser()
        with open(inp, 'rb') as F:
            self._D = OrderedDict(P.parse(F))

    def emit_subst(self, outdir):
        efile = {}
        eall = []
        eout = defaultdict(list)

        stype = self._D['sim_type']
        outnames = set()

        nmarkers = 0

        warnmissing = set()

        # collect and group be element type
        for elem in map(dict, self._D['elements']):
            etype = elem['type']
            ename = elem['name']
            rname = self._M(elem)

            if etype=='marker':
                nmarkers += 1

            if rname in outnames:
                N = 1
                while '%s%d'%(rname,N) in outnames:
                    N += 1
                rname = '%s%d'%(rname,N)
                ename = '%s[%d]'%(ename,N)
            outnames.add(rname)

            _log.debug("Inspect %(name)s : %(type)s"%elem)

            tempfile = os.path.join(self._I, etype+'.template')
            if not os.path.isfile(tempfile):
                if etype not in warnmissing:
                    _log.warn("Unhandled element type %s", etype)
                    warnmissing.add(etype)
                continue # if missing unknown.template, then ignore

            eall.append((elem, ename, rname))

            if etype not in efile:
                _log.debug("emit %s", tempfile)
                efile[etype] = tempfile
            eout[etype].append((elem, ename, rname))

        if len(eout)==0:
            raise RuntimeError("No elements!")

        if not os.path.exists(outdir):
            _log.debug('create directory "%s"'%outdir)
            os.makedirs(outdir)

        shutil.copyfile(self._input, os.path.join(outdir, "input.lat"))
        for fname in ('core.db', 'orbit.db', 'common.template'):
            shutil.copyfile(os.path.join(self._I, fname),
                             os.path.join(outdir, fname))

        with open(os.path.join(outdir, 'lattice.substitutions'), 'w') as F:

            F.write('file "common.template"\n{\n')
            for elem, ename, rname in eall:
                F.write('{P="%s", ELEM="%s"}\n'%(rname, ename))

            F.write('}\n\n')

            for etype, elems in eout.items():
                F.write('file "%s"\n{\n'%os.path.basename(efile[etype]))

                for elem, ename, rname in elems:
                    F.write('{P="%s", ELEM="%s"}\n'%(rname, ename))

                F.write('}\n\n')

        with open(os.path.join(self._I, 'st.cmd'), 'r') as F:
            stcmd = F.read()

        stcmd = stcmd%{'nmarkers':nmarkers, 'latticepath':os.path.join(os.getcwd(),self._input)}

        with open(os.path.join(outdir, 'st.cmd'), 'w') as F:
            F.write(stcmd)
        #os.chmod(os.path.join(outdir, 'st.cmd'), 0755)

        for cfile in set(efile.values()):
            shutil.copyfile(cfile, os.path.join(outdir, os.path.basename(cfile)))

def gendb(args):
    G = Generator(args.input)
    G.emit_subst(args.out)

def getargs():
    import argparse
    P = argparse.ArgumentParser()
    P.add_argument('-d','--debug',action='store_true')
    SP= P.add_subparsers()

    PP= SP.add_parser('gendb', help='Generate .substitions and .template files')
    PP.set_defaults(func=gendb)
    PP.add_argument('input', metavar='FILE', help='Lattice file')
    PP.add_argument('out', metavar='DIR', help='Output directory')

    return P.parse_args()

if __name__=='__main__':
    A = getargs()
    logging.basicConfig(level=logging.DEBUG if A.debug else logging.WARN)
    A.func(A)
