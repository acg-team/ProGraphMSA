#!/usr/bin/python2
import os
import sys
import subprocess
import tempfile

cwd="Align"

if __name__ == "__main__":
        infile = sys.argv[1]
        fin = open(infile,'r')
        fcin = tempfile.NamedTemporaryFile()
        namemap = dict()
        iname = 1

        for l in fin:
                if l.startswith('>'):
                        name = l[1:].strip()
                        newname = 'S%d'%(iname)
                        namemap[newname] = name
                        iname += 1
                        fcin.write(">%s\n"%(newname))
                else:
                        fcin.write(l)

        fin.close()
        fcin.flush()

        out = sys.stdout
        errout = open(os.devnull,'w')
        #errout = open("trust.stderr",'w')
        p = subprocess.Popen(['java', '-Xmx200m', '-Xms200m', '-Xmn50m', '-cp', '.', 'nl.vu.cs.align.SelfSimilarity', '-fasta', fcin.name, '-matrix', 'BLOSUM62', '-noseg', '-gapo', '8', '-gapx', '2', '-force', '-procTotal', '1'], stdout=subprocess.PIPE, stderr=errout, cwd=cwd, bufsize=-1)
        #p = subprocess.Popen(['java', '-Xmx200m', '-Xms200m', '-Xmn50m', '-Xbatch', '-XX:ParallelGCThreads=1', '-XX:+UseSerialGC', '-cp', '.', 'nl.vu.cs.align.SelfSimilarity', '-fasta', fcin.name, '-matrix', 'BLOSUM62', '-noseg', '-gapo', '8', '-gapx', '2', '-force'], stdout=subprocess.PIPE, stderr=errout, cwd=cwd, bufsize=-1)

        lines = iter(p.stdout)

        for l in lines:
                if l[0] == '>':
                        newname = l[1:].strip()
                        name = namemap[newname]
                        out.write(">%s\n"%name)

                elif l.startswith("# START LENGTH"):
                        l = lines.next()
                        l = l.split()
                        start = int(l[0])
                        totlength = int(l[1])
                        msa = []
                        starts = [start]
                        lengths = [totlength]


                        for l in lines:
                                if l.find("# Repeat") < 0:
                                        break
                                l = l.split()
                                starts += [int(l[0])]
                                lengths += [int(l[1])]

                        while l[0] != '>':
                                l = lines.next()


                        l = lines.next().strip().upper()
                        length = len(l)
                        msa += [l]
                        for i in xrange(1,len(starts)+1):
                                if i == len(starts) or starts[i] != starts[i-1]+lengths[i-1]:
                                        end = starts[i-1]+lengths[i-1]-1
                                        gaps = sum(map(lambda x:reduce(lambda y,z:y+(z=='-'),x,0), msa))
                                        totlength = end-start+1-gaps
                                        if len(msa) > 1:
                                                out.write("Length: %d residues - nb: XXX  from  %d to %d - Psim:1.0 region Length:%d\n"%(length,start,end,totlength))
                                                out.write('\n'.join(msa))
                                                out.write("\n**********************\n\n")

                                        if i == len(starts):
                                                break

                                        msa = []
                                        start = starts[i]

                                #assert(starts[i] >= starts[i-1]+lengths[i-1])
                                l = lines.next()
                                assert(l[0] == '>')
                                l = lines.next().strip().upper()
                                assert(length == len(l))
                                msa += [l]

                        #l = lines.next().strip().upper()
                        #length = len(l)
                        #msa += [l]
                        #out.write("Length: %d residues - nb: XXX  from  %d to %d - Psim:1.0 region Length:%d\n"%(length,start,end,totlength))
                        #out.write(l+'\n')
                        #end = int(l[0])+int(l[1])
                        #totlength = end-start+1

                        #l = lines.next()
                        #while l[0] == '>':
                        #        l = lines.next().strip().upper()
                        #        assert(len(l) == length)
                        #        msa += [l]
                        #        out.write(l+'\n')
                        #        l = lines.next()
                        #out.write("**********************\n\n")
        fcin.close()
	p.kill()
