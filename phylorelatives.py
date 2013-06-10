#!/usr/bin/env python
#
# Boris Rebolledo-Jaramillo (boris-at-bx.psu.edu)
#usage: phylorelatives.py [-h] [-i FASTA] [-b INT] [-p] [-r FASTA]
#
#Constructs relatedness of a set of sequences based on the pairwise proportion
#of different sites. It reports the test sequences relatives, tree plots and
#trees in Newick format. One or more test sequences are accepted as long as
#their name includes the strict suffix "_minor" or "_test" (i.e. >seq1_minor).
#IMPORTANT: Sequences must have the same length!
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -i FASTA, --input FASTA
#                        This option can be specified multiple times. Sequences
#                        will be added to "multi-fasta.fa". (e.g. -i major1.fa
#                        -i major2.fa -i minor1.fa)
#  -b INT, --bootstrap INT
#                        Change number of replicas. 0 to deactivate. (Default:
#                        1000)
#  -p, --pairwise        Use pairwise deletion of gaps/missing data. (Default:
#                        Complete deletion)
#  -r FASTA, --root FASTA
#                        Root trees using FASTA sequence as outgroup. (Default:
#                        Display unrooted trees)

import sys
import argparse
import array
import dendropy
import rpy2.rinterface
rpy2.rinterface.set_initoptions(('rpy2','--vanilla','--quiet'))
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


def ape_read_dna(infasta):
    """Read multi-fasta into phylo object"""
    ape = importr('ape')
    data = ape.read_dna(file=infasta,format="fasta")
    return data

def ape_nj(data,missing_info_option=0):
    """Return ape nj tree"""
    ape = importr('ape')
    dist = ape.dist_dna(data,model="raw",pairwise=missing_info_option)
    nj_tree = ape.nj(dist)
    return nj_tree

def ape_consensus(tree,tree_function,data,iters=1000):
    """Return majority rule consensus tree"""
    ape = importr('ape')
    tree_func = robjects.r(tree_function)
    bootstrap = ape.boot_phylo(tree, data, tree_func,
                              B=iters, quiet=1, trees=1)
    bs_trees  = bootstrap.rx('trees')[0]
    cons_tree = ape.consensus(bs_trees,p=0.5)
    return cons_tree

#def parsimony_tree(tree,data):
#    """Return Maximum Parsimony tree from supplied tree"""
#    phangorn = importr('phangorn')
#    mp_tree = phangorn.optim_parsimony(tree,phangorn.phyDat(data), trace=0)
#    return mp_tree


def dendro_relatives(tree,minor):
    """Return minor allele sequence relatives in tree"""
    ape = importr('ape')
    newick = list(ape.write_tree(tree))[0]
    t = dendropy.Tree.get_from_string(newick,"newick")
    minor_leaf = [node for node in t.leaf_nodes()
                  if node.get_node_str() == minor][0]
    parent = minor_leaf.parent_node
    relatives = []
    while len(relatives) == 0:
        output = [relative.get_node_str() for relative in parent.leaf_nodes()]
        relatives = [relative for relative in output if not (relative.endswith('minor') or relative.endswith('test'))]
        parent = parent.parent_node
    return output

def plot_tree(outfile, tree1, tree2=None, root=False):
    """Generate tree(s) plot"""
    ape = importr('ape')
    graphics  = importr('graphics')
    grdevices = importr('grDevices')
    if tree2 is None:
        grdevices.png(file=outfile, width=1024, height=768)
        if root:
            ape.plot_phylo(ape.root(tree1,root),edge_width=2,cex=1,underscore=1)
        else:
            ape.plot_phylo(tree1,edge_width=2,cex=1,underscore=1)
        graphics.title(main='Neighbor Joining')
        grdevices.dev_off()
    elif tree2 is not None:
        grdevices.png(file=outfile, width=1024, height=768)
        graphics.par(mfcol=array.array('i',[1,2]))
        if root:
            ape.plot_phylo(ape.root(tree1,root),edge_width=2,cex=1,underscore=1)
        else:
            ape.plot_phylo(tree1,edge_width=2,cex=1,underscore=1)
        graphics.title(main='Neighbor Joining', cex=1.5, font=2)
        if root:
            ape.plot_phylo(ape.root(tree2,root),edge_width=2,cex=1,underscore=1)
        else:
            ape.plot_phylo(tree2,edge_width=2,cex=1,underscore=1)
        graphics.title(main='Maximum Parsimony',cex=1.5, font=2)
        grdevices.dev_off()
    return

def write_nwk(tree):
    "Write proper Newick string"
    ape = importr('ape')
    nwk = list(ape.write_tree(tree))[0]
    return nwk

def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Constructs relatedness of a set of sequences based on the pairwise proportion of different sites. It reports the test sequences relatives, tree plots and trees in Newick format. One or more test sequences are accepted as long as their name includes the strict suffix "_minor" or "_test" (i.e. >seq1_minor). IMPORTANT: Sequences must have the same length!', epilog='Boris Rebolledo-Jaramillo (boris-at-bx.psu.edu)')
    parser.add_argument('-i', '--input', metavar='FASTA', action='append', type=str, help='This option can be specified multiple times. Sequences will be added to "multi-fasta.fa". (e.g. -i major1.fa -i major2.fa -i minor1.fa)')
    parser.add_argument('-b', '--bootstrap', type=int, metavar='INT',default=1000, help='Change number of replicas. 0 to deactivate. (Default: 1000)')
    parser.add_argument('-p', '--pairwise', action='store_true', help='Use pairwise deletion of gaps/missing data. (Default: Complete deletion)')
    parser.add_argument('-r', '--root', type=str, metavar='FASTA', default=False, help='Root trees using FASTA sequence as outgroup. (Default: Display unrooted trees)')
    parser.add_argument('--relatives-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--newick-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--trees-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('-m', '--multi-fasta', type=str, default='multi-fasta.fa', help=argparse.SUPPRESS)
    args = parser.parse_args()
    
    
    if args.input:
        for fasta in args.input:
            try:
                open(fasta)
            except:
                sys.exit("\nERROR: Could not open %s\n" % fasta)
    try:     
        multi = open(args.multi_fasta, 'w+')
    except:
       sys.exit("\nERROR: Could not create %s\n" % args.multi_fasta)

    for fasta in args.input:
        for line in list(open(fasta)):
            multi.write(line.replace('-', '_'))

    
    if args.root:
        try:
            root = list(open(args.root))
            root_id = [line.strip()[1:].replace('-', '_') for line in root if line.strip().startswith(">")][0]
            for line in root:
                multi.write(line.replace('-', '_'))
        except:
            sys.exit("\nERROR: Could not open %s\n" % args.root)
    else:
        root_id = args.root
    multi.close()
 
    try:
        data = ape_read_dna(args.multi_fasta)
    except:
         sys.exit("\nERROR: Check existence or proper format of %s\n" % args.multi_fasta)
    
    # Get sequence ids in alignment and identify test sequences
    fasta_ids = [seqname.strip()[1:] for seqname in list(open(args.multi_fasta)) if seqname.startswith('>')]
    minor_ids = [seqname for seqname in fasta_ids if seqname.endswith('minor') or seqname.endswith('test')]

    if len(minor_ids) == 0:
        sys.exit("\nERROR: No test sequences found. _minor or _test suffixes are required in the sequence name!\n")
    else:
        pass

    if args.pairwise:
        nj_tree = ape_nj(data,1)
        nj_func = 'function (xx) nj(dist.dna(xx, model="raw", pair=1))'
    else:
        nj_tree = ape_nj(data)
        nj_func = 'function (xx) nj(dist.dna(xx, model="raw"))'
    
    if args.bootstrap == 0:
        cons_tree = nj_tree
    elif args.bootstrap !=1000:
        cons_tree = ape_consensus(nj_tree,nj_func,data,iters=args.bootstrap)
    else:
        cons_tree = ape_consensus(nj_tree,nj_func,data)
    
    # Generate report, trees, and Newick strings
    if args.relatives_out is not None:
        relatives = open(args.relatives_out,'w+')
    else:
        relatives = open(args.multi_fasta+'-relatives.tab','w+')
    if args.newick_out is not None:
        newick = open(args.newick_out,'w+')
    else:
        newick = open(args.multi_fasta+'-newick.txt','w+')
    if args.trees_out is not None:
        tree_plot_file = args.trees_out
    else:
        tree_plot_file = args.multi_fasta+'-tree.png'

    newick.write('>Neighbor_Joining\n%s\n' % (write_nwk(cons_tree)))
    relatives.write('#source\tsample\trelatives\n')
    for node in minor_ids:
        nj_relatives = [relative for relative in dendro_relatives(cons_tree,node) if relative != node]
        relatives.write( 'Neighbor_Joining_tree\t%s\t%s\n' % (node,','.join(sorted(nj_relatives))) )
    plot_tree(tree_plot_file,cons_tree,root=root_id)


# Parsimony requires X11 related dependencies. Tricky
#    if args.max_parsimony:
#        mp_tree = parsimony_tree(cons_tree,data)
#        newick.write('>Maximum_Parsimony\n%s\n' % (write_nwk(mp_tree)))
#        for node in minor_ids:
#            mp_relatives = [relative for relative in dendro_leaves(mp_tree,node) if relative != node]
#            relatives.write( 'Maximum_Parsimony_tree\t%s\t%s\n' % (node,','.join(sorted(mp_relatives))) )
#        plot_tree(tree_plot_file, tree1=cons_tree, tree2=mp_tree, root=args.root)
#    else:
#        plot_tree(tree_plot_file,cons_tree,root=args.root)

    newick.close()
    relatives.close()

if __name__ == "__main__": main()
