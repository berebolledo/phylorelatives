#!/usr/bin/env python
#
# Boris Rebolledo-Jaramillo (boris-at-bx.psu.edu)
#usage: phylorelatives.py [-h] [-i FASTA] [-b INT] [-p] [-r FASTA]
#
#Constructs relatedness of a set of sequences based on the pairwise proportion
#of different sites. It reports the test sequences relatives, tree plot and
#tree in Newick format. One or more test sequences are accepted as long as
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
    robjects.r('''
                add_bs_value = function(tree, boot) {
                tree$node.label = boot$BP
                return(tree)
                }''')
    add_bs_value = robjects.r['add_bs_value']
    cons_tree =  add_bs_value(tree, bootstrap)
    #bs_trees  = bootstrap.rx('trees')[0]
    #cons_tree = ape.consensus(bs_trees,p=0.5)
    return cons_tree


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


def dendro_plot(tree, root=False ):
    """Plot tree to file in ascii format"""
    ape = importr('ape')
    if root:
        newick = list(ape.write_tree(ape.root(tree,root)))[0]
    else:
        newick = list(ape.write_tree(tree))[0]
    t = dendropy.Tree.get_from_string(newick,"newick")
    ascii_tree = t.as_ascii_plot()
    return ascii_tree


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
    parser.add_argument('-j', '--major-only', action='store_true', help='In major-only mode no minor allele sequence is required and only a NJ major alleles tree is generated (Default: Require minor allele sequences)')
    parser.add_argument('--relatives-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--newick-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--trees-out', type=str, metavar='FILE', default=None, help=argparse.SUPPRESS)
    parser.add_argument('-m', '--multi-fasta', type=str, default='multi-fasta.fa', help=argparse.SUPPRESS)
    args = parser.parse_args()
    #parser.print_help()
    
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

    if len(minor_ids) == 0 and not args.major_only:
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
        newick = open(args.multi_fasta+'-newick.nwk','w+')
    if args.trees_out is not None:
        tree_plot_file = open(args.trees_out, 'w+')
    else:
        tree_plot_file = open(args.multi_fasta+'-tree-ascii.txt', 'w+')

    newick.write('%s\n' % (write_nwk(cons_tree)))
    tree_plot_file.write(dendro_plot(cons_tree,root=root_id))
    
    if args.major_only:
        relatives.write('Major allele only mode cannot generate a report')
    else:
        relatives.write('#source\tsample\trelatives\n')
        for node in minor_ids:
            nj_relatives = [relative for relative in dendro_relatives(cons_tree,node) if relative != node]
            relatives.write( 'Neighbor_Joining_tree\t%s\t%s\n' % (node,','.join(sorted(nj_relatives))) )

    newick.close()
    relatives.close()

if __name__ == "__main__": main()
