#! /usr/bin/env python

import sys, re, argparse
from zopy.utils import path2name, warn
from zopy.dictation import polymap

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_proppattern = r"^([^\s]*?): (.*)"
c_gopattern = r"(GO:[0-9]+)"
c_relationships = [
    "part_of",
    "regulates",
    "positively_regulates",
    "negatively_regulates",
]
c_namespace_convert = {
    "biological_process":"BP", 
    "molecular_function":"MF",
    "cellular_component":"CC",
}

# ---------------------------------------------------------------
# begin Term class
# ---------------------------------------------------------------

class Term:

    """Node representing a term in the go ontology"""

    def __init__ ( self, props ):
        
        # set from props dict
        self.id = props["id"][0]
        self.name = props["name"][0]
        self.namespace = props["namespace"][0]
        self.namespace_short = c_namespace_convert[self.namespace]
        self.is_obsolete = True if "is_obsolete" in props else False
        self.parent_ids = self.extract_parent_ids( props )
        # set by the ontology constructor
        self.parents = []
        self.children = []
        self.genes = set()
        self.depth = None
        self.is_root = False
        self.is_leaf = False
        # set by optional methods from the ontology
        self.progeny = None
        self.progeny_genes = None
        self.is_informative = False
        self.is_acceptable = True

    def __repr__ ( self ):
        return "|".join( [
            self.id, 
            self.namespace_short,
            "%02d" % ( self.depth ),
            self.name,
        ] )

    def extract_parent_ids ( self, props ):
        parent_ids = []
        # extract is_a relationships
        # format: list of "goid ! definition"
        for line in props.get( "is_a", [] ):
            parent_ids.append( line.split( )[0] )
        # extract "relationship" relationships
        # format: list of "reltype goid ! definition"
        for line in props.get( "relationship", [] ):
            items = line.strip()
            if items[0] in c_relationships:
                parent_ids.append( items[1] )
        return parent_ids

    def add_genes ( self, gene ):
        gene = [gene] if type( gene ) is not set else gene
        self.genes.update( gene )

    def get_progeny ( self ):
        """ memoized via instance variable """
        if self.progeny is None:
            self.progeny = {self}
            for cnode in self.children:
                self.progeny.update( cnode.get_progeny() )
        return self.progeny

    def get_progeny_genes ( self ):
        """ memoized via instance variable """
        if self.progeny_genes is None:
            self.progeny_genes = set()
            self.progeny_genes.update( self.genes )
            for cnode in self.children:
                self.progeny_genes.update( cnode.get_progeny_genes() )
        return self.progeny_genes
        
# ---------------------------------------------------------------
# end of Term class
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# begin Ontology class
# ---------------------------------------------------------------

class Ontology:

    """Parses an OBO file to build an ontology"""

    def __init__ ( self, path ):

        # get nodes
        all_props = map( self.chunk_to_dict, self.chunk_obo( path ) )
        self.nodes = {}
        for props in all_props:
            node = Term( props )
            if not node.is_obsolete and re.search( c_gopattern, node.id ):
                self.nodes[node.id] = node

        # populate parent/child relationships
        for cnode in self.nodes.values():
            for parent_id in cnode.parent_ids:
                pnode = self.nodes[parent_id]
                cnode.parents.append( pnode )
                pnode.children.append( cnode )

        # identify roots, leaves
        self.roots = []
        self.leaves = []
        for node in self.nodes.values():
            if len( node.parents ) == 0:
                node.is_root = True
                self.roots.append( node )
            if len( node.children ) == 0:
                node.is_leaf = True
                self.leaves.append( node )

        # add depth information
        def recurse_depth ( node ):
            for cnode in node.children:
                if ( cnode.depth is None ) or \
                   ( cnode.depth > node.depth + 1 ):
                    cnode.depth = node.depth + 1
                    recurse_depth( cnode )
        for root in self.roots:
            root.depth = 0
            recurse_depth( root )

        # track genes added to the tree
        self.attached_genes = {}

    # ---------------------------------------------------------------
    # end of init method
    # ---------------------------------------------------------------

    def chunk_obo ( self, path ):
        """__init__ helper function"""
        in_terms = False
        chunks = []
        with open( path ) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    chunks.append( [] )
                    in_terms = True
                elif in_terms and line == "[Typedef]":
                    # come at the end after Terms
                    in_terms = False
                elif in_terms and line != "":
                    chunks[-1].append( line )
        return chunks

    def chunk_to_dict ( self, chunk ):
        """__init__ helper function"""
        props = {}
        current = None
        for line in chunk:
            match = re.search( c_proppattern, line )
            if match:
                current, text = match.groups()
                props.setdefault( current, [] ).append( text )
            elif current is not None:
                props.setdefault( current, [] ).append( line )
        return props

    def attach_genes ( self, path ):
        """assigns genes from a polymap to the tree"""
        for gene, terms_dict in polymap( path ).items():
            self.attached_genes[gene] = 1
            for term in terms_dict:
                self.nodes[term].add_genes( gene )
   
    def set_informative ( self, count ):
        def recurse_set_informative( node ):
            if len( node.get_progeny_genes() ) >= count:
                node.is_informative = True
                # check that no children will be called informative
                for cnode in node.children:
                    if len( cnode.get_progeny_genes() ) >= count:
                        node.is_informative = False
                        break
            else:
                for pnode in node.parents:
                    recurse_set_informative( pnode )
        for leaf in self.leaves:
            recurse_set_informative( leaf )
 
    def iter_nodes ( self ):
        for nodeid, node in self.nodes.items():
            yield node

# ---------------------------------------------------------------
# end of Ontology class
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

if __name__ == "__main__":

    # argument parsing (python argparse)
    parser = argparse.ArgumentParser()
    parser.add_argument( '--obo', default=None, help='' )
    parser.add_argument( '--map', default=None, help='' )
    parser.add_argument( '--depth', default=None, help='' )
    parser.add_argument( '--grep', default=None, nargs="+", help='' )
    parser.add_argument( '--namespace', default=None, \
                         choices=["BP", "MF", "CC"], nargs="+", help='' )
    parser.add_argument( '--informative', default=None, help='' )
    parser.add_argument( '--ignore_progeny', action="store_true", help='' )
    parser.add_argument( '--terms_only', action="store_true", help='' )
    parser.add_argument( '--output', default=None, help='' )
    args = parser.parse_args()

    # warnings
    if args.ignore_progeny:
        warn( "Only considering direct annotations. Annotations will not rise through the DAG." )

    # load data (attach genes?)
    obo = Ontology( args.obo )
    if args.map is not None:
        obo.attach_genes( args.map )

    # set informative?
    if args.informative is not None:
        cutoff = float( args.informative )
        if cutoff < 1:
            print >>sys.stderr, "Intepretting informative cutoff as fraction of annotated genes"
            cutoff *= len( obo.attached_genes )
        cutoff = int( cutoff )
        obo.set_informative( cutoff )

    # trim tree
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    for node in obo.iter_nodes():
        # depth cut
        if args.depth is not None and node.depth != int( args.depth ):
            node.is_acceptable = False
        # grep cut
        if args.grep is not None:
            for pattern in args.grep:
                if not re.search( pattern, node.name ):
                    node.is_acceptable = False
                    break
        # namespaces cut
        if args.namespace is not None:
            if node.namespace_short not in args.namespace:
                node.is_acceptable = False
        # informative cut
        if args.informative and not node.is_informative:
            node.is_acceptable = False
        # output node?
        if node.is_acceptable:
            outline = [str( node )]
            if not args.terms_only:
                outline += list( node.get_progeny_genes() if not args.ignore_progeny else node.genes )
            print >>fh, "\t".join( outline )
