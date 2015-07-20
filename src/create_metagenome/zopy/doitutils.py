import sys, os, re
from doit.tools import config_changed

# ---------------------------------------------------------------
# convert formatted string to doit dict
# ---------------------------------------------------------------

def enact ( task_name, action, serial=False ):
	"""
	input: 
	  "python :myscript.py -i :input.txt -o ::output.txt -c 5 # :database.txt ::logging.txt"
	output (doit dict):
	  action = "python myscript.py -i input.txt -o output.txt -c 5"
	  file_dep = myscript.py, input.txt, database.txt
	  targets = output.txt, logging.txt
	"""
	if type( action ) is list:
		action = " ".join( formatted_action )
	action_original = action
	# process targets
	targets = []
	for match in re.finditer( "(^| )::(.*?)( |;|$)", action ):
		targets.append( match.group( 2 ) )
	action = re.sub( " ::", " ", action )
	action = re.sub( "^::", "", action )
	# process file deps (always include command file!)
	file_dep = []
	for match in re.finditer( "(^| ):(.*?)( |;|$)", action ):
		file_dep.append( match.group( 2 ) )
	action = re.sub( " :", " ", action )
	action = re.sub( "^:", "", action )
	# remove commented items
	action = re.sub( " +#.*", "", action )
	# expected task dictionary
	doitdict = { 
		"targets":targets, 
		"file_dep":file_dep,
		"actions":[action],
		"uptodate":[config_changed( action_original )],
		"verbosity":2,
	}
	# insert task name if a serial task (yield vs return)
	if serial:
		doitdict["name"] = task_name
	# done
	return doitdict

# ---------------------------------------------------------------
# utilities for helping with serial commands
# ---------------------------------------------------------------

def multiplex ( aaSets, aItems ):
	if aaSets == []:
		return [[k] for k in aItems]
	else:
		aaSets2 = []
		for aSet in aaSets:
			for item in aItems:
				aaSets2.append( aSet+[item] )
		return aaSets2

def multipick ( *args ):
	aaSets = []
	for aItems in args:
		aaSets = multiplex( aaSets, aItems )
	return aaSets

def pyzippick ( *args ):
	assert len( set( [len( k ) for k in args] ) ) == 1, \
		"inconsistent list lengths for zipping"
	return [[args[i][j] for i in range( len( args ) )] \
			for j in range( len( args[0] ) )]

def multireplace ( string, swaps ):
	for i, v in enumerate( swaps ):
		original = "{%d}" % ( i+1 )
		string = string.replace( original, v )
	return string

# ---------------------------------------------------------------
# decorators for single and multiple ops
# ---------------------------------------------------------------

"""
Example:

@op
def mono_task ( ):
	return "grep hello :file.txt > ::found.txt"
"""

def op ( funcTask ):
	def wrapped ( ):
		return enact( funcTask.__name__, funcTask( ) )
	wrapped.create_doit_tasks = wrapped
	return wrapped

"""
Example:

@expop
def serial_task( a1=[1, 2], a2=["A", "B"] ):
	return "grep {1} :file{2}.txt > ::found{1}{2}.txt"

Will evaluate for 1A, 1B, 2A, 2B
"""

def expop ( funcTask ):
	def wrapped ( ):
		aaParams = multipick( *list( funcTask.__defaults__ ) )
		for aParams in aaParams:
			aParams = [str( k ) for k in aParams]
			name = "-".join( [funcTask.__name__] + aParams )
			command = multireplace( funcTask( ), aParams )
			yield enact( name, command, serial=True )
	wrapped.create_doit_tasks = wrapped
	return wrapped

"""
Example:

@zipop
def serial_task( a1=[1, 2], a2=["A", "B"] ):
	return "grep {1} :file{2}.txt > ::found{1}{2}.txt"

Will evalute for 1A, 2B
"""

def zipop ( funcTask ):
	def wrapped ( ):
		aaParams = pyzippick( *list( funcTask.__defaults__ ) )
		for aParams in aaParams:
			aParams = [str( k ) for k in aParams]
			name = "-".join( [funcTask.__name__] + aParams )
			command = multireplace( funcTask( ), aParams )
			yield enact( name, command, serial=True )
	wrapped.create_doit_tasks = wrapped
	return wrapped

# ---------------------------------------------------------------
# classes 
# ---------------------------------------------------------------

class workspace:
    def __init__ ( self ):
        self.root = os.getcwd()
        for path in ["input", "output", "temp"]:
            if not os.path.exists( os.path.join( self.root, path ) ):
                os.mkdir( os.path.join( self.root, path ) )
        self.paths = {}      
    def get ( self, f, tarcheck=False, fincheck=False ):
        if fincheck:
            self.paths[f] = os.path.join( "output", f )
        elif tarcheck:
            self.paths[f] = os.path.join( "temp", f )
        else:
            self.paths[f] = os.path.join( "input", f )
        return self.paths[f]
