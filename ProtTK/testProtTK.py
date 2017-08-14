# Testing for the whole ISgeneTK package

import mods.ProtQuant as pq

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("Initializing ProtTK Test...")
 
def teardown_module(module):
    print ("Finalizing GeneTK Test...")
 
def is_setup_function():
    print ("Setting up data tests")
 
def is_teardown_function():
    print ("Tearing down data tests")
 
@with_setup(is_setup_function, is_teardown_function)
def test_data():
    print ('Checking test data')
#    assert multiply(3,4) == 12  <--- Do some assertions with all test data
 

class TestProt:

	def setup(self):
		print ("Setting up next test")

	def teardown(self):
		print ("Tearing down test")

	@classmethod
	def setup(self):
		print ("Setting up test")

	@classmethod
	def teardown(self):
		print ("Tearing down test")

	def test_ProtQuant(self):
		print ("Testing ProtQuant...")
		# Testing readStandards()
		standard_concs, standard_fluors = pq.readStandards(FILENAME)
		assert standard_concs[0] == VAL, "Error reading standard concentrations."
		assert standard_fluors[0][0] == VAL, "Error reading standard data."
		# Testing readProteins()
		prot_fluors = pq.readProteins(FILENAME)
		assert prot_fluors[0][0] == VAL, "Error reading protein data."
		# Testing averageRows()
		assert pq.averageRows(standard_fluors) == VAL, "Error with averaging standard_fluors. Check the test data."
		assert pq.averageRows(prot_fluors) == VAL, "Error with averaging prot_fluors. Check the test data."
		# Testing regressor()
		prot_concs = pq.regressor(standard_concs, standard_fluors, prot_fluors, 'Linear')
		assert prot_concs[0] == VAL, "Error: Regression value is off for protein #1"
		assert prot_concs[1] == VAL, "Error: Regression value is off for protein #2"
		assert prot_concs[2] == VAL, "Error: Regression value is off for protein #3"



		