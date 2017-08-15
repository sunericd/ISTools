# Testing for the whole ISProtTK package

import mods.ProtQuant as pq
import numpy.testing as npt

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
		standard_concs, standard_fluors = pq.readStandards('Test_Data/test_standards.csv')
		assert standard_concs[0] == 2, "Error reading standard concentrations."
		assert standard_fluors[0][0] == 4, "Error reading standard data."
		# Testing readProteins()
		prot_fluors = pq.readProteins('Test_Data/test_proteins.csv')
		assert prot_fluors[0][0] == 3, "Error reading protein data."
		# Testing averageRows()
		assert pq.averageRows(standard_fluors) == [4, 2, 1, 0.5, 0.25], "Error with averaging standard_fluors. Check the test data."
		assert pq.averageRows(prot_fluors) == [4, 2, 1, 0.5], "Error with averaging prot_fluors. Check the test data."
		# Testing regressor()
		prot_concs = pq.regressor(standard_concs, standard_fluors, prot_fluors, 'Linear')
		npt.assert_approx_equal(prot_concs[0], 2), "Error: Regression value is off for protein #1"
		npt.assert_approx_equal(prot_concs[0], 1), "Error: Regression value is off for protein #2"
		npt.assert_approx_equal(prot_concs[0], 0.5), "Error: Regression value is off for protein #3"
		npt.assert_approx_equal(prot_concs[0], 0.25), "Error: Regression value is off for protein #4"



		