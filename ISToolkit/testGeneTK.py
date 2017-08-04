# Testing for the whole ISgeneTK package

import mods.gel_visualizer as gv
import mods.SeqProp as sp
import mods.plasmid_builder as pb
import mods.MSM as msm

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("Initializing GeneTK Test...")
 
def teardown_module(module):
    print ("Finalizing GeneTK Test...")
 
def is_setup_function():
    print ("Setting up data tests")
 
def is_teardown_function():
    print ("Tearing down data tests")
 
@with_setup(is_setup_function, is_teardown_function)
def test_data():
    print 'Checking test data'
#    assert multiply(3,4) == 12  <--- Do some assertions with all test data
 

class TestGene:

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

	def test_GelViz(self):
		print ("Testing Gel.Viz...")
		# Testing that digestSeq() works
		max_lengths, lengths_lists = gv.digestSeq() # Fill in with inputs
		assert len(max_lengths) =< len(length_lists), "Please make sure that your inputs are structured correctly."
		for length in max_length:
			assert isinstance (length, int), "Please make sure input is formatted correctly"
		for lengths in lengths_list:
			for lenght in lengths:
				assert isinstance (length, int), "Please make sure input is formatted correctly"

		# Testing that bigDraw() works

		# Testing that smallDraw() works

	def test_SeqProp(self):
		print ("Testing SeqProp")
	# break up into smaller functions and then assert for each
	# if no re sites, start, stop, or if tm range bad print warning

	def test_BUILDR(self):
		print ("Testing BUILDR")
	# Wait until BUILDR is finished

	def test_MSM(self):
		print ("Testing MSM")
	# Wait until MSM is finished
