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
		max_lengths, lengths_lists = gv.digestSeq() # NEED SAMPLE SEQUENCE
		assert len(max_lengths) ==< len(length_lists), "Please make sure that your inputs are structured correctly."
		for length in max_length:
			assert isinstance (length, int), "Please make sure input is formatted correctly."
		for lengths in lengths_list:
			for length in lengths:
				assert isinstance (length, int), "Please make sure input is formatted correctly."
		# TEST THAT LENGTHS CORRECT FOR TEST DATA
		assert lengths_lists == [[val, val], [val, val], [val, val]], "Restriction lengths do match correct lengths."

		# Testing that bigDraw() works

		# Testing that smallDraw() works

	def test_SeqProp(self):
		print ("Testing SeqProp...")
		# Testing GC, Tm outputs
		# NEED SAMPLE SEQUENCES X 2 (one small < 14, one big > 15)
		gcc_small, gcn_small = sp.getGC(small_seq) # SMALL SEQUENCE HERE
		assert gcc_small == Ns, "Small oligomer GC content does not match." # KNOWN GCC
		gcc_big, gcn_big = sp.getGC(big_seq) # BIG SEQ HERE
		assert gcc_big == Nb, "Large oligomer GC content does not match." # KNOWN GCC
		Tm = sp.getTm(big_seq, gcc_big)
		assert Tm == val, "Melting temperature does not match." 
		# Reverse Complement
		assert len(big_seq) == len(sp.getRevComp(big_seq)), "Reverse complement length does not match."
		# Testing Start, Stop, Exon, and Re Site Indices and Warnings
		starts, stops = sp.getStartStop(big_seq) 
		assert len(starts) == val, "Number of start codon indices do not match."
		assert len(stops) == val, "Number of stop codon indices do not match."
		exons = sp.getExons(big_seq, starts, stops)
		assert len(exons) == min([len(starts), len(stops)]), "Error with exon detection; number of exons does not match min(codons)"
		# READ RE SITES HERE
		re_sites = sp.getRe(sequence, re_sites)
		assert len(re_sites) == val, "Number of restriction sites do not match."


	def test_BUILDR(self):
		print ("Testing BUILDR...")
	# Wait until BUILDR is finished

	def test_MSM(self):
		print ("Testing MSM...")
	# Wait until MSM is finished
