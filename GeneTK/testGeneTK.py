# Testing for the whole ISgeneTK package

import mods.gel_visualizer as gv
import mods.seqprop as sp
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
    print ('Checking test data')
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
		reData = pd.read_csv('C:\\Users\\edsun\\Desktop\\IntegratedSciences\\ISTools\\ISToolkit\\data\\restriction_sites2.csv', sep=',')
		plasmid_seqs = 'GTAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA GGATCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGATCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
		re_list = 'AccI BamHI KasI'
		# Testing that digestSeq() works
		max_lengths, lengths_lists = gv.digestSeq(plasmid_seqs, re_list, reData) # NEED SAMPLE SEQUENCE
		assert len(max_lengths) <= len(length_lists), "Please make sure that your inputs are structured correctly."
		for length in max_length:
			assert isinstance (length, int), "Please make sure input is formatted correctly."
		for lengths in lengths_list:
			for length in lengths:
				assert isinstance (length, int), "Please make sure input is formatted correctly."
		assert max_lengths[0] == 355 and max_lengths[1] == 222 and max_lengths[2] == 277, "Restriction lengths are incorrect."
		assert len(lengths_lists[0]) == 3 and len(lengths_lists[1]) == 1 and len(lengths_lists[2]) == 2, "Number of restriction fragments do not match"
		# TEST THAT LENGTHS CORRECT FOR TEST DATA
		assert lengths_lists == [[230, 355, 133], [222], [277, 60]], "Restriction lengths do match correct lengths."

		# Testing that bigDraw() works

		# Testing that smallDraw() works

	def test_SeqProp(self):
		print ("Testing SeqProp...")
		print ("reading restriction site data...")
		reData = pd.read_csv('C:\\Users\\edsun\\Desktop\\IntegratedSciences\\ISTools\\ISToolkit\\data\\restriction_sites2.csv', sep=',')
		# Testing GC, Tm outputs
		# NEED SAMPLE SEQUENCES X 2 (one small < 14, one big > 15)
		gcc_small, gcn_small = sp.getGC('atgcatgcATGC') # SMALL SEQUENCE HERE
		assert gcc_small == 16.67, "Small oligomer GC content does not match." # KNOWN GCC
		assert sp.getTm('atgcatgcATGC') == 28, "Small oligomer Tm does not match"
		gcc_big, gcn_big = sp.getGC('AAAAATGAAAAGGATCCAAAATGAAAAA') # BIG SEQ HERE
		assert gcc_big == 21.43, "Large oligomer GC content does not match." # KNOWN GCC
		Tm = sp.getTm('AAAAATGAAAAGGATCCAAAATGAAAAA', gcc_big)
		assert Tm == 49.67, "Big oligomer Tm does not match." 
		# Reverse Complement
		assert len('AAAAATGAAAAGGATCCAAAATGAAAAA') == len(sp.getRevComp('AAAAATGAAAAGGATCCAAAATGAAAAA')), "Reverse complement length does not match."
		# Testing Start, Stop, Exon, and Re Site Indices and Warnings
		starts, stops = sp.getStartStop('AAAAATGAAAAGGATCCAAAATGAAAAA') 
		assert len(starts) == 2, "Number of start codon indices do not match."
		assert len(stops) == 2, "Number of stop codon indices do not match."
		exons = sp.getExons('AAAAATGAAAAGGATCCAAAATGAAAAA', starts, stops)
		assert len(exons) == min([len(starts), len(stops)]), "Error with exon detection; number of exons does not match min(codons)"
		# READ RE SITES HERE
		re_sites = sp.getRe('AAAAATGAAAAGGATCCAAAATGAAAAA', reData)
		assert len(re_sites) == 2, "Number of restriction sites do not match."


	def test_BUILDR(self):
		print ("Testing BUILDR...")
	# Wait until BUILDR is finished

	def test_MSM(self):
		print ("Testing MSM...")
		try:
			msm.msm('Test_Data/msm_test.tsv')
		except:
			print ('MSM tsv error')
		try:
			msm.msm('Test_Data/msm_test.csv')
		except:
			print('MSM csv error')
	# Wait until MSM is finished
