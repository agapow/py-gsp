"""
Tests for gsp.motifsearch.py (the meme-related implementation).
"""

DATA_FILE = 'tests/meme.xml'

class test_memesearch (object):
	# overhead
	def setup (self):
		print ("settingup test_memeparse class")
		
		import memeparse
		meme_results = memeparse.parse_meme_results (DATA_FILE)
		from gsp import motifsearch
		# search a chunk of seqs otherwise too big to test
		self.gsp = motifsearch.MotifSearch (meme_results, stop=20)
		
	def teardown (self):
		print ("tearing down test_memeparse class")
		
	# vanilla motif search, not filtering motifs
	def test_10percent (self):
		res = self.gsp.search (.1)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 4, "unexpected number of results returned"
		assert len (res[0][0]) == 9, "pattern of unexpected length"
	
	def test_20percent (self):
		res = self.gsp.search (.2)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 8, "pattern of unexpected length"
		
	def test_30percent (self):
		res = self.gsp.search (.3)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 6, "pattern of unexpected length"
		
	def test_40percent (self):
		res = self.gsp.search (.4)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 3, "pattern of unexpected length"
	
	def test_alphabet (self):
		"""
		Check an alphabet of unique symbols is generated correctly.
		"""
		alpha = self.gsp.alpha
		# we know there are 8 motifs and they should be forward & back
		motif_symbols = []
		for x in range (1, 9):
			motif_symbols.append ('%s+' % x)
			motif_symbols.append ('%s-' % x)
		assert sorted (alpha) == sorted (motif_symbols), "alphabets differ"


class test_memesearch_filtered (object):
	# overhead
	def setup (self):
		print ("settingup test_memeparse class")
		
		import memeparse
		meme_results = memeparse.parse_meme_results (DATA_FILE)
		from gsp import motifsearch
		# search a chunk of seqs otherwise too big to test
		self.gsp = motifsearch.MotifSearch (meme_results, 1.0e-05)
		
	def teardown (self):
		print ("tearing down test_memeparse class")
		
	# vanilla motif search, not filtering motifs
	def test_10percent (self):
		res = self.gsp.search (.1)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 9, "pattern of unexpected length"
	
	def test_20percent (self):
		res = self.gsp.search (.2)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 7, "pattern of unexpected length"
		
	def test_30percent (self):
		res = self.gsp.search (.3)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 6, "pattern of unexpected length"
		
	def test_40percent (self):
		res = self.gsp.search (.4)
		print (res)
		print ("The number of results returned is: ", len (res))
		assert len (res) == 2, "unexpected number of results returned"
		assert len (res[0][0]) == 3, "pattern of unexpected length"
	
	def test_alphabet (self):
		"""
		Check an alphabet of unique symbols is generated correctly.
		"""
		alpha = self.gsp.alpha
		# we know there are 8 motifs and they should be forward & back
		motif_symbols = []
		for x in range (1, 9):
			motif_symbols.append ('%s+' % x)
			motif_symbols.append ('%s-' % x)
		assert sorted (alpha) == sorted (motif_symbols), "alphabets differ"



