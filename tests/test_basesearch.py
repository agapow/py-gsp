"""
Tests for gsp.__init__.py (the standard gsp implementation).

Remeber that GSP returns the longest sequence that occurs in at least the
threshold of transactions.
"""


class test_std_gsp (object):
	def setup (self):
		print ("settingup test_std_gsp class")
		import gsp
		self.trans = [
			'abcdefgh',
			'abcqqfgh',
			'efghabcd',
			'123abc45',
			'babababa',
			'gggggggg',
			'34121234',
			'mnbcdefn',
			'12345678',
			'ABCDEFGH',
		]
		self.gsp = gsp.GspSearch (self.trans)
		
	def test_10percent (self):
		res = self.gsp.search (.1)
		print ("The number of results returned is: ", res)
		assert len (res) == 10, "unexpected number of results returned"
	
	def test_20percent (self):
		res = self.gsp.search (.2)
		print ("The results are: ", res)
		assert len (res) == 1, "unexpected number of results returned"
		
	def test_30percent (self):
		res = self.gsp.search (.3)
		print ("The results are: ", res)
		assert len (res) == 4, "unexpected number of results returned"
		
	def test_40percent (self):
		res = self.gsp.search (.4)
		print ("The results are: ", res)
		assert len (res) == 1, "unexpected number of results returned"
	
	def test_alphabet (self):
		"""
		Check an alphabet of unique symbols is generated correctly.
		"""
		alpha = self.gsp.alpha
		chars = []
		for t in self.trans:
			for ch in t:
				if ch not in chars:
					chars.append (ch)
		assert len (chars) == len (alpha), "alphabet lengths differ"
		assert sorted(chars) == sorted (alpha), "alphabets differ"
		
	def teardown (self):
		print ("tearing down test_std_gsp class")



