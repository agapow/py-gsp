"""
GSP modified to use motifs as extracted from a MEME search.
"""

__version__ = '0.1'


### IMPORTS

import memeparse

from gsp import basesearch


### CONSTANTS & DEFINES

NULL_SYMBOL = 'X'


### CODE ###

def find (needle, haystack):
	"""
	Return the index of the needle in the haystack
	
	Parameters:
		needle: any iterable
		haystack: any other iterable
		
	Returns:
		the index of the start of needle or -1 if it is not found.
	
	Looking for a sub-list of a list is actually a tricky thing. This
	approach uses the Boyer-Moore-Horspool algorithm. Needle and haystack
	should be any iterable, as long as their elements are hashable.

	Example:
	
		>>> find ([1, 2], [1, 1, 2])
		1
		>>> find ((1, 2, 3), range (10))
		1
		>>> find ('gh', 'abcdefghi')
		6
		>>> find ([2, 3], [7, 8, 9])
		-1

	"""
	h = len (haystack)
	n = len (needle)
	skip = {needle[i]: n - i - 1 for i in range(n - 1)}
	i = n - 1
	while i < h:
		for j in range(n):
			if haystack[i - j] != needle[-j - 1]:
				i += skip.get(haystack[i], n)
				break
		else:
			return i - n + 1
	return -1


class MotifSearch (basesearch.GspSearch):
	
	def __init__ (self, raw_transactions, mtf_threshold=0.0, start=None, stop=None):
		"""
		C'tor, simply shaping the raw transactions into a useful form.
		"""
		self.mtf_threshold = mtf_threshold
		self.start = start
		self.stop = stop
		super (MotifSearch, self).__init__ (raw_transactions)
		
	def process_transactions (self, raw_transactions):
		"""
		Create the alphabet & (normalized) transactions.
		"""
		def motif_to_item (m, threshold, rev=False):
			"""
			Given motif information, produce an item (symbol)
	
			Parameters:
				m: a motif information tuple
				threshold: is p-value below this threshold, replace with null symbol
				rev: if True, produce the symbol for the opposite strand, i.e. flip
					the sign
			"""
			if threshold and (threshold < m.pvalue):
				return NULL_SYMBOL
			else:
				strand = m.strand
				if rev:
					if strand == 'plus': strand = 'minus'
					elif strand == 'minus': strand = 'plus'
			
				item_name = '%s_%s' % (m.motif_id, strand)
				item_name = item_name.replace ('motif_', '')
				item_name = item_name.replace ('_plus', '+').replace ('_minus', '-')
			
				return item_name
	
		# extract the sequence of motifs for each sequence
		seq_motifs = []
		for seq in raw_transactions.scanned_sites_summary[self.start:self.stop]:
			seq_id = seq.sequence_id
			seq_motifs.append ({
				'id': seq_id,
				'p_value': seq.pvalue,
				'name': raw_transactions.seq_id_to_name (seq_id),
				'motifs': seq.scanned_sites,
			})
			assert (seq.scanned_sites) == (sorted (seq.scanned_sites, key=lambda x: x.position))
			#print (seq_motifs)
			
		# now create transactions for each sequence:
		# each transactions is a pair of forward and reverse
		self.transactions = []
		for s in seq_motifs:
			self.transactions.append (
				{
					'forward': [motif_to_item (m, self.mtf_threshold) for m in s['motifs']],
					'reverse': [motif_to_item (m, self.mtf_threshold, rev=True) for m in s['motifs'][::-1]],
				}
			)
			
		# build alphabet
		raw_alpha = []
		for v in self.transactions:
			raw_alpha.extend (v['forward'])
			raw_alpha.extend (v['reverse'])
		self.alpha = set ([a for a in raw_alpha if a != NULL_SYMBOL])
		assert self.alpha, "alphabet should have non-zero length"
		print ("There are %s members of the initial alphabet." % len (self.alpha))
	
	def generate_init_candidates (self):
		return [[x] for x in self.alpha]
	
	def filter_candidates (self, trans_min):
		"""
		Return a list of the candidates that occur in at least the given number of transactions.
		"""
		filtered_candidates = []
		for c in self.candidates:
			curr_cand_hits = self.single_candidate_freq (c)
			if trans_min <= curr_cand_hits:
				filtered_candidates.append ((c, curr_cand_hits))
		return filtered_candidates
		
	def single_candidate_freq (self, c):
		"""
		Return true if a candidate is found in the transactions.
		"""
		hits = 0
		for v in self.transactions:
			if self.search_transaction (v['forward'], c) != -1:
				hits += 1
			elif self.search_transaction (v['reverse'], c) != -1:
				hits += 1
		return hits
		
	def search_transaction (self, t, c):
		"""
		Does this candidate appear in this transaction?
		"""
		# print ("T-C: %s %s %s!" % (t, c, find (c, t)))
		return find (c, t)
	

### END ###

