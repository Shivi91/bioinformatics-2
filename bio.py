import collections
import random
import sys

### replication origin ###

def reverse_complement( i, is_rna=False ):
  l = [ c for c in i ]
  l.reverse()
  i = ''.join(l)
  tx = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
  o = ''
  for c in i:
    o += tx[c]
  if is_rna:
    #tx = { 'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C' }
    return rna( o )
  else:
    return o

def dna( i ):
  '''
>>> dna( 'UAC' )
'TAC'
  '''
  tx = { 
    'U': 'T', 
    'A': 'A', 
    'C': 'C', 
    'G': 'G' }
  o = ''
  for c in i:
    o += tx[c]
  return o

def rna( i ):
  tx = { 'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G' }
  o = ''
  for c in i:
    o += tx[c]
  return o

### sequencing antibiotics ###

def spectrum_weights( spectrum ):
  result = collections.Counter()
  for c in spectrum:
    result[c] += 1
  return result

def matches_weights( candidate, spectrum ):
  result = 1 # include 0 weight
  double = candidate + candidate
  spec = spectrum_weights( spectrum )
  for start in xrange( 0, len(candidate) ): # if abcd -> start = [0, 1, 2, 3]
    for length in xrange( 1, len(candidate) ): # l = [1, 2, 3]
      candidate_sum = sum( double[start:start+length] )
      if spec[ candidate_sum ] > 0:
        result += 1
        spec[ candidate_sum ] -= 1
        #print candidate_sum, " now", spec[candidate_sum]
        
  if spec[ sum( candidate ) ] > 0:
    result += 1
  return result

def matches_weights_str( candidate, spectrum ):
  '''
#>>> matches_weights_str( '57-57-71-99-129-137-170-186-194-208-228-265-285-299-307-323-356-364-394-422-493', '99 71 137 57 72 57' )
#21
>>> matches_weights_str( '156-71-113-114-131-156-113-101-129-128-128-114-128-103-97-131-131-113-131-113-128-115-128-113', '0 71 97 101 103 113 113 113 113 114 114 115 128 128 128 128 129 131 131 131 156 156 184 186 186 200 214 227 227 228 230 231 241 242 242 243 244 244 256 257 262 269 270 287 298 299 301 328 331 340 340 343 345 345 356 358 359 370 370 372 375 383 385 397 400 401 429 430 442 453 454 454 459 462 468 471 472 473 474 485 486 487 498 499 501 512 514 514 542 561 567 570 573 575 581 583 585 590 599 600 600 601 602 610 615 615 616 627 627 630 658 695 696 698 698 698 701 703 704 713 723 728 728 728 728 730 730 731 741 744 747 758 761 769 799 810 817 827 829 831 832 841 841 844 844 851 854 854 857 859 862 872 882 884 886 889 928 928 944 945 947 955 955 958 959 960 966 967 972 972 982 985 990 996 997 1000 1000 1003 1041 1056 1059 1062 1068 1068 1068 1073 1075 1075 1084 1087 1089 1095 1097 1103 1113 1114 1128 1128 1131 1152 1172 1172 1181 1182 1184 1189 1190 1190 1196 1197 1199 1200 1202 1210 1212 1227 1231 1242 1259 1259 1283 1295 1298 1303 1303 1303 1303 1304 1311 1312 1317 1318 1325 1325 1328 1330 1338 1340 1345 1355 1356 1388 1396 1416 1426 1426 1427 1431 1432 1432 1434 1440 1442 1443 1445 1451 1453 1453 1454 1458 1459 1459 1469 1489 1497 1529 1530 1540 1545 1547 1555 1557 1560 1560 1567 1568 1573 1574 1581 1582 1582 1582 1582 1587 1590 1602 1626 1626 1643 1654 1658 1673 1675 1683 1685 1686 1688 1689 1695 1695 1695 1696 1701 1703 1704 1713 1713 1733 1754 1757 1757 1771 1772 1782 1788 1790 1796 1798 1801 1810 1810 1812 1817 1817 1817 1823 1826 1829 1844 1882 1885 1885 1888 1889 1895 1900 1903 1913 1913 1918 1919 1925 1926 1927 1930 1930 1938 1940 1941 1957 1957 1996 1999 2001 2003 2013 2023 2026 2028 2031 2031 2034 2041 2041 2044 2044 2053 2054 2056 2058 2068 2075 2086 2116 2124 2127 2138 2141 2144 2154 2155 2155 2157 2157 2157 2157 2162 2172 2181 2182 2184 2187 2187 2187 2189 2190 2227 2255 2258 2258 2269 2270 2270 2275 2283 2284 2285 2285 2286 2295 2300 2302 2304 2310 2312 2315 2318 2324 2343 2371 2371 2373 2384 2386 2387 2398 2399 2400 2411 2412 2413 2414 2417 2423 2426 2431 2431 2432 2443 2455 2456 2484 2485 2488 2500 2502 2510 2513 2515 2515 2526 2527 2529 2540 2540 2542 2545 2545 2554 2557 2584 2586 2587 2598 2615 2616 2623 2628 2629 2641 2641 2642 2643 2643 2644 2654 2655 2657 2658 2658 2671 2685 2699 2699 2701 2729 2729 2754 2754 2754 2756 2757 2757 2757 2757 2770 2771 2771 2772 2772 2772 2772 2782 2784 2788 2814 2885' )
435
  '''
  candidate_list = [ int(x) for x in candidate.split( '-' ) ]
  spectrum_list = [ int(x) for x in spectrum.split( ' ' ) ]
  return matches_weights( candidate_list, spectrum_list )

def convolution_counter( spectrum, return_amino_counter=False ):
  counter = collections.Counter()
  amino_counter = collections.Counter()
  for x in xrange( 0, len(spectrum) ):
    for y in xrange( x+1, len(spectrum) ):
      difference = abs( spectrum[x] - spectrum[y] )
      if difference > 0:
        counter[difference] += 1
        if difference >= 57 and difference <= 200:
          amino_counter[difference] += 1
  if return_amino_counter:
    return amino_counter
  return counter
 
def convolve( spectrum ):
  '''
>>> convolve( [ 0, 137, 186, 323 ] )
[137, 137, 186, 186, 323, 49]
  '''
  counter = convolution_counter( spectrum )
  result = []
  for v in counter:
    for c in xrange(0, counter[v]):
      result.append( v )
  return result

def add_match( leaders, peptide, fitness, n ):
  # position in leaders
  added = False
  for idx in xrange( 0, len(leaders) ):
    if fitness > leaders[idx]['fitness']:
      leaders.insert( idx, { 'fitness': fitness, 'result': peptide } )
      added = True
      break
  if not added:
    leaders.append( { 'fitness': fitness, 'result': peptide } )
  # remove duds
  if len(leaders) > n:
    required_fitness = leaders[n-1]['fitness']
    if leaders[-1]['fitness'] < required_fitness:
      # find cutoff
      for idx in xrange( n, len( leaders ) ):
        if leaders[idx]['fitness'] < required_fitness:
          del leaders[idx:]
          break

def eval_amino( leaders, peptide, best, n, target, spectrum ):
  #print "trying", current
  #added = False
  total = sum( peptide )
  if total > target: # bust
    pass
  else:
    fitness = matches_weights( peptide, spectrum )
    if total == target and fitness >= best['fitness']:
      best['fitness'] = fitness 
      best['result'] = peptide
      #print "best", best
    #if len(leaders) < n or fitness >= leaders[-1]['fitness']:
    add_match( leaders, peptide, fitness, n )

def expand_leaders( leaders, spectrum, aminos, best, n ):
  found = False
  new_leaders = []
  target = max( spectrum )
  for leader in leaders:
    for amino in aminos:
      eval_amino( new_leaders, leader['result'] + [ amino ], best, n, target, spectrum )
  return new_leaders

def find_best( spectrum, aminos, n ):
  '''
>>> find_best( [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], [129, 137, 71, 99, 57, 194, 170, 186, 79, 91, 58, 95, 113, 115, 128, 136, 148, 151, 156, 157, 162, 166, 171, 200, 178, 65, 66, 72, 80, 87, 109, 123, 153, 77, 105, 121], 60 )  
{'result': [57, 137, 71, 99, 57, 72], 'fitness': 22}
  '''
  leaders = []
  add_match( leaders, [], 0, n )
  best = { 'result': [], 'fitness': 0 }
  r = 0
  while len(leaders) > 0:
    #print "pre board size:", len(leaders), " best score:", leaders[0]['fitness'], " cut score:", leaders[min(len(leaders)-1, n-1)]['fitness']
    leaders = expand_leaders( leaders, spectrum, aminos, best, n )
    r += 1
  return best

### regulatory motifs ###
def expand_motif_recursive( head, tail, errors_remaining, result ):
  alphabet = ( 'G', 'T', 'C', 'A' )
  if len(tail) == 0:
    result.add( ''.join(head) )
  else:
    expand_motif_recursive( head + [ tail[0] ], tail[1:], errors_remaining, result )
    if errors_remaining > 0:
      for a in alphabet:
        expand_motif_recursive( head + [ a ], tail[1:], errors_remaining - 1, result )

def expand_motif( motif, mismatches, candidates ):
  '''
>>> expand_motif( 'AA', 1, set() )
set(['AA', 'AC', 'AG', 'CA', 'AT', 'GA', 'TA'])
>>> expand_motif( 'AA', 2, set() )
set(['AA', 'AC', 'GT', 'AG', 'CC', 'TT', 'CG', 'GG', 'GC', 'AT', 'GA', 'TG', 'CT', 'CA', 'TC', 'TA'])
  '''
  expand_motif_recursive( [], list( motif ), mismatches, candidates )
  return candidates

def find_with_mismatches( what, where, mismatches ):
  '''
>>> find_with_mismatches( 'abc', 'dbbce', 1 )
True
>>> find_with_mismatches( 'abc', 'dbbde', 1 )
False
  '''
  for i in xrange( 0, len(where) - len(what) + 1 ):
    remaining = mismatches
    failed = False
    for j in xrange( 0, len(what) ):
      if where[i+j] != what[j]:
        remaining -= 1
        if remaining < 0:
          failed = True
          break
    if not failed:
      return True
  return False

def expand_all_recursive( head, remaining, result ):
  alphabet = ( 'G', 'T', 'C', 'A' )
  if remaining == 0:
    result.add( ''.join( head ) )
  else:
    for a in alphabet:
      expand_all_recursive( head + [a], remaining - 1, result )

def expand_all( k ):
  '''
>>> expand_all( 2 )
set(['AA', 'AC', 'GT', 'AG', 'CC', 'TT', 'CG', 'GG', 'GC', 'AT', 'GA', 'TG', 'CT', 'CA', 'TC', 'TA'])
  '''
  result = set() 
  expand_all_recursive( [], k, result )
  return result

def find_motifs( dnas, k, mismatches ):
  '''
>>> find_motifs( ( 'ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT' ), 3, 1 )
set(['ATT', 'TTT', 'GTT', 'ATA'])
  '''
  # find all k-mers from first dna
  if len( dnas ) == 0:
    return list()
  
  result = set()
  mutations = set()
  for dna in dnas:
    for i in xrange( 0, len(dna) - k + 1 ): # all kmers in dna
      candidate = dnas[0][i:i+k]
      expand_motif( candidate, mismatches, mutations )

  #print "mutations", mutations
  for mutation in mutations:
    found = True
    for dna in dnas: # look in other dnas for mutations
      #print "looking in", dna
      if not find_with_mismatches( mutation, dna, mismatches ):
        #print "looking in", dna, ": not found!"
        found = False
        break
    if found:
      #print "found", mutation, " everywhere"
      result.add( mutation )
  return result

def distance( dnas, candidate, max_score=None ):
  '''
>>> distance( [ 'TTACCTTAAC' ], 'AAA' )
1
  '''
  total = 0
  for dna in dnas: # each dna
    dna_best = None
    for p in xrange( 0, len(dna) - len(candidate) + 1 ): # each position in dna
      target = dna[p: p+len(candidate)]
      target_score = 0
      for x in xrange(0, len(candidate)):
        if candidate[x] != target[x]:
          target_score += 1
          if target_score > dna_best:
            break
      if dna_best is None or target_score < dna_best:
        dna_best = target_score
    total += dna_best
    if max_score is not None and total > max_score:
      break
  return total
  
def find_median( dnas, k ):
  '''
>>> find_median( [ 'AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTACGGGACAG' ], 3 )
(['ACG', 'GAC'], 2)
  '''
  candidates = expand_all( k )
  #print candidates
  best = None
  best_score = None
  for candidate in candidates:
    score = distance( dnas, candidate, best_score )
    if best_score is None or score < best_score:
      best_score = score
      best = [ candidate ]
      #print "new best: ", best_score, "", best
    elif score == best_score:
      best.append( candidate )
  return ( best, best_score )

def profile_most_probable( text, k, probabilities, prob_map = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 } ):
  '''
>>> profile_most_probable( 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, ( ( 0.2, 0.4, 0.3, 0.1 ), ( 0.2, 0.3, 0.3, 0.2 ), (0.3, 0.1, 0.5, 0.1 ), (0.2, 0.5, 0.2, 0.1 ), (0.3, 0.1, 0.4, 0.2) ) )
(['CCGAG'], 0.0048000000000000004)
>>> profile_most_probable( 'ATCCGACAGAGGCAG', 15, [[0.375, 0.16666666666666666, 0.16666666666666666, 0.2916666666666667], [0.4166666666666667, 0.20833333333333334, 0.125, 0.25], [0.2916666666666667, 0.25, 0.20833333333333334, 0.25], [0.25, 0.20833333333333334, 0.20833333333333334, 0.3333333333333333], [0.20833333333333334, 0.3333333333333333, 0.25, 0.20833333333333334], [0.2916666666666667, 0.16666666666666666, 0.25, 0.2916666666666667], [0.3333333333333333, 0.3333333333333333, 0.16666666666666666, 0.16666666666666666], [0.25, 0.25, 0.16666666666666666, 0.3333333333333333], [0.2916666666666667, 0.20833333333333334, 0.25, 0.25], [0.25, 0.375, 0.25, 0.125], [0.20833333333333334, 0.20833333333333334, 0.2916666666666667, 0.2916666666666667], [0.20833333333333334, 0.25, 0.20833333333333334, 0.3333333333333333], [0.4166666666666667, 0.25, 0.08333333333333333, 0.25], [0.25, 0.375, 0.25, 0.125], [0.16666666666666666, 0.2916666666666667, 0.2916666666666667, 0.25]] )
(['ATCCGACAGAGGCAG'], 2.0540357709176262e-09)
  '''
  best = None
  best_score = None
  for p in xrange( 0, len(text) - k +1 ):
    candidate = text[p:p+k]
    current = 1.0
    for x in xrange( 0, k ):
      next_probability = probabilities[x][prob_map[candidate[x]]]
      current *= next_probability
      #print "prob", next_probability, " acc", current
    #print "before rounding", current
    #current = round( current, 6 )
    #print candidate, "", current
    if best_score is None or current > best_score:
      best = [ candidate ]
      best_score = current
      #print "set best", best, "", best_score
    elif current == best_score:
      best.append( candidate )
      #print "appended", candidate, " best", best, "", best_score
  return ( best, best_score )

def build_probability_profile( dnas, prob_map = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }, smooth_laplace=False ):
  '''
>>> build_probability_profile( [ 'AA', 'AG' ] )
[[1.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.5, 0.0]]
>>> build_probability_profile( [ 'AA', 'AG' ], smooth_laplace=True )
[[0.5, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666], [0.3333333333333333, 0.16666666666666666, 0.3333333333333333, 0.16666666666666666]]
  '''
  result = []
  for i in xrange( 0, len(dnas[0])):
    freqs = collections.Counter()
    for dna in dnas:
      freqs[prob_map[dna[i]]] += 1
    row = []
    for j in xrange( 0, len(prob_map) ):
      if smooth_laplace:
        row.append( ( freqs[j] + 1 ) * 1. / ( len(dnas) + len( prob_map ) ) )
      else:
        row.append( freqs[j] * 1. / len(dnas) )
    result.append( row )
  return result

def probability( text, probabilities, prob_map = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 } ):
  '''
>>> probability( 'ATA', ( ( 0.5, 0.2, 0.1, 0 ), ( 0.1, 0.2, 0.4, 0.3 ), ( 0.1, 0.2, 0.3, 0.4 ) ) )
0.015
  '''
  result = 1.0
  for i in xrange(0, len(text)):
    result *= probabilities[i][prob_map[text[i]]]
  return result

def consensus( dnas ):
  '''
>>> consensus( ( 'AAA', 'ATT', 'AAT' ) )
'AAT'
  '''
  result = []
  for i in xrange( len( dnas[0] ) ):
    freqs = collections.Counter()
    for dna in dnas:
      freqs[dna[i]] += 1
    result.append( freqs.most_common(1)[0][0] )
  return ''.join(result)

def hamming( motifs ):
  '''
>>> hamming( ( 'AAA', 'ATT', 'AAT' ) )
2
  '''
  motif_consensus = consensus( motifs )
  motif_hamming = 0
  for i in xrange(0, len(motifs[0]) ):
    for motif in motifs:
      if motif[i] != motif_consensus[i]:
        motif_hamming += 1
  return motif_hamming

def profile_greedy( dnas, k, smooth_laplace=False ):
  '''
>>> profile_greedy( ( 'GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG' ), 3 )
(['CAG', 'CAG', 'CAA', 'CAA', 'CAA'], 2)
>>> profile_greedy( ( 'GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG' ), 3, smooth_laplace=True )
(['TTC', 'ATC', 'TTC', 'ATC', 'TTC'], 2)
  '''
  # initial best
  #print "starting!"
  best = []
  for dna in dnas:
    best.append( dna[0:k] )
  best_score = hamming( best )
  #print "initial best", best, best_score

  for i0 in xrange( 0, len( dnas[0] ) - k + 1 ): # each k-mer in dnas[0]
    motifs = [ dnas[0][i0: i0+k] ]
    for dna_i in xrange(1, len(dnas)): # dnas from 1..t
      # build profile from motifs 0..i-1
      #print "profiling", motifs
      profile = build_probability_profile( motifs, smooth_laplace=smooth_laplace )
      # pick most probable motif from dna[i]
      best_motif = None
      best_motif_prob = None
      for j in xrange(0, len(dnas[dna_i]) - k + 1): # each kmer in dna_i
        motif = dnas[dna_i][j:j+k]
        motif_probability = probability( motif, profile )
        #print motif, " has prob", motif_probability
        if best_motif_prob is None or motif_probability > best_motif_prob:
          best_motif = motif
          best_motif_prob = motif_probability
          #print "best is", best_motif, "", best_motif_prob
      if best_motif is None:
        print "oh no!"
      else:
        motifs.append( best_motif )
      #print "motifs", motifs
    current_score = hamming( motifs )
    if current_score < best_score:
      best_score = current_score
      best = motifs
      #print "best", motifs, best_score
  return ( best, best_score )      

def profile_randomized_repeated( dnas, k, iterations=1000, smooth_laplace=False ):
  '''
>>> profile_randomized_repeated( ( 'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA'), 8, iterations=1000, smooth_laplace=True )[1]
9
  '''
  best = None
  for i in xrange(0, iterations):
    cand = profile_randomized( dnas, k, smooth_laplace=smooth_laplace )
    if best is None or cand[1] < best[1]:
      best = cand
  return best

def profile_randomized( dnas, k, initial_motifs=None, smooth_laplace=False ):
  # choose a random initial set of kmers as the best
  if initial_motifs is None:
    best = []
    for dna in dnas:
      pos = random.randint( 0, len(dna) - k )
      motif = dna[pos:pos+k]
      best.append( motif )
  else:
    best = initial_motifs
  best_score = hamming( best )
  #print "initial best is", best_score, "", best 
  initial_best_score = best_score
  while True:
    new_motifs = []
    profile = build_probability_profile( best, smooth_laplace=smooth_laplace )
    #print profile
    for dna_i in xrange(0, len(dnas)):
      #best_excluding_current = best[0:dna_i] + best[dna_i+1:len(dnas)]
      #profile = build_probability_profile( best, smooth_laplace=smooth_laplace )
      dna = dnas[dna_i]
      most_probable = profile_most_probable( dna, k, profile )
      #print "most_probable for", dna, " given", profile, " is", most_probable
      new_motifs.append( most_probable[0][0] )
    new_motifs_score = hamming( new_motifs )
    #print "new score", new_motifs_score, " new motifs", new_motifs
    if new_motifs_score < best_score: # improvement
      best_score = new_motifs_score
      best = new_motifs
      #print ( best, best_score )
    else: # no improvement
      if best_score < initial_best_score:
        #print "improvement from", initial_best_score, " to", best_score
        pass
      return ( best, best_score )

def motif_gibbs_repeated( dnas, k, iterations=100, starts=50, smooth_laplace=True ):
  '''
>>> motif_gibbs_repeated( ('CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA'), 8, smooth_laplace=True )[1]
9
  '''
  best = None
  for i in xrange(0, starts):
    cand = motif_gibbs( dnas, k, iterations=iterations, smooth_laplace=smooth_laplace )
    if best is None or cand[1] < best[1]:
      best = cand
      print "iteration", i, "", best
  return best
  
def calculate_kmer_probabilities( text, k, profile, prob_map = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 } ):
  '''
>>> calculate_kmer_probabilities( 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, ( ( 0.2, 0.4, 0.3, 0.1 ), ( 0.2, 0.3, 0.3, 0.2 ), (0.3, 0.1, 0.5, 0.1 ), (0.2, 0.5, 0.2, 0.1 ), (0.3, 0.1, 0.4, 0.2) ) )
[0.00024000000000000003, 0.00048000000000000007, 0.0008000000000000003, 6.000000000000001e-05, 0.00018, 8.000000000000003e-05, 0.00012000000000000004, 8.000000000000003e-05, 8.000000000000003e-05, 0.0005000000000000001, 0.00030000000000000003, 0.00027, 0.00072, 0.0019200000000000007, 0.0002400000000000001, 0.00040000000000000013, 6.000000000000001e-05, 0.00030000000000000003, 0.00040000000000000013, 0.00018, 0.0036, 0.00072, 0.0026999999999999997, 0.00024000000000000006, 0.00108, 0.0004800000000000002, 0.0006000000000000002, 0.00020000000000000006, 0.0009, 0.00072, 0.00144, 0.0007200000000000002, 0.00016000000000000007, 0.0003600000000000001, 0.00016000000000000007, 0.0002400000000000001, 0.0005000000000000001, 0.00030000000000000003, 0.0018, 0.00072, 0.0048000000000000004, 0.00288, 0.0024000000000000002, 0.0006000000000000001, 0.00225, 0.0009]
  '''
  probs = []
  for p in xrange( 0, len(text) - k +1 ):
    candidate = text[p:p+k]
    current = 1.0
    for x in xrange( 0, k ):
      next_probability = profile[x][prob_map[candidate[x]]]
      current *= next_probability
    probs.append( current )
  return probs

def choose_random( probabilities ):
  r = random.random()
  s = sum( probabilities )
  t = 0.0
  for idx in xrange(0, len(probabilities)):
    nxt = t + probabilities[idx] / s
    if t <= r < nxt:
      return idx
    t = nxt

def motif_gibbs( dnas, k, iterations=100, smooth_laplace=False ):
  # random initial motifs
  best_motifs = []
  for dna in dnas:
    pos = random.randint( 0, len(dna) - k )
    motif = dna[pos:pos+k]
    best_motifs.append( motif )
  best_score = hamming( best_motifs )
  for iteration in xrange(0, iterations):
    dna_i = random.randint( 0, len(dnas) - 1 ) # any dna
    best_excluding_current = best_motifs[0:dna_i] + best_motifs[dna_i+1:len(dnas)]
    profile = build_probability_profile( best_excluding_current, smooth_laplace=smooth_laplace )
    candidate_motifs = best_motifs[:]
    kmer_probabilities = calculate_kmer_probabilities( dnas[dna_i], k, profile )
    kmer_position = choose_random( kmer_probabilities )
    candidate_motifs[dna_i] = dnas[dna_i][kmer_position: kmer_position + k]
    candidate_score = hamming( candidate_motifs )
    if candidate_score < best_score:
      best_score = candidate_score
      best_motifs = candidate_motifs
  # done
  return ( best_motifs, best_score )
    
if __name__ == "__main__":
  import doctest
  doctest.testmod()
