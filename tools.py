# -*- coding: utf-8 -*-

# from collections import Counter

import copy
import itertools

# given a population N (list of elements) + an integer b
# returns the powerset of the population starting with coalitions of size b; with emptyset otherwise
def generate_powerset(N, b) :
	if b < 0 :
		print("Error in generate_powerset : var b must be >= 0 but is here "+str(b))
		return []
	coal = []
	for i in range(b,len(N)+1) :
		coal += list(itertools.combinations(N,i))
	return coal


# given two coalitions (tuples) X and Y and a preorder (list of equivalence classes)
# returns True if X is before Y in the preorder ; False if after ; -1 if they are equivalent
def is_before(X,Y,order) :
	is_x = False
	is_y = False
	for eq in order :
		for i in eq :
			if i == X :
				is_x=True
			if i == Y :
				is_y=True
		if is_x :
			# if both were present in the same equivalence class, they are equivalent
			if is_y :
				return -1
			return True
		if is_y :
			return False


# given an element el and a ranking prefs (list of equivalence classes)
# returns True if el is in the ranking; False otherwise 
def is_present(el,prefs) :
	for eq in prefs :
	 	if el in eq :
	 		return True
	return False


# given a list of pairwise preferences (lists [x,y] signifying that x > y) and the size N of the studied population
# returns order over individuals based on expressed pairwise preferences
def join_prefs_ind(prefs,N) :
	res = []

	cands = [i+1 for i in range(N)]
	while len(cands) > 1 :
		# retrieve elements dominated by each element
		tmp = [[y for y in cands if [x,y] in prefs] for x in cands]
		# get nb of dominated elements by each element
		lens = [len(l) for l in tmp]
		if max(lens) == min(lens) :
			res.append(cands)
			return res
		# retrieve index of element with highest number of dominated elements
		i = lens.index(max(lens))
		# this element will be in its own equivalence class
		eq = [cands[i]]
		# we add to the equivalence class any other element to which x is equivalent (i.e. prefs contains both [x,y] and [y,x])
		for el in tmp[i] :
			if cands[i] in tmp[cands.index(el)] :
				eq.append(el)
		res.append(eq)
		cands = [i for i in cands if i not in eq]
	if len(cands) > 0 :
		res.append(cands)
	return res


# given a list of pairwise preferences (lists [x,y] signifying that x > y) and a preorder og
# returns the KT distance between preorder expressed via prefs and og, counting equivalence as wrong if not in og
def KT_CP(prefs,og) :
	p2 = copy.deepcopy(prefs)
	dist = 0
	while p2 :
		# in case of an equivalence
		if p2[0][::-1] in p2 :
			tmp = [x for x in og if p2[0][0] in x and p2[0][1] in x]
			# if this equivalence is not present in the preorder og
			if not tmp :
				# then there is a contradiction, and the distance is increased by 1
				dist += 1
			p2 = [x for x in p2 if x != p2[0] and x != p2[0][::-1]] 
		# if there is no equivalence
		else :
			corr = is_before(p2[0][0], p2[0][1], og)
			# if order is reversed in og (or equivalence), we increase the distance
			if not corr or corr == -1 :
				dist += 1
			p2.remove(p2[0])
	return dist


# given a list of pairwise preferences (lists [x,y] signifying that x > y) and a preorder og
# returns the number of pairs of elements for which the preference in one preorder is opposite to that in the other (meaning equivalences not counted)
def count_inverse(prefs,og) :
	p2 = copy.deepcopy(prefs)
	dist = 0
	while p2 :
		# in case of an equivalence, it cannot lead to a contradiction: we simply remove it
		if p2[0][::-1] in p2 :
			p2 = [x for x in p2 if x != p2[0] and x != p2[0][::-1]] 
		# if there is no equivalence
		else :
			corr = is_before(p2[0][0], p2[0][1], og)
			# if order is reversed in og, we increase the distance
			if not corr :
				dist += 1
			p2.remove(p2[0])
	return dist


# given a list of pairwise preferences WITHOUT EQUIVALENCES (lists [A,B] signifying that A > B) and the studied set of coalitions ps
# returns preorder over coalitions based on expressed pairwise preferences
def join_prefs_coal(prefs,ps) :
	res = []
	for i in range(len(ps)) :
		tmp = [x for x in ps if len([c for c in prefs if c[1]==x])== i]
		res.append(tmp)
	return res


# given two binary vectors v_A and v_B describing sets A and B
# proceeds to lexicographical comparison for leximin (i.e. worst to best rank), returning 1 if A > B, 0 if B > A
def lex_comp_min(v_A,v_B) :
	if v_A == v_B :
		return -1
	for i in range(min(len(v_A), len(v_B))) :
		if v_A[i] < v_B[i] :
			return 1
		if v_B[i] < v_A[i] :
			return 0
	# if identical, the larger set is preferred
	if len(v_A) < len(v_B) :
		return 0
	# otherwise, since v_A != v_B, it must hold that len(v_B) < len(v_A), and A is preferred
	return 1


# given two binary vectors v_A and v_B describing sets A and B
# proceeds to lexicographical comparison for leximax (i.e. best to worst rank), returning 1 if A > B, 0 if B > A
def lex_comp_max(v_A,v_B) :
	if v_A == v_B :
		return -1
	for i in range(min(len(v_A), len(v_B))) :
		if v_A[i] < v_B[i] :
			return 1
		if v_B[i] < v_A[i] :
			return 0
	# if identical, the smaller set is preferred
	if len(v_A) < len(v_B) :
		return 1
	# otherwise, since v_A != v_B, it must hold that len(v_B) < len(v_A), and A is preferred
	return 0


# given a list of pairwise prefs representing the result of CP and a preorder (list of equivalence classes) representing the result of another social ranking method
# returns the result of CP "corrected" by the sr method: if a~b according to CP but the sr method expresses a strict preference, then we use that strict preference in CP
def corrected_CP(cp, sr) :
	res = copy.deepcopy(cp)
	ind = [x for x in cp if x[::-1] in cp]
	while ind :
		# find positions of ind[0][0] and ind[0][1] in sr
		for eq in sr :
			if ind[0][0] in eq :
				# if ind[0][0] > ind[0][1]
				if not ind[0][1] in eq :
					res.remove([ind[0][1],ind[0][0]])
				# otherwise if they're equivalent, then we do not change the prefs of CP over them
				break
			# if ind[0][1] placed before ind[0][0]
			elif ind[0][1] in eq :
				res.remove([ind[0][0], ind[0][1]])
				break
		ind.remove(ind[0][::-1])
		ind.remove(ind[0])
	return res


# given two rankings (list of equivalence classes) v1 and v2 over a population pop
# returns the Kendall-Tau distance between v1 and v2 (i.e. number of pairs for which there is a disagreement in preferences)
def Kendall_Tau(v1,v2,pop) :
	res = 0
	for i in range(len(pop)) :
		for j in range(i+1, len(pop)) :
			if is_before(pop[i],pop[j],v1) != is_before(pop[i],pop[j],v2) :
				res += 1
	return res