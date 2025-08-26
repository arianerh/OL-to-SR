# -*- coding: utf-8 -*-


import copy
import itertools
import time

import tools

# given an order over elements of a population pop
# returns the order over the powerset of pop using minmax
def minmax(order, pop) :
	res = []
	to_extend = tools.generate_powerset(pop, 1)
	for i in range(len(order)-1,0,-1) :
		tmp = []
		worst = [C for C in to_extend if order[i][0] in C]
		to_extend = [C for C in to_extend if C not in worst]
		k= 0
		while k in range(i) and worst :
			best = [C for C in worst if order[k][0] in C]
			tmp.append(copy.deepcopy(best))
			worst = [C for C in worst if C not in best]
			k += 1
		res = tmp + [[order[i]]] + res
	return [[order[0]]] + res


# given an order over elements of a population pop
# returns the order over the powerset of pop using maxmin
def maxmin(order, pop) :
	res = []
	to_extend = tools.generate_powerset(pop, 1)
	for i in range(len(order)) :
		tmp = []
		best = [C for C in to_extend if order[i][0] in C]
		to_extend = [C for C in to_extend if C not in best]
		k=len(order)-1
		while k > i and best :
			worst = [C for C in best if order[k][0] in C]
			tmp.insert(0,copy.deepcopy(worst))
			best = [C for C in best if C not in worst]
			k -= 1
		res = res + [[order[i]]] + tmp
	return res


# given two vectors of different sizes s and s+1, each describing a coalition (A or B) by the position of its elements in the ranking (ordered from either worst to best [leximin] or best to worst [leximax])
# using the Gärdenfors principle, determines whether A > B (returns 1) or B > A (returns 0) or indeterminate (returns -1)
def diff_size_comp(v_A,v_B) :
	if len(v_A) > len(v_B) :
		tmp = [x for x in v_A if x not in v_B]
		# reminder that we have access to the position of an element in the initial ranking: if its value is higher, it means its performance is worse
		if tmp[0] > max(v_B) :
			# print("B preferred to larger A by better lowest elem")
			return 0
		if tmp[0] < min(v_B) :
			# print("larger A preferred to B by better highest elem")
			return 1
		return -1
	tmp = [x for x in v_B if x not in v_A]
	if tmp[0] > max(v_A) :
		# print("A preferred to larger B by better lowest elem")
		return 1
	if tmp[0] < min(v_A) :
		# print("larger B preferred to A by better highest elem")
		return 0
	return -1


# given an order over elements of a population pop
# returns the order over the powerset of pop using leximin
def leximin(order, pop) :
	res = []
	to_extend = tools.generate_powerset(pop,0)
	for i in range(len(to_extend)) :
		v_A = [order.index((x,)) for x in to_extend[i]]
		v_A.sort(reverse=True)
		for j in range(i+1,len(to_extend)) :
			# print("comparing "+str(to_extend[i])+" and "+str(to_extend[j]))
			v_B = [order.index((x,)) for x in to_extend[j]]
			v_B.sort(reverse=True)
			# by definition of leximin, the emptyset is the worst possible set
			if len(v_A) == 0 :
				res += [[to_extend[j],to_extend[i]]]
				continue
			# if one is the superset of the other with only one more element
			if abs(len(v_A)-len(v_B))== 1 and (len([x for x in v_A if x not in v_B])==1 or len([x for x in v_B if x not in v_A])==1):
				r = diff_size_comp(v_A,v_B)
				# if A > B
				if r == 1 :
					res += [[to_extend[i],to_extend[j]]]
					continue
				# if B > A
				elif r == 0 :
					res += [[to_extend[j],to_extend[i]]]
					continue
			# lexicographical comparison of v_A and v_B
			rel = tools.lex_comp_min(v_A,v_B)
			# if A > B
			if rel == 1 :
				# print(str(to_extend[i])+" of vector "+str(v_A)+" preferred to "+str(to_extend[j])+" of vector "+str(v_B)+"\n")
				res += [[to_extend[i],to_extend[j]]]
			# if B > A
			elif rel == 0 :
				# print(str(to_extend[j])+" of vector "+str(v_B)+" preferred to "+str(to_extend[i])+" of vector "+str(v_A)+"\n")
				res += [[to_extend[j],to_extend[i]]]
			# we do not deal with equivalences, as only identical coalitions are equivalent
	# print("out of loop at "+str(time.time()))
	res = tools.join_prefs_coal(res, to_extend)
	# print("res finalized at "+str(time.time()))
	return res


# given a list of coalitions (tuples) and the order over individuals
# returns preferences over the given coalitions using recursion
def join_min(coals, order) : 
	if not coals :
		return coals
	v_A = [order.index((x,)) for x in coals[0]]
	v_A.sort(reverse=True)
	better = []
	worse = []
	for i in range(1, len(coals)) :
		v_B = [order.index((x,)) for x in coals[i]]
		v_B.sort(reverse=True)
		# by definition of leximin, the emptyset is the worst possible set
		if len(v_B) == 0 :
			worse.append(coals[i])
			continue
		# lexicographical comparison of v_A and v_B
		rel = tools.lex_comp_min(v_A,v_B)
		# if A > B
		if rel == 1 :
			worse.append(coals[i])
		# if B > A
		elif rel == 0 :
			better.append(coals[i])
		# we do not deal with equivalences, as only identical coalitions are equivalent
	return join_min(better, order)+[[coals[0]]]+join_min(worse, order)


# given an order over elements of a population pop
# returns the order over the powerset of pop using leximin
def leximin2(order, pop) :
	to_extend = tools.generate_powerset(pop,0)
	to_extend = to_extend[::-1]
	return join_min(to_extend, order)


# given an order over elements of a population pop
# returns the order over the powerset of pop using leximax (aka anti-lexcel?)
def leximax(order, pop) :
	res = []
	to_extend = tools.generate_powerset(pop,0)
	for i in range(len(to_extend)) :
		v_A = [order.index((x,)) for x in to_extend[i]]
		v_A.sort()
		for j in range(i+1,len(to_extend)) :
			v_B = [order.index((x,)) for x in to_extend[j]]
			v_B.sort()
			# by definition of leximax, the emptyset is the best possible set
			if len(v_A) == 0 :
				res += [[to_extend[i],to_extend[j]]]
				continue
			# if one is the superset of the other with only one more element
			if abs(len(v_A)-len(v_B))== 1 and (len([x for x in v_A if x not in v_B])==1 or len([x for x in v_B if x not in v_A])==1):
				r = diff_size_comp(v_A,v_B)
				# if A > B
				if r == 1 :
					res += [[to_extend[i],to_extend[j]]]
					continue
				# if B > A
				elif r == 0 :
					res += [[to_extend[j],to_extend[i]]]
					continue
			# lexicographical comparison of v_A and v_B
			rel = tools.lex_comp_max(v_A,v_B)
			# if A > B
			if rel == 1 :
				# print(str(to_extend[i])+" of vector "+str(v_A)+" preferred to "+str(to_extend[j])+" of vector "+str(v_B)+"\n")
				res += [[to_extend[i],to_extend[j]]]
			# if B > A
			elif rel == 0 :
				# print(str(to_extend[j])+" of vector "+str(v_B)+" preferred to "+str(to_extend[i])+" of vector "+str(v_A)+"\n")
				res += [[to_extend[j],to_extend[i]]]
	res = tools.join_prefs_coal(res, to_extend)
	return res


# given a list of coalitions (tuples) and the order over individuals
# returns preferences over the given coalitions using recursion
def join_max(coals, order) : 
	if not coals :
		return coals
	v_A = [order.index((x,)) for x in coals[0]]
	v_A.sort()
	better = []
	worse = []
	for i in range(1, len(coals)) :
		v_B = [order.index((x,)) for x in coals[i]]
		v_B.sort()
		# by definition of leximax, the emptyset is the best possible set
		if len(v_B) == 0 :
			better.append(coals[i])
			continue
		# lexicographical comparison of v_A and v_B
		rel = tools.lex_comp_max(v_A,v_B)
		# if A > B
		if rel == 1 :
			worse.append(coals[i])
		# if B > A
		elif rel == 0 :
			better.append(coals[i])
		# we do not deal with equivalences, as only identical coalitions are equivalent
	return join_max(better, order)+[[coals[0]]]+join_max(worse, order)


# given an order over elements of a population pop
# returns the order over the powerset of pop using leximax
def leximax2(order, pop) :
	to_extend = tools.generate_powerset(pop,0)
	to_extend = to_extend[::-1]
	return join_max(to_extend, order)


# given an order over elements of a population pop
# returns the order over the powerset of pop using indirect-utility ranking
def indirect_utility(order, pop) :
	res = []
	to_extend = tools.generate_powerset(pop,1)
	for i in range(len(order)) :
		best = [C for C in to_extend if order[i][0] in C]
		to_extend = [C for C in to_extend if C not in best]
		res.append(best)
	return res


# given an order over elements of a population pop
# returns the order over the powerset of pop using cardinality-based ordering
def cardinality(order, pop) :
	res = []
	for i in range(len(pop)+1, 0, -1) :
		res.append(list(itertools.combinations([i+1 for i in range(N)],i)))
	return res


# given an order over elements of a population pop
# returns the order over the powerset of pop using indirect utility and cardinality as a tie-breaker 
def indirect_u_card(order, pop) :
	res = indirect_utility(order, pop)
	res2 = []
	for i in range(len(res)) :
		new = []
		for j in range(len(pop), -1, -1) :
			tmp = [C for C in res[i] if len(C)==j]
			if tmp :
				new.append(tmp)
		res2 += new
	return res2


# given an order over elements of a population pop
# returns the order over the powerset of pop using cardinality and indirect utility as a tie-breaker
def card_ind_u(order, pop) :
	res = cardinality(order, pop)
	res2 = []
	for eq in res :
		tmp = [el for el in eq]
		new = []
		for i in range(len(order)) :
			best = [C for C in tmp if order[i][0] in C]
			tmp = [C for C in tmp if C not in best]
			new.append(best)
		res2 += new
	return res2



# given an order over elements of a population pop
# computes ranking over the powerset of pop using borda score of its components
def borda_like(order, pop) :
	to_extend = tools.generate_powerset(pop,0)
	score = [0 for _ in to_extend]
	for i in range(len(to_extend)) :
		for el in to_extend[i] :
			score[i] += len(order)-order.index((el,))
	# print(score)
	res = []
	while to_extend :
		tmp = [x for x in to_extend if score[to_extend.index(x)]==max(score)]
		res.append(tmp)
		to_extend=[x for x in to_extend if x not in tmp]
		score = [s for s in score if s != max(score)]
	# print(res)
	return res