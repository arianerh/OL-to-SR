# -*- coding: utf-8 -*-

import tools

# given an order (list of equivalence classes) of prefs over coalitions and a population N (list of individuals)
# returns the lexcel ranking over individuals
def lexcel(prefs, N) :
	# prefs = [[(2,),(1,),(2,3),(1,2)],[(1,3),(1,2,3),(3,)]]
	res = []
	# cur = [i+1 for i in range(N)]
	cur = [i for i in N]
	cpt = 0
	seen = []

	# while cpt in range(len(prefs)) and len(res) != N :
	while cpt in range(len(prefs)) and len(res) != len(N) :
		score = [0 for _ in cur]
		# print(cpt, prefs[cpt])
		for c in prefs[cpt] :
			# print(c)
			for el in c :
				if el in cur :
					score[cur.index(el)] += 1
		# print("score "+str(score)+" for "+str(cur))
		cur = [el for el in cur if score[cur.index(el)]==max(score)]
		if len(cur) == 1 :
			# print("cur is "+str(cur)+": appending")
			res.append(cur)
			seen += cur
		elif len(cur) < len(N) :
			# print(str(cur)+" of length "+str(len(cur)))
			tmp = lexcel(prefs,cur)
			res += tmp
			seen += [i for i in cur]
			cpt += 1
		else :
			cpt += 1
		# cur = [i+1 for i in range(N) if i+1 not in res]
		cur = [i for i in N if i not in seen]
		# print("updated cur is now "+str(cur))
		#cpt += 1

	# print(cur)
	if cur :
		res.append(cur)
	return res


# given an order (list of equivalence classes) prefs over coalitions and a population N (list of individuals)
# returns a list of pairwise preferences over individuals obtained using CP-majority
# CAREFUL: CP-majority leads to cycles
def CPmaj(prefs, N) :
	# print("prefs are : "+str(prefs))
	res = []
	for i in range(len(N)) :
		for j in range(i+1, len(N)) :
			# retrieve coalitions containing neither N[i] nor N[j]
			t = tools.generate_powerset([e for e in N if e!=N[i] and e!=N[j]],0)
			pro_i = 0
			pro_j = 0
			for c in t :
				c1 = list(c) + [N[i]]
				c1.sort()
				c1 = tuple(c1)
				c2 = list(c) + [N[j]]
				c2.sort()
				c2 = tuple(c2)
				# print("comparing "+str(c1)+" and "+tr(c2))
				if not tools.is_present(c1,prefs) or not tools.is_present(c2,prefs) :
					continue
				before = tools.is_before(c1, c2,prefs) 
				if before == True :
					pro_i += 1
					# print(str(c1)+" in favour of "+str(N[i])+" because preferred to "+str(c2))
				elif before != -1 :
					pro_j += 1
					# print(str(c2)+" in favour of "+str(N[j])+" because preferred to "+str(c1))
			# print(str(pro_i)+" in favour of "+str(N[i])+" and "+str(pro_j)+" for "+str(N[j]))
			if pro_i >= pro_j :
				res.append([N[i],N[j]])
			if pro_j >= pro_i :
				res.append([N[j],N[i]])
	return res


# given an order (list of equivalence classes) prefs over coalitions and a population N (list of individuals)
# returns ranking over individuals using ordinal Banzhaf
def ordinal_banzhaf(prefs, N) :
	res = []
	score = [0 for _ in N]
	for el in N :
		t = tools.generate_powerset([e for e in N if e!=el],0)
		pro_el = 0
		neg_el = 0
		for c in t :
			# print("considering "+str(c))
			c2 = list(c)+[el]
			c2.sort()
			c2 = tuple(c2)
			# print("checking "+str(c2))
			before = tools.is_before(c2, c, prefs)
			# print("before is "+str(before))
			if before == True :
				# print(str(c2)+" preferred to "+str(c))
				pro_el += 1
			elif before != -1 :
				# print(str(c)+" preferred to "+str(c2))
				neg_el += 1
		# print("score of "+str(el)+" is "+str(pro_el - neg_el))
		score[N.index(el)] = pro_el-neg_el
	# print("scores are "+str(score)+"\n")
	N2 = [e for e in N]
	while N2 :
		tmp = [x for x in N2 if score[N2.index(x)]==max(score)]
		res.append(tmp)
		score = [x for x in score if x!=max(score)]
		N2 = [x for x in N2 if x not in tmp]
	return res 