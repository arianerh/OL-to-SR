# -*- coding: utf-8 -*-

import random
import math
import os
import csv
import statistics

import matplotlib.pyplot as plt
import numpy as np

import OL
import tools
import SR


# given an order over a population X and a rule nb
# returns the extension obtained using the specified rule
def get_order_from_rule(order, X, rule) :
	if rule == 1 :
		return OL.minmax(order, X)
	if rule == 2 :
		return OL.maxmin(order, X)
	if rule == 3 :
		# return OL.leximin(order, X)
		return OL.leximin2(order, X)
	if rule == 4 :
		# return OL.leximax(order, X)
		return OL.leximax2(order, X)
	if rule == 5 :
		return OL.borda_like(order, X)
	if rule == 6 :
		return OL.indirect_utility(order, X)
	if rule == 7 :
		return OL.cardinality(order, X)
	if rule == 8 :
		return OL.indirect_u_card(order, X)
	return OL.card_ind_u(order, X)


# given the rule number
# returns the name of the rule (empty string if incorrect rule number)
def get_rule_name(rule) :
	if rule == 1 :
		return "minmax"
	if rule == 2 :
		return "maxmin"
	if rule == 3 :
		return "leximin"
	if rule == 4 :
		return "leximax"
	if rule == 5 :
		return "bordasum"
	return ""


# given the required size of the output order, the complete order, the population X and the number of the rule used
# returns a partial preorder compatible with order of size s and such that at least two singletons are absent
def get_sized_order(s, order, X, rule) :
	res = get_order_from_rule(order, X, rule)
	nb_c = 2**len(X)
	# empty set not considered for maxmin and minmax (rules 1 and 2)
	if rule < 3 :
		nb_c -= 1
	# in our case, any partial preorder containing all singletons makes the problem trivial: we require that at least two singletons be missing
	for _ in range(2) :
		sing = tuple([random.choice(X)])
		for eq in res :
			if sing in eq :
				eq.remove(sing)
				break
		# remove any eventual empty equivalence class
		res = [eq for eq in res if eq]
	nb_c -= 2
	while nb_c > s :
		i = random.randint(0,len(res)-1)
		# if it's the only element in its equivalence class, we remove the entire equivalence class
		if len(res[i]) == 1 :
			res.pop(i)
		# otherwise we remove just one coalition from the equivalence class
		else :
			res[i].remove(random.choice(res[i]))
		nb_c -= 1
	return res


# for partial info: tests which method returns exact initial preorder with only a percentage of coalitions available
def testinfo_exact(rule, N, nb_runs) :
	# sizes of preorders to test (will be computed as percentages of 2^N)
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Testing recovery of exact order with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size:
			continue
		# if it were to contain all coalitions, then it's the complete information case, therefore not relevant to this study
		if math.ceil(size_g) >= 2**N-2 :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)

			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			if o_lex == [[x[0]] for x in order] :
				lex += 1

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			o_CP = tools.join_prefs_ind(pref_CP,N)
			if o_CP == [[x[0]] for x in order] :
				CP += 1

			# ORDINAL BANZHAF
			# returns list of equivalence classes to express ranking over elements
			o_ban = SR.ordinal_banzhaf(res,X)
			if o_ban == [[x[0]] for x in order] :
				o_b += 1
		val.append([lex, CP, o_b])
		# no need to continue if a given percentage of coalitions is sufficient for every method to recover the order correctly: a larger percentage will lead to the same result
		if lex == CP == o_b == nb_runs :
			break

	fig, ax = plt.subplots()
	b = np.arange(len(lab))
	b1 = [x+0.15 for x in b]
	b2 = [x+0.15 for x in b1]
	plt.bar(b,[val[i][0] for i in range(len(val))], label="LEX-CEL", width=0.15)
	plt.bar(b1,[val[i][1] for i in range(len(val))], label="CP-MAJ", width=0.15)
	plt.bar(b2,[val[i][2] for i in range(len(val))], label="ORD BAN", width=0.15)
	plt.axhline(y=nb_runs,color='red',linestyle='-')
	plt.ylim(0,nb_runs+100)
	plt.ylabel("Number of runs")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([i+0.2 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("Nb of runs for which the correct order is recovered from partial preferences")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/exact') :
		os.mkdir('out/plots/exact')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/exact') :
		os.mkdir('out/data/exact')
	plt.savefig("out/plots/exact/N"+str(N)+"_"+get_rule_name(rule)+".png", dpi=200)
	plt.close(fig)

	with open('out/data/exact/N'+str(N)+'_'+get_rule_name(rule)+'.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Correct results [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: tests which method returns exact winner with only a percentage of coalitions available
def testinfo_win(rule, N, nb_runs) :
	# sizes of preorders to test (will be computed as percentages of 2^N)
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Testing recovery of top item with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size :
			continue
		if math.ceil(size_g) == 2**N :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)
			
			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			if len(o_lex[0]) == 1 and o_lex[0][0] == order[0][0] :
				lex += 1

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			o_CP = tools.join_prefs_ind(pref_CP,N)
			if len(o_CP[0]) == 1 and o_CP[0][0] == order[0][0] :
				CP += 1

			# ORDINAL BANZHAF
			# returns list of equivalence classes to express ranking over elements
			o_ban = SR.ordinal_banzhaf(res,X)
			if len(o_ban[0]) == 1 and o_ban[0][0] == order[0][0] :
				o_b += 1
		val.append([lex, CP, o_b])
		if lex == CP == o_b == nb_runs :
			break

	fig, ax = plt.subplots()
	b = np.arange(len(lab))
	b1 = [x+0.15 for x in b]
	b2 = [x+0.15 for x in b1]
	plt.bar(b,[val[i][0] for i in range(len(val))], label="LEX-CEL", width=0.15)
	plt.bar(b1,[val[i][1] for i in range(len(val))], label="CP-MAJ", width=0.15)
	plt.bar(b2,[val[i][2] for i in range(len(val))], label="ORD BAN", width=0.15)
	plt.axhline(y=nb_runs,color='red',linestyle='-')
	plt.ylim(0,nb_runs+100)
	plt.ylabel("Number of runs")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([i+0.2 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("Nb of runs for which the correct top item is recovered from partial preferences")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/win') :
		os.mkdir('out/plots/win')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/win') :
		os.mkdir('out/data/win')
	plt.savefig("out/plots/win/N"+str(N)+"_"+get_rule_name(rule)+"_win.png", dpi=200)
	plt.close(fig)

	with open('out/data/win/N'+str(N)+'_'+get_rule_name(rule)+'_win.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Correct results [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: measures Kendall-Tau distance between initial preorder & that retrieved by each method with only a percentage of coalitions available
def testinfo_KT(rule, N, nb_runs) :
	# sizes of preorders to test (will be computed as percentages of 2^N)
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Measuring (weak) Kendall-Tau distance between recovered order and ground truth with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size :
			continue
		if math.ceil(size_g) == 2**N :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		storage = [[],[],[]]
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)
			
			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			storage[0].append(tools.Kendall_Tau(o_lex, [[x[0]] for x in order], X))

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			storage[1].append(tools.KT_CP(pref_CP, [[x[0]] for x in order]))

			# ORDINAL BANZHAF
			# returns list of equivalence classes to express ranking over elements
			o_ban = SR.ordinal_banzhaf(res,X)
			storage[2].append(tools.Kendall_Tau(o_ban, [[x[0]] for x in order], X))

		val.append([_ for _ in storage])
		# if we reach the point from which we always get a KT distance of 0, we stop
		if max(storage[0]) == max(storage[1]) == max(storage[2]) == 0 :
			break

	fig, ax = plt.subplots()
	bp1 = ax.boxplot([val[i][0] for i in range(len(lab))],positions=[2.2*i for i in range(len(lab))], widths=0.3, boxprops=dict(color="C1"), showfliers=False)
	bp2 = ax.boxplot([val[i][1] for i in range(len(lab))],positions=[2.2*i+0.4 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C2"), showfliers=False)
	bp3 = ax.boxplot([val[i][2] for i in range(len(lab))],positions=[2.2*i+0.8 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C3"), showfliers=False)
	for median in bp1['medians']:
		median.set_color("C1")
	for median in bp2['medians']:
		median.set_color("C2")
	for median in bp3['medians']:
   		median.set_color("C3")
	ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0]], ['LEX-CEL', 'CP-MAJ', 'ORD BANZ'], loc='upper right')

	plt.ylabel("Kendall-Tau distance")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([2.2*i+0.4 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("(Weak) Kendall-Tau distance between recovered order and ground truth")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/KT') :
		os.mkdir('out/plots/KT')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/KT') :
		os.mkdir('out/data/KT')
	
	plt.savefig("out/plots/KT/N"+str(N)+"_"+get_rule_name(rule)+"_KT.png", dpi=200)
	plt.close(fig)

	with open('out/data/KT/N'+str(N)+'_'+get_rule_name(rule)+'_KT.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Results [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: measures nb of incorrect pairwise prefences between initial preorder & that retrieved by each method with only a percentage of coalitions available
def testinfo_errors(rule, N, nb_runs) :
	# sizes of preorders to test (will be computed as percentages of 2^N)
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Counting number of errors between recovered order and ground truth with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size :
			continue
		if math.ceil(size_g) == 2**N :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		storage = [[],[],[]]
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)
			
			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			storage[0].append(tools.Kendall_Tau(o_lex, [[x[0]] for x in order], X))

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			storage[1].append(tools.count_inverse(pref_CP, [[x[0]] for x in order]))

			# ORDINAL BANZHAF
			# returns list of equivalence classes to express ranking over elements
			o_ban = SR.ordinal_banzhaf(res,X)
			storage[2].append(tools.Kendall_Tau(o_ban, [[x[0]] for x in order], X))

		val.append([_ for _ in storage])
		# if at this stage no method makes any error, then they will not make any with higher thresholds either: we can stop here
		if max(storage[0]) == max(storage[1]) == max(storage[2]) == 0 :
			break


	fig, ax = plt.subplots()
	bp1 = ax.boxplot([val[i][0] for i in range(len(lab))],positions=[2.2*i for i in range(len(lab))], widths=0.3, boxprops=dict(color="C1"), showfliers=False)
	bp2 = ax.boxplot([val[i][1] for i in range(len(lab))],positions=[2.2*i+0.4 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C2"), showfliers=False)
	bp3 = ax.boxplot([val[i][2] for i in range(len(lab))],positions=[2.2*i+0.8 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C3"), showfliers=False)
	for median in bp1['medians']:
		median.set_color("C1")
	for median in bp2['medians']:
		median.set_color("C2")
	for median in bp3['medians']:
   		median.set_color("C3")
	ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0]], ['LEX-CEL', 'CP-MAJ', 'ORD BANZ'], loc='upper right')

	plt.ylabel("Nb of errors")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([2.2*i+0.4 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("Nb of errors between recovered order and ground truth")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	# plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/errors') :
		os.mkdir('out/plots/errors')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/errors') :
		os.mkdir('out/data/errors')
		
	plt.savefig("out/plots/errors/N"+str(N)+"_"+get_rule_name(rule)+"_errors.png", dpi=200)
	plt.close(fig)

	with open('out/data/errors/N'+str(N)+'_'+get_rule_name(rule)+'_errors.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Nb of errors [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: measures (weak) Kendall-Tau distance between prefs recovered by combinations of SR methods and the ground truth
def testinfo_KT_sr_CP(rule, N, nb_runs) :
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Measuring (weak) Kendall-Tau distance for combinations of SR methods with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size :
			continue
		if math.ceil(size_g) == 2**N :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		storage = [[],[],[],[]]
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)
			
			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			storage[0].append(tools.Kendall_Tau(o_lex, [[x[0]] for x in order], X))

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			storage[1].append(tools.KT_CP(pref_CP, [[x[0]] for x in order]))

			# LEXCEL CORRECTED WITH CP
			o_comb_lex = tools.corrected_CP(pref_CP, o_lex)
			storage[2].append(tools.KT_CP(o_comb_lex, [[x[0]] for x in order]))

			# ORDINAL BANZHAF CORRECTED WITH CP
			o_ban = SR.ordinal_banzhaf(res,X)
			o_comb_ob = tools.corrected_CP(pref_CP, o_ban)
			storage[3].append(tools.KT_CP(o_comb_ob, [[x[0]] for x in order]))

		val.append([_ for _ in storage])
		# If all methods get a Kendall-Tau distance of 0, no need to explore further
		if max(storage[0]) == max(storage[1]) == max(storage[2]) == 0 :
			break

	fig, ax = plt.subplots()
	bp1 = ax.boxplot([val[i][0] for i in range(len(lab))],positions=[2.2*i for i in range(len(lab))], widths=0.3, boxprops=dict(color="C1"), showfliers=False)
	bp2 = ax.boxplot([val[i][1] for i in range(len(lab))],positions=[2.2*i+0.4 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C2"), showfliers=False)
	bp3 = ax.boxplot([val[i][2] for i in range(len(lab))],positions=[2.2*i+0.8 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C3"), showfliers=False)
	bp4 = ax.boxplot([val[i][3] for i in range(len(lab))],positions=[2.2*i+1.2 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C4"), showfliers=False)
	for median in bp1['medians']:
		median.set_color("C1")
	for median in bp2['medians']:
		median.set_color("C2")
	for median in bp3['medians']:
   		median.set_color("C3")
	for median in bp4['medians']:
   		median.set_color("C4")
	ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0], bp4["boxes"][0]], ['LEX-CEL', 'CP-MAJ', 'CP+LEX', 'CP+BAN'], loc='upper right')

	plt.ylabel("Kendall-Tau distance")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([2.2*i+0.6 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("(Weak) Kendall-Tau distance between recovered order and ground truth")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	# plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/comb') :
		os.mkdir('out/plots/comb')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/comb') :
		os.mkdir('out/data/comb')
		
	plt.savefig("out/plots/comb/N"+str(N)+"_"+get_rule_name(rule)+"_comb_KT.png", dpi=200)
	plt.close(fig)

	with open('out/data/comb/N'+str(N)+'_'+get_rule_name(rule)+'_comb_KT.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Results [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: measures nb of incorrect pairwise preferences between initial preorder & that retrieved by each combination of CP with another SR method, with only a certain percentage of coalitions available
def testinfo_errors_sr_CP(rule, N, nb_runs) :
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Counting number of errors for combinations of methods with rule "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	val = []
	size = 1
	lab = []
	for size_g in [(2**N)*i for i in tests] :
		if size_g <= size :
			continue
		if math.ceil(size_g) == 2**N :
			break
		lab.append(size_g/(2**N))
		size = math.ceil(size_g)
		lex = 0
		CP = 0
		o_b = 0
		storage = [[],[],[],[]]
		for i in range(nb_runs) :
			order = [(x,) for x in X]
			random.shuffle(order)
		
			res = get_sized_order(size, order, X, rule)

			# LEXCEL
			# returns list of equivalence classes to express ranking over elements
			o_lex = SR.lexcel(res, X)
			storage[0].append(tools.Kendall_Tau(o_lex, [[x[0]] for x in order], X))

			# CP-MAJORITY
			# returns list of pairwise preferences (not necessarily transitive)
			pref_CP = SR.CPmaj(res, X)
			storage[1].append(tools.count_inverse(pref_CP, [[x[0]] for x in order]))

			# LEXCEL CORRECTED WITH CP
			o_comb_lex = tools.corrected_CP(pref_CP, o_lex)
			storage[2].append(tools.count_inverse(o_comb_lex, [[x[0]] for x in order]))

			# ORDINAL BANZHAF CORRECTED WITH CP
			o_ban = SR.ordinal_banzhaf(res,X)
			o_comb_ob = tools.corrected_CP(pref_CP, o_ban)
			storage[3].append(tools.count_inverse(o_comb_ob, [[x[0]] for x in order]))

		val.append([_ for _ in storage])
		# if no method makes any error at this point, we can stop
		if max(storage[0]) == max(storage[1]) == max(storage[2]) == 0 :
			break

	fig, ax = plt.subplots()
	bp1 = ax.boxplot([val[i][0] for i in range(len(lab))],positions=[2.2*i for i in range(len(lab))], widths=0.3, boxprops=dict(color="C1"), showfliers=False)
	bp2 = ax.boxplot([val[i][1] for i in range(len(lab))],positions=[2.2*i+0.4 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C2"), showfliers=False)
	bp3 = ax.boxplot([val[i][2] for i in range(len(lab))],positions=[2.2*i+0.8 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C3"), showfliers=False)
	bp4 = ax.boxplot([val[i][3] for i in range(len(lab))],positions=[2.2*i+1.2 for i in range(len(lab))], widths=0.3, boxprops=dict(color="C4"), showfliers=False)
	for median in bp1['medians']:
		median.set_color("C1")
	for median in bp2['medians']:
		median.set_color("C2")
	for median in bp3['medians']:
   		median.set_color("C3")
	for median in bp4['medians']:
   		median.set_color("C4")
	ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0], bp4["boxes"][0]], ['LEX-CEL', 'CP-MAJ', 'CP+LEX', 'CP+BAN'], loc='upper right')

	plt.ylabel("Nb of errors")
	plt.xlabel("Percentage of total coalitions present in the preorder")
	plt.xticks([2.2*i+0.6 for i in range(len(lab))], labels=[int(i*100) for i in lab])
	plt.suptitle("Nb of errors between recovered order and ground truth")
	plt.title( "(N = "+str(N)+", "+str(nb_runs)+" runs, rule = "+get_rule_name(rule)+")")
	# plt.legend()
	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/plots') :
		os.mkdir('out/plots')
	if not os.path.exists('out/plots/comb/') :
		os.mkdir('out/plots/comb')
	if not os.path.exists('out/data') :
		os.mkdir('outy/data')
	if not os.path.exists('out/data/comb') :
		os.mkdir('out/data/comb')
		
	plt.savefig("out/plots/comb/N"+str(N)+"_"+get_rule_name(rule)+"_comb_errors.png", dpi=200)
	plt.close(fig)

	with open('out/data/comb/N'+str(N)+'_'+get_rule_name(rule)+'_comb_errors.csv', 'w') as f:
		pen = csv.writer(f, delimiter=';')
		pen.writerow(['Percentage of total coalitions', 'Nb of errors [Lexcel, CP-maj, Banzhaf]'])
		pen.writerows([[lab[i],val[i]] for i in range(len(lab))])


# for partial info: tests which method returns exact initial preorder with only a percentage of coalitions of size k available
def k_info_exact(rule, N, nb_runs) :
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Testing recovery of exact order from same-sized coalitions with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	stop = N
	# we know that rules 1 and 2 cannot find the correct order when k > n/2, so no need to go further
	if rule < 3 :
		stop = math.ceil(N/2)
	storage = [[],[],[],[]]
	for k in range(2,stop+1) :
		tab_l = []
		tab_cp = []
		tab_comb = []
		size = 0
		# note which percentage we've actually tested 
		pc = []
		nbc = math.comb(N,k)
		for size_g in [nbc*i for i in tests] :
			if size_g <= size or size_g < 2 or math.ceil(size_g) == nbc:
				continue
			pc.append(size_g/nbc)
			size = math.ceil(size_g)
			lex = 0
			CP = 0
			comb = 0
			for i in range(nb_runs) :
				order = [(x,) for x in X]
				random.shuffle(order)
				res = get_order_from_rule(order, X, rule)
				
				res2 = []
				for eq in res :
					tmp = [el for el in eq if len(el)==k]
					if tmp :
						res2.append(tmp)
				while len(res2) > size :
					res2.pop(random.randint(0,len(res2)-1))

				# LEXCEL
				# returns list of equivalence classes to express ranking over elements
				o_lex = SR.lexcel(res2, X)
				if o_lex == [[x[0]] for x in order] :
					lex += 1

				# CP-MAJORITY
				# returns list of pairwise preferences (not necessarily transitive)
				pref_CP = SR.CPmaj(res2, X)
				o_CP = tools.join_prefs_ind(pref_CP,N)
				if o_CP == [[x[0]] for x in order] :
					CP += 1

				# Hybrid CP + lexcel
				pref_hybrid = tools.corrected_CP(pref_CP, o_lex)
				o_comb_lex = tools.join_prefs_ind(pref_hybrid,N)
				if o_comb_lex == [[x[0]] for x in order] :
					comb += 1
			tab_l.append(lex)
			tab_cp.append(CP)
			tab_comb.append(comb)

		# we also keep track of which percentage of information we've actually tested (no redundancy in the final size of population)
		storage[0].append(pc)
		storage[1].append(tab_l)
		storage[2].append(tab_cp)
		storage[3].append(tab_comb)

	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/k-sized') :
		os.mkdir('out/data/k-sized')
	if not os.path.exists('out/data/k-sized/exact') :
		os.mkdir('out/data/k-sized/exact')

	with open('out/data/k-sized/exact/N'+str(N)+'_'+get_rule_name(rule)+'.csv', 'w') as f:
		pen=csv.writer(f, delimiter=';')
		pen.writerow(['k']+[str(i+1) for i in range(len(storage[0]))])
		pen.writerow(['Percentages of available coalitions']+storage[0])
		pen.writerow(['Number of correct retrieved order (Lexcel)']+storage[1])
		pen.writerow(['Number of correct retrieved order (CP-maj)']+storage[2])
		pen.writerow(['Number of correct retrieved order (CP+lex)']+storage[3])


# for partial info: tests which method returns exact winner (i.e. not equivalent to another incorrect winner) with only a percentage of coalitions of size k available
def k_info_win(rule, N, nb_runs) :
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Testing recovery of correct top item from same-sized coalitions with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	storage = [[],[],[],[]]
	for k in range(2,N) :
		tab_l = []
		tab_cp = []
		tab_comb = []
		size = 0
		# note which percentage we've actually tested 
		pc = []
		nbc = math.comb(N,k)
		for size_g in [nbc*i for i in tests] :
			if size_g <= size or size_g < 2 or math.ceil(size_g) == nbc:
				continue
			pc.append(size_g/nbc)
			size = math.ceil(size_g)
			lex = 0
			CP = 0
			comb = 0
			for i in range(nb_runs) :
				order = [(x,) for x in X]
				random.shuffle(order)
				res = get_order_from_rule(order, X, rule)
				
				res2 = []
				for eq in res :
					tmp = [el for el in eq if len(el)==k]
					if tmp :
						res2.append(tmp)
				while len(res2) > size :
					res2.pop(random.randint(0,len(res2)-1))

				# LEXCEL
				# returns list of equivalence classes to express ranking over elements
				o_lex = SR.lexcel(res2, X)
				# we consider that the result is correct only if it returns the correct winner as the UNIQUE correct winner
				if len(o_lex[0])==1 and o_lex[0][0] == order[0][0] :
					lex += 1

				# CP-MAJORITY
				# returns list of pairwise preferences (not necessarily transitive)
				pref_CP = SR.CPmaj(res2, X)
				o_CP = tools.join_prefs_ind(pref_CP,N)
				# we consider that the result is correct only if it returns the correct winner as the UNIQUE correct winner
				if len(o_CP[0]) == 1 and o_CP[0][0] == order[0][0] :
					CP += 1

				# Hybrid CP + lexcel
				pref_hybrid = tools.corrected_CP(pref_CP, o_lex)
				o_comb_lex = tools.join_prefs_ind(pref_hybrid,N)
				if len(o_comb_lex[0]) == 1 and o_comb_lex[0][0] == order[0][0] :
					comb += 1
			tab_l.append(lex)
			tab_cp.append(CP)
			tab_comb.append(comb)

		# we also keep track of which percentage of information we've actually tested (no redundancy in the final size of population)
		storage[0].append(pc)
		storage[1].append(tab_l)
		storage[2].append(tab_cp)
		storage[3].append(tab_comb)

	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/k-sized') :
		os.mkdir('out/data/k-sized')
	if not os.path.exists('out/data/k-sized/win') :
		os.mkdir('out/data/k-sized/win')

	with open('out/data/k-sized/win/N'+str(N)+'_'+get_rule_name(rule)+'_win.csv', 'w') as f:
		pen=csv.writer(f, delimiter=';')
		pen.writerow(['k']+[str(i+1) for i in range(len(storage[0]))])
		pen.writerow(['Percentages of available coalitions']+storage[0])
		pen.writerow(['Number of correct retrieved order (Lexcel)']+storage[1])
		pen.writerow(['Number of correct retrieved order (CP-maj)']+storage[2])
		pen.writerow(['Number of correct retrieved order (CP+lex)']+storage[3])


# for partial info: measures Kendall-Tau distance between initial preorder & that retrieved by each method with only a percentage of coalitions of size k available
def k_info_KT(rule, N, nb_runs) :
	tests = [0.05] + [i/10 for i in range(1,10)]
	tests.append(0.99)
	print("Measuring (weak) Kendall-Tau distance between order recovered and ground truth over same-sized coalitions with "+get_rule_name(rule))
	X = [i+1 for i in range(N)]
	storage = [[],[],[],[]]
	for k in range(2,N) :
		tab_l = []
		tab_cp = []
		tab_comb = []
		size = 0
		# note which percentage we've actually tested 
		pc = []
		nbc = math.comb(N,k)
		for size_g in [nbc*i for i in tests] :
			if size_g <= size or size_g < 2 or math.ceil(size_g) == nbc:
				continue
			pc.append(size_g/nbc)
			size = math.ceil(size_g)
			a1 = []
			a2 = []
			a3 = []
			for i in range(nb_runs) :
				order = [(x,) for x in X]
				random.shuffle(order)
				res = get_order_from_rule(order, X, rule)
				
				res2 = []
				for eq in res :
					tmp = [el for el in eq if len(el)==k]
					if tmp :
						res2.append(tmp)
				while len(res2) > size :
					res2.pop(random.randint(0,len(res2)-1))

				# LEXCEL
				# returns list of equivalence classes to express ranking over elements
				o_lex = SR.lexcel(res2, X)
				a1.append(tools.Kendall_Tau(o_lex, [[x[0]] for x in order], X))

				# CP-MAJORITY
				# returns list of pairwise preferences (not necessarily transitive)
				pref_CP = SR.CPmaj(res2, X)
				o_CP = tools.join_prefs_ind(pref_CP,N)
				a2.append(tools.KT_CP(pref_CP, [[x[0]] for x in order]))
				
				# Hybrid CP + lexcel
				o_comb_lex = tools.corrected_CP(pref_CP, o_lex)
				a3.append(tools.KT_CP(o_comb_lex, [[x[0]] for x in order]))

			tab_l.append([min(a1),statistics.median(a1),statistics.mean(a1),max(a1)])
			tab_cp.append([min(a2),statistics.median(a2),statistics.mean(a2),max(a2)])
			tab_comb.append([min(a3),statistics.median(a3),statistics.mean(a3),max(a3)])

		# we also keep track of which percentage of information we've actually tested (no redundancy in the final size of population)
		storage[0].append(pc)
		storage[1].append(tab_l)
		storage[2].append(tab_cp)
		storage[3].append(tab_comb)

	if not os.path.exists('out') :
		os.mkdir('out')
	if not os.path.exists('out/data') :
		os.mkdir('out/data')
	if not os.path.exists('out/data/k-sized') :
		os.mkdir('out/data/k-sized')
	if not os.path.exists('out/data/k-sized/KT') :
		os.mkdir('out/data/k-sized/KT')

	with open('out/data/k-sized/KT/N'+str(N)+'_'+get_rule_name(rule)+'_KT.csv', 'w') as f:
		pen=csv.writer(f, delimiter=';')
		pen.writerow(['k']+[str(i+1) for i in range(len(storage[0]))])
		pen.writerow(['Percentages of available coalitions']+storage[0])
		pen.writerow(['Results (Lexcel) [min,median, avg,max]']+storage[1])
		pen.writerow(['Results (CP-maj) [min,median, avg,max]']+storage[2])
		pen.writerow(['Results (CP+lex) [min,median, avg,max]']+storage[3])


def start():
	try:
		l = input("By default, 10 000 tests will be ran for each population size from 4 to 9. All measures presented in the paper will be launched and the results saved. Do you wish to change these settings? (Y/N) ")
	except ValueError as ve:
		print("Default settings used")
	if l.lower()=="y" :
		try :
			nb_runs = int(input("Enter number of tests to be ran for each population size: "))
		except ValueError as ve:
			print("Error in value, set by defaut to 10 000")
			nb_runs = 10000
		try :
			min_N = int(input("Enter minimum size of population: "))
		except ValueError as ve:
			print("Error in value, set by default to 4")
			min_N = 4
		try :
			max_N = int(input("Enter maximal size of population: "))
		except ValueError as ve:
			print("Error in value, set by default to 9")
			max_N = 9
		if max_N < min_N:
			print("Max size must be larger than min size: set by default to "+str(min_N+1))
			max_N = min_N + 1
	else :	
		min_N = 4
		max_N = 5
		nb_runs = 10

	unique = False
	try :
		l = input("By default, all 5 order lifting rules will be tested, for each population size. Do you wish to change these settings? (Y/N) ")
	except ValueError as ve:
		print("Default settings used")
	if l.lower()=="y" :
		try :
			rule = int("Select rule number: 1 = minmax\t 2 = maxmin\t 3 = leximin\t 4 = leximax\t 5 = Borda-sum")
		except ValueError as ve :
			print("Error in value: default settings used")
		if rule > 5 or rule < 0 :
			print("Error in value: default settings used")
		else :
			unique = True

	k = False
	try :
		l = input("By default, test will only be ran for preorders over coalitions of all sizes. Do you wish to run tests on k-sized coalitions only instead (no plots)? (Y/N) ")
	except ValueError as ve:
		print("Default settings used")
	if l.lower=="y" :
		k = True

		
	for a in range(min_N,max_N+1):
		if not unique :
			rule = 1
		print("----\nN is "+str(a))
		while rule < 6 :	
			if k :
				k_info_exact(rule, a, nb_runs)
				k_info_win(rule, a, nb_runs)
				k_info_KT(rule, a, nb_runs)
			else :
				testinfo_exact(rule, a, nb_runs)
				testinfo_win(rule, a, nb_runs)
				testinfo_KT(rule, a, nb_runs)
				testinfo_errors(rule, a, nb_runs)
				testinfo_KT_sr_CP(rule, a, nb_runs)
				testinfo_errors_sr_CP(rule, a, nb_runs)
			if unique :
				break
			else :
				rule += 1