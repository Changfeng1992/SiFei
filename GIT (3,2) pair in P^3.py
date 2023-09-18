#!/usr/bin/env python
# coding: utf-8

# In[22]:


# from typing import Counter, DefaultDict
# from math import gcd
# from functools import reduce
# from sympy import *
from itertools import *


# In[4]:


# def constrained_partitions(n, k, min_elem, max_elem):
#     allowed = range(max_elem, min_elem-1, -1)

#     def helper(n, k, t):
#         if k == 0:
#             if n == 0:
#                 yield t
#         elif k == 1:
#             if n in allowed:
#                 yield t + (n,)
#         elif min_elem * k <= n <= max_elem * k:
#             for v in allowed:
#                 yield from helper(n - v, k - 1, t + (v,))

#     return helper(n, k, ())


# In[23]:


# n is degree & k is the number of variables, i.e. dimension of projective space+1
from itertools import combinations

def weakcomps(n: int, k: int):
    """
    Computes weak compositions of n with k elements
    This algorithm uses the following fact:
    The difference between two vectors of size k, whose elements are an ascending k combination of [0, ..., n],
    that is shifted by 1 element, i.e.
      [e1, e2, e3]
      -   [e1, e2, e3] 
    = [e1,  ...  ,-e3]
    will yield a vector of k + 1 elements whose first k elements are a composition the negative of the last element.
    The sums of the t sections of the u, v vectors cancel out,
    so the sum of the elements in (u - v) = m + 1 = n + k.
    Subtracting one per each of the k elements yields a total of n.
    Because `combinations` returns tuples where the elements (i, ..., k-1) are in increasing order,
    we know that each element is positive.
    """
    m: int = n + k - 1

    for t in combinations(range(m), k - 1):
        u = t + (m,)
        v = (-1,) + t
        yield [u[i] - v[i] - 1 for i in range(k)]


# In[24]:


# base degree 3 in P^3
base_mon = list(weakcomps(3,4))
mon = []
for i in base_mon:
    mon.append(tuple(i))

base_mon = mon
base_mon


# In[25]:


fiber_mon = list(weakcomps(2,4))
mon = []
for i in fiber_mon:
    mon.append(tuple(i))

fiber_mon = mon
fiber_mon


# In[26]:


len(fiber_mon)


# In[27]:


def weight(monomial, LP_instance):
    wt = 0
    for i in range(0,len(monomial)):
        wt = wt + monomial[i] * LP_instance[i]
    return(wt)


# In[9]:


time_ops =set([])
potential_wall22 = set([])
potential_wall31 = set([])


# In[10]:


# mon_pair type (3,1) (m1,m2,m3,n) wt(m1)=wt(m2)=wt(m3)=1

from sage.numerical.mip import MIPSolverException
def critical_ops(monomial,monomial2 ,monomial3):
    ops=[]
#     e= MIPSolverException("GLPK: Problem has no feasible solution")
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    w = p.new_variable(integer=True)
    p.add_constraint(weight(monomial,w) == weight(monomial3,w))
    p.add_constraint(weight(monomial2,w) == weight(monomial3,w))
    for i in range(0,len(monomial)-1):
        p.add_constraint(w[i]>= w[i+1])
    s = sum(w[i] for i in range(0,len(monomial)))
    p.add_constraint(s==0)
#     p.add_constraint(w[0]>= w[1])
#     p.add_constraint(w[1]>= w[2])
#     p.add_constraint(w[2]>= w[3])
#     p.add_constraint(w[3]>= w[4])
#     p.add_constraint(w[0]+w[1]+w[2]+w[3]+w[4]==0)
    p.add_constraint(w[0]>=1)
    p.set_objective(None)
    try:
        p.solve()      
    except MIPSolverException:
        return []
    else:
#         print(p.show())
        x = sorted(p.get_values(w, convert=ZZ, tolerance=1e-3).items())
        for i in x:
            ops.append(i[1])
        return(ops)


# In[11]:


for i in range(0,len(base_mon)):
    for j in range(i+1,len(base_mon)):
        for k in range(j+1,len(base_mon)):
            if critical_ops(base_mon[i],base_mon[j],base_mon[k])!=[]:
                temp_ops = critical_ops(base_mon[i],base_mon[j],base_mon[k])
                for l in fiber_mon:
                    if weight(l,temp_ops) !=0:
                        potential_t = -weight(base_mon[i],temp_ops)/weight(l,temp_ops)
                        if 0 < potential_t <1:
                            potential_wall31.add(potential_t)
                            time_ops.add((potential_t, tuple(temp_ops),base_mon[i],base_mon[j],base_mon[k],l))
                        
potential_wall31                   


# In[12]:


# mon_pair type (2,2)
from sage.numerical.mip import MIPSolverException
def critical_ops22(monomial,monomial2 ,monomial3,monomial4):
    ops=[]
#     e= MIPSolverException("GLPK: Problem has no feasible solution")
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    w = p.new_variable(integer=True)
    p.add_constraint(weight(monomial,w) == weight(monomial2,w))
    p.add_constraint(weight(monomial3,w) == weight(monomial4,w))
    for i in range(0,len(monomial)-1):
        p.add_constraint(w[i]>= w[i+1])
    s = sum(w[i] for i in range(0,len(monomial)))
    p.add_constraint(s==0)
#     p.add_constraint(w[0]>= w[1])
#     p.add_constraint(w[1]>= w[2])
#     p.add_constraint(w[2]>= w[3])
#     p.add_constraint(w[3]>= w[4])
#     p.add_constraint(w[0]+w[1]+w[2]+w[3]+w[4]==0)
    p.add_constraint(w[0]>=1)
    p.set_objective(None)
    # p.show()
    try:
        p.solve()      
    except MIPSolverException:
        return []
    else:
        x = sorted(p.get_values(w, convert=ZZ, tolerance=1e-3).items())
        for i in x:
            ops.append(i[1])
        return(ops)


# In[13]:


for i in range(0,len(base_mon)):
    for j in range(i+1,len(base_mon)):
        for k in range(0,len(fiber_mon)):
            for l in range(k+1,len(fiber_mon)):
                if critical_ops22(base_mon[i],base_mon[j],fiber_mon[k],fiber_mon[l])!=[]:
                    temp_ops = critical_ops22(base_mon[i],base_mon[j],fiber_mon[k],fiber_mon[l])
                    if weight(fiber_mon[l],temp_ops) !=0:
                        potential_t = -1/weight(fiber_mon[l],temp_ops)
                        if 0 < potential_t <1:
                            potential_wall22.add(potential_t)
                            time_ops.add((potential_t, tuple(temp_ops),base_mon[i],base_mon[j],fiber_mon[k],fiber_mon[l]))
                        
potential_wall22             


# In[ ]:


time_ops =set([])


# In[14]:


time_ops


# In[28]:


def func_check(vect_mon):
    t = True
    for i in range(0,len(vect_mon)):
        if sum(vect_mon[0:i+1]) < 0:
            t = False
            return t
            break
    return t

def sign_conversion(vect_mon):
    for i in range(0,len(vect_mon)):
        vect_mon[i] = - vect_mon[i]
    return vect_mon

def sign_confirm(vect_mon):
    for i in range(0,len(vect_mon)):
        if vect_mon[i] != 0:
            if vect_mon[i] >0:
                return [True, i]
            else:
                return [False, i]
            break
        else:
            if i == len(vect_mon) -1:
                return ["equal", i]
            continue 

def dot(K, L):
    if len(K) != len(L):
        return 0

    return sum(z[0] * z[1] for z in zip(K, L))


# In[29]:


# This only give mon_1 < mon_2 or not!
def mon_comparison(mon_1, mon_2):
    mon_12 = []
    for i in range(0,len(mon_1)):
        mon_12.append(mon_1[i] - mon_2[i])
    if sign_confirm(mon_12)[0] == True:
        if func_check(mon_12) == True:
            return True
    return False


# In[30]:


baseMonPoset = Poset([base_mon,mon_comparison],facade=True)
baseMonPoset.plot(figsize=10)


# In[ ]:


baseMon


# In[31]:


fiberMonPoset = Poset([fiber_mon,mon_comparison],facade=True)
fiberMonPoset.plot(figsize=10)


# In[32]:


# return true if max(Set_1) >= max(Set_2)
def mon_set_comparison(Set_1,Set_2,Pos):
    Pos_1 = Pos.subposet(Set_1)
    Pos_2 = Pos.subposet(Set_2)
    result = True
    for x in Pos_2.minimal_elements():
        if all(Pos.is_lequal(y,x) ==False for y in Pos_1.minimal_elements()):
            result = False
            break
    return result   

def mon_set_pair_comparison(Set_pair_1,Set_pair_2,Pos_1,Pos_2):
    return (mon_set_comparison(Set_pair_1[0],Set_pair_2[0],Pos_1) and mon_set_comparison(Set_pair_1[1],Set_pair_2[1],Pos_2))


# In[33]:


from sage.numerical.mip import MIPSolverException
def weight_total(monomial, LP_instance):
    wt_1 = weight(monomial[0],LP_instance)
    wt_2 = weight(monomial[1],LP_instance)
    wt = wt_1 + (1/20) * wt_2
    return(wt)


# In[16]:




# poset version

def negative_mon_list(one_ps,Pos):
    M_set = []
    for i in Pos:
        if weight_total(i, one_ps)  < 0:
            M_set.append(i)
    negPos = Pos.subposet(M_set)
    return negPos

def nonnegative_mon_list(one_ps,Pos):
    M_set = []
    for i in Pos:
        if weight_total(i, one_ps)  >= 0:
            M_set.append(i)
    nonnegPos = Pos.subposet(M_set)
    return nonnegPos

def positive_mon_list(one_ps,Pos):
    M_set = []
    for i in Pos:
        if weight_total(i, one_ps)  > 0:
            M_set.append(i)
    plusPos = Pos.subposet(M_set)
    return plusPos

def nonpositive_mon_list(one_ps,Pos):
    M_set = []
    for i in Pos:
        if weight_total(i, one_ps)  <= 0:
            M_set.append(i)
    nonplusPos = Pos.subposet(M_set)
    return nonplusPos

def find_maxmon_nonpositive(one_ps,Pos):
    return nonpositive_mon_list(one_ps,Pos).minimal_elements()

def find_minmon_positive(one_ps,Pos):
    return positive_mon_list(one_ps,Pos).maximal_elements()


# In[37]:


# mon_pair type (4,2) (m1,m2,m3,m4, n1,n2) 

# for those (m1,m2,m3,m4) such that wt(m1,m2,m3,m4)<=0 no solution
#  we consider to clollect those pair such that total weight (m,n)<=0 有解的， 并从中找出极大的pair组合。

from sage.numerical.mip import MIPSolverException
def ops_solution_base(monomial,monomial2 ,monomial3,monomial4):
    ops=[]
#     e= MIPSolverException("GLPK: Problem has no feasible solution")
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    w = p.new_variable(integer=True)
    p.add_constraint(weight(monomial,w) <= 0)
    p.add_constraint(weight(monomial2,w) <=0)
    p.add_constraint(weight(monomial3,w) <=0)
    p.add_constraint(weight(monomial4,w) <=0)
    for i in range(0,len(monomial)-1):
        p.add_constraint(w[i]>= w[i+1])
    s = sum(w[i] for i in range(0,len(monomial)))
    p.add_constraint(s==0)
#     p.add_constraint(w[0]>= w[1])
#     p.add_constraint(w[1]>= w[2])
#     p.add_constraint(w[2]>= w[3])
#     p.add_constraint(w[3]>= w[4])
#     p.add_constraint(w[0]+w[1]+w[2]+w[3]+w[4]==0)
    p.add_constraint(w[0]>=1)
    p.set_objective(None)
    try:
        p.solve()      
    except MIPSolverException:
        return []
    else:
#         print(p.show())
        x = sorted(p.get_values(w, convert=ZZ, tolerance=1e-3).items())
        for i in x:
            ops.append(i[1])
        return(ops)


# In[47]:


def ops_solution_fiber(monomial,monomial2):
    ops=[]
#     e= MIPSolverException("GLPK: Problem has no feasible solution")
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    w = p.new_variable(integer=True)
    p.add_constraint(weight(monomial,w) <= 0)
    p.add_constraint(weight(monomial2,w) <=0)
    for i in range(0,len(monomial)-1):
        p.add_constraint(w[i]>= w[i+1])
    s = sum(w[i] for i in range(0,len(monomial)))
    p.add_constraint(s==0)
#     p.add_constraint(w[0]>= w[1])
#     p.add_constraint(w[1]>= w[2])
#     p.add_constraint(w[2]>= w[3])
#     p.add_constraint(w[3]>= w[4])
#     p.add_constraint(w[0]+w[1]+w[2]+w[3]+w[4]==0)
    p.add_constraint(w[0]>=1)
    p.set_objective(None)
    try:
        p.solve()      
    except MIPSolverException:
        return []
    else:
#         print(p.show())
        x = sorted(p.get_values(w, convert=ZZ, tolerance=1e-3).items())
        for i in x:
            ops.append(i[1])
        return(ops)


# In[57]:


from sage.numerical.mip import MIPSolverException
def ops_solution_total(base_set,fiber_set):
    ops=[]
#     e= MIPSolverException("GLPK: Problem has no feasible solution")
    p = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    w = p.new_variable(integer=True)
    for x in base_set:
        for y in fiber_set:
            p.add_constraint(weight_total((x,y),w) <= 0)
    for i in range(0,len(base_set[0])-1):
        p.add_constraint(w[i]>= w[i+1])
    s = sum(w[i] for i in range(0,len(base_set[0])))
    p.add_constraint(s==0)
#     p.add_constraint(w[0]>= w[1])
#     p.add_constraint(w[1]>= w[2])
#     p.add_constraint(w[2]>= w[3])
#     p.add_constraint(w[0]+w[1]+w[2]+w[3]==0)
    p.add_constraint(w[0]>=1)
    p.set_objective(None)
    try:
        p.solve()      
    except MIPSolverException:
        return []
    else:
#         print(p.show())
        x = sorted(p.get_values(w, convert=ZZ, tolerance=1e-3).items())
        for i in x:
            ops.append(i[1])
        return(ops)


# In[40]:


base_set_pair = []
for i_1 in range(0,len(base_mon)):
    for i_2 in range(i_1 + 1,len(base_mon)):
        for i_3 in range(i_2+1,len(base_mon)):
            for i_4 in range(i_3+1,len(base_mon)):
                if ops_solution_base(base_mon[i_1],base_mon[i_2],base_mon[i_3],base_mon[i_4])==[]:
                    base_set_pair.append([base_mon[i_1],base_mon[i_2],base_mon[i_3],base_mon[i_4]])
                    
base_set_pair 


# In[41]:


len(base_set_pair)


# In[60]:


fiber_set_pair = []
for i_1 in range(0,len(fiber_mon)):
    for i_2 in range(i_1 + 1,len(fiber_mon)):
        if ops_solution_fiber(fiber_mon[i_1],fiber_mon[i_2])!=[]:
            fiber_set_pair.append([fiber_mon[i_1],fiber_mon[i_2]])
                    
len(fiber_set_pair)


# In[61]:


fiber_set_pair


# In[65]:


ops_solution_total( [(0, 0, 0, 3), (0, 0, 1, 2), (0, 0, 2, 1), (1, 1, 1, 0)], [(0, 0, 1, 1)] )


# In[67]:


counter =0
for base in base_set_pair:
    for fiber in fiber_set_pair:
        counter = counter +1 
        print("counter:",counter)
        if ops_solution_total(base,fiber) != []:
            print("warning:bug")
            break


# In[62]:


# the final results:
max_pair_list = []
counter  = 0
for base in base_set_pair:
    for fiber in fiber_set_pair:
        counter = counter +1 
        print("counter:",counter)
        print("length:", len(max_pair_list))
        if ops_solution_total(base,fiber) != []:
            if max_pair_list == []:
                max_pair_list.append((base,fiber))
            else:
                sign = True
                for i in max_pair_list:
                    if mon_set_pair_comparison(i,(base,fiber),baseMonPoset,fiberMonPoset) == True: # i>=(base,fiber)
                        sign = False
                        break
                    else:
                        continue
                if sign ==True:
                    new_max_pair_list = list(filter(lambda y: mon_set_pair_comparison((base,fiber),y,baseMonPoset,fiberMonPoset) ==False,   max_pair_list))
                    new_max_pair_list.append(x)  
                    max_pair_list = max_pair_list
                else:
                    continue
        else:
            continue

                

max_pair_list


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:



# laza & our notation:
# ops_compare(ops_1,ops_2,Pos) == True if (ops1<=0) >= (ops2<=0) i.e. ops_1 is bigger or equal to ops_2
def ops_compare(ops_1, ops_2,Pos):
    nonpositivemon_2 = nonpositive_mon_list(ops_2,Pos)
    max_12 = True
    for x in nonpositivemon_2:
        if weight_total(x ,ops_1) > 0:
            max_12 = False
            break
    return max_12


# In[18]:


def take_maxops(ops_list,Pos):
    maxops_list = []
    counter  = 0
    for x in ops_list:
        counter = counter +1 
        print("counter:",counter)
        print("length:", len(maxops_list), maxops_list)
        if nonpositive_mon_list(x,Pos) != []:
            if maxops_list == []:
                maxops_list.append(x)
            else:
                sign = True
                for i in maxops_list:
                    if ops_compare(i,x,Pos) == True: # i>=x
                        sign = False
                        break
                    else:
                        continue
                if sign ==True:
                    new_maxopslist = list(filter(lambda y: ops_compare(x,y,Pos) ==False,   maxops_list))
                    new_maxopslist.append(x)  
                    maxops_list = new_maxopslist
                else:
                    continue
        else:
            continue
    return maxops_list


# In[14]:


MonPoset = baseMonPoset.product(fiberMonPoset)
MonPoset


# In[26]:


testlist  = []
for x_0 in range(0,40):
    for x_1 in range(-40, x_0 + 1):
        for x_2 in range(-40, x_1 + 1):
            x_3 = - (x_2 +  x_1 + x_0)
            if x_3 <= x_2:
                testlist.append([x_0, x_1, x_2, x_3])
                
testlist.sort()
testlist.remove([0,0,0,0])
len(testlist)


# In[27]:



import time
t = time.perf_counter()
S = take_maxops(testlist,MonPoset)
S
print(len(S),S)
print(f'coast:{time.perf_counter()- t:.8f}s')


# In[19]:


X =[[1, 0, 0, -1], [1, 1, 0, -2], [1, 1, 1, -3], [2, 0, -1, -1], [3, 0, -1, -2], [5, 1, -3, -3], [6, 5, 0, -11], [9, 1, 0, -10], [10, 1, -1, -10], [15, 15, 1, -31], [18, 2, -1, -19], [19, 0, -9, -10], [19, 1, -1, -19], [19, 1, 0, -20], [19, 10, 1, -30], [20, 0, -1, -19], [20, 1, 0, -21], [20, 19, 0, -39], [21, -1, -10, -10], [21, -1, -1, -19], [21, 0, -1, -20], [21, 1, -2, -20], [21, 1, -1, -21], [21, 3, -1, -23], [21, 19, 0, -40], [21, 20, 0, -41], [21, 21, -1, -41], [22, 21, -1, -42], [28, 24, -1, -51], [29, 2, -1, -30], [30, -1, -10, -19], [30, 1, -1, -30], [30, 29, 2, -61], [31, 2, -3, -30], [31, 29, 1, -61], [32, 1, -3, -30], [35, 1, -6, -30], [37, 2, 1, -40], [39, -1, -19, -19], [39, 0, -19, -20], [39, 1, -2, -38], [39, 36, 1, -76]]
len(X)

