#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
import matplotlib.pyplot as plt
import time
import networkx as nx
import sys
import copy
import datetime
import json
import os
from joblib import Parallel, delayed
#%matplotlib inline


# # Main

# # Modules

# ## Accelerated Gradient

# In[211]:


def backtracking(Lips, aux_sol, 
                 obj_val, obj, obj_args,
                 add_fun, add_fun_args, add_fun_to_obj, 
                 grad_val, bt_param):

    i = 0

    objfar = []
    quadar = []

    def judge(Lips, aux_sol, 
              obj_val, obj, obj_args,
              add_fun, add_fun_args, add_fun_to_obj, 
              grad_val):

        aux_sol_plus = aux_sol + (1/Lips*grad_val)

        if add_fun != None and add_fun_to_obj == True:
            add_args_obj, add_args_grad = add_fun(aux_sol_plus, *add_fun_args)

            if add_fun_to_obj == True:
                obj_args_added_inner = *obj_args, *add_args_obj
            else:
                obj_args_added_inner = *obj_args,

        else:
            obj_args_added_inner = *obj_args,

        left = obj(aux_sol_plus, *obj_args_added_inner) - obj_val
        right = (1/(2*Lips))*(grad_val.T@grad_val)
        return (left*(1+1e-7)) >= right

    while judge(Lips, aux_sol, 
                obj_val, obj, obj_args,
                add_fun, add_fun_args, add_fun_to_obj, 
                grad_val) == False:# or convjudge < 1e-2:
        Lips = Lips*bt_param

    return Lips

def add_args(aux_sol, add_fun, add_fun_args, obj_args, grad_args, add_fun_to_obj, add_fun_to_grad):
    
    add_args_obj, add_args_grad = add_fun(aux_sol, *add_fun_args)
            
    if add_fun_to_obj == True:
        obj_args_added = *obj_args, *add_args_obj
    else:
        obj_args_added = *obj_args,

    if add_fun_to_grad == True:
        grad_args_added = *grad_args, *add_args_grad
    else:
        grad_args_added = *grad_args,
        
    return obj_args_added, grad_args_added

def acc_grad(init_sol, obj, grad, *, obj_args=(), grad_args=(), 
             add_fun=None, add_fun_args=(), add_fun_to_obj=False, add_fun_to_grad=False,  
             conv_thres=1e-4, bt_param=1.1, min_iter_for_restart=5, printfreq=100):
    """
    solve an unconstrained minimization problem using the accelerated gradient method (Beck and Teboulle, 2009).
    
    input:
        init_sol(numpy array): initial solution
        obj(callable): A function for culculating the objective function (numpy array + additional argments => float)
                       The first argument must be the current solution of this problem (numpy array).
        grad(callable): A function for culculating the gradient of obj() (numpy array + additional argments => numpy array)
                        The first argument must be the current solution of this problem (numpy array).
                        
        obj_args(tuple): additional arguments for the function obj()
        grad_args(tuple): additional arguments for the function grad()
        
        add_fun(callable): additional function for calculate arguments for obj() and/or grad()
                           The first argument must be the current solution of this problem (numpy array).
                           it must return tuple of two tuples--- arguments for obj and grad.
                           (you should use this if obj() and grad() need the same result of 
                           a heavy calculation using the current solution. 
                           Without this, you have to calculate it twice, in obj() and add().)
        add_fun_args(tuple): additional arguments for the function add_fun()
        add_fun_to_obj(bool): whether the arguments from add_fun() would be passed to obj()
        add_fun_to_grad(bool): whether the arguments from add_fun() would be passed to grad()
        
        conv_thres(float): a threshold for convergence check
                           if all of the squared values of elements of gradient is smaller than conv_thres, stop.
        bt_param(float): a parameter for the backtracking process
        min_iter_for_restart(integer): the minimum iteration number in restart process
        printfreq(integer): determine how frequently output the current situation of the process of this function.
    """
    
    # initial error check =========================
    
    if add_fun == None and (add_fun_to_obj == True or add_fun_to_grad == True):
        print("You want to add arguments to objective or gradient but add_fun is not given.")
        return 1
    
    # initial settings ============================
    
    iteration = 0
    convjudge = 10000000
    sol = init_sol.copy()
    sol_fwd = sol.copy()
    aux_sol = sol.copy()
    Lips = 1 #8は現在の最小num_node.最小num_nodeを変えるときはここも変えること．
    t = 1
    j = 0
    kmin = min_iter_for_restart
    
    # main process =================================
    
    start_time = time.process_time()
    
    x = []
    y_obj = []
    y_conv = []
    
    
    while convjudge > conv_thres:
        
        iteration += 1
        
        #暫定状況の出力
        if iteration%printfreq == 0:
            print("processing...", iteration, "objective function:", obj_val, "conv:",np.max(grad_val**2), "L:", Lips)
         
             #グラフの出力－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－
            
#             x = np.hstack((x, iteration))
#             y_obj = np.hstack((y_obj, obj_val))
#             y_conv = np.hstack((y_conv, (np.max(grad_val**2)**0.5)))
            
#             plt.cla()

#             plt.plot(x, y_obj, label="bt_param_{bt}".format(bt=bt_param))
#             #plt.plot(x2, y2, marker="o", label="label2")
#             #plt.errorbar(x, y, yerr=y_error, fmt='none', capsize=10)

#             plt.xlabel("iteration")
#             plt.ylabel("objective function value")

#             #plt.xscale('log')
#             plt.yscale('log')

#             plt.grid(axis="x")
#             plt.grid(axis="y")
#             #plt.legend()

#             filename = "result/img/free_graph/png"
#             os.makedirs(filename, exist_ok=True)
            
#             plt.savefig(filename + "/iteration_obj_{bt}.png".format(bt=bt_param))


#             plt.cla()
            
#             plt.plot(x, y_conv, label="bt_param_{bt}".format(bt=bt_param))
#             #plt.plot(x2, y2, marker="o", label="label2")
#             #plt.errorbar(x, y, yerr=y_error, fmt='none', capsize=10)

#             plt.xlabel("iteration")
#             plt.ylabel("convergence_condition")

#             #plt.xscale('log')
#             plt.yscale('log')

#             plt.grid(axis="x")
#             plt.grid(axis="y")
#             #plt.legend()

#             filename = "result/img/free_graph/png"
#             os.makedirs(filename, exist_ok=True)


#             plt.savefig(filename + "/iteration_conv_{bt}.png".format(bt=bt_param))
        
        
        
        
        if add_fun_to_obj == True or add_fun_to_grad == True:
            obj_args_added, grad_args_added = add_args(aux_sol, add_fun, add_fun_args, obj_args,
                                                       grad_args, add_fun_to_obj, add_fun_to_grad)
        else:
            obj_args_added = *obj_args,
            grad_args_added = *grad_args,
            
        obj_val = obj(aux_sol, *obj_args_added)
        #print("sol", aux_sol)
        grad_val = grad(aux_sol, *grad_args_added)
        #print("grad", grad_val)
        Lips = backtracking(Lips, aux_sol, 
                            obj_val, obj, obj_args, 
                            add_fun, add_fun_args, add_fun_to_obj, 
                            grad_val, bt_param)
        
        convjudge = np.max(grad_val**2)
        
        sol_fwd = aux_sol + grad_val/Lips
        t_fwd = (1 + np.sqrt(1+4*(t**2)))*0.5
        aux_sol = sol_fwd + (t - 1)/t_fwd*(sol_fwd - sol)
        
        sol = sol_fwd

        t = t_fwd
        
        if grad_val.T@(sol_fwd-sol)>0 and j>=kmin:
            j = 0
            t = 1
        else:
            j += 1
            
    proc_time = time.process_time() - start_time
    
    return sol, obj_val, proc_time


# ## Calculation of choice probability, obj and grad

# In[212]:


def prob(sol, exp_cost_or, cost_rs, exp_cost_sd, num_node, num_depot, logit_param_q, logit_param_n):
    
    price = sol.flatten().reshape([num_depot, num_node])
    exp_cu_rs = np.exp(-logit_param_q*(cost_rs - price))
    exp_od = exp_cost_sd
    
    sumexp_os = exp_cost_or@exp_cu_rs + 1.0e-300 # size: i*s = i*r + r*s
    sumexp_od = sumexp_os@exp_cost_sd + exp_od + 1.0e-300 # size: i*j = i*s + s*j

    Prob_rs = np.zeros([num_node, num_depot, num_node])
    Prob_sd = np.zeros([num_node, num_node, num_node])
    Prob_rs_1 = np.zeros([num_depot, num_node])
    
    cmnt = False

    for o in range(num_node):
        exp_cost_or_o = exp_cost_or[o,:].reshape([1,-1])
        sumexp_os_o = sumexp_os[o,:].reshape([1,-1])
        sumexp_od_o = sumexp_od[o,:].reshape([1,-1])
        exp_ors_o = exp_cost_or_o.T*exp_cu_rs # size: r*s = r*1 times r*s
        exp_osd_o = sumexp_os_o.T*exp_cost_sd # size: s*j = s*1 times s*j
        prob_rs_o = exp_ors_o*(sumexp_os_o**-1) # size: r*s = r*s / 1*s sumが1であることを確認

        if cmnt and i == 6:
            print("== sumprs ======")
            print(np.sum(prs, axis = 0))

        prob_sd_o = exp_osd_o*(sumexp_od_o**-1) # size: s*j = s*j / 1*j sumが1であることを確認

        if cmnt and i == 6:
            print("== sumpsj ======")
            print(np.sum(psj, axis = 0))

        Prob_rs[o] = prob_rs_o
        Prob_sd[o] = prob_sd_o

    exp_cu_rs_0 = np.exp(-logit_param_n*cost_rs)
    exp_cu_rs_1 = np.exp(-logit_param_n*price)
    
    sumexp_rs_k = exp_cu_rs_0 + exp_cu_rs_1
    
    prob_rs_1 = exp_cu_rs_1*(sumexp_rs_k**-1)
    Prob_rs_1 = prob_rs_1
    
    return (sumexp_od.flatten().reshape([-1,1]), sumexp_rs_k.flatten().reshape([-1,1]), ), (Prob_rs, Prob_sd, Prob_rs_1)


def grad(sol, num_drv, num_task, num_node, num_depot, Prob_rs, Prob_sd, Prob_rs_1):
    
    flow_rs = np.zeros([num_depot*num_node, 1])
    
    for o in range(num_node):
        flow_o = num_drv[num_node*o:num_node*(o+1), :].reshape([-1,1])
        flow_s_o = flow_o.T@Prob_sd[o,:,:].T # 1*s
        flow_rs_o = Prob_rs[o,:,:] * flow_s_o # r*s
        flow_rs += flow_rs_o.flatten().reshape([-1,1])
        
    task_rs = np.zeros([num_depot*num_node, 1])
    task_rs = (Prob_rs_1.reshape([-1,1]))*num_task
    
    grad_val = task_rs - flow_rs
    
    return grad_val

def obj(sol, num_drv, num_task, logit_param_q, logit_param_n,sumexp_od, sumexp_rs_k):
    
    S = - (1/logit_param_q) * np.log(sumexp_od)
    SS = - (1/logit_param_n) * np.log(sumexp_rs_k)
    
    obj_val = num_task.T @ SS + num_drv.T @ S
    
    return obj_val[0,0]


# ## Calculation of task allocation

# In[301]:


def allocation_under_price(price, exp_cost_or, cost_rs, exp_cost_sd, num_drv, num_task, num_node, num_depot, logit_param_q, logit_param_n):
    sumV, P = prob(price, exp_cost_or, cost_rs, exp_cost_sd, num_node, num_depot, logit_param_q, logit_param_n)
    p_rs, p_sd, p_rs_1 = P
    sumexp_od, sumexp_rs_k = sumV
    sumV = sumV[0]
    P_rs = np.vstack([np.vstack([p_rs[o].flatten()]*num_node) for o in range(num_node)])
    P_sd = np.hstack([np.vstack([p_sd[o].T for o in range(num_node)])]*num_depot)
    P_rs_1 = p_rs_1.flatten().reshape([-1, 1])
    
    real_allocation_q = P_rs * P_sd * num_drv
    real_allocation_n = P_rs_1 * num_task
    #print(np.sum(num_drv))
    #print(np.sum(num_task))
    #real_allocation_q = np.hstack((real_allocation_q, num_drv - np.sum(real_allocation_q, axis=1).reshape([-1,1])))
    real_allocation_q = np.hstack((real_allocation_q, num_drv / sumexp_od * exp_cost_sd.reshape([-1, 1])))
    real_allocation_n = np.hstack((num_task - np.sum(real_allocation_n, axis=1).reshape([-1, 1]), real_allocation_n))
    #print(np.sum(real_allocation_q))
    #print(np.sum(real_allocation_n))
    
    
    return real_allocation_q, real_allocation_n

def set_priority(pos_key_array, sort_key_array, ascend=False):
    # Return indices of the positive values in "pos_array" 
    #  sorted by the value of "sort_key_array" in descending(ascending if ascend==True) order
        
    include_zero = False
    
    picked_index = np.where(pos_key_array > 0)[0] # pick up indices(od pair) where # of driver is greater(smaller) than that of assigned tasks.
    if picked_index.size == 0:
        print("No elements are positive. zero will be included.")
        picked_index = np.where(pos_key_array >= 0)[0]
        include_zero = True
        
        if picked_index.size == 0:
            print("All elements are negative. Adjust another od pair before.")
            return 1, False

    picked_value = sort_key_array[picked_index] # pick up the value of f_od^rs for all od in pickup_index

    if ascend == True:
        #割り当てられたタスク数の少ない順にrs_indexを並べる
        priority = picked_index[np.argsort(picked_value)] # sort picked up indices in ascending order according to f_od^rs
    else:
        #割り当てられたタスクの多い順にod_indexを並べる
        priority = picked_index[np.argsort(picked_value)][::-1]

    return priority, include_zero

def give_to_other_od(round_f_q, od, resid_driver):
    num_resid = resid_driver[od]
    
    if num_resid < 0: # if od need to reduce tasks,,, #割り当てられたタスクが多い場合(常に成り立つ)
        # otherwise do nothing(wait until tasks are given from other ods.)
        cut_priority = set_priority(round_f_q[od, :], round_f_q[od, :], ascend=True)[0]
        excess = -num_resid #余分に配分されたタスク数
                                    
        for rs in cut_priority: #配分されたタスクが少ないrs_index
            add_od = set_priority(resid_driver, round_f_q[:, rs])[0][0] #タスクが足りない and rsタスクの実行数が多いod_index
                
            if excess >= round_f_q[od, rs]: 
                excess -= round_f_q[od, rs] #割り当ての少ないタスクindexから順にexcess分のタスクを削除
                
                #if round_f_q[od, rs] <= 
                round_f_q[add_od, rs] += round_f_q[od, rs] # give tasks to add_od that task rs is assigned most.
                round_f_q[od, rs] = 0 # cut all tasks assigned to "od".
                #print(resid_driver.shape)
                #resid_driver[od, 0] += round_f_q[od, rs]
                #resid_driver[add_od, 0] -= round_f_q[od, rs]
            
            else:
                round_f_q[add_od, rs] += excess # give tasks to add_od that task rs is assigned most.
                round_f_q[od, rs] -= excess # cut part of tasks assigned to "od".
                #resid_driver[od, 0] += excess
                #resid_driver[add_od, 0] -= excess
                
                break
                
                
    return round_f_q
        
#配列Aの i 列要素の非ゼロ要素から，A_keyの値が大きい順に1ずつ減らし， num_plus 個だけ減らす関数
def minus_driver(A, A_key, i, num_plus):

    reduced_A = np.copy(A)
    num_iteration = 0
    
    while num_plus != 0:
        
        
        
        A_1 = reduced_A[:, i].reshape([-1, 1])
        #print(A_1)
        #num_nonzero = np.count_nonzero(A_1)
        #print(num_nonzero)

        # ind = np.nonzero(A_1)
        # nonzero_index = np.transpose(ind)
        # print(nonzero_index)

        nonzero_index = np.argwhere(A_1 > 0)
        #print(A_1)
        #print(nonzero_index[:,0])

        nonzero_key = A_key[nonzero_index[:, 0], i]
        #print(np.size(nonzero_key))
        

        sorted_index = np.flipud(np.argsort(nonzero_key, axis = 0).flatten())
        #print((sorted_index))

        #順番を取得
        #print(nonzero_index[sorted_index][:, 0])



        for j in nonzero_index[sorted_index][:, 0]:

            if num_plus != 0:
                reduced_A[j, i] -= 1
                num_plus += 1
            #print(j, reduced_A[j, i])
                
        if num_iteration>1000:
            print("Error of driver_rounding\n\n\n")
            break
            
        #print(num_plus)
        num_iteration+=1
                
    return reduced_A
    
    
def all_rs_adjust_rounding(round_f_q, f_q, resid_task):
    
    roundoff_error = round_f_q - f_q #真の配分からのプラス誤差配列　(まるめ　―　真)
    #print(resid_task[-1])
    
    #print(len(resid_task.flatten()))
    #print(len(roundoff_error.flatten()))
    
    for rs, resid_rs in enumerate(resid_task.flatten()): #タスクの余りを一列に並べる
        if resid_rs == 0:
            continue
            
        add_priority = np.argsort(roundoff_error[:, rs])#真の配分からの誤差が小さい順(タスクが足りない順)にソートし，タスク追加の優先順位を取得
        if resid_rs > 0: # need to allocate more
            for i in range(resid_rs):
                round_f_q[add_priority[i], rs] += 1 # プラス誤差が小さい(マイナス誤差が大きい)順に足していく．
        else: # need to decrease allocation
            #for i in range(-resid_rs):
            round_f_q = minus_driver(round_f_q, roundoff_error, rs, resid_rs)# プラス誤差が大きい(=マイナス誤差が小さい)順に引いていく．
        #print(rs, resid_rs)
    return round_f_q
            
            
def roundoff_allocation(allocation_q, allocation_n, num_driver, num_task):
    
    #print(allocation_q.shape, allocation_n.shape)
    
    sum_task = np.array([int(np.round(np.sum(allocation_n)))])
    sum_driver = np.array([int(np.round(np.sum(allocation_q)))])
    #print(sum_task)
    #print(sum_driver)
    
    f_q = allocation_q
    f_n = np.delete(allocation_n, 0, axis=1)
    #print(f_q.shape, f_n.shape)
    
    round_f_q = np.round(allocation_q).astype(int)    #ドライバー配分数を四捨五入
    #round_f_q_last = round_f_q[:, -1]
    #round_f_q = np.delete(round_f_q, -1, axis=1)
    round_f_n = np.round(f_n).astype(int) #実行タスク数を四捨五入
    # print(round_f_n)
    round_f_n = np.vstack((round_f_n, [sum_task - np.sum(round_f_n)]))
    #print(round_f_n, np.sum(round_f_q, axis=0).reshape([-1,1]).astype(int))
    resid_task = round_f_n - np.sum(round_f_q, axis=0).reshape([-1,1]).astype(int) #余るタスク数を計算
    # print(resid_task)
    
    round_f_q = all_rs_adjust_rounding(round_f_q, f_q, resid_task) # residual about task is adjusted to zero.
    
    #round_f_q = np.hstack((round_f_q, round_f_q_last.reshape([-1, 1])))
    #print(round_f_q)
    resid_driver = num_driver.astype(int) - np.sum(round_f_q, axis=1).reshape([-1,1]).astype(int) #ドライバーの余りを計算(真-まるめ)
    #print(resid_driver)
    
    
    while (np.count_nonzero(resid_driver)>0):
        for od in np.where(resid_driver < 0)[0]:
            round_f_q = give_to_other_od(round_f_q, od, resid_driver)
            resid_driver = num_driver.astype(int) - np.sum(round_f_q, axis=1).reshape([-1, 1]).astype(int)
            #print(resid_driver)
            #print(len(np.nonzero(resid_driver)[0]))
            
    
    round_f_n = np.delete(round_f_n, -1, axis=0)
        
    #round_f_q = np.hstack((round_f_q, num_driver - np.sum(round_f_q, axis=1).reshape([-1,1])))
    round_f_n = np.hstack((num_task - np.sum(round_f_n, axis=1).reshape([-1,1]), round_f_n))
    
    round_f_q = round_f_q.astype(int)
    round_f_n = round_f_n.astype(int)
    
    #print(resid_driver)
    
    return round_f_q, round_f_n
    


# ## Data Generation

# In[302]:


# TODO: generator of q, n, epsilon from given the average number of drivers and the number of nodes, depots.

def generate_demand(num_node, num_depot, num_sum_driver):
    
    #q_rand = random.random()*0.2 + 0.9
    
    #whole_driver = numdrv_od*num_node*num_node*q_rand
    #whole_driver = int(whole_driver)
    
    whole_driver = int(num_sum_driver)
        
    qq = np.random.rand(num_node**2)
    q = np.floor((qq/np.sum(qq)) * whole_driver)
    q = np.where(q == 0, 1, q)
    q_resid = whole_driver - int(np.sum(q))
    q_distr = np.random.choice(num_node**2, size = q_resid, replace = False)
    q_add = np.zeros(num_node**2)
    q_add[q_distr] = 1
    q += q_add

    #nn_rand = random.random()*0.2 + 0.9
    
    #whole_task = numdrv_od*num_node*num_node*nn_rand
    #whole_task = int(whole_task)
    
    whole_task = int(num_sum_driver)
    
    nn = np.random.rand(num_depot * num_node)
    n = np.floor((nn/np.sum(nn)) * whole_task)
    n = np.where(n == 0, 1, n)
    n_resid = whole_task - int(np.sum(n))
    n_distr = np.random.choice(num_node*num_depot, size = n_resid, replace = False)
    n_add = np.zeros(num_node*num_depot)
    n_add[n_distr] = 1
    n += n_add
    
    #行列にして返す
    return q.reshape([num_node, num_node]), n.reshape([num_depot, num_node])
        
def create_city_network(num_node):
    
    def calc_weight(n1, n2, pos_dict):
        x1, y1 = pos_dict[n1]
        x2, y2 = pos_dict[n2]
        weight = ((x1-x2)**2 + (y1-y2)**2)**(1/2)
        return weight

    def find_nearest(node, G, pos_dict, num):
        dist_dict = {n:calc_weight(n, node, pos_dict) for n in G.nodes()}
        sort_n_list = sorted(dist_dict.keys(), key = lambda x:dist_dict[x])
        picked = []

        for sn in sort_n_list[1:]:
            if (node,sn) not in G.edges():
                picked.append(sn)
            if len(picked) == num:

                break

        return picked

    def find_net_farthest(node, G, pos_dict, num):
        non_neigh = list(nx.non_neighbors(G, node))
        dist_dict = {n:calc_weight(n, node, pos_dict) for n in non_neigh}
        dijk = nx.single_source_dijkstra_path_length(G, node)
        net_dist_dict = {n:dijk[n] if n in dijk.keys() else 1000000 for n in non_neigh}

        #sort_by_net_farness = sorted(, key = lambda x:net_dist_dict[x])[::-1]
        sort_n_list_master = sorted(non_neigh, key = lambda x:dist_dict[x])
        #pick_priority = sorted(sort_n_list_master[:4], key = lambda x:sort_by_net_farness.index(x))
        pick_priority = sorted(sort_n_list_master[:4], key = lambda x:(dist_dict[x] - net_dist_dict[x]))

        picked = pick_priority[:num]

        return picked

    root = int(num_node**0.5) #都市ネットワークの1辺に並べるノード数
    square_length = root*5 #都市ネットワークの1辺の長さ

    # set coordination of nodes
    X = np.random.rand(root,root)*5 #0.0以上1.0未満の乱数×5(一様分布)(root×rootの個数分)
    Y = np.random.rand(root,root)*5 #上に同じ
    xplace = np.arange(0, square_length, 5) #初項0，公差5，(1辺の長さ未満)の等差数列を作成
    yplace = np.arange(0, square_length, 5) #上に同じ
    Xp, Yp = np.meshgrid(xplace, yplace) #Xpは各格子のX番号，Ypは各格子のY番号
    X += Xp #各ノードのX座標を作成
    Y += Yp #各ノードのY座標を作成
    X = X.flatten() #1次元配列に変換
    Y = Y.flatten() #上に同じ

    # 中心から遠い順に並べ替えたリストを作成
    pos = sorted([(X[i], Y[i]) for i in range(len(X))], 
                 key = lambda x:((x[0] - square_length*0.5)**2 + (x[1] - square_length*0.5)**2),
                 reverse = True)

    nodelist = list(range(num_node)) #ノード数と同じ要素数のリストを作成(要素名は0スタートの昇順整数)
    pos_dict = {node:pos[node] for node in nodelist} #辞書作成

    G = nx.Graph() #無向グラフを作成します
    G.add_nodes_from(nodelist, pos = pos_dict) #ノードリストから頂点リストを作成

    for node in G.nodes():
        number = 1#np.random.randint(1,5)
        picked_node = random.sample(find_nearest(node, G, pos_dict, number), number) #最も近いノードを選択
        for pn in picked_node:
            G.add_edge(node, pn, weight = calc_weight(node, pn, pos_dict)) 
            pass
        pass

    for node in G.nodes():
        number = 1
        picked_node = random.sample(find_net_farthest(node, G, pos_dict, number), number)

        for pn in picked_node:
            G.add_edge(node, pn, weight = calc_weight(node, pn, pos_dict))
            pass
        pass

    # fig, ax = plt.subplots(figsize = [10,10])
    # nx.draw(G, pos = pos_dict, ax = ax, node_size = 100)
    # print(nx.info(G))
    
    return G

def create_distmat(G):
    spl = dict(nx.all_pairs_dijkstra_path_length(G))
    
    distmat = np.array([[spl[r][c] for c in range(nx.number_of_nodes(G))] for r in range(nx.number_of_nodes(G))])
    
    return distmat

def create_costmat(num_node, num_depot, distmat, depot_list):
    def calccost(od_index, rs_index):
        o = od//num_node 
        d = od%num_node
        r = depot_list[rs//num_node]
        s = rs%num_node

        return distmat[o,r] + distmat[r,s] + distmat[s,d]

    costmat = np.zeros([num_node*num_node, (num_node*num_depot+1)])
    for od in range(num_node*num_node):
        o = od//num_node
        d = od%num_node
        costmat[od, num_node*num_depot] = distmat[o,d]
        for rs in range(num_node*num_depot):
            costmat[od, rs] = calccost(od, rs) 

    print("completed making o-r-s-d costs")
    return costmat

def expand_costmat_to_atomic_q(costmat, num_drv):
    # expand costmat(od x rs) to atomic costmat(c_a^rs, a x rs)
    rs = costmat.shape[1]
    for od_index, q_od in enumerate(num_drv.flatten()):
        #タスク実行コストを一律に全員分確保
        atom_od = np.ones([int(q_od), rs]) * (costmat[od_index, :])
        #全てのトリップODペアについて全員分のコストを縦に連結
        if od_index == 0:
            atom_costmat_q = atom_od
        else:
            atom_costmat_q = np.vstack([atom_costmat_q, atom_od])
    
    return atom_costmat_q



def expand_distmat_to_atomic_n(distmat, num_task, depot_list, num_depot, num_node):
    
    for rs_index, n_rs in enumerate(num_task.flatten()):
        r = depot_list[rs_index//num_node]
        s = rs_index%num_node
        atom_rs = np.ones([int(n_rs), 2]) * ([distmat[r,s], 0])
        if rs_index == 0:
            atom_costmat_n = atom_rs
        else:
            atom_costmat_n = np.vstack([atom_costmat_n, atom_rs])
            
    return atom_costmat_n



def create_atomic_costmat(costmat, distmat, num_drv, num_task, logit_param_q, logit_param_n, depot_list, num_depot, num_node):
    #確定項の行列
    fixed_q = expand_costmat_to_atomic_q(costmat, num_drv)
    #誤差項の行列
    error_q = np.random.gumbel(0, 1/logit_param_q, fixed_q.shape)
    
    fixed_n = expand_distmat_to_atomic_n(distmat, num_task, depot_list, num_depot, num_node)
    error_n = np.random.gumbel(0, 1/logit_param_n, fixed_n.shape)
    
    return fixed_q + error_q, fixed_n + error_n


# ## Class for the problem settings of CSD matching

# In[303]:


class CSD:
    # Execute matching under the situation that is given by a DataMaker instance.
    # All input datas are( TODO: SHOULD BE! ) immutable. 
    def __init__(self, num_drv, num_task, dist_matrix, depot_node_list, logit_param_q, logit_param_n):
        self.num_drv = num_drv
        self.num_task = num_task
        self.distmat = dist_matrix
        self.logit_param_q = logit_param_q
        self.logit_param_n = logit_param_n
        self.num_node = dist_matrix.shape[0]
        self.num_depot = len(depot_node_list)
        
        self.cost_or = dist_matrix[:, depot_node_list]
        self.exp_cost_or = np.exp(-logit_param_q * self.cost_or)
        self.cost_rs = self.cost_or.T
        self.cost_sd = dist_matrix
        self.exp_cost_sd = np.exp(-logit_param_q * self.cost_sd)
        
    def solve_dual(self):
        
        price, obj_val, time = acc_grad(np.zeros(self.num_task.shape), obj, grad, 
                                        obj_args=(self.num_drv, self.num_task, self.logit_param_q, self.logit_param_n), 
                                        grad_args=(self.num_drv, self.num_task, self.num_node, self.num_depot), 
                                        add_fun=prob, add_fun_args=(self.exp_cost_or, self.cost_rs, self.exp_cost_sd, 
                                                                    self.num_node, self.num_depot, self.logit_param_q, self.logit_param_n),
                                        add_fun_to_obj=True, add_fun_to_grad=True, 
                                        conv_thres=1.0**2)
        
        self.opt_price = price
        self.opt_welfare = obj_val
        self.time_for_solve_accel = time
    
    def allocate_tasks(self):
        # determine the task allocation for each od by integer
        self.real_allocation_q, self.real_allocation_n = allocation_under_price(self.opt_price, self.exp_cost_or, self.cost_rs, 
                                                      self.exp_cost_sd, self.num_drv, self.num_task, self.num_node, 
                                                      self.num_depot, self.logit_param_q, self.logit_param_n)

        
    def rounded_allocation(self): #タスク配分を整数値に変換
        self.task_allocation_q, self.task_allocation_n = roundoff_allocation(self.real_allocation_q.copy(), self.real_allocation_n.copy(), self.num_drv, self.num_task)
        
    def save_real_allocation(self, path_q, path_n):
        #output the obtained allocation to a .txt file(for MATLAB)
        np.savetxt(path_q, self.real_allocation_q, delimiter=',')
        np.savetxt(path_n, self.real_allocation_n, delimiter=',')
    
    def save_rounded_allocation(self, path_q, path_n):
        np.savetxt(path_q, self.task_allocation_q, delimiter=',')
        np.savetxt(path_n, self.task_allocation_n, delimiter=',')
    

class DataMaker:
    # Generate random data(cost, network structure, demand) or import data from files, and create a CSD instance from them.
    def __init__(self, num_node, num_depot, num_sum_driver, logit_param_q, logit_param_n):
        self.num_node = num_node
        self.num_depot = num_depot
        self.num_sum_driver = num_sum_driver
        self.logit_param_q = logit_param_q
        self.logit_param_n = logit_param_n
        
    def make_demand_supply(self):
        # generate OD demands of driver and task randomly. the generated demand is not 0.
        self.num_drv, self.num_task = generate_demand(self.num_node, self.num_depot, self.num_sum_driver)
        
    def load_demand_supply(self, q, n):
        self.num_drv, self.num_task = q, n
        
    def save_demand_supply(self, path_driver, path_shipper):
        np.savetxt(path_driver, self.num_drv, delimiter = ",")
        np.savetxt(path_shipper, self.num_task, delimiter = ",")
        
    def make_network(self):
        # generate randomly connected? network with link weight based on its euqlid length
        self.network = create_city_network(self.num_node)
        print("completed making a network")
        
    def load_network(self, G):
        self.network = G
        print("completed loading the network")
        
    def save_network(self, path):
        Gj = nx.node_link_data(self.network)
        with open(path, mode = "w") as f:
            json.dump(Gj, f)
    
    def make_distmat(self):
        self.distmat = create_distmat(self.network)
        print("completed making the distance matrix")
        
    def load_distmat(self, distmat):
        self.distmat = distmat
        print("completed loading the distance matrix")
        
    def save_distmat(self, path):
        np.savetxt(path, self.distmat, delimiter=",")
        
    def make_depot_list(self):
        deplist = np.random.choice(self.num_node, self.num_depot, replace=False)
        self.depot_list = deplist
        
    def load_depot_list(self, deplist):
        self.depot_list = deplist
        
    def save_depot_list(self, path):
        np.savetxt(path, self.depot_list, delimiter = ",")
    
    def make_costmat(self):
        # calculate the travel cost for delivery, o-r-s-d.
        self.costmat = create_costmat(self.num_node, self.num_depot, self.distmat, self.depot_list)
        
    def load_costmat(self, costmat):
        self.costmat = costmat.reshape([num_node*num_node, num_node*num_depot])
        
    def save_costmat(self, path):
        np.savetxt(path, self.costmat, delimiter=",")
    
    def make_atomic_costmat(self):
        if self.num_node > 64:
            print("network may be so large that it takes time. Consider get the size down.")
        self.atomic_costmat_q, self.atomic_costmat_n = create_atomic_costmat(self.costmat, self.distmat, self.num_drv, self.num_task, self.logit_param_q, self.logit_param_n, self.depot_list, self.num_depot, self.num_node)
        
    def load_atomic_costmat(self, atom):
        self.atomic_costmat = atom
        
    def save_atomic_costmat(self, path_q, path_n):
        np.savetxt(path_q, self.atomic_costmat_q, delimiter = ',')
        np.savetxt(path_n, self.atomic_costmat_n, delimiter = ',')
        
    def create_CSD(self):
        return CSD(self.num_drv, self.num_task, self.distmat, self.depot_list, self.logit_param_q, self.logit_param_n)
        
        
