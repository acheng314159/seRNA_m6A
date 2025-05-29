# -*- coding: utf-8 -*-
# @Time    : 2018/4/21 15:24
# @Author  : ZENG Yanru
# @Email   : 595438103@qq.com
# @File    : pADJ_FDR.py
# @Software: PyCharm

# import rpy2.robjects as ro
from statsmodels.sandbox.stats.multicomp import multipletests

def FDR(p_list):
    #希望返回一个FDR值，也就是qvalue。
    #下面是我自己写的，明明公式都一样不知道R里面是不是有什么黑科技，决定还是调用R好了
    # length_p = len(p_list)
    # each_p_name = [str(i) for i in range(length_p)]
    # new_for_sort = []
    # for i,j in enumerate(each_p_name):
    #     new_for_sort.append((j,p_list[i])) #这里用tuple如[（“0”,“0.123”）]记录了所有P值的信息，前面仅仅是名字，用于后续改成输入顺序
    #
    # new_for_sort = sorted(new_for_sort,key=lambda p:p[1]) #使用数值从小到大排序，并且tuple中的第一个数，索引为0，就是对应的原本的索引
    #
    # q_value_dict = {}
    # for i,j in enumerate(new_for_sort):
    #     rank = i + 1
    #     p_i = j[1]
    #     q_i = p_i * length_p / rank
    #     q_value_dict[j[0]] = q_i
    #
    # newlist = []
    # for i in range(length_p):
    #     newlist.append(q_value_dict[str(i)])

    #下面是调用R，但是发现服务器上没法用，应该是底层出了问题！
    # plist = ro.FloatVector(p_list)
    # newRlist = ro.r["p.adjust"](plist,method="fdr")
    # newlist = []
    # for i in newRlist:
    #     newlist.append(i)
    # return newlist

    #使用statsmodel包来做
    padjlist = multipletests(p_list,alpha=0.05,method="fdr_bh")[1].tolist()

    return padjlist

if __name__ == '__main__':
    plist = [0.01, 0.02, 0.10, 0.88, 0.5, 0.12, 0.9]
    print(FDR(plist))


