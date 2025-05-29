# -*- coding: utf-8 -*-
# @Time    : 2018/4/23 9:42
# @Author  : ZENG Yanru
# @Email   : 595438103@qq.com
# @File    : count_p_value.py
# @Software: PyCharm
import numpy as np
import pandas as pd
from pADJ_FDR import FDR


class count_p_value:
    __doc__ = "用于计算出p值"

    def __init__(self):
        print("return p value")
        print("输入贡献度的结果列表，也就是类似于：[1,2,3,4]，每个特征一个贡献度。真正预测的和随机打乱的分别输入")
        print("贡献度大于和小于的都要比较，最终返回的p值表示是否差异")

    def cal_p_morethan(self,true_contrib,random_contrib,isFDRAdj=True):
        #true_contrib 一维向量，random_contrib多维向量，有k个样本，k是随机打乱的次数。
        #true_contrib是列表，random_contrib是array或者列表都行，反正要记录两个维度的信息；

        count_feature = len(true_contrib)
        p_value_list = [0 for i in range(count_feature)]
        k = len(random_contrib)
        random_contrib = np.abs(np.array(random_contrib)) #直接取绝对值再算吧，后续应该也不需要random contrib的原始数据了
        df_rd_contrib = pd.DataFrame(random_contrib)
        for i,j in enumerate(true_contrib):
            true_feature_contrib = abs(j)
            df_need = df_rd_contrib[df_rd_contrib[i] > true_feature_contrib] #最终使用了大于号，而不是大于等于，避免所有数据都相等,找到所有符合条件的值
            array_need = np.array(df_need) #转换类型好进行求和
            #这里需要的是统计数目，而不是计算求和！！！
            # sum_need_list = np.sum(array_need,axis=0) #按列求和
            # sum_need = sum_need_list[i] #选取特定的那个特征
            sum_need = len(array_need) #大胆地使用len吧，因为前面已经是挑选出差异的样本了。
            p_value = sum_need/k
            p_value_list[i] = p_value

        if isFDRAdj == True:
            p_value_list = FDR(p_value_list)

        self.morethan_p_value_list = p_value_list
        return p_value_list


    def cal_p_lessthan(self,true_contrib,random_contrib,isFDRAdj=True):
        # true_contrib 一维向量，random_contrib多维向量，有k个样本，k是随机打乱的次数。
        # true_contrib是列表，random_contrib是array或者列表都行，反正要记录两个维度的信息；

        count_feature = len(true_contrib)
        p_value_list = [0 for i in range(count_feature)]
        k = len(random_contrib)
        random_contrib = np.abs(np.array(random_contrib))
        df_rd_contrib = pd.DataFrame(random_contrib)
        # print(k)
        for i, j in enumerate(true_contrib):
            true_feature_contrib = abs(j)
            df_need = df_rd_contrib[df_rd_contrib[i] < true_feature_contrib]  # 最终使用了小于号而不是小于等于，避免所有数据都相等,找到所有符合条件的值
            array_need = np.array(df_need)  # 转换类型好进行求和

            # 这里需要的是统计数目，而不是计算求和！！！
            # sum_need_list = np.sum(array_need, axis=0)  # 按列求和
            # sum_need = sum_need_list[i]  # 选取特定的那个特征
            sum_need = len(array_need)
            p_value = sum_need / k
            p_value_list[i] = p_value

        if isFDRAdj == True:
            # print(p_value_list)
            p_value_list = FDR(p_value_list)

        self.lessthan_p_value_list = p_value_list
        return p_value_list

    def compare_p_value(self,morethan_pvalue_list,lessthan_pvalue_list,cutoff=0.05):
        #more或者less的pvalue list都是前面两个函数得到的形式，真的就是一维的list

        feature_is_important = [0 for i in morethan_pvalue_list]

        for i,j in enumerate(morethan_pvalue_list):
            #首先判断哪一端是显著的
            if j<=cutoff and lessthan_pvalue_list[i]>cutoff:
                feature_is_important[i] = ("Significantly_more",j)
            elif lessthan_pvalue_list[i]<=cutoff and j>=cutoff:
                feature_is_important[i] = ("Significantly_less",lessthan_pvalue_list[i])
            else:
                m = min([j,lessthan_pvalue_list[i]])
                feature_is_important[i] = ("Not_sig",m)

        return feature_is_important

    def return_is_significant(self,true_contrib,random_contrib):
        more = self.cal_p_morethan(true_contrib,random_contrib)
        less = self.cal_p_lessthan(true_contrib,random_contrib)
        feature_is_important = self.compare_p_value(more,less)
        return feature_is_important

