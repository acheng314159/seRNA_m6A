# -*- coding: utf-8 -*-
# @Time    : 2018/4/23 15:51
# @Author  : ZENG Yanru
# @Email   : 595438103@qq.com
# @File    : Construct_model_and_output.py
# @Software: PyCharm

#这里直接使用面向过程即可。
import Get_Contribution_python
import count_p_value
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn import tree
# import pydotplus
import os
import shutil
from sklearn.datasets import load_boston #测试一下脚本是否能用

def export_dgdata_to_storehouse(rf_model,storehouse,feature_list):
    #用于输出多个dg graph的文件
    file_number = len(rf_model.estimators_)
    for estimator_idx in range(file_number):
        dg_filename = storehouse + "dg" + str(estimator_idx) + ".txt"
        dg_data = tree.export_graphviz(rf_model.estimators_[estimator_idx],out_file=None,
                                       feature_names=feature_list,filled=True,rounded=True,
                                       special_characters=True)
        fw = open(dg_filename,"w")
        print(dg_data,file=fw)
        fw.close()

def get_contrib(data,target,storehouse,feature_list):
    # 再构建模型,参数在此处调整
    rf = RandomForestRegressor(n_estimators=100)  # TODO :记得修改参数！
    rf.fit(data, target)

    if os.path.exists(storehouse):
        # 为了避免存储过多，删除storehouse下的所有文件,同时也保证了文件不会重叠
        shutil.rmtree(storehouse)
    os.mkdir(storehouse)

    # 生成多张tree的图片
    export_dgdata_to_storehouse(rf, storehouse, feature_list)

    # 对照片进行处理,得到贡献度
    list_file = os.listdir(storehouse)
    input_file_list = [storehouse + i for i in list_file]
    get_weight = Get_Contribution_python.Get_Contribution_From_RF(input_file_list, feature_list)
    final_contrib = get_weight.cal_weight()
    return final_contrib

def from_dict_to_list(final_contrib,feature_list):
    final_contrib_list = []
    for i in feature_list:
        final_contrib_list.append(final_contrib[i])
    return final_contrib_list

if __name__ == '__main__':
    #首先构造数据 TODO:需要对data，target，featurelist，storehouse，k的次数进行更改，其余的则可以自动弄出来了
    boston = load_boston()
    #data = [[0, 1, 0, 1, 1, 0, 1, 0], [0, 0, 1, 1, 1, 1, 0, 0], [0, 1, 0, 1, 1, 1, 0, 0], [1, 0, 0, 1, 1, 0, 1, 0], [1, 0, 1, 1, 0, 0, 1, 0], [0, 1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 1, 0, 0, 1, 0], [0, 1, 1, 1, 0, 1, 0, 0], [1, 1, 0, 1, 0, 1, 0, 0], [1, 0, 1, 0, 1, 0, 1, 0], [0, 1, 1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0, 0, 1], [1, 0, 1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0, 1, 1], [1, 1, 0, 0, 1, 1, 0, 0], [1, 0, 0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1, 1, 0], [1, 1, 0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 0, 0, 1], [0, 0, 0, 1, 1, 0, 1, 1], [0, 1, 1, 0, 0, 0, 1, 1], [0, 0, 0, 1, 1, 0, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1], [1, 1, 1, 0, 1, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 0, 1, 1], [1, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0, 1, 1], [0, 1, 0, 0, 1, 1, 1, 0], [1, 0, 0, 1, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1, 1, 0], [1, 0, 1, 0, 0, 0, 1, 1], [0, 0, 0, 1, 1, 0, 1, 1], [1, 1, 1, 0, 1, 0, 0, 0], [1, 0, 1, 0, 0, 1, 1, 0], [1, 1, 1, 0, 0, 0, 0, 1]]
    data = xingyang_data
    #target = [ 0.23591076,  0.5298729 ,  0.36944585,  0.12324134,  0.44470259,
    #        0.79279407,  0.85760234,  0.14961072,  0.23067701,  0.06713132,
    #        0.73508112,  0.98167773,  0.66908095,  0.10248019,  0.79662236,
    #        0.62357206,  0.51868634,  0.24064646,  0.08187642,  0.30016405,
    #        0.2939207 ,  0.08140641,  0.39170247,  0.65619811,  0.78948049,
    #        0.12569084,  0.3011816 ,  0.66487847,  0.93774296,  0.81882579,
    #        0.5836517 ,  0.3094207 ,  0.13428659,  0.50406169,  0.14692778,
    #        0.62775501,  0.20250563,  0.43905002,  0.62856858,  0.45497091]
    target = xingyang_target
    #将featurelist做出来，应该是40+个feature
    #feature_list = ["a","b","c","d","e","f","g","h"] #最终输出的结果可以根据这个feature_list来构造新的contribution list
    feature_list = feature_list_xingyang
    # 将模型输出到某个路径下
    storehouse = "temp_dir_xingyang/"

    # 对照片进行处理,得到贡献度
    final_contrib = get_contrib(data,target,storehouse,feature_list)

    #TODO：打印结果，是一个dict，根据需要进行输出即可
    print(final_contrib)
    true_contrib = from_dict_to_list(final_contrib,feature_list)

    #下面进入循环进行k次的随机森林训练，其中因变量和自变量的匹配关系调换一下.
    k = 20 #TODO:重复试验的次数。可以修改这个
    random_contrib = []
    for i in range(k):
        np.random.shuffle(data)
        rd_contrib = get_contrib(data,target,storehouse,feature_list)
        random_contrib.append(from_dict_to_list(rd_contrib,feature_list))

    isSignificant = count_p_value.count_p_value().return_is_significant(true_contrib,random_contrib)
    print(feature_list)
    print(isSignificant)
