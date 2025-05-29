# -*- coding: utf-8 -*-
# @Time    : 2018/4/20 15:53
# @Author  : ZENG Yanru
# @Email   : 595438103@qq.com
# @File    : Get_Contribution.py
# @Software: PyCharm

from pygraph.classes.digraph import digraph
import re
import os

class Get_Contribution_From_RF:
    #整个类最终希望返回一个字典，特征以及贡献度。
    __doc__ = "用于计算整个RF的模型中所有特征的贡献度的。不包含模型。默认输入为从tree.export_graphviz中得到的。而且不同的图存在不同的文件中吧" \
              "毕竟树的数目应该不会太多。"

    def __init__(self,input_file_list,feature_list):
        #也即输入input_file_list 和 feature_list就可以得出最终的feature对应的贡献度的值的dict了
        #file list 是digraph的file值
        self.input_file_list = input_file_list
        self.feature_list = feature_list
        self.init_contribution = eval(str(dict.fromkeys(self.feature_list,0))) #初始化全是0的权值
        # print("This class is for calculate the average contribution of decision tree.")


    def init_relation(self,relation_list):
        # incidents（某个节点）指指向这个节点的所有点；这个节点接收。
        # neighbors（某个节点）指这个节点指向的所有节点。这个节点发出。
        #希望返回一个digraph的类，到时候可以直接调用。
        #relation_list,干脆就直接弄成[[a,b],[c,d]...]就好了，里面是tuple也行
        #主要是为了能够对每棵树进行操作，不然读取文件有点麻烦。
        all_nodes = []
        for i in relation_list:
            all_nodes.append(i[0])
            all_nodes.append(i[1])
        all_nodes = list(set(all_nodes))
        dg = digraph()
        dg.add_nodes(all_nodes)
        for i in relation_list:
            dg.add_edge((i[0],i[1]))

        return dg

    def relation_and_value(self,input_file):
        #返回relation的列表和value的值，这个value的值是不用变化的。value对应的key是节点的id
        #dict value key是node的id，value是这个node的特征名字，以及这个节点的value
        #relation就是为了上一个函数构造的，用于做成一张图，好索引
        #end_nodes，也可以代表是树的最终分叉。
        dict_value = {}
        relation = []
        end_nodes = {}
        f = open(input_file,"r")
        for line in f:
            if line.startswith("digraph") or line.startswith("node") or line.startswith("edge") or line.startswith("}"):
                pass
            else:
                if " -> " in line:
                    nodes = line.strip("\n").split(" ")
                    relation.append([nodes[0],nodes[2]])
                else:
                    node = line.split(" ")[0]
                    feature_search = re.search(r"label=<(.*) &le",line)
                    if feature_search == None: #说明是末端的节点
                        value = float(re.search(r"value = (.*)>",line).group(1))
                        end_nodes[node] = value
                    else:
                        value = float(re.search(r"value = (.*)>",line).group(1))
                        feature_name = feature_search.group(1)
                        dict_value[node] = (feature_name,value)
        self.dict_value = dict_value
        self.relation = relation
        self.end_nodes = end_nodes
        return dict_value,relation,end_nodes

    def relation_and_value_R(self, input_file):
        #为了适应R的结果而整理的，为了代码调用方便，也一样生成和relation and value一样的三个结果
        #除了value变成数字以外，其他的全部用字符串表示就好
        # dict value key是node的id，value是这个node的特征名字，以及这个节点的value.
        # relation就是为了上一个函数构造的，用于做成一张图，好索引
        # end_nodes，也可以代表是树的最终分叉。
        dict_value = {}
        relation = []
        end_nodes = {}
        f = open(input_file,"r")
        for line in f:
            if line.startswith("left"):
                pass
            else: #将一行的需要的值弄出来先
                info = line.strip("\n").split("\t")
                this_node = str(int(info[0])-1) #为了和以前构建出来的脚本的索引一致
                left = str(int(info[1])-1)
                right = str(int(info[2])-1)
                featurename = info[3]
                value = float(info[-1])
                if left == "-1":
                    end_nodes[this_node] = value
                else:
                    dict_value[this_node] = (featurename,value)
                    relation.append([this_node,left])
                    relation.append([this_node,right])
        self.dict_value = dict_value
        self.relation = relation
        self.end_nodes = end_nodes
        return dict_value, relation, end_nodes


    def cal_weight(self):
        count = 0
        for input_file in self.input_file_list:
            count += 1
            dict_value,relation,end_nodes = self.relation_and_value(input_file) #TODO:观察使用什么吧，如果使用R则用这个，如果使用python来做RF则把_R去掉
            dg = self.init_relation(relation)
            for end_node in end_nodes:
                end_node_value = end_nodes[end_node] #最后一个节点的操作，这里是取值
                one_lv_parent_node = dg.incidents(end_node)[0] #倒数第二级的节点
                feature_name = dict_value[one_lv_parent_node][0] #通过id找到名字
                feature_value = dict_value[one_lv_parent_node][1] #通过id找到value
                contrib = end_node_value - feature_value
                try:
                    self.init_contribution[feature_name] += contrib
                except:
                    feature_name = feature_name.replace(".", "-")
                    self.init_contribution[feature_name] += contrib
                while one_lv_parent_node != "0": #在根节点之前
                    feature_value_before = feature_value #是唯一的，不用担心
                    one_lv_parent_node = dg.incidents(one_lv_parent_node)[0] #然后继续找前一级的节点
                    feature_name = dict_value[one_lv_parent_node][0] #找名字
                    feature_value = dict_value[one_lv_parent_node][1] #找值
                    contrib = feature_value_before - feature_value #靠子分支的节点的值减去靠母节点的值，才是母节点的特征的贡献度！
                    try:
                        self.init_contribution[feature_name] += contrib
                    except:
                        feature_name = feature_name.replace(".","-")
                        self.init_contribution[feature_name] += contrib
                assert one_lv_parent_node == "0"

        keylist = []
        for key in self.init_contribution:
            keylist.append(key) #避免因为修改dict而导致顺序改变而报错
        for key in keylist:
            self.init_contribution[key] = self.init_contribution[key]/count

        final_contribution = self.init_contribution
        return final_contribution


if __name__ == '__main__':
    filedir = "./temp"
    list_file = os.listdir(filedir)
    input_file_list = [filedir+i for i in list_file]
    feature_list = ['sepal length (cm)', 'sepal width (cm)', 'petal length (cm)', 'petal width (cm)']
    get_weight = Get_Contribution_From_RF(input_file_list,feature_list)
    final_contrib = get_weight.cal_weight()
    print(final_contrib)








