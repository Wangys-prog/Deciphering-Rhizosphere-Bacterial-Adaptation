import pandas as pd
import numpy as np
#
# #####################bacteria_g_Caulobacter_rho_merge.csv##############
# data = pd.read_csv('./Niche_differentiation/bacteria_niche/bacteria_g_Caulobacter_rho_merge.csv')
# dict_data ={}
# vals = np.unique(data[['Partner', 'Pair']])  # 同时取出两列,作为节点
# df1 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f = df1.index.get_indexer #获得data中对应索引位置
# df1.values[f(data.Partner), f(data.Pair)] =data["lrv"][f(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df1.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Caulobacter_rho_max.csv")
#
# df2 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f2 = df2.index.get_indexer #获得data中对应索引位置
# df2.values[f2(data.Partner), f2(data.Pair)] =data["nucl_dist"][f2(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df2.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Caulobacter_dnalist_max.csv")
#
# ##############bacteria_g_Chryseobacterium_rho_dnadist_merge.csv
#
# data = pd.read_csv('./Niche_differentiation/bacteria_niche/bacteria_g_Chryseobacterium_rho_dnadist_merge.csv')
# dict_data ={}
# vals = np.unique(data[['Partner', 'Pair']])  # 同时取出两列,作为节点
# df1 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f = df1.index.get_indexer #获得data中对应索引位置
# df1.values[f(data.Partner), f(data.Pair)] =data["lrv"][f(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df1.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Chryseobacterium_rho_max.csv")
#
# df2 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f2 = df2.index.get_indexer #获得data中对应索引位置
# df2.values[f2(data.Partner), f2(data.Pair)] =data["val1"][f2(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df2.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Chryseobacterium_dnalist_max.csv")
#
# ##############bacteria_g_Dyadobacter_rho_dnadist_merge.csv
# data = pd.read_csv('./Niche_differentiation/bacteria_niche/bacteria_g_Dyadobacter_rho_dnadist_merge.csv')
# dict_data ={}
# vals = np.unique(data[['Partner', 'Pair']])  # 同时取出两列,作为节点
# df1 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f = df1.index.get_indexer #获得data中对应索引位置
# df1.values[f(data.Partner), f(data.Pair)] =data["lrv"][f(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df1.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Dyadobacter_rho_max.csv")
#
# df2 = pd.DataFrame(0.00000, index=vals, columns=vals)
# f2 = df2.index.get_indexer #获得data中对应索引位置
# df2.values[f2(data.Partner), f2(data.Pair)] =data["val1"][f2(data.Partner)]
# # for i in range(len(f(data.Partner))):
# #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
# df2.to_csv("./Niche_differentiation/bacteria_niche/bacteria_g_Dyadobacter_dnalist_max.csv")
############################g_Flavobacterium
#导入你的数据
def col2max_(input1,output1,output2):
    data = pd.read_csv(input1)
    dict_data ={}
    vals = np.unique(data[['Partner', 'Pair']])  # 同时取出两列,作为节点
    df1 = pd.DataFrame(0.00000, index=vals, columns=vals)
    f = df1.index.get_indexer #获得data中对应索引位置
    df1.values[f(data.Partner), f(data.Pair)] =data["lrv"][f(data.Partner)]
    # for i in range(len(f(data.Partner))):
    #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
    df1.to_csv(output1)

    # df2 = pd.DataFrame(0.00000, index=vals, columns=vals)
    # f2 = df2.index.get_indexer #获得data中对应索引位置
    # df2.values[f2(data.Partner), f2(data.Pair)] =data["val1"][f2(data.Partner)]
    # for i in range(len(f(data.Partner))):
    #   df.values[f(data.Partner)[i], f(data.Pair)[i]] =float(data["lrv"][i])
    # df2.to_csv(output2)

if __name__ == '__main__':
    # feature_labels()
    import argparse
    parser = argparse.ArgumentParser(
        'Script for HYENADNA model')
    parser.add_argument('-i', type=str, help='input CSV')
    parser.add_argument('-o1', type=str, help='path to saved CSV file')
    # parser.add_argument('-o2', type=str, help='path to saved CSV file')
    args = parser.parse_args()
    input = args.i
    outfile1 = args.o1
    # outfile2 = args.o2
    col2max_(input, outfile1)




