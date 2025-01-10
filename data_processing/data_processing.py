# -*- coding: utf-8 -*- 
import ROOT
import os 
import numpy as np
import tensorflow as tf
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
from tensorflow import keras
from sklearn.preprocessing import OneHotEncoder
import sys
import random
from tqdm import tqdm
##
# print ('creating class MyPyClass ... ')
# tf.random.set_seed(1234) 

# ROOT.gInterpreter.ProcessLine('#include "AMPTdata.h"')
NUM_POINTS_LOW = int(sys.argv[1])
NUM_POINTS_HIGH = int(sys.argv[2])
NUM_DIMENSION = int(sys.argv[3])
PTLOWINP=float(sys.argv[4])
PTHIGHINP=float(sys.argv[5])

PbPbsigma=sys.argv[6]
Pbpsigma=sys.argv[7]
os.environ['CUDA_VISIBLE_DEVICES']=sys.argv[8]
EPOCHS = 50 #训练次数
INITIAL_LR = 1e-4 #学习率
MIXEVENT=1 #混合事件数
NUM_CLASSES = 2  #分类类别数
# NUM_CLASSES = 1  #分类类别数
BATCH_SIZE = 32 #BATCH大小
MAXEVENT=100000 #单个类别样本事件数
folders=["ampt_data/PbPb"+PbPbsigma+"5p02TeV", "ampt_data/PbP"+Pbpsigma+"5p02TeV"]

# # print("All the available GPUs:\n",physical_devices)
# physical_devices = tf.config.experimental.list_physical_devices('GPU')#列出所有可见显卡
# if physical_devices:
#     gpu=physical_devices[int(sys.argv[4])]#显示第一块显卡
#     tf.config.experimental.set_memory_growth(gpu, True)#根据需要自动增长显存
#     tf.config.experimental.set_visible_devices(gpu, 'GPU')#只选择第一块
# 使用allow_growth option，刚一开始分配少量的GPU容量，然后按需慢慢的增加，
config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)
ConfigProto(log_device_placement=True) #优化计算资源
##################################################################################
ROOT.gSystem.Load('./analysis_cum_C.so') #加载动态链接库
output_folder="./model_output"
log_dir= os.path.join(output_folder,'NUM_POINTS_PbPb{}_PbP{}_Nch{}_{}_pT{}_{}_{}events_{}D_{}INITIAL_LR_0_2p4'.format(PbPbsigma,Pbpsigma,NUM_POINTS_LOW,NUM_POINTS_HIGH,PTLOWINP,PTHIGHINP,MIXEVENT,NUM_DIMENSION,INITIAL_LR))
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
DATA_DIR = os.path.join(os.path.dirname(__file__), "ampt_data")
# print(os.path.dirname(__file__))
print(DATA_DIR)
enc = OneHotEncoder(sparse=False)
def parse_dataset(NUM_POINTS_LOW,NUM_POINTS_HIGH,NUM_DIMENSION,MIXEVENT,PTLOWINP,PTHIGHINP): #NUM_POINTS,points_data
    train_points = []
    train_labels = []
    validation_points = []
    validation_labels = []
    eventindx = []        
    nevent = [[] for i in range(2)]
    # muti_Buffer = [[] for i in range(2)]
    # eventmuti = []
    # evaluate_points = []
    # evaluate_labels = []
    # nlabel = []
    # vn=[]
    # psi=[]
    class_map = {}
    # NchNevent = ROOT.std.vector[int]()
    # ROOT.gInterpreter.Declare(""" int NchNevent[50];""")
    VAL_SPLIT=0.7
    # folders = glob.glob(os.path.join(DATA_DIR, "*5p02TeV"))

    class_map = [folder.split("/")[-1] for folder in folders]
    print(class_map)
    # start = time.perf_counter()
    for i, folder in enumerate(folders):
        train_file = os.path.join(folder,"train")
        if not train_file:
            os._exit ("can not open list file")
        myana = ROOT.AMPTdata(train_file)
        nevent[i]=myana.GetNevent(NUM_POINTS_LOW,NUM_POINTS_HIGH)  #得到每个系统每个多重数下的事件数
    nsum=0
    NchNevent=[]
    for i in range(NUM_POINTS_HIGH-NUM_POINTS_LOW):
    #     ROOT.NchNevent[i]=(min(nevent[0][i],nevent[1][i]))
    #     nsum=nsum+ROOT.NchNevent[i]
    # if nsum>MAXEVENT:
    #     ratio=MAXEVENT/nsum
    #     for i in range(NUM_POINTS_HIGH-NUM_POINTS_LOW):
    #         ROOT.NchNevent[i]=int(ROOT.NchNevent[i]*ratio)
        NchNevent.append(min(nevent[0][i],nevent[1][i]))
        nsum=nsum+NchNevent[i]
    if nsum>MAXEVENT:
        ratio=MAXEVENT/nsum
        for i in range(NUM_POINTS_HIGH-NUM_POINTS_LOW):
           NchNevent[i]=int(NchNevent[i]*ratio) #得到每个系统每个多重数下的实际输入的事件数

    for i, folder in enumerate(folders):
        train_points_Buffer = [[] for i in range(2)]
        validation_points_Buffer = [[] for i in range(2)]
        # test_points_Buffer = [[] for i in range(2)]
        muti=[0,0]
        nevent = []
        ptrack = []
        eventindx_Buffer = [] 
        itrain=0
        ivalidation=0
        # itest=0
        train_file = os.path.join(folder,"train")#
        if not train_file:
            os._exit ("can not open list file")
        myana = ROOT.AMPTdata(train_file)
        # eventhead=0
        # eventend=0
        print(class_map[i])
        # for NUM_POINTS in tqdm(range(NUM_POINTS_HIGH-NUM_POINTS_LOW)):
        #     eventend+=ROOT.NchNevent[NUM_POINTS]
        #     train_index = int(ROOT.NchNevent[NUM_POINTS] * VAL_SPLIT)
        #     totalset=set(range(eventhead,eventend))
        #     trainset=set(random.sample(totalset, train_index))        
        #     validationset=totalset-trainset
        #     eventhead+=ROOT.NchNevent[NUM_POINTS]
        for NUM_POINTS in tqdm(range(NUM_POINTS_LOW,NUM_POINTS_HIGH)):
            NchNevents=NchNevent[NUM_POINTS-NUM_POINTS_LOW]
            tracks=myana.read(NUM_POINTS,NUM_DIMENSION,PTLOWINP,PTHIGHINP,NchNevents) #读取数据
            nevent=tracks.Idxevent #事件的index
            ptrack=tracks.mp #粒子动量信息列表
            # rtrack=tracks.mr #粒子位置信息列表
            # Psiv=tracks.Psi  #事件平面角
            # for num in tqdm(range(ROOT.NchNevent[NUM_POINTS])):
                
            # eventend+=NchNevent[NUM_POINTS]
            train_index = int(NchNevents * VAL_SPLIT)
            totalset=set(range(0,NchNevents))
            trainset=set(random.sample(totalset, train_index))        
            validationset=totalset-trainset
            for ievent in trainset:
                if len(ptrack[ievent])>4:
                    norm_point_cloud = np.array(ptrack[ievent],dtype='float32').reshape(-1,NUM_DIMENSION)
                    # norm_point_cloud += np.random.uniform( -0.001, 0.001,norm_point_cloud.shape)
                    # angle=Psiv[ievent]
                    # rotate_matrix = np.array([[math2.cos(angle),-math2.sin(angle)],[math2.sin(angle),math2.cos(angle)]],dtype='float32')
                    # norm_point_cloud=np.dot(norm_point_cloud,rotate_matrix)
                    norm_point_cloud = norm_point_cloud - np.mean(norm_point_cloud, axis=0)
                    # norm_point_cloud /= np.mean(np.linalg.norm(norm_point_cloud, axis=1))
                    # norm_point_cloud /= np.linalg.norm(norm_point_cloud, axis=1).reshape(-1,1)
                    train_points_Buffer[i].append(norm_point_cloud.tolist())
                    muti[i]=muti[i]+len(norm_point_cloud)
                    # psi_Buffer.append(Psiv)
                    itrain+=1
                    if itrain==MIXEVENT:
                        train_points.extend(train_points_Buffer[i].copy())
                        train_points_Buffer[i]=[]
                        train_labels.append(i)
                        itrain=0

            for ievent in validationset:
                if len(ptrack[ievent])>4:
                    eventindx_Buffer.append(nevent[ievent])
                    # muti_Buffer.append(ptrack[ievent].size()/NUM_DIMENSION)
                    norm_point_cloud = np.array(ptrack[ievent],dtype='float32').reshape(-1,NUM_DIMENSION)
                    # angle=Psiv[ievent]
                    # rotate_matrix = np.array([[math2.cos(angle),-math2.sin(angle)],[math2.sin(angle),math2.cos(angle)]],dtype='float32')
                    # norm_point_cloud=np.dot(norm_point_cloud,rotate_matrix)
                    norm_point_cloud = norm_point_cloud - np.mean(norm_point_cloud, axis=0)
                    # norm_point_cloud /= np.mean(np.linalg.norm(norm_point_cloud, axis=1))
                    # norm_point_cloud /= np.linalg.norm(norm_point_cloud, axis=1).reshape(-1,1)
                    validation_points_Buffer[i].append(norm_point_cloud.tolist())
                    ivalidation+=1
                    
                    muti[i]=muti[i]+len(norm_point_cloud)
                    if ivalidation==MIXEVENT:
                        validation_points.extend(validation_points_Buffer[i].copy())
                        validation_points_Buffer[i]=[]
                        validation_labels.append(i)
                        eventindx.append(eventindx_Buffer.copy())
                        eventindx_Buffer=[]
                        # eventmuti.append(muti_Buffer)
                        # muti_Buffer=[]
                        ivalidation=0
################################################################################################################################## 185-220

    
    if nsum>MAXEVENT:
        muti[i]=[muti[i]/MAXEVENT for i in range(2)]
    else:
        muti[i]=[muti[i]/nsum for i in range(2)]
    max_len = max((len(l) for l in train_points+validation_points))
    train_points=list(map(lambda l:l + [[0. for col in range(NUM_DIMENSION)]]*(max_len - len(l)), train_points)) #用0补齐
    # ivalidation=int(len(train_points)*(1-VAL_SPLIT))
    # validation_points = train_points[:ivalidation]
    # validation_labels = train_labels[:ivalidation]
    validation_points=list(map(lambda l:l + [[0. for col in range(NUM_DIMENSION)]]*(max_len - len(l)), validation_points))
    train_labels=keras.utils.to_categorical(train_labels, num_classes=NUM_CLASSES) #one-hot
    validation_labels=keras.utils.to_categorical(validation_labels, num_classes=NUM_CLASSES) #one-hot

    print("Num train point clouds:", len(train_points))
    print("Num train point cloud labels:", len(train_labels))
    print("Num val point clouds:", len(validation_points))
    print("Num val point cloud labels:", len(validation_labels))
    
    return (
        int(max_len),
        np.array(train_points,dtype='float32').reshape(-1,MIXEVENT,max_len,NUM_DIMENSION),
        np.array(validation_points,dtype='float32').reshape(-1,MIXEVENT,max_len,NUM_DIMENSION),
        np.array(train_labels,dtype='float32'),
        np.array(validation_labels,dtype='float32'),
        # np.array(evaluate_points),
        # enc.fit_transform(evaluate_labels_array.reshape((evaluate_labels_array.shape[0]),1)),
        # nlabel,
        class_map,
        np.array(eventindx),
        # np.array(eventmuti),
        # np.array(vn)
        muti
    )

"""
Set the number of points to sample and batch size and parse the dataset. This can take
~5minutes to complete.
"""
max_len, train_points, validation_points, train_labels, validation_labels, CLASS_MAP, eventindx, muti= parse_dataset(
    NUM_POINTS_LOW,NUM_POINTS_HIGH,NUM_DIMENSION,MIXEVENT,PTLOWINP,PTHIGHINP
)
tf.debugging.check_numerics(train_points, "train_points is producing nans!")
tf.debugging.check_numerics(validation_points, "validation_points is producing nans!")

print("Saving dataset !!!!!!!!!!!!!!!!!!!!!!!!!")

# The max_len is calculated by finding the maximum length of all the point clouds 
# in the training and validation datasets. 
# This maximum length is then used to pad all the point clouds to this length. 
# Padding is done by adding extra points with coordinates of zeros to the end of the point clouds 
# until they reach the maximum length.
save_maxlen= os.path.join(log_dir, "max_len")
np.save(save_maxlen,max_len)

save_train_points= os.path.join(log_dir, "train_points")
np.save(save_train_points,train_points)

save_train_labels= os.path.join(log_dir, "train_labels")
np.save(save_train_labels,train_labels)

save_validation_points= os.path.join(log_dir, "validation_points")
np.save(save_validation_points,validation_points)

save_validation_labels= os.path.join(log_dir, "validation_labels")
np.save(save_validation_labels,validation_labels)

save_class_map= os.path.join(log_dir, "class_map")
np.save(save_class_map,CLASS_MAP)


save_eventindx= os.path.join(log_dir, "eventindx")
np.save(save_eventindx,eventindx)

save_muti= os.path.join(log_dir, "muti")
np.save(save_muti,muti)

print("Saving dataset !!!!!!!!!!!!!!!!!!!!!!!!! done")

#print(evaluate_labels)evaluate_points, evaluate_labels, 
# train_labels_integers = np.argmax(train_labels, axis=1)
# class_weights = class_weight.compute_class_weight('balanced', classes=np.unique(train_labels_integers), y=train_labels_integers)
# d_class_weights = dict(enumerate(class_weights))
"""
Our data can now be read into a `tf.data.Dataset()` object. We set the shuffle buffer
size to the entire size of the dataset as prior to this the data is ordered by class.
Data augmentation is important when working with point cloud data. We create a
augmentation function to jitter and shuffle the train dataset.
"""


def augment(points, label):
    # jitter points
    # points += tf.random.uniform(points.shape, -0.00001, 0.00001, dtype=tf.float32)
    # shuffle points
    # print(points[0])
    pointlist=[tf.random.shuffle(points[i]) for i in range(MIXEVENT)]
    return pointlist, label

train_dataset = tf.data.Dataset.from_tensor_slices((train_points, train_labels))
validation_dataset = tf.data.Dataset.from_tensor_slices((validation_points, validation_labels))

train_dataset = train_dataset.shuffle(len(train_points)).map(augment).batch(BATCH_SIZE)
validation_dataset = validation_dataset.shuffle(len(validation_points)).batch(BATCH_SIZE)

"""
### Build a model

Each convolution and fully-connected layer (with exception for end layers) consits of
Convolution / Dense -> Batch Normalization -> sigmoid Activation.
"""

