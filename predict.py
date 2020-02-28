import functools
import itertools
import os
import random
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import csv
import matplotlib.pyplot as plt
import numpy as np
from flask import Flask, make_response, send_from_directory, request, render_template, url_for, redirect
from sklearn import metrics
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, KFold, cross_val_score

from keras.layers import Dense, Activation, Flatten, Dropout, Reshape
from keras.layers import Conv1D,Conv2D, MaxPooling2D
from keras.models import Sequential,Model
from keras.utils.np_utils import to_categorical
from keras import optimizers
from keras.optimizers import Adam,SGD
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l2
import copy

def predict_for_deepphos(file_name,sites,predictFrame = 'general',
                         hierarchy=None, kinase=None):
    '''

    :param train_file_name: input of your prdict file
                            it must be a .csv file and theinput format  is proteinName, postion,sites, shortseq
    :param sites: the sites predict: site = 'S','T' OR 'Y'
    :param predictFrame: 'general' or 'kinase'
    :param hierarchy: if predictFrame is kinse: you must input the hierarchy:
            group,family,subfamily,kinase to choose corresponding model
    :param kinase: kinase name
    :return:
     a file with the score
    '''
    print("----------------------------")
    print("in predict_for_deepphos")
    win1 = 51
    win2 = 33
    win3 = 15
    from methods.dataprocess_predict import getMatrixInput

    [X_test1,y_test,ids,position] = getMatrixInput(file_name, sites, win1)
    [X_test2,_,_,_] = getMatrixInput(file_name, sites, win2)
    [X_test3,_,_,_]  = getMatrixInput(file_name, sites, win3)
    print("X_test1 ----------------------------")
    print(X_test1)
    print("X_test2 ----------------------------")
    print(X_test2)
    print("X_test3 ----------------------------")
    print(X_test3)
    print("----------------------------")
    #     print X_test1.shape
#     print len(position)

    from methods.model_n import model_net
    model = model_net(X_test1, X_test2, X_test3, y_test,nb_epoch = 0)

    #load model weight
    if predictFrame == 'general':

        if site == ('S','T'):
            outputfile = 'general_S_T'
            model_weight = './models/model_general_S,T.h5'
        if site == 'Y':
            outputfile = 'general_Y'
            model_weight = './models/model_general_Y.h5'


    if predictFrame == 'kinase':
        outputfile = 'kinase_{:s}_{:s}'.format(hierarchy, kinase)
        model_weight = './models/model_{:s}_{:s}.h5'.format(hierarchy, kinase)
    print(model_weight)
    print("-----------------------------------------------------")
    model.load_weights(model_weight)
    predictions_t = model.predict([X_test1, X_test2, X_test3])
    results_ST = np.column_stack((ids, position,predictions_t[:, 1]))

    result = pd.DataFrame(results_ST)
    result.to_csv(outputfile + "prediction_phosphorylation.txt", index=False, header=None, sep='\t',
                  quoting=csv.QUOTE_NONNUMERIC)

def predict_for_deepphos_from_json(input, organism):
    '''

    :param train_file_name: input of your prdict file
                            it must be a .csv file and theinput format  is proteinName, postion,sites, shortseq
    :param sites: the sites predict: site = 'S','T' OR 'Y'
    :param predictFrame: 'general' or 'kinase'
    :param hierarchy: if predictFrame is kinse: you must input the hierarchy:
            group,family,subfamily,kinase to choose corresponding model
    :param kinase: kinase name
    :return:
     a file with the score
    '''
    print("----------------------------")
    print("in predict_for_deepphos")
    win1 = 51
    win2 = 33
    win3 = 15


    from methods.dataprocess_predict import getMatrixInputFromJson

    print('running X_test1 ----------------------------')
    [X_test1,y_test,ids,position,full_names,names] = getMatrixInputFromJson(input, organism, win1)
    print('running X_test2 ----------------------------')
    [X_test2,_,_,_,_,_] = getMatrixInputFromJson(input, organism, win2)
    print('running X_test3 ----------------------------')
    [X_test3,_,_,_,_,_]  = getMatrixInputFromJson(input, organism, win3)

    print('----------------------------')
    result_json = {}
    #     print X_test1.shape
#     print len(position)

    from methods.model_n import model_net
    model = model_net(X_test1, X_test2, X_test3, y_test,nb_epoch = 0)

    #load model weight
    # if predictFrame == 'general':
    #
    #     if site == ('S','T'):
    #         outputfile = 'general_S_T'
    #         model_weight = './models/model_general_S,T.h5'
    #     if site == 'Y':
    #         outputfile = 'general_Y'
    #         model_weight = './models/model_general_Y.h5'

    #for kin in ['family_CDK']:

    for kin in ['family_CDK','family_CK2','family_MAPK','family_PKC','family_Src',
                   'group_AGC','group_Atypical','group_CAMK','group_CMGC','group_TK',
                   'kinase_CDC2','kinase_CK2a1','kinase_PKACa','kinase_PKCa','kinase_SRC',
                   'subfamily_CDC2','subfamily_CDK2','subfamily_ERK1','subfamily_PKCa']:
    # if predictFrame == 'kinase':
    #     outputfile = 'kinase_{:s}_{:s}'.format(hierarchy, kinase)
    #     model_weight = './models/model_{:s}_{:s}.h5'.format(hierarchy, kinase)
        print(kin)
        outputfile = 'kinase_{:s}'.format(kin)

        model_weight = './models/model_{:s}.h5'.format(kin)
        print("-----------------------------------------------------")
        print(model_weight)

        try:
            model.load_weights(model_weight)
            predictions_t = model.predict([X_test1, X_test2, X_test3])

            results_ST = np.column_stack((ids, names, full_names ,position,predictions_t[:, 1]))

            result = pd.DataFrame(results_ST, columns= ['target', 'name', 'full_name', 'position', 'score'])
            print(result.to_json(orient='records'))
            result_json[kin] = result.to_json(orient='records')
        #print(result_json)





            # result.to_csv(outputfile + "_api_prediction_phosphorylation.txt", index=False, header=None, sep='\t',
            #               quoting=csv.QUOTE_NONNUMERIC)
        except:
            print("there was an error")
            result_json[kin] = json.dumps([{'target':"", 'full_name':"", 'position':"", 'score':""}])

    print(result_json)
    return json.dumps(result_json)

def predict_for_deepphos_from_file(file,sites):
    '''

    :param train_file_name: input of your prdict file
                            it must be a .csv file and theinput format  is proteinName, postion,sites, shortseq
    :param sites: the sites predict: site = 'S','T' OR 'Y'
    :param predictFrame: 'general' or 'kinase'
    :param hierarchy: if predictFrame is kinse: you must input the hierarchy:
            group,family,subfamily,kinase to choose corresponding model
    :param kinase: kinase name
    :return:
     a file with the score
    '''
    print("----------------------------")
    print("in predict_for_deepphos")
    win1 = 51
    win2 = 33
    win3 = 15
    from methods.dataprocess_predict import getMatrixInputFromFile
    print("running X_test1 ----------------------------")
    [X_test1,y_test,ids,position] = getMatrixInputFromFile(file, sites, win1)
    print("running X_test2 ----------------------------")
    [X_test2,_,_,_] = getMatrixInputFromFile(file, sites, win2)
    print("running X_test3 ----------------------------")
    [X_test3,_,_,_]  = getMatrixInputFromFile(file, sites, win3)

    result_json = {}
    #     print X_test1.shape
#     print len(position)

    from methods.model_n import model_net
    try:
        model = model_net(X_test1, X_test2, X_test3, y_test,nb_epoch = 0)

        #load model weight
        # if predictFrame == 'general':
        #
        #     if site == ('S','T'):
        #         outputfile = 'general_S_T'
        #         model_weight = './models/model_general_S,T.h5'
        #     if site == 'Y':
        #         outputfile = 'general_Y'
        #         model_weight = './models/model_general_Y.h5'

        for kin in ["family_CDK"]:
        # for kin in ["family_CDK","family_CK2","family_MAPK","family_PKC","family_Src",
        #                "group_AGC","group_Atypical","group_CAMK","group_CMGC","group_TK",
        #                "kinase_CDC2","kinase_CK2a1","kinase_PKACa","kinase_PKCa","kinase_SRC",
        #                "subfamily_CDC2","subfamily_CDK2","subfamily_ERK1","subfamily_PKCa"]:
        # if predictFrame == 'kinase':
        #     outputfile = 'kinase_{:s}_{:s}'.format(hierarchy, kinase)
        #     model_weight = './models/model_{:s}_{:s}.h5'.format(hierarchy, kinase)
            print(kin)
            outputfile = 'kinase_{:s}'.format(kin)


            model_weight = './models/model_{:s}.h5'.format(kin)
            print("-----------------------------------------------------")
            print(model_weight)

            # try:

            model.load_weights(model_weight)
            predictions_t = model.predict([X_test1, X_test2, X_test3])

            results_ST = np.column_stack((ids, position,predictions_t[:, 1]))

            result = pd.DataFrame(results_ST, columns= ['target', 'position', 'score'])
            print(result)
            #print(result.to_json(orient='records'))
            result_json[kin] = result.to_json(orient='records')
                #print(result_json)
            #keras.backend.tensorflow_backend.clear_session()





                # result.to_csv(outputfile + "_api_prediction_phosphorylation.txt", index=False, header=None, sep='\t',
                #               quoting=csv.QUOTE_NONNUMERIC)
    except:
         print("there was an error")
         result_json[kin] = []

    print(result_json)
    return result_json

def runserver():
    #ensure_root()

    # for pyc_file_path in get_valid_pyc_files():
    #     print 'removing', pyc_file_path
    #     os.remove(pyc_file_path)
    port = int(os.environ.get('PORT', 4000))
    app.run(host='0.0.0.0', port=port)

app = Flask(__name__)


@app.route("/")
@app.route("/home")
@app.route("/about")
@app.route("/clustergram")
@app.route("/bootstrap")
def basic_pages():
    # train_file_name = 'test data.csv'
    #
    # site = 'S', 'T'
    # #
    # predict_for_deepphos(train_file_name, site, predictFrame='kinase',
    #                      hierarchy='group', kinase='AGC')
    return make_response(open('static/index.html').read())




@app.route('/handle_form', methods=['POST', 'GET'])
def handle_form():
    print("Posted file: {}".format(request.files['file']))
    res = {}
    #from methods.dataprocess_predict import convertFileToReader
    file = request.files['file']
    #reader = convertFileToReader(file)
    site = 'S', 'T', 'Y'
    res = predict_for_deepphos_from_file(file, site)
    print("=============")
    print("=============")
    print(res)
    return str(res)

@app.route('/api/predict/organism/<organism>/<input>', methods=['GET'])
def predict_api(organism,input):
    print("input:")
    print(input)
    res = {}
    #from methods.dataprocess_predict import convertFileToReader
    #file = request.files['file']
    #reader = convertFileToReader(file)
    site = 'S', 'T'
    res = predict_for_deepphos_from_json(input, organism)
    print("=============")
    print("=============")
    print(res)
    return str(res)

if __name__ == '__main__':
    train_file_name = 'testdata2.csv'

    site = 'S', 'T'
    #
    # print("second iter ++++++++++++++++////////++++++++++++")
    # predict_for_deepphos(train_file_name, site, predictFrame='kinase',
    #                      hierarchy='group', kinase='CAMK')
    # print("first iter ++++++++++++++++////////++++++++++++")
    # predict_for_deepphos(train_file_name, site, predictFrame='kinase',
    #                      hierarchy='group', kinase='AGC')



    # threaded should be false to reload the weights in
    app.run(debug=False, threaded=False)
    #app.run(debug=True, use_reloader=False)




