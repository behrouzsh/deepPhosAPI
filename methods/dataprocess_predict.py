import csv
import io
import numpy as np
import keras.utils.np_utils as kutils
import re
import urllib
import json

#input format   label,proteinName, postion,sites, shortsequence,
#input must be a .csv file
#positive_position_file_name is an csv file

def getMatrixInputFromJsonPrev(inputText, organism, window_size=51, empty_aa = '*'):
    # input format  proteinName, postion, shortsequence,
    prot = []  # list of protein name
    pos = []  # list of position with protein name
    full_names = []
    names = []
    rawseq = []
    # all_label = []
    protein_info = {}

    short_seqs = []
    half_len = int((window_size - 1) / 2)

    # with open(positive_position_file_name, 'r') as rf:
    #     reader = csv.reader(rf)
    inputItems = inputText.split(",")
    # f.stream.seek(0)
    # stream = io.StringIO(f.stream.read().decode("UTF8"), newline=None)
    # reader = csv.reader(stream)
    print("==========================")
    for row in inputItems:
        print(row)
        modification = re.match(r"[^[]*\[([^]]*)\]", row).groups()[0]
        #print(re.match(r"[^[]*\[([^]]*)\]", row))
        print(modification)
        sseq = ""
        position = 0
        center = 0
        protein = row.split("[")[0]
        position = int(modification.split("@")[1])
        input_center = str(modification.split("@")[0])[-1]

        if protein not in protein_info.keys():
            url2 = "http://eh3.uc.edu/pinet/api/uniprotdb/organism/{}/accession/{}".format(organism, protein)
            try:
                with urllib.request.urlopen(url2) as url:
                    response2 = json.loads(url.read().decode())
                    print(response2)

                    sseq = response2["sequence"]
                    protein_info[protein] = response2
                    center = sseq[position - 1]
                    gene = response2["primary_gene_name"][0]


            except:
                protein_info[protein]["sequence"] = "*"

                protein_info[protein]["gene"] = "*"


        else:
            sseq = protein_info[protein]["sequence"]
            gene = protein_info[protein]["primary_gene_name"][0]
            center = sseq[position - 1]




        # response = urllib.urlopen(url)
        #
        # data = json.loads(response.read())






        #sseq = row[2]
        #position = int(row[1])
        #center = sseq[position-1]

        print("------------------------------")
        print(protein)
        print(position)
        print(sseq)
        print(center)
        print(input_center)



        sites = 'S', 'T', 'Y'
        if center in sites:
            prot.append(protein)
            pos.append(position)
            full_name = protein+ "[p" + center + "@" + str(position) + "](" + gene + ")"
            name = protein + "[ph" + center + "@" + str(position) + "]"
            names.append(name)
            full_names.append(full_name)
            rawseq.append(sseq)
            # print rawseq

            #short seq
            if position - half_len > 0:
                start = position - half_len
                left_seq = sseq[start - 1:position - 1]
            else:
                left_seq = sseq[0:position - 1]

            end = len(sseq)
            if position + half_len < end:
                end = position + half_len
            right_seq = sseq[position:end]

            if len(left_seq) < half_len:
                nb_lack = half_len - len(left_seq)
                left_seq = ''.join([empty_aa for count in range(nb_lack)]) + left_seq

            if len(right_seq) < half_len:
                nb_lack = half_len - len(right_seq)
                right_seq = right_seq + ''.join([empty_aa for count in range(nb_lack)])
            shortseq = left_seq + center + right_seq
            short_seqs.append(shortseq)
            # coding = one_hot_concat(shortseq)
            # all_codings.append(coding)

    all_label = [0] *5 + [1]*(len(short_seqs)-5)
    targetY = kutils.to_categorical(all_label)

    ONE_HOT_SIZE = 21
    # _aminos = 'ACDEFGHIKLMNPQRSTVWY*'
    letterDict = {}
    letterDict["A"] = 0
    letterDict["C"] = 1
    letterDict["D"] = 2
    letterDict["E"] = 3
    letterDict["F"] = 4
    letterDict["G"] = 5
    letterDict["H"] = 6
    letterDict["I"] = 7
    letterDict["K"] = 8
    letterDict["L"] = 9
    letterDict["M"] = 10
    letterDict["N"] = 11
    letterDict["P"] = 12
    letterDict["Q"] = 13
    letterDict["R"] = 14
    letterDict["S"] = 15
    letterDict["T"] = 16
    letterDict["V"] = 17
    letterDict["W"] = 18
    letterDict["Y"] = 19
    letterDict["*"] = 20

    # print len(short_seqs)
    Matr = np.zeros((len(short_seqs), window_size, ONE_HOT_SIZE))
    samplenumber = 0
    for seq in short_seqs:
        AANo = 0
        for AA in seq:
            index = letterDict[AA]
            # print index
            Matr[samplenumber][AANo][index] = 1
            # print samplenumber
            AANo = AANo+1
        samplenumber = samplenumber + 1

    return Matr, targetY, prot, pos, full_names, names


def getMatrixInputFromJson(prot, pos, full_names, names, short_seqs,window_size=51, empty_aa = '*'):
    # input format  proteinName, postion, shortsequence,


    all_label = [0] *5 + [1]*(len(short_seqs)-5)
    targetY = kutils.to_categorical(all_label)

    ONE_HOT_SIZE = 21
    # _aminos = 'ACDEFGHIKLMNPQRSTVWY*'
    letterDict = {}
    letterDict["A"] = 0
    letterDict["C"] = 1
    letterDict["D"] = 2
    letterDict["E"] = 3
    letterDict["F"] = 4
    letterDict["G"] = 5
    letterDict["H"] = 6
    letterDict["I"] = 7
    letterDict["K"] = 8
    letterDict["L"] = 9
    letterDict["M"] = 10
    letterDict["N"] = 11
    letterDict["P"] = 12
    letterDict["Q"] = 13
    letterDict["R"] = 14
    letterDict["S"] = 15
    letterDict["T"] = 16
    letterDict["V"] = 17
    letterDict["W"] = 18
    letterDict["Y"] = 19
    letterDict["*"] = 20

    # print len(short_seqs)
    Matr = np.zeros((len(short_seqs), window_size, ONE_HOT_SIZE))
    samplenumber = 0
    for seq in short_seqs:
        AANo = 0
        for AA in seq:
            index = letterDict[AA]
            # print index
            Matr[samplenumber][AANo][index] = 1
            # print samplenumber
            AANo = AANo+1
        samplenumber = samplenumber + 1

    return Matr, targetY, prot, pos, full_names, names

def getMatrixInputFromFile(f,sites, window_size=51, empty_aa = '*'):
    # input format  proteinName, postion, shortsequence,
    prot = []  # list of protein name
    pos = []  # list of position with protein name
    rawseq = []
    # all_label = []

    short_seqs = []
    half_len = int((window_size - 1) / 2)

    # with open(positive_position_file_name, 'r') as rf:
    #     reader = csv.reader(rf)
    f.stream.seek(0)
    stream = io.StringIO(f.stream.read().decode("UTF8"), newline=None)
    reader = csv.reader(stream)
    print("==========================")
    for row in reader:
        print(row)
        sseq = row[2]
        position = int(row[1])
        center = sseq[position-1]

        if center in sites:
            prot.append(row[0])
            pos.append(row[1])
            rawseq.append(sseq)
            # print rawseq

            #short seq
            if position - half_len > 0:
                start = position - half_len
                left_seq = sseq[start - 1:position - 1]
            else:
                left_seq = sseq[0:position - 1]

            end = len(sseq)
            if position + half_len < end:
                end = position + half_len
            right_seq = sseq[position:end]

            if len(left_seq) < half_len:
                nb_lack = half_len - len(left_seq)
                left_seq = ''.join([empty_aa for count in range(nb_lack)]) + left_seq

            if len(right_seq) < half_len:
                nb_lack = half_len - len(right_seq)
                right_seq = right_seq + ''.join([empty_aa for count in range(nb_lack)])
            shortseq = left_seq + center + right_seq
            short_seqs.append(shortseq)
            # coding = one_hot_concat(shortseq)
            # all_codings.append(coding)

    all_label = [0] *5 + [1]*(len(short_seqs)-5)
    targetY = kutils.to_categorical(all_label)

    ONE_HOT_SIZE = 21
    # _aminos = 'ACDEFGHIKLMNPQRSTVWY*'
    letterDict = {}
    letterDict["A"] = 0
    letterDict["C"] = 1
    letterDict["D"] = 2
    letterDict["E"] = 3
    letterDict["F"] = 4
    letterDict["G"] = 5
    letterDict["H"] = 6
    letterDict["I"] = 7
    letterDict["K"] = 8
    letterDict["L"] = 9
    letterDict["M"] = 10
    letterDict["N"] = 11
    letterDict["P"] = 12
    letterDict["Q"] = 13
    letterDict["R"] = 14
    letterDict["S"] = 15
    letterDict["T"] = 16
    letterDict["V"] = 17
    letterDict["W"] = 18
    letterDict["Y"] = 19
    letterDict["*"] = 20

    # print len(short_seqs)
    Matr = np.zeros((len(short_seqs), window_size, ONE_HOT_SIZE))
    samplenumber = 0
    for seq in short_seqs:
        AANo = 0
        for AA in seq:
            index = letterDict[AA]
            # print index
            Matr[samplenumber][AANo][index] = 1
            # print samplenumber
            AANo = AANo+1
        samplenumber = samplenumber + 1

    return Matr, targetY, prot, pos

def getMatrixInput(positive_position_file_name,sites, window_size=51, empty_aa = '*'):
    # input format  proteinName, postion, shortsequence,
    prot = []  # list of protein name
    pos = []  # list of position with protein name
    rawseq = []
    # all_label = []

    short_seqs = []
    half_len = int((window_size - 1) / 2)

    with open(positive_position_file_name, 'r') as rf:
        reader = csv.reader(rf)
        print("In reading input files")
        for row in reader:
            print(row)
            print("-------------------")
            sseq = row[2]
            position = int(row[1])
            center = sseq[position-1]
            if center in sites:
                prot.append(row[0])
                pos.append(row[1])
                rawseq.append(sseq)
                # print rawseq

                #short seq
                if position - half_len > 0:
                    start = position - half_len
                    left_seq = sseq[start - 1:position - 1]
                else:
                    left_seq = sseq[0:position - 1]

                end = len(sseq)
                if position + half_len < end:
                    end = position + half_len
                right_seq = sseq[position:end]

                if len(left_seq) < half_len:
                    nb_lack = half_len - len(left_seq)
                    left_seq = ''.join([empty_aa for count in range(nb_lack)]) + left_seq

                if len(right_seq) < half_len:
                    nb_lack = half_len - len(right_seq)
                    right_seq = right_seq + ''.join([empty_aa for count in range(nb_lack)])
                shortseq = left_seq + center + right_seq
                short_seqs.append(shortseq)
                # coding = one_hot_concat(shortseq)
                # all_codings.append(coding)

        all_label = [0] *5 + [1]*(len(short_seqs)-5)
        targetY = kutils.to_categorical(all_label)

        ONE_HOT_SIZE = 21
        # _aminos = 'ACDEFGHIKLMNPQRSTVWY*'
        letterDict = {}
        letterDict["A"] = 0
        letterDict["C"] = 1
        letterDict["D"] = 2
        letterDict["E"] = 3
        letterDict["F"] = 4
        letterDict["G"] = 5
        letterDict["H"] = 6
        letterDict["I"] = 7
        letterDict["K"] = 8
        letterDict["L"] = 9
        letterDict["M"] = 10
        letterDict["N"] = 11
        letterDict["P"] = 12
        letterDict["Q"] = 13
        letterDict["R"] = 14
        letterDict["S"] = 15
        letterDict["T"] = 16
        letterDict["V"] = 17
        letterDict["W"] = 18
        letterDict["Y"] = 19
        letterDict["*"] = 20

        # print len(short_seqs)
        Matr = np.zeros((len(short_seqs), window_size, ONE_HOT_SIZE))
        samplenumber = 0
        for seq in short_seqs:
            AANo = 0
            for AA in seq:
                index = letterDict[AA]
                # print index
                Matr[samplenumber][AANo][index] = 1
                # print samplenumber
                AANo = AANo+1
            samplenumber = samplenumber + 1

    return Matr, targetY, prot, pos












