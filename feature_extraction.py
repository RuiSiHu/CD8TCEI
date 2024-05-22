#!/usr/bin/env python
#_*_coding:utf-8_*_

from collections import Counter
import numpy as np
import re
import math
from collections import Counter

def Count_1(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum

def Count_2(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code


def get_features(fastas):
    def AAC():
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        encodings = []
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            count = Counter(sequence)
            for key in count:
                count[key] = count[key] / len(sequence)
            code = []
            for aa in AA:
                code.append(count[aa])
            encodings.append(code)
        return encodings

    def ASDC():
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)

        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            sum = 0
            pair_dict = {}
            for pair in aaPairs:
                 pair_dict[pair] = 0
            for j in range(len(sequence)):
                for k in range(j + 1, len(sequence)):
                     if sequence[j] in AA and sequence[k] in AA:
                        pair_dict[sequence[j] + sequence[k]] += 1
                        sum += 1
            for pair in aaPairs:
                code.append(pair_dict[pair] / sum)
            encodings.append(code)
        return encodings

    def CTDC():
        group1 = {
            'hydrophobicity_PRAM900101': 'RKEDQN',
            'hydrophobicity_ARGP820101': 'QSTNGDE',
            'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
            'hydrophobicity_PONP930101': 'KPDESNQT',
            'hydrophobicity_CASG920101': 'KDEQPSRNTG',
            'hydrophobicity_ENGD860101': 'RDKENQHYP',
            'hydrophobicity_FASG890101': 'KERSQD',
            'normwaalsvolume': 'GASTPDC',
            'polarity': 'LIFWCMVY',
            'polarizability': 'GASDT',
            'charge': 'KR',
            'secondarystruct': 'EALMQKRH',
            'solventaccess': 'ALFCGIVW'
        }
        group2 = {
            'hydrophobicity_PRAM900101': 'GASTPHY',
            'hydrophobicity_ARGP820101': 'RAHCKMV',
            'hydrophobicity_ZIMJ680101': 'HMCKV',
            'hydrophobicity_PONP930101': 'GRHA',
            'hydrophobicity_CASG920101': 'AHYMLV',
            'hydrophobicity_ENGD860101': 'SGTAW',
            'hydrophobicity_FASG890101': 'NTPG',
            'normwaalsvolume': 'NVEQIL',
            'polarity': 'PATGS',
            'polarizability': 'CPNVEQIL',
            'charge': 'ANCQGHILMFPSTWYV',
            'secondarystruct': 'VIYCWFT',
            'solventaccess': 'RKQEND'
        }
        group3 = {
            'hydrophobicity_PRAM900101': 'CLVIMFW',
            'hydrophobicity_ARGP820101': 'LYPFIW',
            'hydrophobicity_ZIMJ680101': 'LPFYI',
            'hydrophobicity_PONP930101': 'YMFWLCVI',
            'hydrophobicity_CASG920101': 'FIWC',
            'hydrophobicity_ENGD860101': 'CVLIMF',
            'hydrophobicity_FASG890101': 'AYHWVMFLIC',
            'normwaalsvolume': 'MHKFRYW',
            'polarity': 'HQRKNED',
            'polarizability': 'KMHFRYW',
            'charge': 'DE',
            'secondarystruct': 'GNPSD',
            'solventaccess': 'MSPTHY'
        }

        groups = [group1, group2, group3]
        property = (
            'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
            'hydrophobicity_PONP930101',
            'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
            'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

        encodings = []
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            for p in property:
                c1 = Count_1(group1[p], sequence) / len(sequence)
                c2 = Count_1(group2[p], sequence) / len(sequence)
                c3 = 1 - c1 - c2
                code = code + [c1, c2, c3]
            encodings.append(code)
        return encodings

    def CTDT():
        group1 = {
            'hydrophobicity_PRAM900101': 'RKEDQN',
            'hydrophobicity_ARGP820101': 'QSTNGDE',
            'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
            'hydrophobicity_PONP930101': 'KPDESNQT',
            'hydrophobicity_CASG920101': 'KDEQPSRNTG',
            'hydrophobicity_ENGD860101': 'RDKENQHYP',
            'hydrophobicity_FASG890101': 'KERSQD',
            'normwaalsvolume': 'GASTPDC',
            'polarity': 'LIFWCMVY',
            'polarizability': 'GASDT',
            'charge': 'KR',
            'secondarystruct': 'EALMQKRH',
            'solventaccess': 'ALFCGIVW'
        }
        group2 = {
            'hydrophobicity_PRAM900101': 'GASTPHY',
            'hydrophobicity_ARGP820101': 'RAHCKMV',
            'hydrophobicity_ZIMJ680101': 'HMCKV',
            'hydrophobicity_PONP930101': 'GRHA',
            'hydrophobicity_CASG920101': 'AHYMLV',
            'hydrophobicity_ENGD860101': 'SGTAW',
            'hydrophobicity_FASG890101': 'NTPG',
            'normwaalsvolume': 'NVEQIL',
            'polarity': 'PATGS',
            'polarizability': 'CPNVEQIL',
            'charge': 'ANCQGHILMFPSTWYV',
            'secondarystruct': 'VIYCWFT',
            'solventaccess': 'RKQEND'
        }
        group3 = {
            'hydrophobicity_PRAM900101': 'CLVIMFW',
            'hydrophobicity_ARGP820101': 'LYPFIW',
            'hydrophobicity_ZIMJ680101': 'LPFYI',
            'hydrophobicity_PONP930101': 'YMFWLCVI',
            'hydrophobicity_CASG920101': 'FIWC',
            'hydrophobicity_ENGD860101': 'CVLIMF',
            'hydrophobicity_FASG890101': 'AYHWVMFLIC',
            'normwaalsvolume': 'MHKFRYW',
            'polarity': 'HQRKNED',
            'polarizability': 'KMHFRYW',
            'charge': 'DE',
            'secondarystruct': 'GNPSD',
            'solventaccess': 'MSPTHY'
        }

        groups = [group1, group2, group3]
        property = (
            'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
            'hydrophobicity_PONP930101',
            'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
            'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

        encodings = []

        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
            for p in property:
                c1221, c1331, c2332 = 0, 0, 0
                for pair in aaPair:
                    if (pair[0] in group1[p] and pair[1] in group2[p]) or (
                            pair[0] in group2[p] and pair[1] in group1[p]):
                        c1221 = c1221 + 1
                        continue
                    if (pair[0] in group1[p] and pair[1] in group3[p]) or (
                            pair[0] in group3[p] and pair[1] in group1[p]):
                        c1331 = c1331 + 1
                        continue
                    if (pair[0] in group2[p] and pair[1] in group3[p]) or (
                            pair[0] in group3[p] and pair[1] in group2[p]):
                        c2332 = c2332 + 1
                code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
            encodings.append(code)
        return encodings

    def CTDD():
        group1 = {
            'hydrophobicity_PRAM900101': 'RKEDQN',
            'hydrophobicity_ARGP820101': 'QSTNGDE',
            'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
            'hydrophobicity_PONP930101': 'KPDESNQT',
            'hydrophobicity_CASG920101': 'KDEQPSRNTG',
            'hydrophobicity_ENGD860101': 'RDKENQHYP',
            'hydrophobicity_FASG890101': 'KERSQD',
            'normwaalsvolume': 'GASTPDC',
            'polarity': 'LIFWCMVY',
            'polarizability': 'GASDT',
            'charge': 'KR',
            'secondarystruct': 'EALMQKRH',
            'solventaccess': 'ALFCGIVW'
        }
        group2 = {
            'hydrophobicity_PRAM900101': 'GASTPHY',
            'hydrophobicity_ARGP820101': 'RAHCKMV',
            'hydrophobicity_ZIMJ680101': 'HMCKV',
            'hydrophobicity_PONP930101': 'GRHA',
            'hydrophobicity_CASG920101': 'AHYMLV',
            'hydrophobicity_ENGD860101': 'SGTAW',
            'hydrophobicity_FASG890101': 'NTPG',
            'normwaalsvolume': 'NVEQIL',
            'polarity': 'PATGS',
            'polarizability': 'CPNVEQIL',
            'charge': 'ANCQGHILMFPSTWYV',
            'secondarystruct': 'VIYCWFT',
            'solventaccess': 'RKQEND'
        }
        group3 = {
            'hydrophobicity_PRAM900101': 'CLVIMFW',
            'hydrophobicity_ARGP820101': 'LYPFIW',
            'hydrophobicity_ZIMJ680101': 'LPFYI',
            'hydrophobicity_PONP930101': 'YMFWLCVI',
            'hydrophobicity_CASG920101': 'FIWC',
            'hydrophobicity_ENGD860101': 'CVLIMF',
            'hydrophobicity_FASG890101': 'AYHWVMFLIC',
            'normwaalsvolume': 'MHKFRYW',
            'polarity': 'HQRKNED',
            'polarizability': 'KMHFRYW',
            'charge': 'DE',
            'secondarystruct': 'GNPSD',
            'solventaccess': 'MSPTHY'
        }

        groups = [group1, group2, group3]
        property = (
            'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
            'hydrophobicity_PONP930101',
            'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
            'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

        encodings = []

        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            for p in property:
                code = code + Count_2(group1[p], sequence) + Count_2(group2[p], sequence) + Count_2(group3[p], sequence)
            encodings.append(code)
        return encodings

    def GDPC():
        group = {
            'alphaticr': 'GAVLMI',
            'aromatic': 'FYW',
            'postivecharger': 'KRH',
            'negativecharger': 'DE',
            'uncharger': 'STCPNQ'
        }

        groupKey = group.keys()
        baseNum = len(groupKey)
        dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

        index = {}
        for key in groupKey:
            for aa in group[key]:
                index[aa] = key

        encodings = []


        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])

            code = []
            myDict = {}
            for t in dipeptide:
                myDict[t] = 0

            sum = 0
            for j in range(len(sequence) - 2 + 1):
                myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                    sequence[j + 1]]] + 1
                sum = sum + 1

            if sum == 0:
                for t in dipeptide:
                    code.append(0)
            else:
                for t in dipeptide:
                    code.append(myDict[t] / sum)
            encodings.append(code)

        return encodings

    def GTPC():
        group = {
            'alphaticr': 'GAVLMI',
            'aromatic': 'FYW',
            'postivecharger': 'KRH',
            'negativecharger': 'DE',
            'uncharger': 'STCPNQ'
        }

        groupKey = group.keys()
        baseNum = len(groupKey)
        triple = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

        index = {}
        for key in groupKey:
            for aa in group[key]:
                index[aa] = key

        encodings = []

        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])

            code = []
            myDict = {}
            for t in triple:
                myDict[t] = 0

            sum = 0
            for j in range(len(sequence) - 3 + 1):
                myDict[index[sequence[j]] + '.' + index[sequence[j + 1]] + '.' + index[sequence[j + 2]]] = myDict[index[
                                                                                                                      sequence[
                                                                                                                          j]] + '.' +
                                                                                                                  index[
                                                                                                                      sequence[
                                                                                                                          j + 1]] + '.' +
                                                                                                                  index[
                                                                                                                      sequence[
                                                                                                                          j + 2]]] + 1
                sum = sum + 1

            if sum == 0:
                for t in triple:
                    code.append(0)
            else:
                for t in triple:
                    code.append(myDict[t] / sum)
            encodings.append(code)

        return encodings

    print('Feature extraction...')
    encoding = []

    encoding.append(AAC())
    encoding.append(ASDC())
    encoding.append(CTDC())
    encoding.append(CTDT())
    encoding.append(CTDD())
    encoding.append(GDPC())
    encoding.append(GTPC())

    return np.column_stack(encoding)

