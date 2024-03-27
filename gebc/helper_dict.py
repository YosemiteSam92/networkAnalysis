from subprocess import list2cmdline


def printDict(dict):
    print()
    for key, value in dict.items():
        print(key, value, "\n")


def dictWithLists(dict, key, value):
    if key not in dict:
        dict[key] = []
    dict[key].append(value)


def sortListsBasedOnValuesFromFirstList(listOfLists, reverse=False):
    numLists = len(listOfLists)
    sortedIterator = sorted(
        zip(*listOfLists),
        key=lambda cartesianProduct: cartesianProduct[0],
        reverse=reverse
        )
    # return list of sorted lists
    # x is an entry of the Cartesian product of the sorted lists
    # the second for expression is evaluated first
    return [[x[index] for x in sortedIterator] for index in range(numLists)]


def sortDictByValue(myDict):
    sortedValues, sortedKeys = sortListsBasedOnValuesFromFirstList(
        [list(myDict.values()), list(myDict.keys())]
    )
    return dict(zip(sortedKeys, sortedValues))


def reverseDict(dict):
    dictToReturn = {}
    for key, value in dict.items():
        dictToReturn[value] = key
    return dictToReturn


def reverseDictWithTupleValues(dict):

    '''input dict:
    {1 : [('ca', 'ca', 'sy', 'hs'), ('ca', 'ca', 'sy', 'o')]}
    output dict:
    {('ca', 'ca', 'sy', 'hs') : 1, ('ca', 'ca', 'sy', 'o') : 2}'''

    dictToReturn = {}
    for key, value in dict.items():
        for tupl in value:
            dictToReturn[tupl] = key
    return dictToReturn


def renumberKeysAndValuesFrom0(dict):

    """ 
    Input: a dictionary with integer keys
    and lists of integers as values.
    Output: a dictionary with consecutive, incremental
    integers starting from 0 as keys, mapped in a one-to-one
    fashion to the original keys of dict. The values are 
    lists of integers remapped in the same way. The map 
    dictionary is also returned.
    """

    renumberedDict = {}
    consecutiveKeysToSparseKeys = {}
    sparseKeysToConsecutiveKeys = {}
    # create conversion maps
    for index, key in enumerate(dict):
        consecutiveKeysToSparseKeys[index] = int(key)
        sparseKeysToConsecutiveKeys[int(key)] = index
    # apply conversion to dict
    for key in dict:
        renumberedDict[sparseKeysToConsecutiveKeys[int(key)]] =\
            [sparseKeysToConsecutiveKeys[int(sparseLabel)]
                for sparseLabel in dict[key]]
    return renumberedDict, consecutiveKeysToSparseKeys


def getFirstKeyFromValue(dict, valueToFind):
    """If there are multiple keys with value,
       return the first one"""
    for key, value in dict.items():
        if value == valueToFind:
            return key


def changeTypeOfDictKeys(dict, ifType, newType):
    listOfKeys = list(dict.keys())
    if type(listOfKeys[0]) == ifType:
        for key in listOfKeys:
            dict[newType(key)] = dict[key]
            del dict[key]
    return dict
