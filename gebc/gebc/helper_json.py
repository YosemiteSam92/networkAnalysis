import json


def loadJson(fileName):
    with open(fileName, "r") as f:
        return json.load(f)


def printJson(dataToDump, fileName):
    with open(fileName, "w") as f:
        json.dump(dataToDump, f, indent=4)
