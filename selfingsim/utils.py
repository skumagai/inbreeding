def getdata(f):
    with open(f, "r") as rf:
        return json.load(rf)

def getn(data):
    return len(data)

def getnloc(data):
    if len(data) > 0:
        return len(data[0][2])
    else:
        return 0
