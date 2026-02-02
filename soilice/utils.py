import dill

def save(filename,output):
    from .src_soil import modelInOut
    data=modelInOut()
    for k,v in vars(output).items():
        setattr(data, k, v)
    with open(filename, "wb") as f:
        dill.dump(data, f)

def load(filename):
    with open(filename, "rb") as f:
        return dill.load(f)

def loadModel(filename):
    from .src_soil import model
    sim=model(opts=None)
    for k,v in vars(load(filename)).items():
        setattr(sim, k, v)
    return sim
        