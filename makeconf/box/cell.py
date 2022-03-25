class Cell:

    def __init__(self, id, rc):
        self.id = id
        self.voroid = None
        self.rc = rc
        self.verts = []
        self.type = "passive"
        self.area = None
        self.perim = None
        self.outer = False
        self.kappa = None
        self.gamma = None
        self.lam = None
        self.k = None
        self.l0 = []
        self.cell_myo = None
        self.myo = []
