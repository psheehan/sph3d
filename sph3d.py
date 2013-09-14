import ctypes
import numpy

def VectorFactory(D):
    class Vector(ctypes.Structure):
        _fields_ = [
            ("v", ctypes.c_double*D),
            ]

        def __getitem__(self, i):
            return self.v[i]

        def __setitem__(self, i, item):
            self.v[i] = item

    return Vector

def ParamsFactory(D):
    class Params(ctypes.Structure):
        _fields_ = [
            ("rmax", ctypes.c_double*D),
            ("time", ctypes.c_double),
            ("dt", ctypes.c_double),
            ("nsteps", ctypes.c_int),
            ("eta", ctypes.c_double),
            ("alpha_min", ctypes.c_double),
            ("beta", ctypes.c_double),
            ("epsilon", ctypes.c_double),
            ("sigma", ctypes.c_double),
            ("kernel_const", ctypes.c_double),
            ("np", ctypes.c_int),
            ("Lboundary", ctypes.c_int*D),
            ("Rboundary", ctypes.c_int*D),
            ]

    return Params

def ParticleFactory(D):
    class Particle(ctypes.Structure):
        _fields_ = [
            ("rho", ctypes.c_double),
            ("m", ctypes.c_double),
            ("r", VectorFactory(D)),
            ("v", VectorFactory(D)),
            ("a", VectorFactory(D)),
            ("g", VectorFactory(D)),
            ("vu", ctypes.c_double),
            ("ve", ctypes.c_double),
            ("adiss", VectorFactory(D)),
            ("vudiss", ctypes.c_double),
            ("vediss", ctypes.c_double),
            ("h", ctypes.c_double),
            ("P", ctypes.c_double),
            ("u", ctypes.c_double),
            ("e", ctypes.c_double),
            ("cs", ctypes.c_double),
            ("alpha", ctypes.c_double),
            ("alphau", ctypes.c_double),
            ("valpha", ctypes.c_double),
            ("valphau", ctypes.c_double),
            ("Omega", ctypes.c_double),
            ]

    return Particle

lib = ctypes.cdll.LoadLibrary('./libsph3d.so')

def SPHFactory(D):
    if (D == 1):
        class SPH(ctypes.Structure):
            def __init__(self, S, P):
                lib.sph_new_1D(ctypes.byref(self),ctypes.byref(S), \
                        ctypes.byref(P))

            def step(self):
                lib.sph_step_1D(ctypes.byref(self))

            def initialize(self):
                lib.sph_initialize_1D(ctypes.byref(self))

    elif (D == 2):
        class SPH(ctypes.Structure):
            def __init__(self, S, P):
                lib.sph_new_2D(ctypes.byref(self),ctypes.byref(S), \
                        ctypes.byref(P))

            def step(self):
                lib.sph_step_2D(ctypes.byref(self))

            def initialize(self):
                lib.sph_initialize_2D(ctypes.byref(self))

    if (D == 3):
        class SPH(ctypes.Structure):
            def __init__(self, S, P):
                lib.sph_new_3D(ctypes.byref(self),ctypes.byref(S), \
                        ctypes.byref(P))

            def step(self):
                lib.sph_step_3D(ctypes.byref(self))

            def initialize(self):
                lib.sph_initialize_3D(ctypes.byref(self))

    return SPH
