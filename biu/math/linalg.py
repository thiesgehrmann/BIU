from .. import utils

np = utils.py.loadExternalModule("numpy")

###########################################################################

def scalar_projection(vector_a, vector_b):
    """ Project vector_a onto vector_b"""
    return np.dot(vector_a, (vector_b / (np.linalg.norm(vector_b))))
#edef

def angle_between_vectors(vector_a, vector_b):
    """
    Determine the angle (in degrees) between the vector_a and vector_b
    Always returns the acute angle, regardless of orientation of the vectors
    """
    return 180 * np.arccos(np.dot(vector_a, vector_b)/(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))) / np.pi
#edef