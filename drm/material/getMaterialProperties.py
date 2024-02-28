import ctypes
# import os
import numpy as np

# NOTE: the c script has to be converted to .so file with the following command:
# cc -fPIC -shared -o material_property_relative_V13.so material_property_relative_V13.c
# Reference: 
# https://medium.com/spikelab/calling-c-functions-from-python-104e609f2804
# https://www.journaldev.com/31907/calling-c-functions-from-python

def getMaterialProperties(points: dict[int, list[float, float, float]]) -> dict[int, list[float]]:
    # NOTE: x, y, z in each point should be in meter, and the directions are NS, EW, and pointing down to the earth (same as Hercules)
    # fileFolder = os.path.dirname(os.path.abspath(__file__))
    cLib = np.ctypeslib.load_library('material_property_relative_V13.so', 'drm/material')
    getMaterialPropertiesInC = cLib.material_property_relative_V13
    getMaterialPropertiesInC.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(7, )) # NOTE: There are 7 returned values
    elementMaterialProperties = {}
    for tag, [x, y, z] in points.items():
        # NOTE: The directions of x, y, z in material_property_relative_V13.c are EW, NS, and pointing up from the earth
        materialProperties = getMaterialPropertiesInC(ctypes.c_double(y), ctypes.c_double(x), ctypes.c_double(-z))
        # NOTE: The 4th to 7th material properties in the list are for nonlinear material. We don't need it for now.
        # NOTE 2: The properties from index 0 to 2 are Vs (m/s), Vp (m/s), and rho (density, ton/m^3)
        materialProperties[2] = 1000.0*materialProperties[2] # Convert density from ton/m^3 to kg/m^3 (SI consistent)
        elementMaterialProperties[tag] = materialProperties[:3]
    # NOTE: The final elements in the return array are Vs (m/s), Vp (m/s), and rho (density, kg/m^3)
    return elementMaterialProperties

if __name__ == '__main__':
    # NOTE: This is just an example.
    materialProperties = getMaterialProperties({1: [8200, 1300, 30]})
    print(materialProperties)