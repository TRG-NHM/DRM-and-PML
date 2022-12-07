import pandas as pd
import numpy as np
from func.getNodes import getNodeCoordinates
from func.getDistance import getDistanceBetweenTwoCoordinates
from func.getElements import getElements
from func.getPMLElementInputLine import getUserElementLines, getParametersLines
from materialDatabase.getMaterialProperties import getMaterialProperties

def writeRawData(materialProperties: dict[int, list[float, float, float]], fileName=None, minVs=None) -> None:
    df = pd.DataFrame(materialProperties)
    df = df.transpose()
    df.columns = ['Vs', 'Vp', 'rho']
    df['poissonsRatio'] = (df['Vp']**2 - 2*df['Vs']**2)/(2*(df['Vp']**2 - df['Vs']**2))
    df['G'] = df['rho']*df['Vs']**2 # Shear Modulus
    df['youngsModulus'] = 2*df['G']*(1+df['poissonsRatio'])
    if minVs is not None:
        df['Vs_corrected'] = df['Vs']
        df['Vs_corrected'][df['Vs'] < minVs] = 250
        df['poissonsRatio_corrected'] = (df['Vp']**2 - 2*df['Vs_corrected']**2)/(2*(df['Vp']**2 - df['Vs_corrected']**2))
        df['G_corrected'] = df['rho']*df['Vs_corrected']**2 # Shear Modulus
        df['youngsModulus_corrected'] = 2*df['G_corrected']*(1+df['poissonsRatio_corrected'])
    if fileName is None:
        fileName = 'rawMaterial.csv'
    df.to_csv(fileName)
    return

# NOTE: Since we moved the model from its original place to another location, we have to get the material properties with given coordinates
def writeMaterialPropertiesForElements(jobName, partName, origin, isCoordinateConverted=False, minVs=None, 
    materialFileName='ElementMaterial.txt', elementSetFileName='ElementSet.txt', sectionFileName='ElementSection.txt',
    istanbulMaterialModelOrigin=(40.9565, 28.675), elementTypeOnPart='C3D8'):
    ''' NOTE: If the directions used in Abaqus model are EW, NS, and pointing up 
    from the earth, then isCoordinateConverted can be set to True to correct the 
    differences between the directions used in Hercules and Abaqus '''
    distance, x_dis, y_dis = getDistanceBetweenTwoCoordinates(istanbulMaterialModelOrigin, origin)
    nodes = getNodeCoordinates(jobName, partName)
    elements = getElements(jobName, partName, elementType=elementTypeOnPart)
    elementCentroids = {}
    for eleLabel, nodeLabels in elements.items():
        nodesOnElement = np.array([nodes[label] for label in nodeLabels])
        centroid = sum(nodesOnElement)/len(nodesOnElement)
        if isCoordinateConverted:
            elementCentroids[eleLabel] = [centroid[1]+x_dis, centroid[0]+y_dis, -centroid[2]]
        else:
            elementCentroids[eleLabel] = [centroid[0]+x_dis, centroid[1]+y_dis, centroid[2]]
    # ===== used for debugging =====
    # df = pd.DataFrame(elementCentroids)
    # df = df.transpose()
    # df.columns = ['x', 'y', 'z']
    # df.to_csv('elementCentroids.csv')
    # return
    # ===== =====
    materialProperties = getMaterialProperties(elementCentroids)
    # ===== used for debugging =====
    writeRawData(materialProperties, fileName='rawMaterialDRM.csv', minVs=minVs)
    # return
    # ===== =====
    # Element Set File Writing
    with open(elementSetFileName, 'w') as f:
        f.writelines(['*Elset, elset=SET-ELEMENT-%d\n%d,\n'%(eleLabel, eleLabel) for eleLabel in elements.keys()])
    # Element Section File Writing
    with open(sectionFileName, 'w') as f:
        f.writelines(['*Solid Section, elset=SET-ELEMENT-%d, material=MATERIAL-%d\n,\n'%(eleLabel, eleLabel) for eleLabel in elements.keys()])
    # Element Material File Writing
    lines = []
    for eleLabel in elements.keys():
        Vs, Vp, rho = materialProperties[eleLabel]
        # NOTE: Applying this may make the material properties weird (even with some negative Poisson's Ratios). Use with cautions.
        if minVs is not None and Vs < minVs:
            Vs = minVs
        poissonsRatio = (Vp**2 - 2*Vs**2)/(2*(Vp**2 - Vs**2)) # Ref: Eq. (5.34) in Geotechnical Earthquake Engineering (Kramer)
        G = rho*Vs**2 # Shear Modulus
        youngsModulus = 2*G*(1+poissonsRatio)
        lines += ['*Material, name=MATERIAL-%d\n'%eleLabel,
            '*Density\n', '%f,\n'%rho,
            '*Elastic\n', '%f, %f\n'%(youngsModulus, poissonsRatio)]
    with open(materialFileName, 'w') as f:
        f.writelines(lines)
    return

def writeMaterialPropertiesForPMLElements(jobName, partName, origin, PML_depth, lengths, alpha=0.0, beta=0.0,
    isCoordinateConverted=False, minVs=None, materialFileName='ElementPML.txt', maxVpFileName='maxVpInPMLDomain.txt',
    istanbulMaterialModelOrigin=(40.9565, 28.675), elementTypeOnPart='C3D8R', userElementType='U3'):
    ''' NOTE: If the directions used in Abaqus model are EW, NS, and pointing up 
    from the earth, then isCoordinateConverted can be set to True to correct the 
    differences between the directions used in Hercules and Abaqus '''
    distance, x_dis, y_dis = getDistanceBetweenTwoCoordinates(istanbulMaterialModelOrigin, origin)
    nodes = getNodeCoordinates(jobName, partName)
    elements = getElements(jobName, partName, elementType=elementTypeOnPart)
    elementCentroids = {}
    for eleLabel, nodeLabels in elements.items():
        nodesOnElement = np.array([nodes[label] for label in nodeLabels])
        centroid = sum(nodesOnElement)/len(nodesOnElement)
        if isCoordinateConverted:
            elementCentroids[eleLabel] = [centroid[1]+x_dis, centroid[0]+y_dis, -centroid[2]]
        else:
            elementCentroids[eleLabel] = [centroid[0]+x_dis, centroid[1]+y_dis, centroid[2]]
    materialProperties = getMaterialProperties(elementCentroids)
    # ===== used for debugging =====
    writeRawData(materialProperties, fileName='rawMaterialPML.csv', minVs=minVs)
    # return
    # ===== =====
    # File Writing
    lines = getUserElementLines(userElementType)
    maxVp = 0
    minVp = 1E9
    for eleLabel, nodeLabels in elements.items():
        Vs, Vp, rho = materialProperties[eleLabel]
        if minVs is not None and Vs < minVs:
            Vs = minVs
        if Vp > maxVp:
            maxVp = Vp
        if Vp < minVp:
            minVp = Vp
        poissonsRatio = (Vp**2 - 2*Vs**2)/(2*(Vp**2 - Vs**2)) # Ref: Eq. (5.34) in Geotechnical Earthquake Engineering (Kramer)
        G = rho*Vs**2 # Shear Modulus
        youngsModulus = 2*G*(1+poissonsRatio)
        lines += ['*Element, type=%s, elset=EL-PML-ELEMENT-%d\n'%(userElementType, eleLabel),
            '%i, '%eleLabel+', '.join([str(nodeLabel) for nodeLabel in nodeLabels])+'\n']
        lines += getParametersLines(youngsModulus, poissonsRatio, rho, PML_depth, lengths, 
            alpha, beta, elementSetName='EL-PML-ELEMENT-%d'%eleLabel)
    with open(materialFileName, 'w') as f:
        f.writelines(lines)
    with open(maxVpFileName, 'w') as f:
        f.writelines(['MaxVp,MinVp\n', str(maxVp)+','+str(minVp)+'\n'])
    return

# if __name__ == '__main__':
#     from func.getDataFromInputFile import getDataFromHerculesInputFile
#     # ===== Abaqus model information =====
#     jobName = 'Model_complete'
#     partName = 'Part-Soil'
#     origin = (41.0318, 28.9417)
#     lengths = {'x': 100, 'y': 100, 'z': 30}
#     DRM_depth = 2 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
#     PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
#     # ===== Hercules model information =====
#     HerculesData = getDataFromHerculesInputFile('Istanbul_sim55/inputfiles/parameters_FullRegion_all_station_topo_final.in')
#     minVs = float(HerculesData['simulation_shear_velocity_min'])
#     # ===== Writing Material Properties =====
#     writeMaterialPropertiesForElements(jobName, partName, origin, isCoordinateConverted=True, minVs=minVs)
#     writeMaterialPropertiesForPMLElements(jobName, partName, origin, PML_depth, 
#         lengths, isCoordinateConverted=True, minVs=minVs)