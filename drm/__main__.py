import argparse
from func import hercules
from func.getEquivalentForces import getEquivalentForces
from func.writeHistoryOutputForDRM import writeDisplacementHistoryForDRM
from func.getNodes import writeNodeTable
from abaqusInputProcessing import getMaterialPropertiesAtCentroid, writeDataForCreatingAbaqusModel, modifyInput
from writeMaterialProperties import writeMaterialPropertiesForElements, writeMaterialPropertiesForPMLElements

def prospectusModel(stepNum: int) -> None:
    ''' This function can be used as an example of homogeneous model without static step.
        step 1: writeDataForCreatingAbaqusModel
        NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare 
        the model (run `createDomain.py`. abaqus cae noGUI=createDomain.py can do the trick).
        After the _pre.inp file generated, move the _pre.inp file back.
        step 2: modifyInput 
        step 3: writeDisplacementHistoryForDRM (MPI can be used in this step)
        step 4: getEquivalentForces 
        NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file '''
    dbPath = '../HerculesDatabaseInquirySystem/database/planedisplacements.hdf5'
    # ===== Abaqus model information =====
    # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
    # total length of the part (including the interested domain, DRM layer, and PML layer).
    lengths = {'x': 100, 'y': 100, 'z': 100}
    DRM_depth = 5 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
    PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    partName = 'Part-Soil'
    materialName = 'soil'
    # ===== Rayleigh Damping =====
    alpha = 0
    beta = 0
    # =====
    jobName = 'DRM_PML_Analysis'
    subroutineFileName = 'PML3dUEL912.for'
    timeIncrement = 0.05 # unit: sec. This should be read from the site response file
    duration = 29.95 # unit: sec. This should be read from the site response file
    origin = (41.022024713821054, 28.88575077152881)
    noDamping = True
    stepApplication = 'moderate dissipation'
    isStaticStepIncluded = False
    if stepNum == 1:
        # ===== Hercules model information =====
        HerculesInputFilePath = '../Mw5p83/inputfiles/parameters_ZDistrictDRM.in'
        HerculesData = hercules.getInputData(HerculesInputFilePath)
        minVs = float(HerculesData['simulation_shear_velocity_min'])
        # ===== Get material properties =====
        centroid = [x/2 for x in lengths.values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(centroid, origin=origin, minVs=minVs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_x': lengths['x'], 'length_y': lengths['y'], 'length_z': lengths['z'], 
            'DRM_depth': DRM_depth, 'PML_depth': PML_depth,
            'partName': partName, 'materialName': materialName, 
            'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density, 
            'alpha': alpha, 'beta': beta, 
            'jobName': jobName, 'subroutineFileName': subroutineFileName,
            'timeIncrement': timeIncrement, 'duration': duration,
            'noDamping': noDamping, 'stepApplication': stepApplication,
            'isStaticStepIncluded': isStaticStepIncluded}
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # For homogeneous model without static step:
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta)
    elif stepNum == 3:
        # ===== Getting the displacement histories for DRM nodes =====
        # NOTE: mpirun works here
        import timeit
        numExec = 1
        print('Averaged Elapsed Time: %.2f secs' % (timeit.timeit(lambda: writeDisplacementHistoryForDRM(dbPath, jobName, partName, origin), number=numExec)/numExec))
        # df = writeDisplacementHistoryForDRM(jobName, partName, origin)
    elif stepNum == 4:
        # ===== Compute the equivalent forces =====
        getEquivalentForces(jobName, partName, isHomogeneous=True, dispHistoryFileName='DispHistory.csv', materialName=materialName)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file
    else:
        raise ValueError('There are only 4 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., DispHistory.csv generated based on HDF5 database).
    return

def IstanbulModel(stepNum: int) -> None:
    ''' This function can be used as an example of heterogeneous model with static step.
        step 1: writeDataForCreatingAbaqusModel
        NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare 
        the model (run `createDomain.py`. abaqus cae noGUI=createDomain.py can do the trick).
        After the _pre.inp files generated, move the _pre.inp files back.
        step 2: modifyInput 
        NOTE: After this step, move the .inp files (for inhomogeneous model, also move ElementSet.txt, 
        ElementSection.txt, ElementMaterial.txt, ElementPML.txt, and PML3dUEL_Inhomogeneous.for manually 
        modified with maxVp in maxVpInPMLDomain.txt) to the Abaqus working directory and run the static 
        analysis and get RF.txt by running getReactionForce.py, then move RF.txt back.
        In addition, the generated nodeTable.csv should be used to specify stations in Hercules' input 
        file and run Hercules analysis to generate station files.
        step 3: getEquivalentForces 
        NOTE: After this step, move Cload.txt to the Abaqus working directory and run the complete analysis '''
    # ===== Abaqus model information =====
    # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
    # total length of the part (including the interested domain, DRM layer, and PML layer).
    lengths = {'x': 100, 'y': 100, 'z': 30}
    DRM_depth = 2 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
    PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    partName = 'Part-Soil'
    materialName = 'soil'
    # ===== Rayleigh Damping =====
    alpha = 0
    beta = 0
    # =====
    jobName = 'Istanbul_model'
    subroutineFileName = 'PML3dUEL912.for'
    timeIncrement = 0.01 # unit: sec. This should be read from the site response file
    duration = 39.99 # unit: sec. This should be read from the site response file
    origin = (41.0318, 28.9417)
    isCoordinateConverted = True
    isOriginAtCenter = True
    # ===== Hercules model information =====
    HerculesInputFilePath = '../Abaqus Steel Building Model from Bulent/Istanbul_sim55/inputfiles/parameters_FullRegion_all_station_topo_final.in'
    HerculesData = hercules.getInputData(HerculesInputFilePath)
    minVs = float(HerculesData['simulation_shear_velocity_min'])
    if stepNum == 1:
        # ===== Get material properties =====
        centroid = [x/2 for x in lengths.values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(centroid, origin=origin, minVs=minVs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_x': lengths['x'], 'length_y': lengths['y'], 'length_z': lengths['z'], 
            'DRM_depth': DRM_depth, 'PML_depth': PML_depth,
            'partName': partName, 'materialName': materialName, 
            'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density, 
            'alpha': alpha, 'beta': beta, 
            'jobName': jobName, 'subroutineFileName': subroutineFileName,
            'timeIncrement': timeIncrement, 'duration': duration,
            'isCoordinateConverted': isCoordinateConverted,
            'isOriginAtCenter': isOriginAtCenter}
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # # For homogeneous model:
        # # ===== Modify the preliminary Abaqus input file =====
        # modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, modelType='static')
        # modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, modelType='complete')

        # For heterogeneous model
        # ===== Writing Material Properties =====
        writeMaterialPropertiesForElements(jobName+'_static_pre', partName, origin, isCoordinateConverted=True, minVs=minVs)
        writeMaterialPropertiesForPMLElements(jobName+'_static_pre', partName, origin, PML_depth, 
            lengths, isCoordinateConverted=True, minVs=minVs)
        # ===== Modify the preliminary Abaqus input file =====
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, isHomogeneous=False, modelType='static')
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, isHomogeneous=False, modelType='complete')
        # ===== Write nodeTable =====
        writeNodeTable(jobName+'_static', partName, coordinateCorrection=True, AbaqusModelOrigin=origin, isCoordinateConverted=True)
        # NOTE: This node table should be used to specify stations in Hercules' input file and run Hercules analysis to generate station files.
    elif stepNum == 3:
        # ===== Compute the equivalent forces =====
        stationFolder = '../Abaqus Steel Building Model from Bulent/Stations_flat/'
        getEquivalentForces(jobName+'_complete', partName, stationFolder=stationFolder, RFFileName='RF.txt', isCoordinateConverted=True)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file.
    else:
        raise ValueError('There are only 3 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., Hercules' station files).
    return

if __name__ == '__main__':
    # /// Parse input argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--stepNum', '-s', type=int, help='The step number for execution.')
    stepNum = parser.parse_args().stepNum
    # /// Pick which model you want to perform
    prospectusModel(stepNum)
    # IstanbulModel(stepNum)