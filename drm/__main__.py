import argparse
from func import hercules
from func.getEquivalentForces import getEquivalentForces
from func.writeHistoryOutputForDRM import writeDisplacementHistoryForDRM
# from func.getNodes import writeNodeTable
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
    # NOTE: Change the parameters in kwargs accordingly
    kwargs = {
        'dbPath': '/Volumes/CORSAIR/Hercules/Cases/Istanbul M6p81/1Hz_Hoffman2/database/planedisplacements.hdf5',
        'planeIndices': list(range(133, 141)), # The plane indices used for writeDisplacementHistoryForDRM.
        'matchMethod': 'interpolated',
        'dispHistoryFileName': 'DispHistory.hdf5',
        # ===== Abaqus model information =====
        # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
        # total length of the part (including the interested domain, DRM layer, and PML layer).
        'lengths': {'x': 100, 'y': 100, 'z': 100},
        'DRM_depth': 5, # The thickness of the DRM layer should be the mesh size used in the Abaqus model
        'partName': 'Part-Soil',
        'materialName': 'soil',
        'jobName': 'DRM_PML_Analysis',
        'subroutineFileName': 'PML3dUEL_homogeneous.for',
        'timeIncrement': 0.1, # unit: sec. This should be read from the site response file
        'duration': 44.9, # unit: sec. This should be read from the site response file
        'origin': (41.021874, 28.885526),
        'noDamping': True,
        'stepApplication': 'moderate dissipation',
        'isStaticStepIncluded': False,
        # ===== Rayleigh Damping =====
        'alpha': 0,
        'beta': 0
    }
    kwargs['PML_depth'] = kwargs['DRM_depth']*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    # =====
    if stepNum == 1:
        # ===== Hercules model information =====
        HerculesInputFilePath = 'inputfiles/parameters.in'
        HerculesData = hercules.getInputData(HerculesInputFilePath)
        kwargs['minVs'] = float(HerculesData['simulation_shear_velocity_min'])
        # ===== Get material properties =====
        kwargs['centroid'] = [x/2 for x in kwargs['lengths'].values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(**kwargs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_'+axis: kwargs['lengths'][axis] for axis in ['x', 'y', 'z']} \
            | {'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density} \
            | kwargs
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # For homogeneous model without static step:
        modifyInput(**kwargs)
    elif stepNum == 3:
        # ===== Getting the displacement histories for DRM nodes =====
        # NOTE: mpirun works here
        import timeit
        numExec = 1
        print('Averaged Elapsed Time: %.2f secs' % (timeit.timeit(lambda: writeDisplacementHistoryForDRM(**kwargs), number=numExec)/numExec))
    elif stepNum == 4:
        # ===== Compute the equivalent forces =====
        kwargs['isHomogeneous'] = True
        getEquivalentForces(**kwargs)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file
    else:
        raise ValueError('There are only 4 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., Hercules' station files).
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
        ('abaqus viewer odb=Istanbul_model_static.odb script=getReactionForce.py' can do the trick.)
        ===== The following is no longer needed with the use of plane files ==== 
        In addition, the generated nodeTable.csv should be used to specify stations in Hercules' input 
        file and run Hercules analysis to generate station files.
        ===== =====
        step 3: writeDisplacementHistoryForDRM (MPI is not supported for matchMethod='nearest' yet)
        step 4: getEquivalentForces 
        NOTE: After this step, move Cload.txt to the Abaqus working directory and run the complete analysis.
        ('abaqus viewer odb=Istanbul_model_complete.odb script=plotAndSave.py' can be used to save the 
        responses at designated points.) '''
    # NOTE: Change the parameters in kwargs accordingly
    kwargs = {
        'dbPath': '/Volumes/CORSAIR/Hercules/Cases/Istanbul M6p81/1Hz_Hoffman2/database/planedisplacements.hdf5',
        # planeIndices are the plane indices used for writeDisplacementHistoryForDRM.
        'planeIndices': list(range(12, 23)), # Case 1, 8-meter element size
        # 'planeIndices': list(range(56, 67)), # Case 2, 8-meter element size
        # 'planeIndices': list(range(100, 111)), # Case 3, 8-meter element size
        # 'planeIndices': list(range(133, 141)), # Case 4, 50-meter element size
        # === Version 3 ===
        # 'planeIndices': list(range(1, 12)), # Case 1, 8-meter element size
        # 'planeIndices': list(range(12, 23)), # Case 2, 8-meter element size
        # 'planeIndices': list(range(23, 34)), # Case 3, 8-meter element size
        # 'planeIndices': list(range(34, 45)), # Case 4, 8-meter element size
        'matchMethod': 'nearest', # can be 'nearest' or 'interpolated'
        'allowance': 0.1, # unit: m.
        'dispHistoryFileName': 'DispHistory.hdf5',
        # ===== Abaqus model information =====
        # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
        # total length of the part (including the interested domain, DRM layer, and PML layer).
        'lengths': {'x': 400, 'y': 400, 'z': 120},
        'DRM_depth': 8, # The thickness of the DRM layer should be the mesh size used in the Abaqus model
        'partName': 'Part-Soil',
        'materialName': 'soil',
        'jobName': 'Istanbul_model',
        'subroutineFileName': 'PML3dUEL_Inhomogeneous.for',
        'timeIncrement': 0.1, # unit: sec. This should be read from the site response file
        'duration': 44.9, # unit: sec. This should be read from the site response file
        'origin': (41.0318, 28.9417), # Case 1
        # 'origin': (41.001970, 29.085491), # Case 2
        # 'origin': (41.031450, 28.784960), # Case 3
        # 'origin': (41.022314, 28.886133), # Case 4
        # === Version 3 ===
        # 'origin': (41.0318, 28.9417), # Case 1
        # 'origin': (41.029468, 28.944425), # Case 2
        # 'origin': (41.040714, 28.921306), # Case 3
        # 'origin': (41.001530, 28.985787), # Case 4
        'isCoordinateConverted': True,
        'isOriginAtCenter': True,
        # ===== Rayleigh Damping =====
        'alpha': 0,
        'beta': 0,
    }
    kwargs['PML_depth'] = kwargs['DRM_depth']*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    # ===== Hercules model information =====
    HerculesInputFilePath = 'inputfiles/parameters.in'
    HerculesData = hercules.getInputData(HerculesInputFilePath)
    kwargs['minVs'] = float(HerculesData['simulation_shear_velocity_min'])
    if stepNum == 1:
        # ===== Get material properties =====
        kwargs['centroid'] = [x/2 for x in kwargs['lengths'].values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(**kwargs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_'+axis: kwargs['lengths'][axis] for axis in ['x', 'y', 'z']} \
            | {'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density} \
            | kwargs
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # For homogeneous model:
        # ===== Modify the preliminary Abaqus input file =====
        # modifyInput(**kwargs, modelType='static')
        # modifyInput(**kwargs, modelType='complete')

        # For heterogeneous model
        prev_jobName = kwargs['jobName']
        # ===== Writing Material Properties =====
        kwargs['jobName'] += '_static_pre'
        writeMaterialPropertiesForElements(**kwargs)
        writeMaterialPropertiesForPMLElements(**kwargs)
        # ===== Modify the preliminary Abaqus input file =====
        kwargs['jobName'] = prev_jobName
        modifyInput(**kwargs, isHomogeneous=False, modelType='static')
        modifyInput(**kwargs, isHomogeneous=False, modelType='complete')
        # ===== Write nodeTable =====
        # NOTE: This is the old method and has been replaced with the use of plane files.
        # writeNodeTable(prev_jobName+'_static', partName, coordinateCorrection=True, AbaqusModelOrigin=origin, isCoordinateConverted=isCoordinateConverted)
        # NOTE: This node table should be used to specify stations in Hercules' input file and run Hercules analysis to generate station files.
    elif stepNum == 3:
        # ===== Getting the displacement histories for DRM nodes =====
        import timeit
        numExec = 1
        kwargs['jobName'] += '_complete'
        print('Averaged Elapsed Time: %.2f secs' % 
            (timeit.timeit(lambda: writeDisplacementHistoryForDRM(**kwargs), number=numExec)/numExec))
    elif stepNum == 4:
        kwargs.update({
            # 'isHomogeneous': True, # NOTE: Uncomment this line for homogeneous model
            'jobName': kwargs['jobName'] + '_complete',
            'RFFileName': 'RF.txt',
        })
        # ===== Compute the equivalent forces =====
        getEquivalentForces(**kwargs)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file.
    else:
        raise ValueError('There are only 4 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., Hercules' station files).
    return

if __name__ == '__main__':
    # /// Parse input argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--stepNum', '-s', type=int, help='The step number for execution.')
    stepNum = parser.parse_args().stepNum
    # /// Pick which model you want to perform
    # prospectusModel(stepNum)
    IstanbulModel(stepNum)