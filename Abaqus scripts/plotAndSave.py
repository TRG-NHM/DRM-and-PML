# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import csv

def openOdb(odbName):
    if odbName is None:
        # ===== Read the .odb file opened in current viewport, so make sure to open the correct .odb file in Abaqus CAE =====
        vp = session.viewports[session.currentViewportName]
        odbName = vp.odbDisplay.name
        odb = session.odbs[odbName]
    else:
        odb = session.openOdb(name=odbName)
    return odb

def cleanXYDataInCurrentSession():
    for xyDataName in session.xyDataObjects.keys():
        del session.xyDataObjects[xyDataName]
    return

def plotResponseAndSave(nodeLabels, odbName=None, stepNames=None, instanceName=None, timeReduction=0):
    # NOTE: Abaqus' built-in xyPlot can get field outputs way faster than getting them frame by frame in Python.
    # So it's recommended to use plotResponseAndSave() over getResponseAndSave.
    # NOTE 2: Unfortunately, for unknown reason, plotResponseAndSave() can only be run with GUI. xyPlot.xyDataListFromField will raise an error with noGUI option.
    # /// Handle input variables
    if type(nodeLabels) is not list:
        nodeLabels = [nodeLabels]
    odb = openOdb(odbName)
    odbName = odb.name # NOTE: In case the input odbName is None
    if stepNames is None:
        stepNames = odb.steps.keys()
    elif stepNames is not None and type(stepNames) is not list:
        stepNames = [stepNames]
    if instanceName is None and len(odb.rootAssembly.instances) == 1:
        instanceName = odb.rootAssembly.instances.keys()[0]
    else:
        raise ValueError('The model contains multiple instances. instanceName must be defined.')
    cleanXYDataInCurrentSession()
    quantities = ['U', 'V', 'A']
    for nodeLabel in nodeLabels:
        xyData = [['Time(s)', 'X|(m)', 'Y-(m)', 'Z.(m)', 'X|(m/s)', 'Y-(m/s)', 'Z.(m/s)', 'X|(m/s2)', 'Y-(m/s2)', 'Z.(m/s2)']]
        suffix = 'PI: %s N: %d'%(instanceName, nodeLabel)
        for stepName in stepNames:
            session.odbData[odbName].setValues(activeFrames=((stepName, ('0:-1', )), ))
            xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(
                ('U', NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), (COMPONENT, 'U3'), )), 
                ('V', NODAL, ((COMPONENT, 'V1'), (COMPONENT, 'V2'), (COMPONENT, 'V3'), )), 
                ('A', NODAL, ((COMPONENT, 'A1'), (COMPONENT, 'A2'), (COMPONENT, 'A3'), )),
                ), nodeLabels=((instanceName, (nodeLabel, )), ))
            times = [x[0] for x in session.xyDataObjects['_U:U1 '+suffix]]
            for i, time in enumerate(times):
                line = ['%.5f'%(time - timeReduction)]
                for quantity in quantities:
                    for component in [1, 2, 3]:
                        line.append(session.xyDataObjects['_%s:%s%d %s'%(quantity, quantity, component, suffix)][i][1])
                xyData.append(line)
        with open('disp_vel_accel.%s.%s.csv'%(instanceName, nodeLabel), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(xyData)
    return

# def getNodeIndex(odb, instanceName, nodeLabel, stepName, fieldOutputName):
#     ''' Check whether instance and node label exist and find the corresponding index for getting field output values '''
#     if len(odb.rootAssembly.instances) > 1:
#         if instanceName is None:
#             raise ValueError('Instance Name must be specified if there are more than 1 instance in the model.')
#         if instanceName not in odb.rootAssembly.instances.keys():
#             raise ValueError('The specified instance %s does not exist in the model. Available names include %s'%(instanceName, ', '.join(odb.rootAssembly.instances.keys())))
#         nodeIndex = None
#         for i, value in enumerate(odb.steps[stepName].frames[0].fieldOutputs[fieldOutputName].values):
#             if value.instance is not None and value.instance.name == instanceName and value.nodeLabel == nodeLabel:
#                 nodeIndex = i
#                 print('Corresponding node index for node label %d is %d\n'%(nodeLabel, nodeIndex))
#                 break
#         if nodeIndex is None:
#             raise ValueError('Node label %d does not exist on instance %s.'%(nodeLabel, instanceName))
#     else:
#         nodeIndex = nodeLabel - 1
#     return nodeIndex

# def getResponseAndSave(nodeLabels, odbName=None, stepNames=None, instanceName=None, quantities=['U', 'V', 'A']):
#     # /// Handle input variables
#     if type(nodeLabels) is not list:
#         nodeLabels = [nodeLabels]
#     odb = openOdb(odbName)
#     if stepNames is None:
#         stepNames = odb.steps.keys()
#     elif stepNames is not None and type(stepNames) is not list:
#         stepNames = [stepNames]
#     # =====
#     components = {'U': ['X|(m)', 'Y-(m)', 'Z.(m)'], 'V': ['X|(m/s)', 'Y-(m/s)', 'Z.(m/s)'], 'A': ['X|(m/s2)', 'Y-(m/s2)', 'Z.(m/s2)']}
#     for nodeLabel in nodeLabels:
#         xyData = [['Time(s)']]
#         for quantity in quantities:
#             xyData[0].extend(components[quantity])
#         for stepName in stepNames:
#             nodeIndex = getNodeIndex(odb, instanceName, nodeLabel, stepName, quantities[0])
#             for frame in odb.steps[stepName].frames:
#                 line = ['%.5f'%frame.frameValue]
#                 for quantity in quantities:
#                     for i, componentLabels in enumerate(frame.fieldOutputs[quantity].componentLabels):
#                         line.append(frame.fieldOutputs[quantity].values[nodeIndex].data[i])
#                         # NOTE: For history output, can use odb.steps[stepName].historyRegions['Node %s.%d'%(instanceName, nodeLabel)].historyOutputs['%s%d'%(quantity, component)].data[i]
#                 xyData.append(line)
#         fileNames = {'U': 'disp', 'V': 'vel', 'A': 'accel'}
#         fileName = '_'.join([fileNames[quantity] for quantity in quantities]) + '.%d.csv'%nodeLabel
#         with open(fileName, 'wb') as f:
#             writer = csv.writer(f)
#             writer.writerows(xyData)

def plotAndSaveHistoryOutput(nodeLabels, odbName=None):
    if type(nodeLabels) is not list:
        nodeLabels = [nodeLabels]
    odb = openOdb(odbName)
    odbName = odb.name # NOTE: In case the input odbName is None
    historyOutputNames = session.odbData[odb.name].historyVariables.keys()
    quantities = ['Spatial displacement: U', 'Spatial velocity: V', 'Spatial acceleration: A']
    for nodeLabel in nodeLabels:
        xyData = []
        for quantity in quantities:
            for direction in [1, 2, 3]:
                prefix = quantity+str(direction)+' at Node '+str(nodeLabel)
                for historyOutputName in historyOutputNames:
                    if historyOutputName.startswith(prefix):
                        xyData.append(session.XYDataFromHistory(name=quantity[-1]+str(direction), 
                            odb=odb, outputVariableName=historyOutputName))
                        break
        session.writeXYReport(
            fileName = 'disp_vel_accel.%s.rpt'%nodeLabel, appendMode=OFF, xyData=tuple(xyData))

if __name__ == '__main__':
    # nodeLabels = [8302, 8922]
    nodeLabels = [29406, 36011]
    # plotResponseAndSave(nodeLabels)
    plotAndSaveHistoryOutput(nodeLabels)