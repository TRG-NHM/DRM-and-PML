# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import visualization
import xyPlot
import displayGroupOdbToolset as dgo

def openOdb(odbName):
    if odbName is None:
        # ===== Read the .odb file opened in current viewport, so make sure to open the correct .odb file in Abaqus CAE =====
        vp = session.viewports[session.currentViewportName]
        odbName = vp.odbDisplay.name
        odb = session.odbs[odbName]
    else:
        odb = session.openOdb(name=odbName)
    return odb

def getNodeIndex(odb, instanceName, nodeLabel, stepName, fieldOutputName):
    ''' Check whether instance and node label exist and find the corresponding index for getting field output values '''
    if len(odb.rootAssembly.instances) > 1:
        if instanceName is None:
            raise ValueError('Instance Name must be specified if there are more than 1 instance in the model.')
        if instanceName not in odb.rootAssembly.instances.keys():
            raise ValueError('The specified instance %s does not exist in the model. Available names include %s'%(instanceName, ', '.join(odb.rootAssembly.instances.keys())))
        nodeIndex = None
        for i, value in enumerate(odb.steps[stepName].frames[0].fieldOutputs[fieldOutputName].values):
            if value.instance is not None and value.instance.name == instanceName and value.nodeLabel == nodeLabel:
                nodeIndex = i
                print('Corresponding node index for node label %d is %d\n'%(nodeLabel, nodeIndex))
                break
        if nodeIndex is None:
            raise ValueError('Node label %d does not exist on instance %s.'%(nodeLabel, instanceName))
    else:
        nodeIndex = nodeLabel - 1
    return nodeIndex

def writeRF(DRMNodeSetName='outDRM', odbName=None, instanceName=None, 
    stepName=None, outputFileName='RF.txt'):
    odb = openOdb(odbName)
    odbName = odb.name # NOTE: In case the input odbName is None
    if stepName is None:
        stepName = odb.steps.keys()[0]
    step = odb.steps[stepName]
    if instanceName is None and len(odb.rootAssembly.instances) == 1:
        instanceName = odb.rootAssembly.instances.keys()[0]
    else:
        raise ValueError('The model contains multiple instances. instanceName must be defined.')
    instance = odb.rootAssembly.instances[instanceName]
    # =====
    DRMNodeLabels = [node.label for node in instance.nodeSets[DRMNodeSetName.upper()].nodes]
    lines = []
    for nodeLabel in DRMNodeLabels:
        nodeIndex = getNodeIndex(odb, instanceName, nodeLabel, stepName, 'RF')
        RF = step.frames[-1].fieldOutputs['RF'].values[nodeIndex].data
        for i, direction in enumerate(['x', 'y', 'z']):
            lines += ['*Cload, amplitude=Amp-RF-%d%s\n'%(nodeLabel, direction),
                '%s.%d, %s, 1\n'%(instanceName, nodeLabel, i+1),
                '*Amplitude, name=Amp-RF-%d%s, time=TOTAL TIME\n'%(nodeLabel, direction)]
            lines += ['0., 0., %f, %f\n'%(step.timePeriod, RF[i])]
    with open(outputFileName, 'w') as f:
        f.writelines(lines)
    return

if __name__ == '__main__':
    odbName = 'Istanbul_model_static.odb'
    writeRF(odbName=odbName)