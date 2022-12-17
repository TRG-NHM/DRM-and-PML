# -*- coding: mbcs -*-
from abaqus import *
from abaqusConstants import *
import __main__
# ===== Abaqus Viewer does not contain these modules =====
# import section
# import regionToolset
# import displayGroupMdbToolset as dgm
# import part
# import material
# import assembly
# import optimization
# import step
# import interaction
# import load
# import mesh
# import job
# import sketch
# import connectorBehavior
# ===== =====
import visualization
import xyPlot
import displayGroupOdbToolset as dgo

def set_frame_to_maximum():
    # ===== Read the .odb file opened in current viewport, so make sure to open the correct .odb file in Abaqus CAE =====
    vp = session.viewports[session.currentViewportName]
    odbName = vp.odbDisplay.name
    # ///note that this script can only handle the .odb file with only one step
    active_frame = session.odbData[odbName].activeFrames[0][1][-1]
    if type(active_frame) is not int:
        stepName = session.odbData[odbName].activeFrames[0][0]
        active_frame = len(session.odbData[odbName].steps[stepName].frames)-1

    maxValue = -float('Inf')
    maxFrame = 0
    for index in xrange(active_frame+1):
        vp.odbDisplay.setFrame(step=0, frame=index )
        autoMaxValue = vp.odbDisplay.contourOptions.autoMaxValue
        if autoMaxValue > maxValue:
            maxValue = autoMaxValue
            maxFrame = index

    # print 'maxValue = '+str(maxValue)+', frame = '+str(maxFrame)
    vp.odbDisplay.setFrame(step=0, frame=maxFrame )

if __name__ == '__main__':
    set_frame_to_maximum()