# -*- coding: mbcs -*-
# Do not delete the following import lines
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
# import step
# import interaction
# import load
# import mesh
# import optimization
# import job
# import sketch
# import connectorBehavior
# ===== =====
import visualization
import xyPlot
import displayGroupOdbToolset as dgo

def getNumberOfActiveFrames():
    # ===== Read the .odb file opened in current viewport, so make sure to open the correct .odb file in Abaqus CAE =====
    vp = session.viewports[session.currentViewportName]
    odbName = vp.odbDisplay.name
    # NOTE: this script can only handle the .odb file with only one step
    num_active_frame = session.odbData[odbName].activeFrames[0][1][-1]
    if type(num_active_frame) is not int:
        stepName = session.odbData[odbName].activeFrames[0][0]
        num_active_frame = len(session.odbData[odbName].steps[stepName].frames)-1
    return num_active_frame

if __name__ == '__main__':
    targetVideoDuration = 7 # sec
    numActiveFrames = getNumberOfActiveFrames()
    suggestedFrameRate = int(numActiveFrames/targetVideoDuration)

    session.linkedViewportCommands.setValues(linkViewports=True)
    for viewport in session.viewports.values():
        viewport.animationController.setValues(animationType=TIME_HISTORY)
    viewport = session.viewports[session.currentViewportName]
    viewport.animationController.play(duration=UNLIMITED)
    viewport.animationController.stop()

    session.imageAnimationOptions.setValues(vpDecorations=OFF, vpBackground=ON, compass=ON, timeScale=1, frameRate=suggestedFrameRate)
    # session.writeImageAnimation(fileName='output', format=QUICKTIME, canvasObjects=(
    #     viewport2, viewport))
    session.writeImageAnimation(fileName='output', format=QUICKTIME, canvasObjects=session.viewports.values())