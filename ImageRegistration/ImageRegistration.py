from __main__ import vtk, qt, ctk, slicer

from slicer.ScriptedLoadableModule import *

#
# ImageRegistration
#

class ImageRegistration:
  def __init__(self, parent):
    parent.title = "Image Registration"
    parent.categories = ["Registration"]
    parent.dependencies = []
    parent.contributors = ["Pietro Nardelli (UCC)",]
    parent.helpText = """
    Register two images using roll rotation
    """
    parent.acknowledgementText = """
Developed by Steve Pieper, Isomics, Inc.,
partially funded by NIH grant 3P41RR013218-12S1 (NAC) and is part of the National Alliance
for Medical Image Computing (NA-MIC), funded by the National Institutes of Health through the
NIH Roadmap for Medical Research, Grant U54 EB005149."""
    self.parent = parent

#
# ImageRegistrationWidget
#

class ImageRegistrationWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent

    self.layout = self.parent.layout()

    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):

    # Collapsible button
    self.selectionCollapsibleButton = ctk.ctkCollapsibleButton()
    self.selectionCollapsibleButton.text = "Selection"
    self.layout.addWidget(self.selectionCollapsibleButton)

    # Layout within the collapsible button
    self.formLayout = qt.QFormLayout(self.selectionCollapsibleButton)

    #
    # the volume selectors
    #
    self.fixedFrame = qt.QFrame(self.selectionCollapsibleButton)
    self.fixedFrame.setLayout(qt.QHBoxLayout())
    self.formLayout.addWidget(self.fixedFrame)
    self.fixedSelector = qt.QLabel("Fixed Image: ", self.fixedFrame)
    self.fixedFrame.layout().addWidget(self.fixedSelector)
    self.fixedSelector = slicer.qMRMLNodeComboBox(self.fixedFrame)
    self.fixedSelector.nodeTypes = ( ("vtkMRMLVectorVolumeNode"), "" )
    self.fixedSelector.addEnabled = False
    self.fixedSelector.removeEnabled = False
    self.fixedSelector.setMRMLScene( slicer.mrmlScene )
    self.fixedFrame.layout().addWidget(self.fixedSelector)

    self.movingFrame = qt.QFrame(self.selectionCollapsibleButton)
    self.movingFrame.setLayout(qt.QHBoxLayout())
    self.formLayout.addWidget(self.movingFrame)
    self.movingSelector = qt.QLabel("Moving Image: ", self.movingFrame)
    self.movingFrame.layout().addWidget(self.movingSelector)
    self.movingSelector = slicer.qMRMLNodeComboBox(self.movingFrame)
    self.movingSelector.nodeTypes = ( ("vtkMRMLVectorVolumeNode"), "" )
    self.movingSelector.setMRMLScene( slicer.mrmlScene )
    self.movingSelector.addEnabled = False
    self.movingSelector.removeEnabled = False
    self.movingFrame.layout().addWidget(self.movingSelector)

    '''self.outputFrame = qt.QFrame(self.selectionCollapsibleButton)
    self.outputFrame.setLayout(qt.QHBoxLayout())
    self.formLayout.addWidget(self.outputFrame)
    self.outputSelector = qt.QLabel("Output Image: ", self.outputFrame)
    self.outputFrame.layout().addWidget(self.outputSelector)
    self.outputSelector = slicer.qMRMLNodeComboBox(self.outputFrame)
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.addEnabled = True
    self.outputSelector.renameEnabled = True
    self.outputFrame.layout().addWidget(self.outputSelector)'''

    self.numberOfAnglesBox = qt.QSpinBox()
    self.numberOfAnglesBox.setMinimum(1)
    self.numberOfAnglesBox.setMaximum(360)
    self.numberOfAnglesBox.setValue(36)

    self.formLayout.layout().addWidget(self.numberOfAnglesBox)

    # Apply button
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Register the two images."
    self.formLayout.addWidget(self.applyButton)
    self.applyButton.connect('clicked(bool)', self.onApply)

    # Add vertical spacer
    self.layout.addStretch(1)

  def onApply(self):
    fixedImage = self.fixedSelector.currentNode()
    movingImage = self.movingSelector.currentNode()
    #outputImage = self.outputSelector.currentNode()

    fixedExtract = vtk.vtkImageExtractComponents()
    fixedExtract.SetComponents(0,1,2)
    fixedLuminance = vtk.vtkImageLuminance()
    fixedExtract.SetInputConnection(fixedImage.GetImageDataConnection())
    fixedLuminance.SetInputConnection(fixedExtract.GetOutputPort())
    fixedLuminance.Update()

    fixedIJKToRAS = vtk.vtkMatrix4x4()
    fixedImage.GetIJKToRASMatrix(fixedIJKToRAS)
    fixedScalarVolume = slicer.vtkMRMLScalarVolumeNode()
    fixedScalarVolume.SetName('fixedScalarImage')
    fixedScalarVolume.SetIJKToRASMatrix(fixedIJKToRAS)
    fixedScalarVolume.SetImageDataConnection(fixedLuminance.GetOutputPort())
    slicer.mrmlScene.AddNode(fixedScalarVolume)

    movingExtract = vtk.vtkImageExtractComponents()
    movingExtract.SetComponents(0,1,2)
    movingLuminance = vtk.vtkImageLuminance()
    movingExtract.SetInputConnection(movingImage.GetImageDataConnection())
    movingLuminance.SetInputConnection(movingExtract.GetOutputPort())
    movingLuminance.Update()

    movingIJKToRAS = vtk.vtkMatrix4x4()
    movingImage.GetIJKToRASMatrix(movingIJKToRAS)
    movingScalarVolume = slicer.vtkMRMLScalarVolumeNode()
    movingScalarVolume.SetName('movingScalarImage')
    movingScalarVolume.SetIJKToRASMatrix(movingIJKToRAS)
    movingScalarVolume.SetImageDataConnection(movingLuminance.GetOutputPort())
    slicer.mrmlScene.AddNode(movingScalarVolume)

    anglesNumber = self.numberOfAnglesBox.value   
    print anglesNumber
    imageRegistration = slicer.modules.imageregistrationcli
    parameters = {
          "fixedImage": fixedScalarVolume,
          "movingImage": movingScalarVolume,
          "anglesNumber": anglesNumber,
          }
    slicer.cli.run( imageRegistration,None,parameters,wait_for_completion=True )
