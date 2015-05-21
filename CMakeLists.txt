cmake_minimum_required(VERSION 2.8.9)

project(ImageRegistrationModule)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.3/Extensions/FlipImage")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Pietro Nardelli (UCC)")
set(EXTENSION_DESCRIPTION "Image registration for virtual bronchoscopy")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/ImageRegistrationModule/ImageRegistrationModule.png?revision=21745&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-ImageRegistrationModule-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-ImageRegistrationModule-screenshot-2.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(ImageRegistrationCLI)
add_subdirectory(ImageRegistration)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
