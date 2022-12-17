def getUserElementLines(userElementType: str) -> list[str]:
    # [NOTE] Info for parameters under *USER ELEMENT:
    #   TYPE: Must be 'Un' where n is a positive integer less than 10000, and it must be the same as the element type key used to identify this element on the *ELEMENT option.
    #   NODE: Number of nodes associated with an element of this type.
    #   COORDINATES: 3 for 3D, 2 for 2D.
    #   PROPERTIES: The number of properties needed as data in UEL to define such an element.
    #   VARIABLES: 360 for 3D problems and 80 for 2D problems. This is determined by Wenyang's UEL subroutine.
    #   [TODO (maybe?)] So how to calculate the number of or variables (VARIABLES) that are needed for *USER ELEMENT?
    #   Numbers in the second line: The active DoF
    # TODO: Should it be Unsymm?
    # lines = ['*User Element, Type=%s, Nodes=8, Coordinates=3, Properties=12, Variables=360, Unsymm\n'%userElementType,
    lines = ['*User Element, Type=%s, Nodes=8, Coordinates=3, Properties=12, Variables=360\n'%userElementType,
        '1, 2, 3, 21, 22, 23, 24, 25, 26\n']
    return lines

def getParametersLines(youngsModulus: float, poissonsRatio: float, density: float, PML_depth: float, lengths: list[float], 
    alpha: 0.0, beta: 0.0, elementSetName='PML', eleType_pos=0.0, afp=2.0, PML_Rcoef=1e-10) -> list[str]:
    # [NOTE] Info for parameters used for UEL
    #   E: Young's Modulus
    #   xnu: Poisson Ratio
    #   rho: Density
    #   EleType_pos: 0.0 (A trivial property)
    #   PML_L: Depth of PML (distance from the boundary of DRM layer to the exterior of PML region)
    #   afp: Described as `m` and is equal to 2.0 in Wenyang's paper. It is the polynomial degree for alpha and beta functions.
    #   PML_Rcoef: Described as `R` and is equal to 10^(-10) in Wenyang's paper. It is a user-tunable reflection coefficient.
    #   RD_half_width_x: Half width of the domain in x direction without PML layers (i.e., the distance in x direction from the center to the boundary of DRM layer)
    #   RD_half_width_y: Half width of the domain in y direction without PML layers (i.e., the distance in y direction from the center to the boundary of DRM layer)
    #   RD_depth: The depth (in z direction) of the domain without PML (i.e., the depth of the interested region + the thickness of DRM layer)
    #   Damp_alpha, Damp_beta: alpha and beta used in Rayleigh Damping    
    RD_half_width_x = lengths['x']/2 - PML_depth
    RD_half_width_y = lengths['y']/2 - PML_depth
    RD_depth = lengths['z'] - PML_depth
    lines = ['*UEL Property, elset=%s\n'%elementSetName,
        '%f, %f, %f, %f, %f, %f, %e, %f,\n'%(youngsModulus, poissonsRatio, density, eleType_pos, PML_depth, afp, PML_Rcoef, RD_half_width_x),
        '%f, %f, %f, %f\n'%(RD_half_width_y, RD_depth, alpha, beta)]
    return lines