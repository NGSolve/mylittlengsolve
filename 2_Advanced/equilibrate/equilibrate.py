from ngsolve import *
from libequilibrate import *

def Equilibrate(flux, source, fes, order=None):
    """
    Computes an equilibrated flux using the method by Braess, Pillwein and Schoeberl.
    input:
         flux  ... vectorial coefficient function
         source .. coefficient function
         fes  .... finite element space providing Dirichlet bc 
         order ... order of reconstructed flux
    output:
         an HDiv finite element function eqflux satisfying
             div equflux = -f  

    flux and source must be piecewise polynomials, and satisfy
    the lowest order Galerkin condition:
       \int flux grad(phi) = \int f phi    forall phi in V_h

    fes provides boundary condition of primal solution
    """
    
    if (order==None): order=fes.globalorder

    mesh = fes.mesh
    Xsigma = HDiv(mesh, order=order, discontinuous=True)
    Xw     = L2(mesh,order=order-1)
    Xwf    = FESpace("facet", mesh, order=order)
    Xlam   = FESpace("number",mesh)
    X = FESpace([Xsigma,Xw, Xwf, Xlam])

    sigma, w, wf, lam = X.TrialFunction()
    tau,   v, vf, mu  = X.TestFunction()

    n = specialcf.normal(mesh.dim)

    aequ = BilinearForm(X)
    aequ += SymbolicBFI(tau*sigma + w*div(tau) + v*div(sigma) + w*mu + lam*v)
    aequ += SymbolicBFI(-sigma*n*vf-tau*n*wf, element_boundary=True)
    aequ.Assemble()

    sigmacorr = GridFunction(Xsigma)
    EquilibratePatches (flux, source, aequ, sigmacorr, fes)
            
    fesflux = HDiv(fes.mesh, order=order, discontinuous=False)
    hdivflux = GridFunction(fesflux, name="equflux")
    hdivflux.Set(flux-sigmacorr)
    return hdivflux
