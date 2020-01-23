
#include <solve.hpp>

#include <nginterface.h>
 

namespace ngsolve
{
  using namespace ngsolve;

#include "mbs.hpp"


  void NumProcDynamic::ComputeU(const FlatVector<double> & point,
				FlatVector<double> & up,
				const S_GridFunction<double>& bu,
				LocalHeap & lh)
  {
    int elnr;
    double lami[3]; //lokale koordinaten
    Array<int> dnums;
    ElementTransformation eltrans;

    elnr = Ng_FindElementOfPoint (const_cast<double*>(&point(0)), lami);

    if (!elnr) return;
    elnr--;

    const S_GridFunction<double> & u =
      dynamic_cast<const S_GridFunction<double>&> (bu);

    const FESpace & fes = u.GetFESpace();
    const FiniteElement & fel = fes.GetFE (elnr, lh);
    ma.GetElementTransformation (elnr, eltrans, lh);
    fes.GetDofNrs (elnr, dnums);

    FlatVector<double> elu(dnums.Size() * fes.GetDimension(),lh);

    u.GetElementVector (dnums, elu);
    //fes.TransformVec (elnr, 0, elu, TRANSFORM_SOL);


    if (fel.ElementType() == ET_TRIG)
      {
	double hlam1 = lami[0];
	double hlam2 = lami[1];
	lami[0] = 1-hlam1-hlam2;
	lami[1] = hlam1;
      }

    if (fel.ElementType() == ET_TET)
      {
	double hlam1 = lami[0];
	double hlam2 = lami[1];
	double hlam3 = lami[2];
	lami[0] = 1-hlam1-hlam2-hlam3;
	lami[1] = hlam1;
	lami[2] = hlam2;
      }

    IntegrationPoint ip(lami, 1);
    FlatVector<double> shape(dnums.Size(),lh);

    const ScalarFiniteElement<3> & nfel =
      dynamic_cast<const ScalarFiniteElement<3> &>(fel);

    nfel.CalcShape(ip,shape);
    //FlatVector<double> sum(D);
    up = 0;

    int i,j;
    for (i=0; i<=DIM; i++)
      {
        for (j=0; j<dnums.Size(); j++)
	  up(i)+=shape(j)*elu(DIM*j+i);
      }
  }




  NumProcDynamic :: NumProcDynamic (PDE & apde, const Flags & flags)
    : NumProc (apde), points(0) // , masslin(0,0), invmasslin(0,0), cnts2bnds(0,0)
  {
    int i;

    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", NULL));
    bfm = pde.GetBilinearForm (flags.GetStringFlag ("bilinearformm", NULL));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunctionu", NULL));
    gfusmall = pde.GetGridFunction (flags.GetStringFlag ("gridfunctionusmall", NULL));
    gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunctionv", NULL));
    dt = flags.GetNumFlag ("dt", 1e-5);
    te = flags.GetNumFlag ("te", 0.1);

    string timeint;
    timeint = flags.GetStringFlag ("TI", "R2");
    tit = TI_R2; //default
    if (timeint == "R2") tit = TI_R2;
    else if (timeint == "R3") tit = TI_R3;
    else if (timeint == "R4") tit = TI_R4;
    else if (timeint == "MP") tit = TI_MP;
    else if (timeint == "IE") tit = TI_IE;
    else if (timeint == "L2") tit = TI_L2;
    else if (timeint == "L3") tit = TI_L3;
    else cout << "time integration method not recognized, using RadauIIA,2stage" << endl;
    cout << "Time-Integration method=" << timeint << ", id=" << tit << endl;

    // connection array  0..N-1
    const Array<double> & pcon = flags.GetNumListFlag ("connections");
    cout << pcon.Size()/2 << " connections:" << endl;
    connections.SetSize (pcon.Size()/2);
    for (i = 0; i < pcon.Size(); i+=2)
      {
	connections[i/2][0] = int(pcon[i])-1;
	connections[i/2][1] = int(pcon[i+1])-1;
	cout << pcon[i] << " - " << pcon[i+1] << endl;
      }

    const Array<double> & pconx = flags.GetNumListFlag ("connectionsx");
    cout << pconx.Size()/2 << " connectionsx:" << endl;
    connectionsx.SetSize (pconx.Size()/2);
    for (i = 0; i < pconx.Size(); i+=2)
      {
	connectionsx[i/2][0] = int(pconx[i])-1;
	connectionsx[i/2][1] = int(pconx[i+1])-1;
	cout << pconx[i] << " - " << pconx[i+1] << endl;
      }

    const Array<double> & pcony = flags.GetNumListFlag ("connectionsy");
    cout << pcony.Size()/2 << " connectionsy:" << endl;
    connectionsy.SetSize (pcony.Size()/2);
    for (i = 0; i < pcony.Size(); i+=2)
      {
	connectionsy[i/2][0] = int(pcony[i])-1;
	connectionsy[i/2][1] = int(pcony[i+1])-1;
	cout << pcony[i] << " - " << pcony[i+1] << endl;
      }
    const Array<double> & pconz = flags.GetNumListFlag ("connectionsz");
    cout << pconz.Size()/2 << " connectionsz:" << endl;
    connectionsz.SetSize (pconz.Size()/2);
    for (i = 0; i < pconz.Size(); i+=2)
      {
	connectionsz[i/2][0] = int(pconz[i])-1;
	connectionsz[i/2][1] = int(pconz[i+1])-1;
	cout << pconz[i] << " - " << pconz[i+1] << endl;
      }
  }

  NumProcDynamic :: ~NumProcDynamic()
  {
    ;
  }

  void NumProcDynamic :: Do(LocalHeap & lh)
  {
    cout.precision(4);
    cout << "call MBS solver" << endl;
    int i, j, k, ii, jj, rk;

    const MeshAccess & ma = pde.GetMeshAccess();

    int nd = pde.GetFESpace ("v") -> GetNDof();
    int nv = ma.GetNV();
    int ne = ma.GetNE();
    int nse = ma.GetNSE();

    int rksteps; //stages of runge-kutta formula
    //int intsteps = 10;
    if (tit == TI_IE) {rksteps = 1;}
    else if (tit == TI_MP || tit == TI_R2 || tit == TI_L2) {rksteps = 2;}
    else if (tit == TI_R3 || tit == TI_L3) {rksteps = 3;}
    else if (tit == TI_R4) {rksteps = 4;}

    cout << "Runge-Kutta stages = " << rksteps << endl;

    Vector<> rktaui(rksteps);
    Vector<> rkdtaui(rksteps);
    Matrix<> rkmata(rksteps);
    Matrix<> rkmataa(rksteps);

    Matrix<> rkprea(rksteps);  // preconditioner for RK
    Matrix<> rkpreaa(rksteps);

    rkmata = 0;

    // Runge-Kutta formula for preconditioner (default: impl Euler)

    if (tit == TI_R2)
      {
        // Radau IIA, order 3
        rktaui(0) = 1.0/3.0;
        rktaui(1) = 1;
        rkmata(0, 0) = 5.0/12.0;
        rkmata(0, 1) = -1.0/12.0;
        rkmata(1, 0) = 0.75;
        rkmata(1, 1) = 0.25;
      }
    else if (tit == TI_R3)
      {
        // Radau IIA, order 5
        double sq6 = sqrt(6.);
        rktaui(0) = (4. - sq6) / 10.;
        rktaui(1) = (4. + sq6) / 10.;
        rktaui(2) = 1.;

        rkmata(0,0) = (88.-7.*sq6) / 360.;
        rkmata(0,1) = (296.-169.*sq6) / 1800.;
        rkmata(0,2) = (-2.+3.*sq6) / 225.;

        rkmata(1,0) = (296.+169.*sq6) / 1800.;
        rkmata(1,1) = (88.+7.*sq6) / 360.;
        rkmata(1,2) = (-2.-3.*sq6) / 225.;

        rkmata(2,0) = (16.-sq6) / 36.;
        rkmata(2,1) = (16.+sq6) / 36.;
        rkmata(2,2) = 1.0 / 9.0;
      }
    else if (tit == TI_R4)
      {
        rktaui(0) = 0.088587959512703947;
        rktaui(1) = 0.40946686444073471; 
        rktaui(2) = 0.78765946176084706; 
        rktaui(3) = 1;

        rkmata(0,0) = 0.11299947932315619;
        rkmata(0,1) = -0.040309220723522206;
        rkmata(0,2) = 0.025802377420336391;
        rkmata(0,3) = -0.0099046765072664239;

        rkmata(1,0) = 0.23438399574740026;
        rkmata(1,1) = 0.2068925739353589;
        rkmata(1,2) = -0.047857128048540719;
        rkmata(1,3) = 0.016047422806516273;

        rkmata(2,0) = 0.21668178462325034;
        rkmata(2,1) = 0.40612326386737331;
        rkmata(2,2) = 0.18903651817005634;
        rkmata(2,3) = -0.02418210489983294;

        rkmata(3,0) = 0.22046221117676838;
        rkmata(3,1) = 0.38819346884317188;
        rkmata(3,2) = 0.32884431998005974;
        rkmata(3,3) = 0.0625;
      }
    else if (tit == TI_IE)
      {
        // impl Euler
        rktaui(0) = 1;
        rkmata(0, 0) = 1;
      }
    else if (tit == TI_MP)
      {
        // midpoint
        rktaui(0) = 0.5;
        rktaui(1) = 1;
        rkmata(0, 0) = 0.5;
        rkmata(0, 1) = 0.;
        rkmata(1, 0) = 0.;
        rkmata(1, 1) = 1.;
      }
    else if (tit == TI_L2)
      {
        // trapez
        rktaui(0) = 0;
        rktaui(1) = 1;
        rkmata(0, 0) = 0;
        rkmata(0, 1) = 0;
        rkmata(1, 0) = 0.5;
        rkmata(1, 1) = 0.5;
      }
    else if (tit == TI_L3)
      {
        // LobattoIIIa, order 4
        rktaui(0) = 0;
        rktaui(1) = 0.5;
        rktaui(2) = 1;
        rkmata(0, 0) = 0;
        rkmata(0, 1) = 0;
        rkmata(0, 2) = 0;
        rkmata(1, 0) = 5.0/24.0;
        rkmata(1, 1) = 1.0/3.0;
        rkmata(1, 2) = -1.0/24.0;
        rkmata(2, 0) = 1.0/6.0;
        rkmata(2, 1) = 2.0/3.0;
        rkmata(2, 2) = 1.0/6.0;
      }

    // R-K formula for preconditioner (impl Euler):
    cout << "rkmata = " << rkmata << endl;

    rkprea = 0;
    for (i = 0; i < rksteps; i++)
      {
	if (i == 0)
	  rkprea(0,0) = rktaui(0);
	else
	  {
	    rkprea(i,i) = rktaui(i) - rktaui(i-1);
	    for (j = 0; j < i; j++)
	      rkprea(i,j) = rkprea(i-1, j);
	  }
      }
    cout << "rkprea = " << rkprea << endl;


    rkdtaui(0) = rktaui(0);
    for (i = 1; i < rksteps; i++)
      rkdtaui(i) = rktaui(i) - rktaui(i-1);


    // check Runge Kutta consistency:
    for (i = 1; i <= 8; i++)
      {
        for (j = 0; j < rksteps; j++)
          {
            double sum = pow (rktaui(j), double(i)) / double(i);
            for (k = 0; k < rksteps; k++)
              sum -= rkmata(j, k) * pow (rktaui(k), double(i-1));
            cout << "order " << i << ", consistent = " << sum << endl;
          }
      }

    rkmataa = rkmata * rkmata;
    rkpreaa = rkprea * rkprea;


    MBS_System mbssystem;
    PrepareSystem(mbssystem);

    VVector<Vec<DIM> > & vecu =
      dynamic_cast<VVector<Vec<DIM> >&> (gfu->GetVector());

    VVector<Vec<DIM> > & vecu0 =
      dynamic_cast<VVector<Vec<DIM> >&> 
      (pde.GetGridFunction("u0")->GetVector());

    VVector<Vec<DIM> > & vecusmall =
      dynamic_cast<VVector<Vec<DIM> >&> (gfusmall->GetVector());

    VVector<Vec<DIM> > & vecv = 
      dynamic_cast<VVector<Vec<DIM> >& > (gfv->GetVector());
    
    VVector<Vec<DIM> > & veca = 
      dynamic_cast<VVector<Vec<DIM> >&>
      (pde.GetGridFunction("veca")->GetVector());
    
    const VVector<Vec<DIM> > & vecf = 
      dynamic_cast<const VVector<Vec<DIM> >&>
      (pde.GetLinearForm("f")->GetVector());


    
    VVector<Vec<DIM> > hv1(nd);
    VVector<Vec<DIM> > hv2(nd);
    VVector<Vec<DIM> > hv3(nd);
    VVector<Vec<DIM> > hv4(nd);

    VVector<Vec<DIM> > vecunew(nd);
    VVector<Vec<DIM> > veculin(nd);
    VVector<Vec<DIM> > vecanew(nd);
    VVector<Vec<DIM> > vecf2(nd);
    VVector<Vec<DIM> > vecfnew(nd);
    VVector<Vec<DIM> > vecfold(nd);

    vecf2 = 0;
    vecfnew = 0;
    vecfold = 0;


    ofstream dynout ("dynamic.out");
    ofstream dynout2 ("dynamic2.out");
    

    const BaseSparseMatrix & mata =
      dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());

    const BaseSparseMatrix & matm = 
      dynamic_cast<const BaseSparseMatrix&> (bfm->GetMatrix());



    BaseSparseMatrix & hmat = const_cast<BaseSparseMatrix&> 
      (dynamic_cast<const BaseSparseMatrix&> 
       (pde.GetBilinearForm ("hmat")->GetMatrix()));


    Array<BaseMatrix*> invhmats(rksteps);
    for (i = 0; i < rksteps; i++)
      {
	hmat = 0;
	hmat.AsVector() += matm.AsVector();
	hmat.AsVector() += (dt * dt * rkpreaa(i,i)) * mata.AsVector();
	invhmats[i] = hmat.InverseMatrix();
      }


    Array<BaseMatrix*> invhmats2(rksteps);
    for (i = 0; i < rksteps; i++)
      {
	hmat = 0;
	hmat.AsVector() += matm.AsVector();
	hmat.AsVector() += (dt * dt * rkdtaui(i) * rkdtaui(i) / (rksteps * rksteps)) * mata.AsVector();
	//      hmat.AsVector() += (dt * dt * rkdtaui(i) * rkdtaui(i) / (intsteps * intsteps)) * mata.AsVector(); //hannes

	invhmats2[i] = hmat.InverseMatrix();
      }


    BaseMatrix & invm = *matm.InverseMatrix();

    /*
      BaseMatrix & hmat = const_cast<BaseMatrix&> 
      (pde.GetBilinearForm ("hmat")->GetMatrix());

      double alpha = 0;
      double beta = sqr (1 - alpha) / 4;
      double gamma = 0.5 - alpha;

      hmat.SetScalar (0);
      hmat.Add (1, matm);
      hmat.Add (dt * dt * beta * (1+alpha), mata);

      BaseMatrix & invhmat = *hmat.InverseMatrix();
      BaseMatrix & invm = *matm.InverseMatrix();
    */







    // constraints:
    // 
    // \int_Gamma_i v ds  =   bnd_integral(i) . v
    // for each boundary part i
    // 
    Array<VVector<double>* > bnd_integral(nbnds);
    Array<double> bcoeffs(nbnds);

    const FESpace & fespace = *pde.GetFESpace ("vscal");
    for (i = 0; i < nbnds; i++)
      {
	bcoeffs = 0;
	bcoeffs[i] = 1;
        
        Flags flags;
	T_LinearForm<double> lf (fespace, "lf", flags);
	lf.AddIntegrator (new NeumannIntegrator<DIM>
			  (new DomainConstantCoefficientFunction (bcoeffs)));
	lf.Assemble(lh);
	bnd_integral[i] = new VVector<double> (nd);
	(*bnd_integral[i]) = 1.0 * lf.GetVector();
      }

    // compute \int_{Gamma_i} 1 ds
    Array<double> bnd_integral1(nbnds);

    VVector<double> shv1 (nd);
    shv1 = 1.0;
    for (i = 0; i < nbnds; i++)
      {
	bnd_integral1[i] = InnerProduct (*bnd_integral[i], shv1);
	cout << "bnd_integral1 (" << i << ") = " << bnd_integral1[i] << endl;
      }



    // compute constraint-functionals (such as ux_\Gamma 3 - ux_\Gamma 4 = 0)

    // number of constraints
    int ncnts = DIM * connections.Size() + connectionsx.Size()
      + connectionsy.Size() + connectionsz.Size();
  
    cnts.SetSize (ncnts);
    //   Array<VVector<Vec<DIM> >*> minvcnts(ncnts); ????

    // boundary integral to constraint
    cnts2bnds.SetSize (DIM*nbnds, ncnts);
    cnts2bnds = 0;

    for (i = 0; i < ncnts; i++)
      {
	cnts[i] = new VVector<Vec<DIM> > (nd);
	*cnts[i] = 0;
	//      minvcnts.Elem(i) = dynamic_cast<VVector<Vec<DIM> >*>(vecu.Copy()); ?????
      }

    for (j = 0; j < connections.Size(); j++)
      {
	int f1 = connections[j][0];
	int f2 = connections[j][1];
	for (i = 0; i < nd; i++)
	  {
	    double val = 0;
	    if (f1 != -1)
	      val =
		(*bnd_integral[f1])(i) / bnd_integral1[f1];
	    if (f2 != -1)
	      val -=
		(*bnd_integral[f2])(i)  / bnd_integral1[f2];
	  
	    if (DIM == 3)
	      (*cnts[3*j])(i)(0) =
		(*cnts[3*j+1])(i)(1) =
		(*cnts[3*j+2])(i)(2) = val;
	    else
	      (*cnts[2*j])(i)(0) =
		(*cnts[2*j+1])(i)(1) = val;
	  }	

	for (jj = 0; jj < DIM; jj++)
	  {
	    if (f1 != -1)
	      cnts2bnds(DIM*f1+jj, DIM*j+jj) = 1/bnd_integral1[f1];
	    if (f2 != -1)
	      cnts2bnds(DIM*f2+jj, DIM*j+jj) = -1/bnd_integral1[f2];
	  }
      }

    int hcntnr = DIM * connections.Size();
    for (k = 0; k < DIM; k++)
      {
	int nconxk;
	switch (k)
	  {
	  case 0: nconxk = connectionsx.Size(); break;
	  case 1: nconxk = connectionsy.Size(); break;
	  case 2: nconxk = connectionsz.Size(); break;
	  }

	for (j = 0; j < nconxk; j++)
	  {
	    int f1, f2;
	    switch (k)
	      {
	      case 0:
		f1 = connectionsx[j][0];
		f2 = connectionsx[j][1];
		break;
	      case 1:
		f1 = connectionsy[j][0];
		f2 = connectionsy[j][1];
		break;
	      case 2:
		f1 = connectionsz[j][0];
		f2 = connectionsz[j][1];
		break;
	      }
	  
	    for (i = 0; i < nd; i++)
	      {
		double val = 0;
		if (f1 != -1)
		  val = (*bnd_integral[f1])(i) / bnd_integral1[f1];
		if (f2 != -1)
		  val -=
		    (*bnd_integral[f2])(i) / bnd_integral1[f2];
	      
		(*cnts[hcntnr])(i)(k) = val;
	      }	
	    if (f1 != -1)
	      cnts2bnds(DIM*f1+k, hcntnr) = 1/bnd_integral1[f1];
	    if (f2 != -1)
	      cnts2bnds(DIM*f2+k, hcntnr) = -1/bnd_integral1[f2];

	    hcntnr++;
	  }
      }

    (*testout) << "constraint vectors:" << endl
	       << cnts << endl;

    // precompute schur complement without rotations:
    Array<Matrix<>*> bnd_int_schurs (rksteps);

    for (rk = 0; rk < rksteps; rk++)
      {
	bnd_int_schurs[rk] = new Matrix<> (DIM*nbnds);
	Matrix<> & bnd_int_schur = *bnd_int_schurs[rk];
	for (i = 0; i < nbnds; i++)
	  for (ii = 0; ii < DIM; ii++)
	    {
	      hv1 = 0;
	      for (k = 0; k < nd; k++)
		hv1(k)(ii) = (*bnd_integral[i])(k);
	      hv2 = (*invhmats[rk]) * hv1;
	    
	      for (j = 0; j < nbnds; j++)
		for (jj = 0; jj < DIM; jj++)
		  {
		    double sum = 0;
		    for (k = 0; k < nd; k++)
		      sum += hv2(k)(jj) * (*bnd_integral[j]) (k);

		    bnd_int_schur(DIM*i+ii, DIM*j+jj) = sum;
		  }
	    }

	(*testout) << "bnd_int_schur[" << rk << "] = " << endl
		   << bnd_int_schur << endl;
      }


    // precompute schur complement without rotations:
    Array<Matrix<>*> bnd_int_schurs2 (rksteps);

    for (rk = 0; rk < rksteps; rk++)
      {
	bnd_int_schurs2[rk] = new Matrix<> (DIM*nbnds);
	Matrix<> & bnd_int_schur = *bnd_int_schurs2[rk];
	for (i = 0; i < nbnds; i++)
	  for (ii = 0; ii < DIM; ii++)
	    {
	      hv1 = 0;
	      for (k = 0; k < nd; k++)
		hv1(k)(ii) = (*bnd_integral[i])(k);
	      hv2 = (*invhmats2[rk]) * hv1;
	    
	      for (j = 0; j < nbnds; j++)
		for (jj = 0; jj < DIM; jj++)
		  {
		    double sum = 0;
		    for (k = 0; k < nd; k++)
		      sum += hv2(k)(jj) * (*bnd_integral[j])(k);
		    bnd_int_schur(DIM*i+ii, DIM*j+jj) = sum;
		  }
	    }

	(*testout) << "bnd_int_schur2[" << rk << "] = " << endl
		   << bnd_int_schur << endl;
      }


    Vector<> lami(ncnts);
    Array<Vector<>*> lamis(rksteps);
    for (rk = 0; rk < rksteps; rk++)
      lamis[rk] = new Vector<> (ncnts);
 
    vecu = 0;
    vecv = 0;
    veca = 0;
    /*
      invm.Mult (vecf, veca);
    */

    CalcRMatrix (vecu, mbssystem);
    /*
      veca = (*invhmats[0]) * vecf;
      Project2Constraints (mbssystem, *invhmats[0],
      *bnd_int_schurs[0], veca, *lamis[0]);
      */
    //  veca *= -0.5;

    veca = 0;

    int cnt = 0;
    double t = 0;

    double energycor = 0;
    double sumkin = 0, sumpot = 0;

    double tpot;

    // Runge Kutta with Euler preconditioner:

    Array<VVector<Vec<DIM> >*> vecanews(rksteps);
    for (rk = 0; rk < rksteps; rk++)
      {
	vecanews[rk] = new VVector<Vec<DIM> > (nd);
	*vecanews[rk] = 0;
      }

    Array<VVector<Vec<DIM> >*> vecress(rksteps);
    for (rk = 0; rk < rksteps; rk++)
      vecress[rk] = new VVector<Vec<DIM> > (nd);

    Array<VVector<Vec<DIM> >*> veccorrs(rksteps);
    for (rk = 0; rk < rksteps; rk++)
      veccorrs[rk] = new VVector<Vec<DIM> > (nd);

    Array<VVector<Vec<DIM> >*> vecunews(rksteps);
    for (rk = 0; rk < rksteps; rk++)
      {
	vecunews[rk] = new VVector<Vec<DIM> > (nd);
	*vecunews[rk] = 0;
      }

    VVector<Vec<DIM> > hvecu(nd);
    VVector<Vec<DIM> > hvecv(nd);

    cout << "rktaui = " << endl << rktaui << endl;
    cout << "rkmat = " << endl << rkmataa << endl;
    cout << "rkpre = " << endl << rkpreaa << endl;


    ofstream def("kurbel_mbs.out");
    ofstream def2("kurbel_mbs.dat");
    ofstream def3("kurbel_mbs.drift");
    LocalHeap lh(50000);

    double t_kraftaus = 10000000.*0.2; //0.1000; // kurbel:0.5

    double acttime = 0;
    double steptime = 0;

    clock_t starttimea;
    starttimea = clock();

    int max_it_RK = 100;
    int max_it_pre = 20;

    double tol_pre = 1e-6; //orig 1e-6
    double tol_RK = 1e-2; //orig: 1e-1

    double rkpre_needed;

    while (t <= te)
      {
	lh.CleanUp();

	cout.precision(4);
	cout << "t = " << t << endl;

        for (rk = 0; rk < rksteps; rk++)
	  (*vecanews[rk]) = veca;

	// nonlinear preconditioning step
	for (rk = 0; rk < rksteps; rk++)
	  {
	    (*vecanews[rk]) = veca;
	    double res_max = 1e10;
	    k = 0;

	    while (k < max_it_pre && res_max>tol_pre)
	      {
		k++;
		vecunew = vecu;
		vecunew += (dt * rktaui(rk)) * vecv;


		res_max = 0;
		for (j = 0; j < rksteps; j++)
		  vecunew += (dt*dt*rkpreaa(rk,j)) * (*vecanews[j]);

		*vecunews[rk] = vecunew;

		// compute f(vecunew)

		CalcLinearApproximation (vecunew, veculin);
		CalcRMatrix (veculin, mbssystem);

		ComputePotentialPrime (mbssystem, vecunew, veculin, vecfnew, 1); // f=Ku

		if (t < t_kraftaus)
		  vecfnew -= vecf;
		else if (fabs(t - t_kraftaus) < 1e-2) cout << "kraft aus" << endl;
		// -b + K u


		hv2 = vecfnew + matm * (*vecanews[rk]);

		ApplyRotation (hv2, hv1, mbssystem, 1);
		hv2 = (*invhmats[rk]) * hv1;
		ApplyRotation (hv2, hv1, mbssystem, 0);

		Project2Constraints (mbssystem, *invhmats[rk],
				     *bnd_int_schurs[rk],
				     hv1, lami);

		//if (rktaui(rk) > 0) ...Hannes, wegen L2/L3
		(*vecanews[rk]) -= hv1;
		//} Hannes: check residuum each iteration!


		veca = *vecanews[rk];

		// check residuum:
		vecunew = vecu;
		vecunew += (dt * rktaui(rk)) * vecv;

		for (j = 0; j < rksteps; j++)
		  vecunew += (dt*dt*rkpreaa(rk, j)) * (*vecanews[j]);

		CalcLinearApproximation (vecunew, veculin);
		CalcRMatrix (veculin, mbssystem);

		ComputePotentialPrime (mbssystem, vecunew, veculin, vecfnew, 1); // f=Ku
		vecfnew -= vecf;
		// -b + K u
		vecfnew += matm * (*vecanews[rk]);
		for (i = 0; i < ncnts; i++)
		  vecfnew -= lami(i) * (*cnts[i]);

		// 	  cout << " lami = " << lami << endl;
		// 	  cout << " cnts:";
		// 	  for (i = 1; i <= ncnts; i++)
		// 	    cout <<  *cnts.Get(i) * vecunew << " ";

		res_max = max(fabs(vecfnew.L2Norm()),res_max);

		cout << "\r pre steps = " << k << flush;
	      }
	    cout << "pre, rk = " << rk << ", |vecu| = " << vecunew.L2Norm();
	    cout << " |veca| = " << vecanews[rk]->L2Norm();
	    cout << " |res| = " << vecfnew.L2Norm();
	    cout << endl;
	    rkpre_needed = rk;
	  }



	k = 1;
	double max_res = 1e10;
	while ( k <= max_it_RK && max_res > tol_RK)
	  {
	    k++;
	    max_res = 0;
	    for (rk = 0; rk < rksteps; rk++)
	      {
		vecunew = vecu;
		vecunew += (dt * rktaui(rk)) * vecv;
		for (j = 0; j < rksteps; j++)
		  vecunew += (dt*dt*rkmataa(rk, j)) * (*vecanews[j]);

		*vecunews[rk] = vecunew;

		// compute f(vecunew)

		CalcLinearApproximation (vecunew, veculin);
		CalcRMatrix (veculin, mbssystem);

		// vecfnew = M a + K u - b
		tpot = ComputePotentialPrime (mbssystem, vecunew, veculin, vecfnew, 0); // f=Ku
		tpot -= InnerProduct (vecunew, vecf);

		if (t < t_kraftaus)
		  vecfnew -= vecf;

		vecfnew += matm * (*vecanews[rk]);
		*vecress[rk] = vecfnew;
	      }

	    for (rk = 0; rk < rksteps; rk++)
	      {
		vecunew = *vecunews[rk];

		CalcLinearApproximation (vecunew, veculin);
		CalcRMatrix (veculin, mbssystem);

		hv1 = 0;
		for (j = 0; j < rk; j++)
		  hv1 += (dt * dt * rkpreaa(rk, j)) * (*veccorrs[j]);

		ApplyRotation (hv1, hv2, mbssystem, 1);
		mata.Mult (hv2, hv1);
		ApplyRotation (hv1, hv2, mbssystem, 0);

		hv2 *= -1;
		hv2 += *vecress[rk];

		ApplyRotation (hv2, hv1, mbssystem, 1);
		hv2 = (*invhmats[rk]) * hv1;

		//Rotation ist neu:  ?????
		ApplyRotation (hv2, hv1, mbssystem, 0);

		Project2Constraints (mbssystem, *invhmats[rk],
				     *bnd_int_schurs[rk],
				     hv1, *lamis[rk]);

		if (rktaui(rk) > 0)
		  {
		    *veccorrs[rk] = hv1;
		    *vecanews[rk] -= hv1;
		    max_res = max(max_res, hv1.L2Norm());
		  }
		else
		  *veccorrs[rk] = 0;

	      }
 	    cout << "\rk = " << k << ", |w| = " << hv1.L2Norm() << flush;
	  }


	// check equations:
	for (rk = 0; rk < rksteps; rk++)
	  {
	    vecunew = vecu;
	    vecunew += (dt * rktaui(rk)) * vecv;

	    for (j = 0; j < rksteps; j++)
	      vecunew += (dt*dt*rkmataa(rk, j)) * (*vecanews[j]);

	    CalcLinearApproximation (vecunew, veculin);
	    CalcRMatrix (veculin, mbssystem);

	    // M a + K u - b
	    ComputePotentialPrime (mbssystem, vecunew, veculin, vecfnew, 0); // f=Ku
	    if (t < 0.1)
	      vecfnew -= vecf;
	    vecfnew += matm * (*vecanews[rk]);

	    for (i = 0; i < ncnts; i++)
	      vecfnew -= (*lamis[rk])(i) * (*cnts[i]);

	    cout << "|vecu| = " << L2Norm (*vecunews[rk]);
	    cout << ", |veca| = " << L2Norm (*vecanews[rk]);
	    //	  cout << ", lami = " << *lamis[rk];
	    //	  cout << " cnts:";
	    //	  for (i = 0; i < ncnts; i++)
	    //	    cout <<  InnerProduct (*cnts[i], *vecunews[rk]) << " ";

	    cout << ", |res| = " << L2Norm (vecfnew) << ", t_cpu=" << acttime  << ", stept_cpu=" << steptime << endl;
	    cout << endl;
	  }


	// update

	vecu = *vecunews[rksteps-1];

	for (rk = 0; rk < rksteps; rk++)
	  vecv += (dt * rkmata(rksteps-1, rk)) * (*vecanews[rk]);

	veca = *vecanews[rksteps-1];


	CalcLinearApproximation (vecu, veculin);
	CalcRMatrix (veculin, mbssystem);
	SetLargeDeformation (veculin, mbssystem);
	hvecu = vecu - veculin;
	ApplyRotation (hvecu, vecusmall, mbssystem, 1);



	t += dt;

	hv1 = matm * vecf;
	double tkin = 0.5 * InnerProduct (vecv, hv1);

	//      dynout2 << t << " " << vecu.VElem(4, 1) ;
	//<< " " << vecu.VElem(7, 2) << endl;

	/*
	  if (DIM == 2)
	  {
	  SystemVector<SysVector2d> & u2d =
	  dynamic_cast<SystemVector<SysVector2d> &> (vecu);

	  SysVector2d v1 = u2d.Elem(39);
	  v1.Add (1, SysVector2d(1,0));

	  SysVector2d v2 = u2d.Elem(73);
	  v2.Add (1, SysVector2d(0.5,0));

	  SysVector2d v1n(-v1.Elem(2), v2.Elem(1));

	  v2 *= -1;
	  v2.Add (0.5, v1);
	  double d = (v2 * v1n) / v1n.L2Norm();

	  //	  double d = u2d.Elem(2,2);
	  dynout2 << " " << d << endl;
	  }

	  if (DIM == 2)
	  {
	  u3d.SetScalar(0);
	  for (i = 1; i <= np; i++)
	  {
	  u3d.VElem(i, 1) = vecu.VElem(i,1);
	  u3d.VElem(i, 2) = vecu.VElem(i,2);
	  }
	  }
	*/

	Ng_Redraw ();

        steptime = -acttime;
        acttime = double(clock() - starttimea) / CLOCKS_PER_SEC;
	steptime += acttime;


	dynout << t << " "
	       << tpot << " "
	       << tkin << " "
	       << (tpot + tkin) << " "
	       << rkpre_needed << " "
	       << k << " "
	       << " " << steptime
	       << " " << acttime << endl;


#ifdef notused
	/*
	//slidercrank:
	FlatVector<double> pl(DIM,lh); //left support
	FlatVector<double> pr(DIM,lh); //right support
	FlatVector<double> pc(DIM,lh); //connect-point
	FlatVector<double> pc1(DIM,lh); //connect-point
	FlatVector<double> pc2(DIM,lh); //connect-point
	FlatVector<double> upl(DIM,lh);
	FlatVector<double> upr(DIM,lh);
	FlatVector<double> upc(DIM,lh);
	FlatVector<double> upc1(DIM,lh);
	FlatVector<double> upc2(DIM,lh);
	FlatVector<double> upcd(DIM,lh);

	FlatVector<double> point1(DIM,lh);
	FlatVector<double> point2(DIM,lh);
	FlatVector<double> pointm(DIM,lh);
	FlatVector<double> up1(DIM,lh);
	FlatVector<double> up2(DIM,lh);
	FlatVector<double> upm(DIM,lh);
	FlatVector<double> vm(DIM,lh);
	FlatVector<double> v1(DIM,lh);
	FlatVector<double> v2(DIM,lh);

	double off = 0.5;
	double l1 = 1.;
	double h2 = 0.015;

	pl(0)=0; pl(1)=0.025; pl(2)=0.005;
	pr(0)=1.5; pr(1)=0.025; pr(2)=0.001;
	pc(0)=0.5; pc(1)=0.025; pc(2)=0.001;
	pc1(0)=0.5; pc1(1)=0.010; pc1(2)=0.001;
	pc2(0)=0.5; pc2(1)=0.040; pc2(2)=0.001;
	ComputeU(pl, upl, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pr, upr, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pc, upc, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pc1, upc1, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pc2, upc2, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);

	upcd = 0.5*(upc1+upc2)-upc;

	point1(0)=off; point1(1)=h2; point1(2)=0.00001;
	point2(0)=l1+off; point2(1)=h2; point2(2)=0.00001;
	pointm(0)=0.5*l1+off; pointm(1)=h2; pointm(2)=0.00001;
	ComputeU(point1, up1, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(point2, up2, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pointm, upm, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);

	vm(0) = 0.5*(up2(0)-up1(0)+l1);
	vm(1) = 0.5*(up2(1)-up1(1));
	v1(0) = upm(0)-up1(0)+0.5*l1;
	v1(1) = upm(1)-up1(1);
	v2(0) = v1(0)-vm(0);
	v2(1) = v1(1)-vm(1);
	double l = sqrt(vm(0)*vm(0)+vm(1)*vm(1));
	if (l == 0) {l=1;}
	vm(0) /= l;vm(1) /= l;
	double vmx = vm(0);
	vm(0) = vm(1);
	vm(1) =-vmx;
	double wm = vm(0)*v2(0) + vm(1)*v2(1);

	double ang = atan2(up2(1)-up1(1),up2(0)-up1(0)+l1);
	if (ang > 1.) {ang-=2.*M_PI;}

	def.precision(8);
	def << t
        << " " << wm
        << " " << ang << endl;

	def2 << t << " "  << up1(0) << " " << up2(0) << " " << upm(0)
        << "  " << up1(1) << " " << up2(1) << " " << upm(1)
        << "  " << up1(2) << " " << up2(2) << " " << upm(2) << endl;

	double upcdl = sqrt(sqr(upcd(0))+sqr(upcd(1))+0*sqr(upcd(2)));
	double upll = sqrt(sqr(upl(0))+sqr(upl(1))+0*sqr(upl(2)));
	double uprl = sqrt(sqr(upr(0))+sqr(upr(1))+0*sqr(upr(2)));

	def3 << t << " "  << upcdl << " " << upll << " " << uprl
        << "  " << upcd(0) << " " << upcd(1) << " " << upcd(2)
        << "  " << upl(0) << " " << upl(1) << " " << upl(2)
        << "  " << upr(0) << " " << upr(1) << " " << upr(2)
        << endl;

	*/

	//2d-pendel:
	FlatVector<double> point1(DIM,lh);
	FlatVector<double> point2(DIM,lh);
	FlatVector<double> pointm(DIM,lh);
	FlatVector<double> up1(DIM,lh);
	FlatVector<double> up2(DIM,lh);
	FlatVector<double> upm(DIM,lh);
	FlatVector<double> vm(DIM,lh);
	FlatVector<double> v1(DIM,lh);
	FlatVector<double> v2(DIM,lh);

	//pendeld:
	double off = 0.0;
	double l1 = 3.;
	double h2 = 0.5;





	point1(0)=off; point1(1)=h2; point1(2)=0.00001;
	point2(0)=l1+off; point2(1)=h2; point2(2)=0.00001;
	pointm(0)=0.5*l1+off; pointm(1)=h2; pointm(2)=0.00001;
	ComputeU(point1, up1, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(point2, up2, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	ComputeU(pointm, upm, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);


	FlatVector<> stressp1(6,lh);
	CalcPointFlux<double> (ma, dynamic_cast<const S_GridFunction<double>&> (*gfusmall), 
			       pointm, stressp1, *bfa->GetIntegrator(0), 1, lh);

	cout << "p = " << pointm << ", stress = " << stressp1 << endl;




	vm(0) = 0.5*(up2(0)+l1);
	vm(1) = 0.5*(up2(1));
	v1(0) = upm(0)+0.5*l1;
	v1(1) = upm(1);
	v2(0) = v1(0)-vm(0);
	v2(1) = v1(1)-vm(1);
	double l = sqrt(vm(0)*vm(0)+vm(1)*vm(1));
	if (l == 0) {l=1;}
	vm(0) /= l;vm(1) /= l;
	double vmx = vm(0);
	vm(0) = vm(1);
	vm(1) =-vmx;
	double wm = vm(0)*v2(0) + vm(1)*v2(1);

	double ang = atan2(up2(1),up2(0)+l1);
	if (ang > 1.) {ang-=2.*M_PI;}

        steptime = -acttime;
        acttime = double(clock() - starttimea) / CLOCKS_PER_SEC;
	steptime += acttime;

	def.precision(8);
	def << t
	    << " " << wm
	    << " " << ang
	    << " " << steptime
	    << " " << acttime << endl;

	def2 << t << " "  << up1(0) << " " << up2(0) << " " << upm(0)
	     << "  " << up1(1) << " " << up2(1) << " " << upm(1)
	     << "  " << up1(2) << " " << up2(2) << " " << upm(2) << endl;

	

	/* kurbelwelle:

	FlatVector<double> point(DIM, lh);
	FlatVector<double> up(DIM, lh);

	point(0)=15;
	point(1)=0;
	point(2)=20;
	ComputeU(point, up, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	//cout << "x=" << x << ", u=" << up(1) << endl;

	def << t << " " << (sqrt(up(1)*up(1)+(up(2)+20)*(up(2)+20))-20)
	<< " " << atan2(up(1),up(2)+20) << endl;


	def2 << t << " " << up(0)  << " " << up(1)  << " " << up(2);

	point(0)=-14.99; point(1)=0; point(2)=0;
	ComputeU(point, up, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	def2 << " " << up(0)  << " " << up(1)  << " " << up(2);
      
	point(0)=44.99; point(1)=0; point(2)=0;
	ComputeU(point, up, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	def2 << " " << up(0)  << " " << up(1)  << " " << up(2);
      
	point(0)=15.; point(1)=0; point(2)=80; //Kolben
	ComputeU(point, up, dynamic_cast<const S_GridFunction<double>&> (*gfu),lh);
	def2 << " " << up(0)  << " " << up(1)  << " " << up(2);
      
	def2 << endl;
	*/

#endif
      }
  
    cout << "stored" << endl;
  }
  


  void NumProcDynamic :: PrintReport (ostream & ost)
  {
    ost << "NumProc Dynamic:" << endl;
  }



  double NumProcDynamic :: 
  ComputePotentialPrime (MBS_System & mbssystem,
			 const VVector<Vec<DIM> > & u, 
			 const VVector<Vec<DIM> > & ulin, 
			 VVector<Vec<DIM> > & prime,
			 int simplified)
  {
    int i, j, k;

    const MeshAccess & ma = pde.GetMeshAccess();
    const BaseMatrix & mata = bfa->GetMatrix();

    int nd = pde.GetFESpace ("v") -> GetNDof();
    int nv = ma.GetNV();
    int ne = ma.GetNE();
    int nbody = mbssystem.Size();

    VVector<Vec<DIM> > hv1(nd);
    VVector<Vec<DIM> > u0(nd);
    VVector<Vec<DIM> > usmall(nd);
    VVector<Vec<DIM> > hv2(nd);
    VVector<Vec<DIM> > hv3(nd);
    VVector<Vec<DIM> > hv4(nd);

    double pot;

    Vector<> upol1(pol1.Size()), hupol1(pol1.Size());

    CalcRMatrix (ulin, mbssystem);
    SetLargeDeformation (u0, mbssystem);
  
    hv1 = u - u0;

    ApplyRotation (hv1, usmall, mbssystem, 1);
    hv2 = mata * usmall;
    pot = 0.5 * InnerProduct(usmall, hv2);
    ApplyRotation (hv2, prime, mbssystem, 0);


    if (!simplified)
      {
	ComputeFCorrection (mbssystem, ulin, hv2, prime, hv1, hv3);

	// dual project to linear:
	for (i = 0; i < pol1.Size(); i++)
	  upol1(i) = InnerProduct (*pol1[i], hv3);
	hupol1 = invmasslin * upol1;
	hv3 = 0;
	for (i = 0; i < pol1.Size(); i++)
	  hv3 += hupol1(i) * (*mpol1[i]);
      
	prime += hv3;
      }

    return pot;
  }

  void NumProcDynamic :: 
  Project2Constraints (MBS_System & mbssystem,
		       const BaseMatrix & invhmat,
		       Matrix<> & bnd_int_schur,
		       VVector<Vec<DIM> > & v,
		       Vector<> & lami)
  {
    // compute Schur complement matrix:
    int i, ii, jj;

    const MeshAccess & ma = pde.GetMeshAccess();

    int nd = pde.GetFESpace ("v") -> GetNDof();
    int nv = ma.GetNV();
    int ne = ma.GetNE();
    int nbody = mbssystem.Size();

    int ncnts = cnts2bnds.Width();

    VVector<Vec<DIM> > hv1(nd);
    VVector<Vec<DIM> > hv2(nd);

    Vector<> hlami(ncnts);

    Matrix<> rotate_bnd(DIM*nbnds);
    rotate_bnd = 0;

    Matrix<> rot_cnts2bnds (DIM*nbnds, ncnts);
    Matrix<> hmat1(DIM*nbnds, ncnts);
    Matrix<> hmat2(ncnts, ncnts);
    Matrix<> invschur(ncnts);


    for (i = 0; i < nbnds; i++)
      {
	Mat<3,3> rot;
	int body = bnd2domain[i];
	mbssystem[body]->GetRot (rot);
	for (ii = 0; ii < DIM; ii++)
	  for (jj = 0; jj < DIM; jj++)
	    rotate_bnd(DIM*i+ii, DIM*i+jj) = rot(jj, ii);
      }

    rot_cnts2bnds = rotate_bnd * cnts2bnds;
    hmat1 = bnd_int_schur * rot_cnts2bnds;
    hmat2 = Trans (rot_cnts2bnds) * hmat1;
    /*
      Mult (rotate_bnd, cnts2bnds, rot_cnts2bnds);
      Mult (bnd_int_schur, rot_cnts2bnds, hmat1);
      CalcAtB (rot_cnts2bnds, hmat1, hmat2);
    */

    //	  (*testout) << "schur, new = " << hmat2 << endl;

    CalcInverse (hmat2, invschur);
    // compute lagragne parameters and solution:
  
    for (i = 0; i < ncnts; i++)
      hlami(i) = InnerProduct (*cnts[i], v);
    lami = invschur * hlami;

    hv1 = 0;
    for (i = 0; i < ncnts; i++)
      hv1 += lami(i) * (*cnts[i]);
  
    ApplyRotation (hv1, hv2, mbssystem, 1);
    hv1 = invhmat * hv2;
    ApplyRotation (hv1, hv2, mbssystem, 0);

    v -= hv2;
  }
				


  void NumProcDynamic :: 
  CalcLinearApproximation (const VVector<Vec<DIM> > & u, 
			   VVector<Vec<DIM> > & ulin)
  {
    Vector<> upol1(pol1.Size()), hupol1(pol1.Size());

    for (int i = 0; i < pol1.Size(); i++)
      upol1(i) = InnerProduct (*mpol1[i], u);
    hupol1 = invmasslin * upol1;
    ulin = 0;
    for (int i = 0; i < pol1.Size(); i++)
      ulin += hupol1(i) * (*pol1[i]);
  }


  void NumProcDynamic :: 
  CalcRMatrix (const VVector<Vec<DIM> > & u, 
	       MBS_System & mbssystem)
  {
    if (DIM == 3)
      {
	Mat<3,3> rot;
	Vec<3> p1, p2, p3;
      
	for (int b = 0; b < mbssystem.Size(); b++)
	  {
	  
	    MBS_Body* body = mbssystem[b];
	  
	    int pi1 = body->GetRefPoint(0);
	    int pi2 = body->GetRefPoint(1);
	    int pi3 = body->GetRefPoint(2);
      
	    pde.GetMeshAccess().GetPoint (pi1, p1);
	    pde.GetMeshAccess().GetPoint (pi2, p2);
	    pde.GetMeshAccess().GetPoint (pi3, p3);

	  
	    //Referenzvektoren berechnen
	    Vec<3> vref1, vref2, vref3;
	    vref1 = p2 - p1;
            Vec<3> v13 = p3-p1;
	    vref3 = Cross (vref1, v13);
	    vref2 = Cross (vref3, vref1);

	  
	    //Vektoren in verdrehter Lage berechnen
	    Vec<3> v1, v2, v3;
	    v1 = vref1 + u(pi2) - u(pi1);
	    Vec<3> h = (p3-p1)+u(pi3)-u(pi1);
	    v3 = Cross (v1, h);
	    v2 = Cross (v3, v1);

	    vref1 /= L2Norm (vref1);
	    vref2 /= L2Norm (vref2);
	    vref3 /= L2Norm (vref3);
	    v1 /= L2Norm (v1);
	    v2 /= L2Norm (v2);
	    v3 /= L2Norm (v3);

	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		rot (i, j) = 
		  v1(i) * vref1(j) +
		  v2(i) * vref2(j) +
		  v3(i) * vref3(j);
	  
	    body->SetRot(rot);
	  
	    body->SetTrans (u(pi1));
	    //  cout  << " phi = " << atan2 (v1.Elem(3), v1.Elem(1)) / M_PI * 180
	    //	<< ", rot=" << rot << endl;

	    (*testout) << "Body " << b << ", rot = " << rot << endl;
	  }
      }
    else
      {
	/*
	  const SystemVector<SysVector2d> & u = 
	  dynamic_cast<const SystemVector<SysVector2d> &> (hu);
      
	  SysMatrix3d rot; // just rot(1,1)
	  SysVector2d p1, p2;

	  for (b = 1; b <= mbssystem.Size(); b++)
	  {
	  
	  MBS_Body* body = mbssystem.Get(b);
	  
	  int pi1 = body->GetRefPoint(1);
	  int pi2 = body->GetRefPoint(2);
      
	  pde.GetMeshAccess().GetPoint (pi1, &p1.Elem(1));
	  pde.GetMeshAccess().GetPoint (pi2, &p2.Elem(1));

	  //Referenzvektoren berechnen
	  SysVector2d vref1 = p2;
	  vref1 -= p1;

	  //Vektoren in verdrehter Lage berechnen
	  SysVector2d v1 = vref1;
	  v1 += u.Get(pi2);
	  v1 -= u.Get(pi1);
      
	  vref1 *= 1/vref1.L2Norm();
	  v1 *= 1/v1.L2Norm();
	  
	  
	  SysVector2d v2, vref2;
	  vref2.Elem(1) = -vref1.Elem(2);
	  vref2.Elem(2) = vref1.Elem(1);
	  v2.Elem(1) = -v1.Elem(2);
	  v2.Elem(2) = v1.Elem(1);
	  
	  for (i = 1; i <= 3; i++)
	  for (j = 1; j <= 3; j++)
	  rot.Elem (i, j) = 0;
	  rot.Elem(3,3) = 1;

	  for (i = 1; i <= 2; i++)
	  for (j = 1; j <= 2; j++)
	  rot.Elem (i, j) = 
	  v1.Get(i) * vref1.Get(j) +
	  v2.Get(i) * vref2.Get(j);

	  body->SetRot(rot);
	  body->SetTrans (u.Get (pi1));
	  //  cout  << " phi = " << atan2 (v1.Elem(3), v1.Elem(1)) / M_PI * 180
	  //	<< ", rot=" << rot << endl;
	  }
	*/
      }
  }


 
  void NumProcDynamic :: PrepareSystem (MBS_System & mbssystem)
  {
    //find min-x point
    const MeshAccess & ma = pde.GetMeshAccess();

    const FESpace & fes = *pde.GetFESpace ("v");
    int nd = fes.GetNDof();
    int nv = ma.GetNV();
    int ne = ma.GetNE();
    int nse = ma.GetNSE();
    int i, j, k, l;
    int maxind = 0;
    int minind = 1000000;

    for (i = 0; i < ne; i++)
      {
	if (ma.GetElIndex(i) > maxind) {maxind = ma.GetElIndex(i);}
	if (ma.GetElIndex(i) < minind) {minind = ma.GetElIndex(i);}
      }

    int no_bodys = maxind - minind + 1;

    /*
      int nosubbodies = 1;
      no_bodys += nosubbodies;
    */

    cout << "NO vertices=" << nv << endl;
    cout << "NO DOF=" << nd << endl;
    cout << "NO Bodys=" << no_bodys << endl;

    points.SetSize (nv);

    for (i = 0; i < nv; i++)
      ma.GetPoint (i, points(i));


    Vec<3> point;
    double minpoint[3];

    Array <int> pnums;
    Array <int> usedpnums(nd);
    usedpnums = 0;

    point2domain.SetSize (nd);

    for (k = 0; k < no_bodys; k++)
      {
	MBS_Body* body = new MBS_Body();
	mbssystem.Append(body);

	cout << "body = " << k << endl;

	body->SetDomainNr(k);

	//alle Punkte, die zum Körper gehören herausfinden
	for (i = 0; i < ne; i++)
	  {
	    if (ma.GetElIndex(i) == k)
	      {
		pnums.SetSize(0);
		fes.GetDofNrs (i, pnums);

		for (j = 0; j < pnums.Size(); j++)
		  if (!usedpnums[pnums[j]])
		    {
		      usedpnums[pnums[j]] = 1;
		      body->AddNode(pnums[j]);
		      point2domain[pnums[j]] = k;
		    }
	      }
	  }
	cout << "set points" << endl;
      
	//3 geeignete Punkte für ein Referenzkoordinatensystem finden:
      
	minpoint[0] = 1E100;
	minpoint[1] = 1E100;
	minpoint[2] = 1E100;
	int pi1 = 0;
	int pi2 = 0;
	int pi3 = 0;
      
	//find point with minimal X
	for (l = 0; l < body->NONodes(); l++)
	  {
	    i = body->GetNode(l);
	    if (i >= nv) continue;

	    ma.GetPoint (i, point);
	    if (point[0] < minpoint[0] || (point[0] == minpoint[0] && 
					   (point[1] < minpoint[1] || point[2] < minpoint[2])))
	      {
		pi1 = i;
		minpoint[0] = point[0];
		minpoint[1] = point[1];
		minpoint[2] = point[2];
	      }
	  }

	cout << "p1 = (" << minpoint[0] << ", " << minpoint[1] << ", " << minpoint[2] << ")" << endl ;
      
	minpoint[0] = -1E100;
	minpoint[1] = 1E100;
	minpoint[2] = 1E100;
	//find point with maximal X
	for (l = 0; l < body->NONodes(); l++)
	  {
	    i = body->GetNode(l);
	    if (i >= nv) continue;

	    ma.GetPoint (i, point);
	    if (point[0] > minpoint[0] || (point[0] == minpoint[0] && 
					   (point[1] < minpoint[1] || point[2] < minpoint[2])))
	      {
		pi2 = i;
		minpoint[0] = point[0];
		minpoint[1] = point[1];
		minpoint[2] = point[2];
	      }
	  }

	cout << "p2 = (" << minpoint[0] << ", " << minpoint[1] << ", " << minpoint[2] << ")" << endl ;

      
	minpoint[0] = 1E100;
	minpoint[1] = 1E100;
	minpoint[2] = -1E100;
	//find point with maximal Z
	for (l = 0; l < body->NONodes(); l++)
	  {
	    i = body->GetNode(l);
	    if (i >= nv) continue;

	    ma.GetPoint (i, point);
	    if (point[2] > minpoint[2] || (point[2] == minpoint[2] && 
					   (point[0] < minpoint[0] || point[1] < minpoint[1])))
	      {
		pi3 = i;
		minpoint[0] = point[0];
		minpoint[1] = point[1];
		minpoint[2] = point[2];
	      }
	  }

	cout << "p3 = (" << minpoint[0] << ", " << minpoint[1] << ", " << minpoint[2] << ")" << endl ;

      
	body->SetRefPoint(0, pi1);  
	body->SetRefPoint(1, pi2);  
	body->SetRefPoint(2, pi3);
      
	cout << "Body Nr. " << k << endl;
      
	cout << "pi1=" << pi1 << ", pi2=" << pi2 << ", pi3=" << pi3 << endl;
	
      }


    nbnds = -1;
    for (i = 0; i < ma.GetNSE(); i++)
      if (ma.GetSElIndex(i) > nbnds)
	nbnds = ma.GetSElIndex(i);
    nbnds++;

    bnd2domain.SetSize (nbnds);
    bnd2domain = 0;

    for (i = 0; i < nse; i++)
      {
	int index = ma.GetSElIndex(i);
	ma.GetSElPNums (i, pnums);

	bnd2domain[index] = point2domain[pnums[0]];
      }

    for (i = 0; i < nbnds; i++)
      cout << "bnd2domain(" << i << ") = " << bnd2domain[i] << endl;



    // linear functions on individual bodies:  

    const BaseMatrix & matm = bfm->GetMatrix();

    pol1.SetSize (DIM*(DIM+1) * mbssystem.Size());
    mpol1.SetSize (pol1.Size());
    masslin.SetSize (pol1.Size(), pol1.Size());
    invmasslin.SetSize (pol1.Size(), pol1.Size());
	  
    for (i = 0; i < pol1.Size(); i++)
      {
	pol1[i] = new VVector<Vec<DIM> > (nd); 
	mpol1[i] = new VVector<Vec<DIM> > (nd); 
	(*pol1[i]) = 0;
      }

    //double point[3];
    if (DIM == 3)
      {
	for (k = 0; k < mbssystem.Size(); k++)
	  for (i = 0; i < nv; i++)
	    if (point2domain[i] == k)
	      {
		ma.GetPoint (i, point);
		for (j = 0; j < 3; j++)
		  {
		    (*pol1[12*k+j])   (i)(j) = 1;
		    (*pol1[12*k+3+j]) (i)(j) = point[0];
		    (*pol1[12*k+6+j]) (i)(j) = point[1];
		    (*pol1[12*k+9+j]) (i)(j) = point[2];
		  }
	      }
      }
    else
      {
	for (k = 0; k < mbssystem.Size(); k++)
	  for (i = 0; i < nv; i++)
	    if (point2domain[i] == k)
	      {
		ma.GetPoint (i, point);
		for (j = 0; j < 2; j++)
		  {
		    (*pol1[6*k-6+j]) (i)(j) = 1;
		    (*pol1[6*k-4+j]) (i)(j) = point[0];
		    (*pol1[6*k-2+j]) (i)(j) = point[1];
		  }
	      }
      }

    for (i = 0; i < pol1.Size(); i++)
      *mpol1[i] = matm * (*pol1[i]);

    for (i = 0; i < pol1.Size(); i++)
      for (j = 0; j < pol1.Size(); j++)
	masslin(i,j) = InnerProduct (*pol1[i], *mpol1[j]);
    (*testout) << "masslin = " << endl << masslin << endl;

    CalcInverse (masslin, invmasslin);

    (*testout) << "invmasslin = " << endl << invmasslin << endl;

    cout << "prepare done" << endl;
  }
 

  void NumProcDynamic :: SetLargeDeformation (VVector<Vec<DIM> > & u, 
					      const MBS_System & mbssystem)
  {
    int b, i, j;
    int nv = ma.GetNV();

    if (DIM == 3)
      {
	u = 0;
	for (b = 0; b < mbssystem.Size(); b++)
	  {
	  
	    MBS_Body * body = mbssystem[b];
	  
	    Vec<3> point;
	    Mat<3,3> rot;
	    body->GetRot(rot);
	  
	    Vec<3> p_off;
	    pde.GetMeshAccess().GetPoint(body->GetRefPoint(0), p_off);
	  
	    //  (*testout) << "rotmat = " << rot << endl;
	    for (j = 0; j < body->NONodes(); j++)
	      {
		i = body->GetNode(j);
		if (i >= nv) continue;
		
		point = points(i) - p_off;
		u(i) = rot * point - point + body->GetTrans();
	      }
	  }
      }
    else
      {
	/*
	  SystemVector<SysVector2d> & u =
	  dynamic_cast<SystemVector<SysVector2d> &> (hu);

	  for (b = 1; b <= mbssystem.Size(); b++)
	  {
	  
	  MBS_Body* body = mbssystem.Get(b);
	  
	  SysVector3d point, hui;
	  SysMatrix3d rot;
	  body->GetRot(rot);
	  
	  SysVector3d p_off;
	  pde.GetMeshAccess().GetPoint(body->GetRefPoint(1), &p_off.Elem(1));
	  
	  //  (*testout) << "rotmat = " << rot << endl;
	  for (j = 1; j <= body->NONodes(); j++)
	  {
	  i = body->GetNode(j);
	      
	  //	  pde.GetMeshAccess().GetPoint (i, &point.Elem(1));
	  point = points.Get(i);
	  point -= p_off;
	  rot.Mult (point, hui); 
	  hui -= point;	  
	  hui += body->GetTrans();
	      
	  u.VElem(i, 1) = hui.Get(1);
	  u.VElem(i, 2) = hui.Get(2);
	  }
	  }
	*/
      }    
  }

  void NumProcDynamic :: ApplyRotation (const VVector<Vec<DIM> > & x,
					VVector<Vec<DIM> > & y, 
					const MBS_System & mbssystem,
					int inverse)
  {
    int b, i, j;
  
    for (b = 0; b < mbssystem.Size(); b++)
      {
	MBS_Body* body = mbssystem[b];
      
	Mat<DIM> rot;
	body->GetRot(rot);
      
	if (!inverse)
	  {
	    for (j = 0; j < body->NONodes(); j++)
	      {
		i = body->GetNode(j);
		y(i) = rot * x(i);
	      }
	  }
	else
	  {
	    Mat<DIM> irot;
	    CalcInverse (rot, irot);
	  
	    for (j = 0; j < body->NONodes(); j++)
	      {
		i = body->GetNode(j);
		y(i) = irot * x(i);
	      }
	  }
      }

    /*
      else
      {
      const SystemVector<SysVector2d> & x =
      dynamic_cast<const SystemVector<SysVector2d> &> (hx);
      SystemVector<SysVector2d> & y =
      dynamic_cast<SystemVector<SysVector2d> &> (hy);      
      
      for (b = 1; b <= mbssystem.Size(); b++)
      {
	  
      MBS_Body* body = mbssystem.Get(b);
	  
      SysMatrix3d rot;
      body->GetRot(rot);

      if (!inverse)
      {
      for (j = 1; j <= body->NONodes(); j++)
      {
      i = body->GetNode(j);
      SysVector3d xi(x.Get(i,1), x.Get(i,2), 0);
      SysVector3d yi;
      rot.Mult (xi, yi);
      y.Elem(i,1) = yi.Get(1);
      y.Elem(i,2) = yi.Get(2);
      }
      }
      else
      {
      for (j = 1; j <= body->NONodes(); j++)
      {
      i = body->GetNode(j);
      SysVector3d xi(x.Get(i,1), x.Get(i,2), 0);
      SysVector3d yi;
      rot.Solve (xi, yi);
      y.Elem(i,1) = yi.Get(1);
      y.Elem(i,2) = yi.Get(2);
      }
      }
      }
      }
    */
  }


  void NumProcDynamic :: 
  ComputeFCorrection (MBS_System & mbssystem,
		      const VVector<Vec<DIM> > & ularge,
		      const VVector<Vec<DIM> > & kusmall,
		      const VVector<Vec<DIM> > & rkusmall,
		      const VVector<Vec<DIM> > & usmallrot,
		      VVector<Vec<DIM> > & f2)
  {
    const MeshAccess & ma = pde.GetMeshAccess();

    int i, j, k, ii;

    const FESpace & fes = *pde.GetFESpace ("v");
    int nd = fes.GetNDof();
    //    int np = ma.GetNP();
    int nbody = mbssystem.Size();

    if (DIM == 3)
      {
	static Matrix<> dvdr(nbody, DIM*DIM);
	static VVector<Vec<DIM> > hv(nd);
	static VVector<Vec<DIM> > u0l(nd);
	static VVector<Vec<DIM> > u0r(nd);
	hv.SetSize (nd);
	u0l.SetSize (nd);
	u0r.SetSize (nd);
      
	Mat<DIM,DIM> drot;
	drot = 0;
      
	for (j = 0; j < nbody; j++)
	  mbssystem[j] -> SetRot (drot);
    
	for (j = 0; j < nbody; j++)
	  {
	    for (i = 0; i < 9; i++)
	      {
		drot = 0;
		drot(i) = 1;
		mbssystem[j] -> SetRot (drot);
		ApplyRotation (kusmall, hv, mbssystem);
		dvdr(j, i) = InnerProduct (hv, usmallrot);
	      }
	    drot = 0;
	    mbssystem[j] -> SetRot (drot);
	  }
      
      
	Mat<DIM,DIM> rotr, rotl;
	f2 = 0;
	hv = ularge;
	double eps = 1e-4;
      
	for (i = 0; i < nd; i++)
	  {
	    int pointofbody = -1;
	  
	    for (j = 0; j < nbody; j++)
	      for (ii = 0; ii < 3; ii++)
		if (mbssystem[j]->GetRefPoint(ii) == i)
		  pointofbody = j;

	    if (pointofbody == -1) continue;
	  
	  
	    for (ii = 0; ii < 3; ii++)
	      {
		hv(i)(ii) = ularge(i)(ii) + eps;
		CalcRMatrix (hv, mbssystem);
		SetLargeDeformation (u0r, mbssystem);
		mbssystem[pointofbody] -> GetRot (rotr);
	      
		hv(i)(ii) = ularge(i)(ii) - eps;
		CalcRMatrix (hv, mbssystem);
		SetLargeDeformation (u0l, mbssystem);
		mbssystem[pointofbody] -> GetRot (rotl);
	      
		hv(i)(ii) = ularge(i)(ii);
	      
		rotr -= rotl;
		//	  (*testout) << "diff(rot) = " << rotr << endl;
		double sum = 0;
		for (j = 0; j < 9; j++)
		  sum += dvdr(pointofbody, j) * rotr(j) / (2*eps);
	      
		/*
		  u0r.Add (-1, u0l);
		  sum -= (u0r * rkusmall) / (2*eps);
		*/
		f2 (i)(ii) = sum;
	      }
	  }
      
	CalcRMatrix (ularge, mbssystem);
	//  (*testout) << "f2 = " << f2 << endl;
      }
    else
      {
#ifdef NONE
	const SystemVector<SysVector2d> & ularge =
	  dynamic_cast<const SystemVector<SysVector2d> &> (hularge);
	const SystemVector<SysVector2d> & kusmall =
	  dynamic_cast<const SystemVector<SysVector2d> &> (hkusmall);
	const SystemVector<SysVector2d> & rkusmall =
	  dynamic_cast<const SystemVector<SysVector2d> &> (hrkusmall);
	const SystemVector<SysVector2d> & usmallrot =
	  dynamic_cast<const SystemVector<SysVector2d> &> (husmallrot);
	SystemVector<SysVector2d> & f2 =
	  dynamic_cast<SystemVector<SysVector2d> &> (hf2);

	static DenseMatrix dvdr(nbody, 9);
	static SystemVector<SysVector2d> hv(np);
	static SystemVector<SysVector2d> u0l(np);
	static SystemVector<SysVector2d> u0r(np);
	hv.SetSize (np);
      
	SysMatrix3d drot;
	drot.SetScalar (0);
      
	for (j = 1; j <= nbody; j++)
	  mbssystem.Elem(j) -> SetRot (drot);
    
	for (j = 1; j <= nbody; j++)
	  {
	    for (i = 1; i <= 9; i++)
	      {
		drot.SetScalar (0);
		drot.Elem(i) = 1;
		mbssystem.Elem(j) -> SetRot (drot);
		ApplyRotation (kusmall, hv, mbssystem);
		dvdr.Elem(j, i) = hv * usmallrot;
	      }
	    drot.SetScalar (0);
	    mbssystem.Elem(j) -> SetRot (drot);
	  }
      
      
	SysMatrix3d rotr, rotl;
	f2.SetScalar (0);
	hv.Set (1, ularge);
	double eps = 1e-4;
      
	for (i = 1; i <= np; i++)
	  {
	    int pointofbody = 0;
	  
	    for (j = 1; j <= nbody; j++)
	      for (ii = 1; ii <= 2; ii++)
		if (mbssystem.Elem(j)->GetRefPoint(ii) == i)
		  pointofbody = j;
	    if (!pointofbody) continue;
	  
	  
	    for (ii = 1; ii <= 2; ii++)
	      {
		hv.Elem(i, ii) = ularge.Get(i, ii) + eps;
		CalcRMatrix (hv, mbssystem);
		SetLargeDeformation (u0r, mbssystem);
		mbssystem.Elem(pointofbody) -> GetRot (rotr);
	      
		hv.Elem(i, ii) = ularge.Get(i, ii) - eps;
		CalcRMatrix (hv, mbssystem);
		SetLargeDeformation (u0l, mbssystem);
		mbssystem.Elem(pointofbody) -> GetRot (rotl);
	      
		hv.Elem(i, ii) = ularge.Get(i, ii);
	      
		rotr.Add (-1, rotl);
		//	  (*testout) << "diff(rot) = " << rotr << endl;
		double sum = 0;
		for (j = 1; j <= 9; j++)
		  sum += dvdr.Get(pointofbody, j) * rotr.Get(j) / (2*eps);
	      
		/*
		  u0r.Add (-1, u0l);
		  sum -= (u0r * rkusmall) / (2*eps);
		*/
		f2.VElem (i, ii) = sum;
	      }
	  }
      
	CalcRMatrix (ularge, mbssystem);
	//  (*testout) << "f2 = " << f2 << endl;

	//      f2.SetScalar (0);
#endif
      }
  }




  namespace IMBS {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetNumProcs().AddNumProc ("mbs", NumProcDynamic::Create);
    }

    Init init;
  }
}
