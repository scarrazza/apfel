************************************************************************
*
*     TabulationFrag.f:
*
*     Example program used for fragmentation function evolution.
*     The fragmentation functions are hardcoded and for the following
*     options are available:
*                          
*     DSS                  HKNS:
*                          
*     - "dss_pip_lo"       - "hknsff07_pip_lo"
*     - "dss_pip_nlo"      - "hknsff07_pip_nlo"
*     - "dss_pim_lo"       - "hknsff07_pim_lo"
*     - "dss_pim_nlo"      - "hknsff07_pim_nlo"
*     - "dss_Kp_lo"        - "hknsff07_Kp_lo"
*     - "dss_Kp_nlo"       - "hknsff07_Kp_nlo"
*     - "dss_Km_lo"        - "hknsff07_Km_lo"
*     - "dss_Km_nlo"       - "hknsff07_Km_nlo"
*     - "dss_p_lo"         - "hknsff07_p_lo"
*     - "dss_p_nlo"        - "hknsff07_p_nlo"
*     - "dss_pb_lo"        - "hknsff07_pb_lo"
*     - "dss_pb_nlo"       - "hknsff07_pb_nlo"
*     - "dss_h_lo"    
*     - "dss_h_nlo"   
*     - "dss_hb_lo"   
*     - "dss_hb_nlo"  
*     
************************************************************************
      program TabulationFrag
*
      implicit none
*
      integer ipt,iq,ilha
      double precision Q0,Q
      double precision Q02,Q2(0:4)
      double precision AlphaQCD
      double precision xPDF
      double precision xlha(11),xSigA(11),xSigB(11),xgA(11),xgB(11)
      double precision xT3A(11),xT3B(11),xV3A(11),xV3B(11)
      double precision srel,grel,t3rel,v3rel
      data xlha / 1d-2, 5d-2, 1d-1, 2d-1, 3d-1, 4d-1,
     1            5d-1, 6d-1, 7d-1, 8d-1, 9d-1 /
      data Q2 / 1d0, 1d1, 1d2, 1d3, 1d4 /
c      data Q2 / 1d0, 2.0448d0, 2.045d0, 18.489d0, 18.491d0 /
c      data Q2 / 18.5d0, 20d0, 1d2, 1d3, 1d4 /
*
*     Settings
*
*     HKNS
*
c      call SetPDFSet("hknsff07_pip_lo")
c      call SetPerturbativeOrder(0)
c      call SetMaxFlavourPDFs(6)
c      call SetMaxFlavourAlpha(6)
c      call SetPoleMasses(1.43d0,4.3d0,175d0)
c      call SetAlphaEvolution("lambda")
c      call SetLambdaQCDRef(0.220d0,4)
cc      call SetLambdaQCDRef(0.323d0,4)
c      call SetNumberOfGrids(2)
c      call SetGridParameters(1,50,3,1d-2)
c      call SetGridParameters(2,40,3,7d-1)
*
*     DSS
*
c      call SetPDFSet("dss_p_lo")
c      call SetPerturbativeOrder(0)
c      call SetMaxFlavourPDFs(5)
c      call SetMaxFlavourAlpha(5)
c      call SetPoleMasses(1.43d0,4.3d0,175d0)
c      call SetAlphaEvolution("lambda")
c      call SetLambdaQCDRef(0.220d0,4)
cc      call SetLambdaQCDRef(0.334d0,4)
c      call SetNumberOfGrids(2)
c      call SetGridParameters(1,50,3,5d-2)
c      call SetGridParameters(2,40,3,7d-1)
*
*     Evolve PDFs on the grids
*
      Q02 = 1d0
      Q0  = dsqrt(Q02)
*
      do ipt=0,1
*
*     DSS
*
         write(6,*) "======================= DSS ======================"
         call SetTimeLikeEvolution(.true.)
         if(ipt.eq.0)then
            call SetPDFSet("dss_pip_lo")
            call SetLambdaQCDRef(0.220d0,4)
         elseif(ipt.eq.1)then
            call SetPDFSet("dss_pip_nlo")
            call SetLambdaQCDRef(0.334d0,4)
         endif
         call SetPerturbativeOrder(ipt)
         call SetPoleMasses(1.43d0,4.3d0,175d0)
         call SetAlphaEvolution("lambda")
         call SetNumberOfGrids(2)
         call SetGridParameters(1,150,3,5d-2)
         call SetGridParameters(2,100,3,7d-1)
         call InitializeAPFEL
*
         do iq=1,4
            Q = dsqrt(Q2(iq))
*
            call EvolveAPFEL(Q0,Q)
*
            write(6,"(a,i1,a)") " Perturbative order = N",ipt,"LO"
            write(6,"(a,f9.3,a)") " Final energy = ",Q2(iq)," GeV^2"
            write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
            write(6,*) "  "
*
            do ilha=2,11
                    xSigA(ilha) = 
     2              xPDF(1,xlha(ilha)) + xPDF(-1,xlha(ilha)) +
     3              xPDF(2,xlha(ilha)) + xPDF(-2,xlha(ilha)) +
     4              xPDF(3,xlha(ilha)) + xPDF(-3,xlha(ilha)) +
     5              xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)) +
     6              xPDF(5,xlha(ilha)) + xPDF(-5,xlha(ilha)) +
     7              xPDF(6,xlha(ilha)) + xPDF(-6,xlha(ilha))
                    xT3A(ilha) = xPDF(2,xlha(ilha))
     1                         + xPDF(-2,xlha(ilha))
     2                         - xPDF(1,xlha(ilha))
     3                         - xPDF(-1,xlha(ilha))
                    xV3A(ilha) = xPDF(2,xlha(ilha))
     1                         - xPDF(-2,xlha(ilha))
     2                         - xPDF(1,xlha(ilha))
     3                         + xPDF(-1,xlha(ilha))
                    xgA(ilha)  = xPDF(0,xlha(ilha))
            enddo
            write(*,*) "  "
*
            call EvolveAPFEL(Q,Q)
*
c            write(6,"(a5,10a12)") "x","SigmaA","SigmaB","dSig",
c     1                               "gluonA","gluonB","dglu"
            write(6,"(a5,10a12)") "x","T3A","T3B","dT3",
     1                                "V3A","V3B","dV3"
            do ilha=2,11
                    xSigB(ilha) =
     2              xPDF(1,xlha(ilha)) + xPDF(-1,xlha(ilha)) +
     3              xPDF(2,xlha(ilha)) + xPDF(-2,xlha(ilha)) +
     4              xPDF(3,xlha(ilha)) + xPDF(-3,xlha(ilha)) +
     5              xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)) +
     6              xPDF(5,xlha(ilha)) + xPDF(-5,xlha(ilha)) +
     7              xPDF(6,xlha(ilha)) + xPDF(-6,xlha(ilha))
                    xT3B(ilha) = xPDF(2,xlha(ilha))
     1                         + xPDF(-2,xlha(ilha))
     2                         - xPDF(1,xlha(ilha))
     3                         - xPDF(-1,xlha(ilha))
                    xV3B(ilha) = xPDF(2,xlha(ilha))
     1                         - xPDF(-2,xlha(ilha))
     2                         - xPDF(1,xlha(ilha))
     3                         + xPDF(-1,xlha(ilha))
                    xgB(ilha)  = xPDF(0,xlha(ilha))
*
               srel  = 1d2 * ( xSigA(ilha) - xSigB(ilha) ) / xSigA(ilha)
               t3rel = 1d2 * ( xT3A(ilha) - xT3B(ilha) ) / xT3A(ilha)
               v3rel = 1d2 * ( xV3A(ilha) - xV3B(ilha) ) / xV3A(ilha)
               grel  = 1d2 * ( xgA(ilha) - xgB(ilha) ) / xgA(ilha)

               write(6,"(es7.1,6es12.4)") xlha(ilha),
c     1             xSigA(ilha),xSigB(ilha),srel,xgA(ilha),xgB(ilha),grel
     1         xT3A(ilha),xT3B(ilha),t3rel,xV3A(ilha),xV3B(ilha),v3rel
               write(30+ipt,"(es7.1,6es12.4)") xlha(ilha),
c     1             xSigA(ilha),xSigB(ilha),srel,xgA(ilha),xgB(ilha),grel
     1         xT3A(ilha),xT3B(ilha),t3rel,xV3A(ilha),xV3B(ilha),v3rel
            enddo
            write(30+ipt,*) "  "
            write(30+ipt,*) "  "
            write(*,*) "  "
         enddo
         call CleanUp
c$$$*
c$$$*     HKNS
c$$$*
c$$$         write(6,*) "====================== HKNS ======================"
c$$$         call SetTimeLikeEvolution(.true.)
c$$$         if(ipt.eq.0)then
c$$$            call SetPDFSet("hknsff07_pip_nlo")
c$$$            call SetLambdaQCDRef(0.220d0,4)
c$$$         elseif(ipt.eq.1)then
c$$$            call SetPDFSet("hknsff07_pip_nlo")
c$$$            call SetLambdaQCDRef(0.323d0,4)
c$$$         endif
c$$$         call SetPerturbativeOrder(ipt)
c$$$         call SetPoleMasses(1.43d0,4.3d0,175d0)
c$$$         call SetAlphaEvolution("lambda")
c$$$         call SetNumberOfGrids(2)
c$$$         call SetGridParameters(1,50,3,5d-2)
c$$$         call SetGridParameters(2,40,3,7d-1)
c$$$         call InitializeAPFEL
c$$$*
c$$$         do iq=0,4
c$$$            Q = dsqrt(Q2(iq))
c$$$*
c$$$            call EvolveAPFEL(Q0,Q)
c$$$*
c$$$            write(6,"(a,i1,a)") " Perturbative order = N",ipt,"LO"
c$$$            write(6,"(a,f8.3,a)") " Final energy = ",Q2(iq)," GeV^2"
c$$$            write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
c$$$            write(6,*) "  "
c$$$*
c$$$            do ilha=2,11
c$$$                    xSigA(ilha) = 
c$$$     2              xPDF(1,xlha(ilha)) + xPDF(-1,xlha(ilha)) +
c$$$     3              xPDF(2,xlha(ilha)) + xPDF(-2,xlha(ilha)) +
c$$$     4              xPDF(3,xlha(ilha)) + xPDF(-3,xlha(ilha)) +
c$$$     5              xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)) +
c$$$     6              xPDF(5,xlha(ilha)) + xPDF(-5,xlha(ilha)) +
c$$$     7              xPDF(6,xlha(ilha)) + xPDF(-6,xlha(ilha))
c$$$                    xT3A(ilha) = xPDF(2,xlha(ilha))
c$$$     1                         + xPDF(-2,xlha(ilha))
c$$$     2                         - xPDF(1,xlha(ilha))
c$$$     3                         - xPDF(-1,xlha(ilha))
c$$$                    xV3A(ilha) = xPDF(2,xlha(ilha))
c$$$     1                         - xPDF(-2,xlha(ilha))
c$$$     2                         - xPDF(1,xlha(ilha))
c$$$     3                         + xPDF(-1,xlha(ilha))
c$$$                    xgA(ilha)  = xPDF(0,xlha(ilha))
c$$$            enddo
c$$$            write(*,*) "  "
c$$$*
c$$$            call EvolveAPFEL(Q,Q)
c$$$*
c$$$c            write(6,"(a5,6a12)") "x","SigmaA","SigmaB","dSig",
c$$$c     1                               "gluonA","gluonB","dglu"
c$$$            write(6,"(a5,10a12)") "x","T3A","T3B","dT3",
c$$$     1                                "V3A","V3B","dV3"
c$$$            do ilha=2,11
c$$$                    xSigB(ilha) =
c$$$     2              xPDF(1,xlha(ilha)) + xPDF(-1,xlha(ilha)) +
c$$$     3              xPDF(2,xlha(ilha)) + xPDF(-2,xlha(ilha)) +
c$$$     4              xPDF(3,xlha(ilha)) + xPDF(-3,xlha(ilha)) +
c$$$     5              xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)) +
c$$$     6              xPDF(5,xlha(ilha)) + xPDF(-5,xlha(ilha)) +
c$$$     7              xPDF(6,xlha(ilha)) + xPDF(-6,xlha(ilha))
c$$$                    xT3B(ilha) = xPDF(2,xlha(ilha))
c$$$     1                         + xPDF(-2,xlha(ilha))
c$$$     2                         - xPDF(1,xlha(ilha))
c$$$     3                         - xPDF(-1,xlha(ilha))
c$$$                    xV3B(ilha) = xPDF(2,xlha(ilha))
c$$$     1                         - xPDF(-2,xlha(ilha))
c$$$     2                         - xPDF(1,xlha(ilha))
c$$$     3                         + xPDF(-1,xlha(ilha))
c$$$                    xgB(ilha)  = xPDF(0,xlha(ilha))
c$$$*
c$$$               srel  = 1d2 * ( xSigA(ilha) - xSigB(ilha) ) / xSigA(ilha)
c$$$               t3rel = 1d2 * ( xT3A(ilha) - xT3B(ilha) ) / xT3A(ilha)
c$$$               v3rel = 1d2 * ( xV3A(ilha) - xV3B(ilha) ) / xV3A(ilha)
c$$$               grel  = 1d2 * ( xgA(ilha) - xgB(ilha) ) / xgA(ilha)
c$$$               write(6,"(es7.1,6es12.4)") xlha(ilha),
c$$$c     1             xSigA(ilha),xSigB(ilha),srel,xgA(ilha),xgB(ilha),grel
c$$$     1         xT3A(ilha),xT3B(ilha),t3rel,xV3A(ilha),xV3B(ilha),v3rel
c$$$               write(33+ipt,"(es7.1,6es12.4)") xlha(ilha),
c$$$c     1             xSigA(ilha),xSigB(ilha),srel,xgA(ilha),xgB(ilha),grel
c$$$     1         xT3A(ilha),xT3B(ilha),t3rel,xV3A(ilha),xV3B(ilha),v3rel
c$$$            enddo
c$$$            write(33+ipt,*) "  "
c$$$            write(33+ipt,*) "  "
c$$$            write(*,*) "  "
c$$$         enddo
c$$$         call CleanUp
      enddo
*
      end
