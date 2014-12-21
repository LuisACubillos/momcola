TOP_OF_MAIN_SECTION
  arrmblsize=500000; // 
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
  gradient_structure::set_MAX_NVAR_OFFSET(1000);

DATA_SECTION
	int iseed
	!!long int lseed=iseed;
	!!CLASS random_number_generator rng(iseed);

  init_int nyear  
  init_int nages
  init_int nlength
  init_matrix year_matrix(1,nyear,1,13)
  //!!cout<<"data:"<<year_matrix<<endl;exit(1);
  
  init_vector ages(1,nages)
  init_vector length2(1,nlength)
  init_matrix Cps(1,nyear,1,nages)
  init_matrix Cpsl(1,nyear,1,nlength)
  init_matrix Cst(1,nyear,1,nages)
  init_matrix Cct(1,nyear,1,nages)
  init_matrix Chb(1,nyear,1,nages)
  init_matrix Wmed(1,nyear,1,nages)
  //!!cout<<"Wmed:"<<Wmed<<endl;exit(1);

  init_matrix Win(1,nyear,1,nages)
  init_matrix Wps(1,nyear,1,nages)

  //!!cout<<"Wps:"<<Wps<<endl;exit(1);
 
  !! ad_comm::change_datafile_name("esmeco_op.ctl");
	init_vector matur(1,nages);
	init_number Bsa;
	init_vector Csal(1,nlength);
	init_vector nmus(1,6);
	init_vector dt(1,3);
	init_vector growth(1,2);
	init_vector indSt(1,2);
	init_vector indCt(1,2);
	init_vector indSrv(1,2);
	init_vector indq(1,2);
	init_number Nyeq_ini; //  numbers of initial years on equilibrium condition 
	init_int    phase1;
	init_int    phase2;
	init_int    phase3;
	init_int    phase4;
	init_int    phase5;
	init_int    phase6;
	init_int    phase7;
	init_int    phase8;
	init_int    phase9;
	init_int    phase10;
	init_int    phase11;
	init_int    phase12;
	init_vector Mandh_pars(1,2);
	init_vector cv_penal(1,6);
	//!!cout<<"Penlaties:"<<cv_penal<<endl;exit(1);

	!! ad_comm::change_datafile_name("proyectaFs.ctl");
    int yr_sim
    init_int nFt
    init_vector mf(1,nFt)
	!! yr_sim=1;


INITIALIZATION_SECTION

// Some initial values

//  log_Bo 8.45

  log_Bo 8.45
  log_qacus 0
  log_A50ps 1.4
  log_A50st 1.6
  log_Dst 0
  log_A50ct 1.6
  log_Dct 0
  Lo        29.8 // 
  cv        0.09
  log_A50sa 0
  log_Ssa 0
  log_qsa -0.3567 // by default q=0.7
  log_M -1.0498
  log_h -0.2877


PARAMETER_SECTION

// Selectivity: acoustic and swept-area
 init_bounded_vector log_A50hb(1,3,0,2.5,phase1)// acoustic logistic 
 init_bounded_vector log_Dhb(1,3,-2.3,1.3,phase1)

 init_bounded_number log_A50sa(-2.3,0,phase1) // Swept-area  dome-shape
 init_bounded_vector log_Ssa(1,2,-2.3,0.705,phase1)// 

// Selectivity: purse-seine
 init_bounded_number log_A50ps(0,1.38,phase1)  // dome-shape
 init_bounded_vector log_Sps(1,2,-4.5,1.7,phase1)// 

// Selectivity: trawl
 init_bounded_vector log_A50st(1,3,1,2.5,phase1)// 
 init_bounded_vector log_Dst(1,4,-2.3,10.3,phase1)//
 init_bounded_vector log_A50ct(1,3,1,2.5,phase1)//
 init_bounded_vector log_Dct(1,3,-2.3,1.3,phase1)//

// Recruits 
 init_bounded_number log_Bo(4,12,phase2)
 init_bounded_vector log_dev_Rt(1,nyear-1,-20,20,phase3)
 init_bounded_vector log_dev_No(1,nages,-10,10,phase4)
 init_bounded_number log_avgR(2,20,phase2)

// catchabilities
 init_vector log_qst(1,3,phase6) // trawl
 init_number log_qps(phase6) // purse-seine
 init_number log_qacus(phase7) // acoustic
 init_number log_qsa(phase8) // swept-area


// growth parameters
 init_bounded_number Lo(1,100,phase9) // Mean length for first age-group
 init_bounded_number cv(0.01,0.5,phase10) // cv


// M and h estimation
 init_number log_M(phase11) // 
 init_number log_h(phase12) // 

 //mu40 en fase 3
 //init_bounded_number mu40(0.01,0.5,3)
 //number phi
 //number spr40
 //number sprpen
 //matrix npr(1,2,1,nages)
 //-------------------------------------------------
 // Array definition

  vector years(1,nyear)
  vector Acous(1,nyear)
  vector CPUEst(1,nyear)
  vector CPUEps(1,nyear)
  vector Yps(1,nyear)
  vector Yst(1,nyear)
  vector Yct(1,nyear)
  vector Yps_pred(1,nyear)
  vector Yst_pred(1,nyear)
  vector Yct_pred(1,nyear)
  vector min_N1(1,nyear)
  vector min_N2(1,nyear)
  vector min_N3(1,nyear)

  vector sigma_age(1,nages)
  vector mu_age(1,nages)
  vector N1(1,nages)
  vector N2(1,nages)
  vector N3(1,nages)
  vector N4(1,nages)

  vector Neq(1,nages)

  vector Ones_age(1,nages)
  vector Ones_length(1,nlength)
  vector Ones_yr(1,nyear)
  vector ppred_sa(1,nlength)
  vector pobs_sa(1,nlength)

  sdreport_vector SSB(1,nyear)
  sdreport_vector Biomass(1,nyear)
  vector Bio_hb(1,nyear)
  vector Bio_st(1,nyear)
  vector Bio_ps(1,nyear)
  vector Bio_ct(1,nyear)
  sdreport_vector Recruits(1,nyear)
  vector Rpred(1,nyear)

  vector pred_Acous(1,nyear)
  vector pred_CPUEps(1,nyear)
  vector pred_CPUEst(1,nyear)
  vector likeval(1,15);
  vector penalty(1,10)
  vector S_sa(1,nages)
  vector N_sa(1,nages)
  vector Bio_trw(1,nyear)

 
  number So
  number suma1
  number suma2
  number suma3
  number suma4
  number suma5
  number alfa
  number beta
  number Bio_sa
  sdreport_number log_Ro
  number M;
  number h;
  number BPRo
  number nn
  number log_Rmed
  

  matrix S_ps(1,nyear,1,nages)
  matrix S_st(1,nyear,1,nages)
  matrix S_ct(1,nyear,1,nages)
  matrix S_hb(1,nyear,1,nages)

  matrix N(1,nyear,1,nages)
  matrix N_st(1,nyear,1,nages)
  matrix N_hb(1,nyear,1,nages)
  matrix N_ps(1,nyear,1,nages)
  matrix N_ssb(1,nyear,1,nages)
  matrix mu_ct(1,nyear,1,nages)
  matrix mu_st(1,nyear,1,nages)
  matrix mu_ps(1,nyear,1,nages)
  matrix Nj(1,nyear,1,nages)
  matrix S_tot(1,nyear,1,nages)
  matrix mu_tot(1,nyear,1,nages)
  vector mu_lastyr(1,4)

  matrix Cps_pred(1,nyear,1,nages)
  matrix Cpsl_pred(1,nyear,1,nlength)
  matrix Cst_pred(1,nyear,1,nages)
  matrix Cct_pred(1,nyear,1,nages)
  matrix pobs_ps(1,nyear,1,nages)
  matrix pobs_st(1,nyear,1,nages)
  matrix pobs_ct(1,nyear,1,nages)
  matrix ppred_ps(1,nyear,1,nages)
  matrix ppred_st(1,nyear,1,nages)
  matrix ppred_ct(1,nyear,1,nages)
  matrix pobs_hb(1,nyear,1,nages)
  matrix ppred_hb(1,nyear,1,nages)
  matrix pobs_psl(1,nyear,1,nlength)
  matrix ppred_psl(1,nyear,1,nlength)
  matrix P1(1,nages,1,nlength)
  matrix P2(1,nages,1,nlength)
  matrix P3(1,nages,1,nlength)
  matrix p_age2length(1,nages,1,nlength)

  vector cv1(1,nyear); // 
  vector cv2(1,nyear); // 
  vector cv3(1,nyear); // 
  vector cv4(1,nyear); //
  vector cv5(1,nyear); // 
  vector cv6(1,nyear); //
  vector sigma1(1,nyear);
  vector sigma2(1,nyear);
  vector sigma3(1,nyear);
  
  //Para la proyeccion de un yr calculo de CTP
  number Rp //Reclutamiento promedio
  number Nplus //Abundancia 1 yr despues
  number Yp    
  number mean_BDop

  vector Np(1,nages)
  vector Nvp(1,nages)  //N sin explotacion
  vector wp(1,nages)  //peso del ultimo yr
  vector Selp(1,nages) //Selectividad del ultimo yr
  vector mu_ref(1,nages)
  vector NDp(1,nages)  //numero desovantes
  vector Ctp(1,nages)  //captura en numero
  vector NDvp(1,nages); //Numero desovante en ausencia de pesca

  matrix Yproy(nyear+1,nyear+yr_sim,1,nFt);
  matrix BDp(nyear+1,nyear+yr_sim,1,nFt);
  matrix BDop(nyear+1,nyear+yr_sim,1,nFt);
  matrix RPDdin(nyear+1,nyear+yr_sim,1,nFt);
  matrix RPDmed(nyear+1,nyear+yr_sim,1,nFt);
  matrix RPDeqp(nyear+1,nyear+yr_sim,1,nFt);
  matrix RPDeqm(nyear+1,nyear+yr_sim,1,nFt);

  objective_function_value f

PRELIMINARY_CALCS_SECTION

  years=column(year_matrix,1);// years
  Acous=column(year_matrix,2);// acoustic
  CPUEst=column(year_matrix,4); // CPUE trawl
  CPUEps=column(year_matrix,6); // CPUE purse-seine
  Yps=column(year_matrix,8); // Landings purse-seine
  Yst=column(year_matrix,10);// Landings trawl
  Yct=column(year_matrix,12);// Landings trawl, central zone

  cv1=column(year_matrix,3);
  cv2=column(year_matrix,5);
  cv3=column(year_matrix,7);
  cv4=column(year_matrix,9);
  cv5=column(year_matrix,11);
  cv6=column(year_matrix,13);

  Ones_age=1;// ones vector (age size)
  Ones_yr=1;// ones vector (year size)
  Ones_length=1;// ones vector (length size)
  M=Mandh_pars(1);
  h=Mandh_pars(2);

PROCEDURE_SECTION

  f=0;
  Eval_selectivity();
  Eval_mortality();
  Eval_abundance();
  Eval_age2length();
  Eval_C_predicted();
  Eval_stock_recruit();
  Eval_indexes();
  Eval_likelihood();
  
  //if(active(mu40)){
	  //compute_spr40();
	  //}
  
  //if(mceval_phase())
  if(last_phase())
  {
	  CTP_cal();
  }
  /*
  if (mceval_phase())
    {
       MCWrite();   
    } 
 */


  /*
	FUNCTION compute_spr40

	//Proxy de mu_msy es mu40
	//dvariable phi;
	//dvar_matrix npr(1,2,1,nages);
    //dvariable spr40;
	
	//phi=BPRo*exp(-0.67*M);
	phi=0; //*exp(-0.67*M);
	spr40=0;	
	npr(1,1)=1;
	npr(2,1)=1;
	//spr40=npr(1)*exp(-0.67*M)*matur(1)*Wmed(nyear,1);
	
	for(int j=2;j<=nages;j++)
	{
		npr(1,j)=npr(1,j-1)*exp(-1.0*M);
		npr(2,j)=npr(2,j-1)*(1.-mu40*S_tot(nyear,j-1))*exp(-1.0*M);
		//spr40  +=npr(j)*exp(-0.67*M)*matur(j)*Wmed(nyear,j);
	}
	npr(1,nages)=npr(1,nages)+npr(1,nages)*(1-mu40*S_tot(nyear,nages))*exp(-1.0*M);
	//spr40+= npr(nages)*exp(-0.67*M)*matur(nages)*Wmed(nyear,nages);
	for(int j=1;j<=nages;j++)
	{
		phi    =phi+npr(1,j)*matur(j)*Wmed(nyear,j);
		spr40  =spr40+npr(2,j)*matur(j)*Wmed(nyear,j);		
	}

	sprpen = 300.*square(spr40/phi-0.40);

 */
  
FUNCTION CTP_cal

	ofstream ctp_out("ctpmcol.dat");
	//ofstream byr_out("byrmcol.dat");
	//ofstream ryr_out("ryrmcol.dat");
	//ofstream muyr_out("muyrmcol.dat");
	//ofstream syr_out("syrmcol.dat");
	ofstream ind_out("indices.dat");
	//ofstream ssb_out("ssbmcol.dat");
	//ofstream rpr1_out("rpr1mcol.dat");
  	for (int j=1; j<=nFt; j++) // ciclo de multiplicadores de F
  	{
  		Np = N(nyear); //Abundancia en el ultimo anio
  		Nvp = N(nyear); //Abundancia en el ultimo anio en ausencia de pesca
  		Rp = mean(Recruits(nyear-5,nyear)); // Supuesto: reclutamiento constante igual al promedio de los ultimos 5 años
  		wp = Wmed(nyear); //Pesos medio ultimo anio
  		Selp = S_tot(nyear); //Selectividad total en el ultimo anio con datos
  		mu_ref = 0.18*mf(j)*Selp; //Tasa de explotacion por edad Hay q multiplicar por mu de referencia (mu40)!!

  		for (int i=nyear+1; i<=nyear+yr_sim; i++)
  		{
  			Nplus = Np(nages)*(1-mu_ref(nages))*exp(-1.0*M); //a utilizar en grupo plus
  			Np(2,nages) = ++elem_prod(Np(1,nages-1),(1-mu_ref(1,nages-1)))*exp(-1.0*M);
  			Np(nages) += Np(nages) + Nplus;// Grupo plus
  			Np(1) = Rp;
  			NDp = elem_prod(elem_prod(Np*exp(-0.67*M),(1-mu_ref)),matur); //Abundancia desovante proyectada
  			BDp(i,j) = sum(elem_prod(NDp,wp));
  			Ctp = elem_prod(mu_ref,Np)*exp(-0.5*M); //captura numero
  			Yp = sum(elem_prod(Ctp,wp)); //Captura en peso
  			Yproy(i,j) = Yp;
  			Nplus = Nvp(nages)*exp(-1.0*M); //Nplus es re-utilizado en el grupo plus de abundancia en ausencia de pesca
  			Nvp(2,nages) = ++Nvp(1,nages-1)*exp(-1.0*M);
  			Nvp(nages) = Nvp(nages) + Nplus;// Grupo plus en ausencia de pesca
  			Nvp(1) = Rp;
  			NDvp = elem_prod(Nvp*exp(-0.67*M),matur); //Abundancia desovante proyectada en ausencia de pesca
  			BDop(i,j) = sum(elem_prod(NDvp,wp));
  		}
  	}

  	//biomasa desovante y razones (varias) de potencial desovante
  	mean_BDop = mean(BDop);
  	RPDdin = elem_div(BDp,BDop); //dinamica
  	RPDmed = BDp/mean_BDop; //media
  	//RPDeqp = BDp/So; //de BDo F40% 
  	//RPDeqm = BDp/; //de equilibrio modelo
    ctp_out << Yproy << endl; ctp_out.close();
	//Biomasa total
	//byr_out<< Biomass(nyear) <<endl; byr_out.close();
	//Reclutas
	//ryr_out<<Recruits(nyear)<<endl; ryr_out.close();
	//Tasa de explotacion
	//muyr_out<<(mu_ct(nyear)+mu_ps(nyear)+mu_st(nyear))<< endl; muyr_out.close();
	//Stock desovante
	//syr_out<< SSB(nyear)<<endl;syr_out.close();
	ind_out<< Biomass(nyear)  << " "
	       << SSB(nyear)      << " "
		   << Recruits(nyear) << " "
		   << mu_lastyr       << endl;
		   ind_out.close();

 //FUNCTION MCWrite

  //report3 << Biomass << endl;
  //report4 << SSB << endl;
  //report6 << Recruits << endl;

FUNCTION Eval_selectivity
  int i;

  for (i=1;i<=nyear;i++)

   {
// Purse_seine
    S_ps(i)=mfexp(-0.5/square(exp(log_Sps(1)))*square(-1.*(ages-exp(log_A50ps))));
      for (int j=1;j<=nages;j++){
      if(ages(j)>mfexp(log_A50ps)){
      S_ps(i,j)=mfexp(-0.5/square(exp(log_Sps(2)))*square(-1.*(ages(j)-exp(log_A50ps))));}
    }

// Trawling southern area

//   domme or asymptotic at the begining, depending on the penalties
 
     S_st(i)=mfexp(-0.5/square(exp(log_Dst(1)))*square(-1.*(ages-exp(log_A50st(1)))));
     for (int j=1;j<=nages;j++){
     if(ages(j)>mfexp(log_A50st(1))){
     S_st(i,j)=mfexp(-0.5/square(exp(log_Dst(2)))*square(-1.*(ages(j)-exp(log_A50st(1)))));}}


// by year-block (3)
    if (years(i)>=indSt(1))
    {S_st(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50st(2)))/exp(log_Dst(3))))));}
    if (years(i)>=indSt(2))
    {S_st(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50st(3)))/exp(log_Dst(4))))));}



// Trawling central area
    S_ct(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50ct(1)))/exp(log_Dct(1))))));
 // by year-block (3)
    if (years(i)>=indCt(1))
    {S_ct(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50ct(2)))/exp(log_Dct(2))))));}
    if (years(i)>=indCt(2))
    {S_ct(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50ct(3)))/exp(log_Dct(3))))));}
 

// Acoustic
     S_hb(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50hb(1)))/exp(log_Dhb(1))))));
// by year-block (3)
     if (years(i)>=indSrv(1))
     {S_hb(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50hb(2)))/exp(log_Dhb(2))))));}
     if (years(i)>=indSrv(2))
     {S_hb(i)=(elem_div(Ones_age,(1+exp(-1.0*log(19)*(ages-exp(log_A50hb(3)))/exp(log_Dhb(3))))));}
  }


  // Swept-area
      S_sa=mfexp(-0.5/square(exp(log_Ssa(1)))*square(-1.*(ages-exp(log_A50sa))));

      for (int j=1;j<=nages;j++){
      if(ages(j)>mfexp(log_A50sa)){
      S_sa(j)=mfexp(-0.5/square(exp(log_Ssa(2)))*square(-1.*(ages(j)-exp(log_A50sa))));}}


FUNCTION Eval_mortality

  if(phase11>0){M=mfexp(log_M);}
  if(phase12>0){h=mfexp(log_h);}


FUNCTION Eval_abundance
 int i, j;

  Neq(1)=1;
 // BPRo estimation

  for (j=2;j<=nages;j++)
   {Neq(j)=Neq(j-1)*mfexp(-1*M);} // steady-state condition per recruit
    Neq(nages)=Neq(nages)/(1-exp(-1*M)); // plus group
    BPRo=sum(elem_prod(Neq,Win(1)));  

 // Ro estimation
  log_Ro=log_Bo-log(BPRo);
  log_Rmed=log_avgR;

// Abundance at the first year by age
  Neq=Neq*mfexp(log_Ro);// steady-state as reference

  if (phase4<0){
  N(1)=Neq;}
  else{
  N(1)=mfexp(log(Neq)+log_dev_No-0.5*square(cv_penal(2)));}

 // Recruitments by year starting in t=1
  for (i=2;i<=nyear;i++){ 
     
    N(i,1)=mfexp(log_Rmed+log_dev_Rt(i-1)-0.5*square(cv_penal(3)));
    }


   N_sa=elem_prod(S_sa,N(1));

// Abundance by seson within year
   dvariable fpen=0.0;

   Biomass(1)=sum(elem_prod(N(1),Win(1)));
   N2=N(1)*exp(-dt(1)*M);// central

   Bio_ct(1)=sum(elem_prod(elem_prod(N2,S_ct(1)),Wmed(1))+1E-10);//central
   mu_ct(1)=Yct(1)/Bio_ct(1)*S_ct(1);
   Cct_pred(1)=elem_prod(N2,mu_ct(1));
   N3=posfun(elem_prod(N2,1-mu_ct(1)),0.5,fpen);// central
   f+=10000*fpen; 

   fpen=0.0;

   Nj(1)=N3*exp(-0.25*M);// It will be use for total CPUE
   N2=N3*exp(-dt(2)*M);// southern
   Bio_st(1)=sum(elem_prod(elem_prod(N2,S_st(1)),Wmed(1))+1E-10);
   mu_st(1)=Yst(1)/Bio_st(1)*S_st(1);
   Cst_pred(1)=elem_prod(N2,mu_st(1));
   N3=posfun(elem_prod(N2,1-mu_st(1)),0.5,fpen);// 
   N_hb(1)=elem_prod(N3,S_hb(1));// 
   Bio_hb(1)=sum(elem_prod(N_hb(1),Wmed(1))+1E-10);
   SSB(1)=sum(elem_prod(elem_prod(N3,matur),Wmed(1))+1E-10);
   f+=10000*fpen; 

   fpen=0.0;

   N2=N3*exp(-dt(3)*M);// purse-seine
   Bio_ps(1)=sum(elem_prod(elem_prod(N2,S_ps(1)),Wps(1))+1E-10);
   mu_ps(1)=Yps(1)/Bio_ps(1)*S_ps(1);
   Cps_pred(1)=elem_prod(N2,mu_ps(1));
   N3=posfun(elem_prod(N2,1-mu_ps(1)),0.5,fpen);// central

   N4=N3*exp(-(1.-sum(dt))*M);// final
   N4(nages)=N4(nages)/(1-exp(-1.*M));   

   f+=10000*fpen; 

// Survival by year

  for (i=1;i<nyear;i++){

   dvariable fpen=0.0;

   N(i+1)(2,nages)=++N4(1,nages-1);//
//   N(i+1)(nages)+=N4(nages);

   Biomass(i+1)=sum(elem_prod(N(i+1),Win(i+1)));
   N2=N(i+1)*exp(-dt(1)*M);// central
   Bio_ct(i+1)=sum(elem_prod(elem_prod(N2,S_ct(i+1)),Wmed(i+1))+1E-10);//central
   mu_ct(i+1)=Yct(i+1)/Bio_ct(i+1)*S_ct(i+1);
  
   Cct_pred(i+1)=elem_prod(N2,mu_ct(i+1));
   N3=posfun(elem_prod(N2,1-mu_ct(i+1)),0.5,fpen);// central
   f+=10000*fpen; 

   fpen=0.0;

   Nj(i+1)=N3*exp(-0.25*M);// It will be use for total CPUE (in July)
   N2=N3*exp(-dt(2)*M);// southern
   Bio_st(i+1)=sum(elem_prod(elem_prod(N2,S_st(i+1)),Wmed(i+1))+1E-10);
   mu_st(i+1)=Yst(i+1)/Bio_st(i+1)*S_st(i+1);
   Cst_pred(i+1)=elem_prod(N2,mu_st(i+1));
   N3=posfun(elem_prod(N2,1-mu_st(i+1)),0.5,fpen);// central

   N_hb(i+1)=elem_prod(N3,S_hb(i+1));// 
   Bio_hb(i+1)=sum(elem_prod(N_hb(i+1),Wmed(i+1))+1E-10);
   SSB(i+1)=sum(elem_prod(elem_prod(N3,matur),Wmed(i+1))+1E-10);
   f+=10000*fpen; 

   fpen=0.0;

   N2=N3*exp(-dt(3)*M);// purse-seine
   Bio_ps(i+1)=sum(elem_prod(elem_prod(N2,S_ps(i+1)),Wps(i+1))+1E-10);
   mu_ps(i+1)=Yps(i+1)/Bio_ps(i+1)*S_ps(i+1);
   Cps_pred(i+1)=elem_prod(N2,mu_ps(i+1));
   N3=posfun(elem_prod(N2,1-mu_ps(i+1)),0.5,fpen);// central
   f+=10000*fpen; 

   N4=N3*exp(-(1.-sum(dt))*M);// final
   N4(nages)=N4(nages)/(1-exp(-1.*M));   


   f+=10000*fpen; 
  }

  S_tot=mu_ct+mu_ps+mu_st+1e-10;
  for (i=1;i<=nyear;i++){
  S_tot(i)=S_tot(i)/max(S_tot(i));}
  mu_lastyr(1)=max(mu_ps(nyear));
  mu_lastyr(2)=max(mu_ct(nyear));
  mu_lastyr(3)=max(mu_st(nyear));
  mu_lastyr(4)=mu_lastyr(1)+mu_lastyr(2)+mu_lastyr(3);

  Bio_trw=rowsum(elem_prod(elem_prod(Nj,Wmed),S_tot));


FUNCTION Eval_age2length

// 
 int i, j;
 
  mu_age(1)=Lo;
  for (i=2;i<=nages;i++){
   mu_age(i)=growth(1)*(1-exp(-growth(2)))+exp(-growth(2))*mu_age(i-1);
   }
   sigma_age=cv*mu_age;

  for (i=1;i<=nages;i++){
    P1(i)=(length2-mu_age(i))/sigma_age(i);

     for (j=2;j<=nlength;j++){
       P2(i,j)=cumd_norm(P1(i,j));}}
   
  for (i=1;i<=nages;i++){
     for (j=2;j<=nlength;j++){
       P3(i,j)=P2(i,j)-P2(i,j-1);}}

  p_age2length=elem_div(P3+1e-16,outer_prod(rowsum(P3+1e-16),Ones_length));



FUNCTION Eval_C_predicted

// pred=predicted, obs=observed

// Catch-at-age and length
  Cpsl_pred=Cps_pred*p_age2length;

// landings
  Yps_pred=rowsum(elem_prod(Cps_pred,Wps));
  Yst_pred=rowsum(elem_prod(Cst_pred,Wmed));
  Yct_pred=rowsum(elem_prod(Cct_pred,Wmed));


// proprotions
  pobs_ps=elem_div(Cps,outer_prod(rowsum(Cps+1e-10),Ones_age));
  pobs_st=elem_div(Cst,outer_prod(rowsum(Cst+1e-10),Ones_age));
  pobs_ct=elem_div(Cct,outer_prod(rowsum(Cct+1e-10),Ones_age));
  pobs_psl=elem_div(Cpsl,outer_prod(rowsum(Cpsl+1e-10),Ones_length));

  ppred_ps=elem_div(Cps_pred,outer_prod(rowsum(Cps_pred+1e-10),Ones_age));
  ppred_st=elem_div(Cst_pred,outer_prod(rowsum(Cst_pred+1e-10),Ones_age));
  ppred_ct=elem_div(Cct_pred,outer_prod(rowsum(Cct_pred+1e-10),Ones_age));
  ppred_psl=elem_div(Cpsl_pred,outer_prod(rowsum(Cpsl_pred+1e-10),Ones_length));

  pobs_hb=elem_div(Chb,outer_prod(rowsum(Chb+1e-10),Ones_age));
  ppred_hb=elem_div(N_hb,outer_prod(rowsum(N_hb+1e-10),Ones_age));

  ppred_sa=(N_sa*p_age2length)/sum(N_sa);
  pobs_sa=Csal/sum(Csal);

FUNCTION Eval_stock_recruit

  Recruits=column(N,1);
  So=sum(elem_prod(elem_prod(Neq*exp(-0.67*M),matur),Wmed(1)));
  alfa=4*h*mfexp(log_Ro)/(5*h-1);

  beta=So*(1-h)/(5*h-1);
  Rpred=elem_div(alfa*SSB,beta+SSB);

FUNCTION Eval_indexes
 
 int i;

   pred_Acous=mfexp(log_qacus)*Bio_hb;
   pred_CPUEps=mfexp(log_qps)*Bio_ps;
   pred_CPUEst=mfexp(log_qst(1))*Bio_trw; 
   Bio_sa=exp(log_qsa)*sum(elem_prod(N_sa,Win(1)));// Swept-area

   for (i=1;i<=nyear;i++)
  {
     if (years(i)>=indq(1)){
     pred_CPUEst(i)=mfexp(log_qst(2))*Bio_trw(i);}

     if (years(i)>=indq(2)){
     pred_CPUEst(i)=mfexp(log_qst(3))*Bio_trw(i);}
  } 



FUNCTION Eval_likelihood
  int i;

  suma1=0; suma2=0; suma3=0; suma4=0;

  for (i=2;i<=nyear;i++)
  {
   if (Acous(i)>0){
    suma1+=square((log(Acous(i))-log(pred_Acous(i)))/cv1(i));}
   if (CPUEst(i)>0){
    suma2+=square((log(CPUEst(i))-log(pred_CPUEst(i)))/cv2(i));}
   if (CPUEps(i)>0){
    suma3+=square((log(CPUEps(i))-log(pred_CPUEps(i)))/cv3(i));}
   if(i>=Nyeq_ini){
    suma4+=square(log(Rpred(i-1))-log(Recruits(i)));} // Recruits
  }


   likeval(1)=0.5*suma1;//survey
   likeval(2)=0.5*suma2;//cpue_st
   likeval(3)=0.5*suma3;//cpue_ps

// multinomial
   likeval(7)=-nmus(1)*sum(elem_prod(pobs_ps,log(ppred_ps+1E-10)));
   likeval(8)=-nmus(2)*sum(elem_prod(pobs_st,log(ppred_st+1E-10)));
   likeval(9)=-nmus(3)*sum(elem_prod(pobs_ct,log(ppred_ct+1E-10)));
   likeval(10)=-nmus(4)*sum(elem_prod(pobs_hb,log(ppred_hb+1E-10)));
   likeval(11)=-nmus(5)*sum(elem_prod(pobs_psl,log(ppred_psl+1E-10)));
   likeval(12)=-nmus(6)*sum(elem_prod(pobs_sa,log(ppred_sa+1E-10)));

// penalties
// lognormal Ninicial y Reclutas
   penalty(1)=0.5*norm2(log_dev_No)/square(cv_penal(1));
   penalty(2)=0.5*suma4/square(cv_penal(2)); // Rpred
   penalty(3)=0.5*norm2(log_dev_Rt(1,Nyeq_ini))/square(cv_penal(3));
   penalty(4)=0.5*square((log_M-log(Mandh_pars(1)))/cv_penal(4));
   penalty(5)=0.5*square((log_h-log(Mandh_pars(2)))/cv_penal(5));
   penalty(6)=0.5*square((log_Dst(2)-log(nages)-log_A50st(1))/cv_penal(6));
 
   f+=sum(likeval)+sum(penalty); //+sprpen;


REPORT_SECTION

  rep(years)
  rep(Yst)
  rep(Yst_pred)
  rep(Yct)
  rep(Yct_pred)
  rep(Yps)
  rep(Yps_pred)
  rep(Acous)
  rep(pred_Acous)
  rep(CPUEst)
  rep(pred_CPUEst)
  rep(CPUEps)
  rep(pred_CPUEps)
  rep(SSB)
  rep(Biomass)
  rep(Recruits)
  rep(Rpred)
  rep(N)
  rep(mu_ps)
  rep(mu_st)
  rep(mu_ct)
  rep(pobs_st)
  rep(ppred_st)
  rep(pobs_ct)
  rep(ppred_ct)
  rep(pobs_ps)
  rep(ppred_ps)
  rep(pobs_hb)
  rep(ppred_hb)
  rep(Bio_sa)
  rep(p_age2length)
  rep(pobs_psl)
  rep(ppred_psl)
  rep(S_st)
  rep(S_ct)
  rep(S_ps)
  rep(S_hb)
  rep(S_sa)
  rep(pobs_sa)
  rep(ppred_sa)
  rep(likeval)
  rep(penalty)
  rep(BPRo)
  rep(h)
  rep(M)
  rep(exp(log_qacus))
  rep(sum(likeval))
  //rep(mu40)

GLOBALS_SECTION
  //const double pi = 3.141592654;
  #include <admodel.h> 
  #include <time.h>
  #include <string.h>
  //ofstream report3("rep_byr_mcmc",ios::app);
  //ofstream report4("rep_syr_mcmc",ios::app);
  //ofstream report5("rep_fyr_mcmc",ios::app);
  //ofstream report6("rep_ryr_mcmc",ios::app);

  #undef rep
  #define rep(object) report << #object "\n" << object << endl;

  #undef depur
  #define depur(object) cout << #object "\n" << object << endl; exit(1);
