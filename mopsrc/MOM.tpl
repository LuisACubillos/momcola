TOP_OF_MAIN_SECTION
  arrmblsize=300000; // 
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
  init_vector ages(1,nages)
  init_vector length2(1,nlength)
  init_matrix Cps(1,nyear,1,nages)
  init_matrix Cpsl(1,nyear,1,nlength)
  init_matrix Cst(1,nyear,1,nages)
  init_matrix Cct(1,nyear,1,nages)
  init_matrix Chb(1,nyear,1,nages)
  init_matrix Wmed(1,nyear,1,nages)
  init_matrix Win(1,nyear,1,nages)
  init_matrix Wps(1,nyear,1,nages)

  !! ad_comm::change_datafile_name("MOM_op.ctl");
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
	init_int SrType;
	init_int ErrObsType;
	init_int SrErrProc;
	init_int nsimtmp;
	//init_int numFst;
	!! ad_comm::change_datafile_name("proyectaFs.ctl");
    //init_int yr_sim
    init_int nFt
    init_vector mf(1,nFt)
	number pi;
	int opm; //activa el modelo operativo

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
 init_bounded_number log_avgR1(2,20,phase2)  //reclutamiento promedio del primer periodo (1981-1998)
 init_bounded_number log_PR2(-4,0,phase2)   //reclutamiento promedio del primer periodo (1999-2010)


// catchabilities
 init_vector log_qst(1,3,phase6) // trawl
 init_number log_qps(phase6) // purse-seine
 init_number log_qacus(phase7) // acoustic
 init_number log_qsa(phase8) // swept-area
 vector qst(1,3)
 number qps
 number qacus
 number qsa 

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
  number beta1
  number beta2
  number Bio_sa
  sdreport_number log_Ro
  number M;
  number h;
  number BPRo
  number nn
  number log_Rmed1
  number log_Rmed2
  number Prec  
  
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
  
  
  // ####### ACA SE ESPECIFICA LAS MATRICES, VECTORES Y NUMEROS DEL MOP en la fase de proyeccion
  vector Yps_fut(nyear+1,nyear+nsimtmp);
  vector Yst_fut(nyear+1,nyear+nsimtmp);
  vector Yct_fut(nyear+1,nyear+nsimtmp);
  vector N2_fut(1,nages);
  vector N3_fut(1,nages);
  vector N4_fut(1,nages);
  vector Biomass_fut(nyear+1,nyear+nsimtmp);
  matrix N_fut(nyear+1,nyear+nsimtmp,1,nages);
  vector SSB_fut(nyear+1,nyear+nsimtmp);
  vector Bio_hb_fut(nyear+1,nyear+nsimtmp);
  vector Bio_st_fut(nyear+1,nyear+nsimtmp);
  vector Bio_ps_fut(nyear+1,nyear+nsimtmp);
  vector Bio_ct_fut(nyear+1,nyear+nsimtmp);
  vector Bio_trw_fut(nyear+1,nyear+nsimtmp);
  vector acustica_epsilon(nyear+1,nyear+nsimtmp);
  vector cpuest_epsilon(nyear+1,nyear+nsimtmp);
  vector cpueps_epsilon(nyear+1,nyear+nsimtmp);
  vector rec_epsilon(nyear+1,nyear+nsimtmp);
  vector sim_Acous(nyear+1,nyear+nsimtmp);
  vector sim_CPUEst(nyear+1,nyear+nsimtmp);
  vector sim_CPUEps(nyear+1,nyear+nsimtmp);
  vector yrs(nyear+1,nyear+nsimtmp);
  
  matrix S_ps_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix S_st_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix S_ct_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix S_hb_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix N_hb_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix mu_ct_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix mu_st_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix mu_ps_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Cps_pred_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Cst_pred_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Cct_pred_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Cpsl_pred_fut(nyear+1,nyear+nsimtmp,1,nlength);
  matrix Win_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Wmed_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Wps_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix Nj_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix S_tot_fut(nyear+1,nyear+nsimtmp,1,nages);
  matrix catch_fut(1,nFt,nyear+1,nyear+nsimtmp);
  
  //Para la funcion de desempeño
  matrix fut_Btot(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_SSBt(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_Rt(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_mups(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_muct(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_must(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_mutot(1,nFt,nyear+1,nyear+nsimtmp);
  //Para guardar las estimaciones futuras del estimador
  matrix keep_Btot(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_SSBt(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_Rt(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_Yps(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_Yct(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_Yst(1,nFt,nyear+1,nyear+nsimtmp);
  matrix fut_Ytot(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_mups(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_muct(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_must(1,nFt,nyear+1,nyear+nsimtmp);
  matrix keep_mutot(1,nFt,nyear+1,nyear+nsimtmp);
  
  //proyeccion de un yr calculo de CTP
  number Nplus; //Abundancia 1 yr despues
  number mean_BDvp;
  vector Nvp(1,nages)  //N sin explotacion
  vector NDvp(1,nages); //Numero desovante en ausencia de pesca
  matrix BDvp(1,nFt,nyear+1,nyear+nsimtmp);
  matrix RPDdin(1,nFt,nyear+1,nyear+nsimtmp);
  matrix RPDeqp(1,nFt,nyear+1,nyear+nsimtmp);
  //matrix RPDeqm(nyear+1,nyear+yr_sim,1,nFt);

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
  pi = 3.141592653590;

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
  
  // AQUI LLAMA AL MODELO OPERATIVO ya sea en la ultima fase o despues del mcmc con -mceval 
  //if(last_phase()){Oper_model();}
  if(mceval_phase())
  {
	  Oper_model();
  }

FUNCTION Oper_model

     
	 //Vectores de error de proceso y observacion para los indices:
	 dvector ran_rec(nyear+1,nyear+nsimtmp); //random para reclutamiento
	 dvector ran_acustica(nyear+1,nyear+nsimtmp); //random para cruceros acusticos
	 dvector ran_cpuest(nyear+1,nyear+nsimtmp);
	 dvector ran_cpueps(nyear+1,nyear+nsimtmp);
	 
	 int yrfav; //year con un nuevo regimen favorable (elegido entre 2018 y 2022)
	 int upk; //actualiza los años
	 int l;
	 int k;
	 int numyear;
	 dvariable mupen=0.01;
	 
	 simname = "esmeco.dat"; //nombre de los datos del estimador de ifop
	 dvector CatchNow(1,3); //para recibir las capturas del estimador
     dvector eval_now(1,7); //son 7 los indicadores del estimador
	 //dvector p(1,nages); //vector proporcion de edades
	 //dvector freq(1,nages); //muestra de comp x edad
 	 ofstream SaveOM("OPM_Out.csv",ios::app);
	 
	  	 
	 //Lectura del algoritmo
	 ran_rec.fill_randn(rng);
	 ran_acustica.fill_randn(rng);
	 ran_cpuest.fill_randn(rng);
	 ran_cpueps.fill_randn(rng);	 
	 yrfav=(int)round((46-38)*randu(rng)+38); //busca un salto aleatorio para un regimen favorable entre el yr 2025 y 2015 (5 years)
     //Sigmas
	 dvariable sigma_acustica=cv1(nyear);
	 dvariable sigma_cpuest=cv2(nyear);
	 dvariable sigma_cpueps=cv3(nyear);
	 dvariable sigma_rec=cv_penal(2);

	 //corre el modelo esmeco para estimar las ctps iniciales q se leeran aqui como capturas
	 for(l=1;l<=nFt;l++) //por cada CTP
	 {
		 int rv=system("cp esmeco_real.dat esmeco.dat");
		 rv=system("./esmeco -nox -nohess");
		 //comienza el ciclo por cada una de las nFt CTPS
		 for(int i=nyear+1;i<=nyear+nsimtmp;i++)
		 {
			 //Error de observacionn de los indices
			 rec_epsilon(i)=(mfexp(ran_rec(i)*sigma_rec)); //desv. anual reclutamiento
			 acustica_epsilon(i)=(mfexp(ran_acustica(i)*sigma_acustica));
			 cpuest_epsilon(i)=(mfexp(ran_cpuest(i)*sigma_cpuest));
			 cpueps_epsilon(i)=(mfexp(ran_cpueps(i)*sigma_cpueps));
			 
			 //Lee las Ctps 
			 ifstream CTP_tmp("ctpmcol.dat");
			 CTP_tmp >> CatchNow;
			 CTP_tmp.close();
			 //Se usa 95% ya que se deja el 5% para investigación
			 catch_fut(l,i) = 0.95*CatchNow(l);
			 //Reparte la cuota entre pesquerias (70% V-X Región, y 30% X-XII)
			 Yps_fut(i)=0.01*catch_fut(l,i);
			 Yct_fut(i)=0.59*catch_fut(l,i);
			 Yst_fut(i)=0.4*catch_fut(l,i);
			 //Guarda las capturas por flota
			 fut_Yps(l,i)=Yps_fut(i);
			 fut_Yct(l,i)=Yct_fut(i);
			 fut_Yst(l,i)=Yst_fut(i);
			 fut_Ytot(l,i)=Yps_fut(i)+Yct_fut(i)+Yst_fut(i);
			 
			 //cout << "catch_fut " << catch_fut(l,i)<<endl;exit(1);
			 //DINAMICA EN NUMERO DEL FUTURO (yrfut 1)
			 if(i==nyear+1)
			 {
				 //copia el peso inicial del último año en el peso inicial del futuro
				 //declarar la matriz Win_fut y Wmed_fut
				 Win_fut(i)=Win(nyear);
				 //cout<< "Win_fut "<< Win_fut <<endl;exit(1);
				 
				 Wmed_fut(i)=Wmed(nyear);
				 Wps_fut(i)=Wps(nyear);
				 //cout<< "Wps_fut "<< Wps_fut <<endl;exit(1);
				 
				 //copia selectividades del uñtimo año al futuro
				 S_ct_fut(i)=S_ct(nyear);
				 S_st_fut(i)=S_st(nyear);
				 S_ps_fut(i)=S_ps(nyear);
				 S_hb_fut(i)=S_hb(nyear);
				 S_tot_fut(i)=S_tot(nyear);
				 //cout<< "S_tot_fut "<< S_tot_fut <<endl;exit(1);
				 
				 //Aqui considerar la sobrevivencia de 2010 al 2011
				 N_fut(i)(2,nages)=++N4(1,nages-1); //N4 desde la edad 1 hasta la n-1
				 //En ausencia de pesca
				 Nvp = N(nyear); //Copia la abundancia del ultimo anio
				 Nplus = Nvp(nages)*exp(-1.0*M);
	   			 Nvp(2,nages) = ++Nvp(1,nages-1)*exp(-1.0*M);
	   			 Nvp(nages) = Nvp(nages) + Nplus;// Grupo plus en ausencia de pesca
				 if(SrType==1)
				 {
					 //Relacion S-R regimen 1er con error de proceso (rec_epsilon)
					 if(SrErrProc==1)
					 {
						N_fut(i,1)=SSB(nyear)/(alfa+beta2*SSB(nyear))*rec_epsilon(i); 					 	
 					 	Nvp(1)=N_fut(i,1);
					 }
					 if(SrErrProc==2)
					 {
					 	N_fut(i,1)=SSB(nyear)/(alfa+beta2*SSB(nyear));
					 	Nvp(1)=N_fut(i,1);
					 }
 //Regimen bajo para el 1er year				
				 }
				 else
				 {
					 //Trigonometrico
					 if(SrErrProc==1)
					 {
	 				 	N_fut(i,1)=(mean(Recruits)+exp(7.150)*cos(pi*0.0662*double(i)-2.4838))*rec_epsilon(i);
					 	Nvp(1)=N_fut(i,1);
					 }
					 if(SrErrProc==2)
					 {
 	 				 	N_fut(i,1)=mean(Recruits)+exp(7.150)*cos(pi*0.0662*double(i)-2.4838);
					 	Nvp(1)=N_fut(i,1);
					 }
				 }
				 
				 //Dinamica de la abundancia
				 Biomass_fut(i)=sum(elem_prod(N_fut(i),Win_fut(i)));
				 //Biomasa desovante sin pesca
				 NDvp=elem_prod(Nvp*exp(-0.67*M),matur);
				 BDvp(l,i)=sum(elem_prod(NDvp,Wmed_fut(i)));

				 //cout<< "N_fut "<< N_fut(i) <<endl;
				 //cout<< "B_fut "<< Biomass_fut(i) <<endl;exit(1);
				 
				 N2_fut=N_fut(i)*exp(-dt(1)*M); 
				 Bio_ct_fut(i)=sum(elem_prod(N2_fut,S_ct_fut(i)));
				 mu_ct_fut(i)=Yct_fut(i)/Bio_ct(i)*S_ct_fut(i);
				 Cct_pred_fut(i)=elem_prod(N2_fut,mu_ct_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_ct_fut(i)),0.5,mupen);
				 
				 Nj_fut(i)=N3_fut*exp(-0.25*M); //abundancia al 1 de julio
				 N2_fut=N3_fut*exp(-dt(2)*M); //sur
				 Bio_st_fut(i)=sum(elem_prod(elem_prod(N2_fut,S_st_fut(i)),Wmed_fut(i))+1E-10);
				 mu_st_fut(i)=Yst_fut(i)/Bio_st_fut(i)*S_st_fut(i);
				 Cst_pred_fut(i)=elem_prod(N2_fut,mu_st_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_st_fut(i)),0.5,mupen);
				 N_hb_fut(i)=elem_prod(N3_fut,S_hb_fut(i));
				 Bio_hb_fut(i)=sum(elem_prod(N_hb_fut(i),Wmed_fut(i))+1E-10);
				 SSB_fut(i)=sum(elem_prod(elem_prod(N3_fut,matur),Wmed_fut(i))+1E-10);
				 
				 N2_fut=N3_fut*exp(-dt(3)*M);
				 Bio_ps_fut(i)=sum(elem_prod(elem_prod(N2_fut,S_ps_fut(i)),Wps_fut(i))+1E-10);
				 mu_ps_fut(i)=Yps_fut(i)/Bio_ps_fut(i)*S_ps_fut(i);
				 Cps_pred_fut(i)=elem_prod(N2_fut,mu_ps_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_ps_fut(i)),0.5,mupen);
				 
				 N4_fut=N3_fut*exp(-(1.-sum(dt))*M);
				 N4_fut(nages)=N4_fut(nages)/(1-exp(-1.*M));
				 				 
				 Bio_trw_fut(i)=sum(elem_prod(elem_prod(Nj_fut(i),Wmed_fut(i)),S_tot_fut(i)));
				 
				 //Catch-at-length
				 Cpsl_pred_fut(i)=Cps_pred_fut(i)*p_age2length;
				 //cout<< "Cpsl_pred_fut "<< Bio_trw_fut(i) <<endl;exit(1);
				 
				 //Simula los indices
				 if(ErrObsType==1)
				 {
					 //sin error de observacion
				 	sim_Acous(i)=qacus*Bio_hb_fut(i);
					if(Yps_fut(i)>0)
					{
						sim_CPUEps(i)=qps*Bio_ps_fut(i);
					}
					else
					{
						sim_CPUEps(i)=0;
					} 				 	
				 	sim_CPUEst(i)=qst(3)*Bio_trw_fut(i);
				 }
				 else
				 {
					 //Con error de observacion
 				 	sim_Acous(i)=qacus*Bio_hb_fut(i)*acustica_epsilon(i);
					if(Yps_fut(i)>0)
					{
						sim_CPUEps(i)=qps*Bio_ps_fut(i)*cpueps_epsilon(i);
					}
					else
					{
						sim_CPUEps(i)=0;
					} 				 	
 				 	sim_CPUEst(i)=qst(3)*Bio_trw_fut(i)*cpuest_epsilon(i);
				 }
				 //cout<<"SIM cpuest"<< sim_CPUEst(i)<<endl;exit(1);
				 yrs(i)=years(1)+double(i)-1;
				 upk=i;
				 numyear = i;
				 //AHORA ESCRIBE EN EL ARCHIVO *.dat DEL ESTIMADOR
				 ofstream simdata(simname);
				 simdata << "#num_yr" << endl;
				 simdata << numyear << endl;
				 simdata << "#num_edades"<<endl;
				 simdata << nages <<endl;
				 simdata << "#num_tallas"<<endl;
				 simdata << nlength << endl;
				 simdata << "#Data historica"<<endl;
				 simdata <<"#Yr	HB cv CPUEst cv CPUEps cv Yps cv Yst cv Yct cv "<<endl;
				 				 simdata << year_matrix <<endl;
				 simdata <<"#Data simulada futura"<<endl;
				 simdata <<"#Yr	HB cv CPUEst cv CPUEps cv Yps cv Yst cv Yct cv "<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << yrs(k) << " " << sim_Acous(k) << " " << sigma_acustica << " " << sim_CPUEst(k) << " " << sigma_cpuest << " " << sim_CPUEps(k) << " " << sigma_cpueps << " " << Yps_fut(k) << " " << 0.05 << " " << Yst_fut(k) << " " << 0.05 << " " << Yct_fut(k) << " " << 0.05 << endl;
				 }
				 simdata <<"#edades"<<endl;
				 simdata << ages << endl;
				 simdata << "#tallas" << endl;
				 simdata << length2 << endl;
				 simdata << "#Captura a la edad Cerco"<<endl;
				 simdata << Cps << endl;
				 simdata << "#Captura a la edad cerco Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cps_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la talla cerco "<<endl;
				 simdata << Cpsl <<endl;
				 simdata << "#Captura a la talla simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cpsl_pred_fut(k) << endl;
				 }
				 simdata << "#Captura a la edad arrastre sur"<<endl;
				 simdata << Cst << endl;
				 simdata << "#Captura a la edad arrastre sur Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cst_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la edad arrastre centro"<<endl;
				 simdata << Cct << endl;
				 simdata << "#Captura a la edad arrastre centro Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cct_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la edad crucero"<<endl;
				 simdata << Chb << endl;
				 simdata << "#Captura a la edad crucero Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << qacus*N_hb_fut(k)<<endl; 
				 }
				 simdata << "#Peso medio a la edad"<<endl;
				 simdata << Wmed << endl;
				 simdata << "#Peso medio sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Wmed_fut(k) << endl;
				 }
				 simdata << "#Peso inicial a la edad"<<endl;
				 simdata << Win << endl;
				 simdata << "#Peso inicial sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Win_fut(k) << endl;
				 }
				 simdata << "#Peso cerco a la edad"<<endl;
				 simdata << Wps << endl;
				 simdata << "#Peso cerco sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Wps_fut(k) << endl;
				 }
				 //cout<<"N:"<<N(nyear)<<endl;
				 rv=system("./esmeco -nox -nohess");
				 //lee los indices de evaluacion actuales y los guarda para
				 //evaluar el desempeño
				 ifstream eval_indices("indices.dat");
				 eval_indices >> eval_now;
				 eval_indices.close();
				 //Guarda los indices del estimador
				 keep_Btot(l,i)=eval_now(1);
				 keep_SSBt(l,i)=eval_now(2);
				 keep_Rt(l,i)=eval_now(3);
				 keep_mups(l,i)=eval_now(4);
				 keep_muct(l,i)=eval_now(5);
				 keep_must(l,i)=eval_now(6);
				 keep_mutot(l,i)=eval_now(7);
				 //Guarda los indices del Mop
				 fut_Btot(l,i)=Biomass_fut(i);
				 fut_SSBt(l,i)=SSB_fut(i);
				 fut_Rt(l,i)=N_fut(i,1);
				 fut_mups(l,i)=max(mu_ps_fut(i));
				 fut_muct(l,i)=max(mu_ct_fut(i));
				 fut_must(l,i)=max(mu_st_fut(i));
				 fut_mutot(l,i)=fut_mups(l,i)+fut_muct(l,i)+fut_must(l,i);
			 }
			 else
			 {
				 //copia el peso inicial del último año en el peso inicial del futuro
				 //declarar la matriz Win_fut y Wmed_fut
				 Win_fut(i)=Win(nyear);
				 Wmed_fut(i)=Wmed(nyear);
				 Wps_fut(i)=Wps(nyear);
				 //copia selectividades del uñtimo año al futuro
				 S_ct_fut(i)=S_ct(nyear);
				 S_st_fut(i)=S_st(nyear);
				 S_ps_fut(i)=S_ps(nyear);
				 S_hb_fut(i)=S_hb(nyear);
				 S_tot_fut(i)=S_tot(nyear);
				 
				 //Aqui considerar la sobrevivencia de 2010 al 2011
				 N_fut(i)(2,nages)=++N4_fut(1,nages-1); //N4 desde la edad 1 hasta la n-1
				 //En ausencia de pesca
				 Nvp = N_fut(i-1); //Copia la abundancia del ultimo anio proyectado
				 Nplus = Nvp(nages)*exp(-1.0*M);
	   			 Nvp(2,nages) = ++Nvp(1,nages-1)*exp(-1.0*M);
	   			 Nvp(nages) = Nvp(nages) + Nplus;// Grupo plus en ausencia de pesca		 
				 if(SrType==1)
				 {
					 //Relacion S-R regimen 1er con error de proceso (rec_epsilon)
					 if(SrErrProc==1)
					 {
			 			if(i<=yrfav)
			 			{
							N_fut(i,1)=SSB_fut(i-1)/(alfa+beta2*SSB_fut(i-1))*rec_epsilon(i);	//desfavorable
							Nvp(1)=N_fut(i,1);			
			 			}
			 			else
			 			{
							N_fut(i,1)=SSB_fut(i-1)/(alfa+beta1*SSB_fut(i-1))*rec_epsilon(i);  //favorable
							Nvp(1)=N_fut(i,1);			
			 			}
					 }
					 if(SrErrProc==2)
					 {
 			 			if(i<=yrfav)
 			 			{
							N_fut(i,1)=SSB_fut(i-1)/(alfa+beta2*SSB_fut(i-1));	//desfavorable
							Nvp(1)=N_fut(i,1);			
 			 			}
 			 			else
 			 			{
							N_fut(i,1)=SSB_fut(i-1)/(alfa+beta1*SSB_fut(i-1));  //favorable
							Nvp(1)=N_fut(i,1);			
 			 			}
					 }
//Regimen bajo para el 1er year				
				 }
				 else
				 {
					 //Trigonometrico
					 if(SrErrProc==1)
					 {
	 				 	N_fut(i,1)=(mean(Recruits)+exp(7.150)*cos(pi*0.0662*double(i)-2.4838))*rec_epsilon(i);
						Nvp(1)=N_fut(i,1);			
					 }
					 if(SrErrProc==2)
					 {
 	 				 	N_fut(i,1)=mean(Recruits)+exp(7.150)*cos(pi*0.0662*double(i)-2.4838);
						Nvp(1)=N_fut(i,1);			
					 	
					 }
				 }
				 			 
				 //Dinamica de la abundancia
				 Biomass_fut(i)=sum(elem_prod(N_fut(i),Win_fut(i)));
				 //Biomasa desovante sin pesca
				 NDvp=elem_prod(Nvp*exp(-0.67*M),matur);
				 BDvp(l,i)=sum(elem_prod(NDvp,Wmed_fut(i)));
				 
				 
				 N2_fut=N_fut(i)*exp(-dt(1)*M);
				 Bio_ct_fut(i)=sum(elem_prod(N2_fut,S_ct_fut(i)));
				 mu_ct_fut(i)=Yct_fut(i)/Bio_ct(i)*S_ct_fut(i);
				 Cct_pred_fut(i)=elem_prod(N2_fut,mu_ct_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_ct_fut(i)),0.5,mupen);
				 
				 Nj_fut(i)=N3_fut*exp(-0.25*M);
				 N2_fut=N3_fut*exp(-dt(2)*M); //sur
				 Bio_st_fut(i)=sum(elem_prod(elem_prod(N2_fut,S_st_fut(i)),Wmed_fut(i))+1E-10);
				 mu_st_fut(i)=Yst_fut(i)/Bio_st_fut(i)*S_st_fut(i);
				 Cst_pred_fut(i)=elem_prod(N2_fut,mu_st_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_st_fut(i)),0.5,mupen);
				 N_hb_fut(i)=elem_prod(N3_fut,S_hb_fut(i));
				 Bio_hb_fut(i)=sum(elem_prod(N_hb_fut(i),Wmed_fut(i))+1E-10);
				 SSB_fut(i)=sum(elem_prod(elem_prod(N3_fut,matur),Wmed_fut(i))+1E-10);
				 
				 N2_fut=N3_fut*exp(-dt(3)*M);
				 Bio_ps_fut(i)=sum(elem_prod(elem_prod(N2_fut,S_ps_fut(i)),Wps_fut(i))+1E-10);
				 mu_ps_fut(i)=Yps_fut(i)/Bio_ps_fut(i)*S_ps_fut(i);
				 Cps_pred_fut(i)=elem_prod(N2_fut,mu_ps_fut(i));
				 N3_fut=posfun(elem_prod(N2_fut,1-mu_ps_fut(i)),0.5,mupen);
				 
				 N4_fut=N3_fut*exp(-(1.-sum(dt))*M);
				 N4_fut(nages)=N4_fut(nages)/(1-exp(-1.*M));
				 
				 Bio_trw_fut(i)=sum(elem_prod(elem_prod(Nj_fut(i),Wmed_fut(i)),S_tot_fut(i)));
				 //Catch-at-length
				 Cpsl_pred_fut(i)=Cps_pred_fut(i)*p_age2length;
				 //Simula los indices
				 if(ErrObsType==1)
				 {
					 //sin error de observacion
				 	sim_Acous(i)=qacus*Bio_hb_fut(i);
					if(Yps_fut(i)>0)
					{
						sim_CPUEps(i)=qps*Bio_ps_fut(i);
					}
					else
					{
						sim_CPUEps(i)=0;
					} 				 	
				 	sim_CPUEst(i)=qst(3)*Bio_trw_fut(i);
				 }
				 else
				 {
					 //Con error de observacion
 				 	sim_Acous(i)=qacus*Bio_hb_fut(i)*acustica_epsilon(i);
					if(Yps_fut(i)>0)
					{
						sim_CPUEps(i)=qps*Bio_ps_fut(i)*cpueps_epsilon(i);
					}
					else
					{
						sim_CPUEps(i)=0;
					} 				 	
 				 	sim_CPUEst(i)=qst(3)*Bio_trw_fut(i)*cpuest_epsilon(i);
				 }
				 yrs(i)=years(1)+i-1;
				 upk=i;
				 numyear = i;
				 //AHORA ESCRIBE EN EL ARCHIVO *.dat DEL ESTIMADOR
				 ofstream simdata(simname);
				 simdata << "#num_yr" << endl;
				 simdata << numyear << endl;
				 simdata << "#num_edades"<<endl;
				 simdata << nages <<endl;
				 simdata << "#num_tallas"<<endl;
				 simdata << nlength << endl;
				 simdata << "#Data historica"<<endl;
				 simdata <<"#Yr	HB cv CPUEst cv CPUEps cv Yps cv Yst cv Yct cv "<<endl;
				 simdata << year_matrix <<endl;
				 simdata <<"#Data simulada futura"<<endl;
				 simdata <<"#Yr	HB cv CPUEst cv CPUEps cv Yps cv Yst cv Yct cv "<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << yrs(k) << " " << sim_Acous(k) << " " << sigma_acustica << " " << sim_CPUEst(k) << " " << sigma_cpuest << " " << sim_CPUEps(k) << " " << sigma_cpueps << " " << Yps_fut(k) << " " << 0.05 << " " << Yst_fut(k) << " " << 0.05 << " " << Yct_fut(k) << " " << 0.05 << endl;
				 }
				 simdata <<"#edades"<<endl;
				 simdata << ages << endl;
				 simdata << "#tallas" << endl;
				 simdata << length2 << endl;
				 simdata << "#Captura a la edad Cerco"<<endl;
				 simdata << Cps << endl;
				 simdata << "#Captura a la edad cerco Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cps_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la talla cerco "<<endl;
				 simdata << Cpsl <<endl;
				 simdata << "#Captura a la talla simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cpsl_pred_fut(k) << endl;
				 }
				 simdata << "#Captura a la edad arrastre sur"<<endl;
				 simdata << Cst << endl;
				 simdata << "#Captura a la edad arrastre sur Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cst_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la edad arrastre centro"<<endl;
				 simdata << Cct << endl;
				 simdata << "#Captura a la edad arrastre centro Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << Cct_pred_fut(k)<<endl; 
				 }
				 simdata << "#Captura a la edad crucero"<<endl;
				 simdata << Chb << endl;
				 simdata << "#Captura a la edad crucero Simulado"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
				 	 simdata << qacus*N_hb_fut(k)<<endl; 
				 }
				 simdata << "#Peso medio a la edad"<<endl;
				 simdata << Wmed << endl;
				 simdata << "#Peso medio sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Wmed_fut(k) << endl;
				 }
				 simdata << "#Peso inicial a la edad"<<endl;
				 simdata << Win << endl;
				 simdata << "#Peso inicial sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Win_fut(k) << endl;
				 }
				 simdata << "#Peso cerco a la edad"<<endl;
				 simdata << Wps << endl;
				 simdata << "#Peso cerco sim"<<endl;
				 for(k=nyear+1;k<=upk;k++)
				 {
					 simdata << Wps_fut(k) << endl;
				 }
				 rv=system("./esmeco -nox -nohess");
				 ifstream eval_indices("indices.dat");
				 eval_indices >> eval_now;
				 eval_indices.close();
				 //Guarda los indices del estimador
				 keep_Btot(l,i)=eval_now(1);
				 keep_SSBt(l,i)=eval_now(2);
				 keep_Rt(l,i)=eval_now(3);
				 keep_mups(l,i)=eval_now(4);
				 keep_muct(l,i)=eval_now(5);
				 keep_must(l,i)=eval_now(6);
				 keep_mutot(l,i)=eval_now(7);
				 //Guarda los indices del Mop
				 fut_Btot(l,i)=Biomass_fut(i);
				 fut_SSBt(l,i)=SSB_fut(i);
				 fut_Rt(l,i)=N_fut(i,1);
				 fut_mups(l,i)=max(mu_ps_fut(i));
				 fut_muct(l,i)=max(mu_ct_fut(i));
				 fut_must(l,i)=max(mu_st_fut(i));
				 fut_mutot(l,i)=fut_mups(l,i)+fut_muct(l,i)+fut_must(l,i);
				 //cout << "==========fut_Btot:======"<< " "
				 //     << fut_Btot<<endl;
				 //cout << "keep_Btot:"<<keep_Btot<<endl;				 
			 }
			 //Variables de desempeño
			 //razon desovante potencial dinamica
			 //RPDdin(i,l)=fut_SSBt(i,l)/BDvp(i,l);	 		 
		 }
	     //cout << "fut_Btot:"<<fut_Btot<<endl;
	     //cout << "keep_Btot:"<<keep_Btot<<endl;exit(1);

	 }
	 
	 
	 RPDdin=elem_div(fut_SSBt,BDvp);
	 //Biomasa total futura modelo operante
	 byr_out<< fut_Btot <<endl;
	 syr_out<< fut_SSBt <<endl;
	 ryr_out<< fut_Rt << endl;
	 mps_out<< fut_mups<<endl;
	 mct_out<< fut_muct<<endl;
	 mst_out<< fut_must<<endl;
	 mt_out<< fut_mutot<<endl;//mt_out.close();
	 yps_out << fut_Yps<<endl;//yps_out.close();
	 yct_out << fut_Yct<<endl;//yct_out.close();
	 yst_out << fut_Yst<<endl;//yst_out.close();
	 ytot_out << fut_Ytot<<endl;//ytot_out.close();

	 //Biomasa total futura estimador
	 ebyr_out<< keep_Btot <<endl;//ebyr_out.close();
	 esyr_out<< keep_SSBt <<endl;//esyr_out.close();
	 eryr_out<< keep_Rt << endl;//eryr_out.close();
	 emps_out<< keep_mups<<endl;//emps_out.close();
	 emct_out<< keep_muct<<endl;//emct_out.close();
	 emst_out<< keep_must<<endl;//emst_out.close();
	 emt_out<< keep_mutot<<endl;//emt_out.close();
	 eytot_out << catch_fut <<endl;//eytot_out.close();



	 
	 //Razon de potencial reproductivo dinamico
	 rpr_out << RPDdin <<endl;//rpr_out.close();
	 //Razon de potencial reproductivo estatico	 
	 rpr2_out << fut_SSBt/So <<endl;//rpr2_out.close();
     yrfav_out << yrfav << endl;


 /*  
   FUNCTION CTP_cal

	ofstream ctp_out("ctpmcol.dat");
	ofstream ssb_out("ssbmcol.dat");
	ofstream rpr1_out("rpr1mcol.dat");
  	for (int j=1; j<=nFt; j++) // ciclo de multiplicadores de F
  	{
  		Np = N(nyear); //Abundancia en el ultimo anio
  		Nvp = N(nyear); //Abundancia en el ultimo anio en ausencia de pesca
  		Rp = mean(Recruits(nyear-5,nyear)); // Supuesto: reclutamiento constante igual al promedio de los ultimos 5 a√±os
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
	ssb_out << BDp << endl;ssb_out.close();
	rpr1_out<< RPDdin << endl; rpr1_out.close();

  */
 
  //FUNCTION MCWrite
  //
  //report3 << Biomass_fut << endl;
  //report4 << SSB_fut << endl;
  //report6 << Yst_fut << endl;

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
  log_Rmed1=log_avgR1;
  Prec=exp(log_PR2);
  log_Rmed2=log_PR2+log_avgR1;

// Abundance at the first year by age
  Neq=Neq*mfexp(log_Ro);// steady-state as reference

  if (phase4<0){
  N(1)=Neq;}
  else{
  N(1)=mfexp(log(Neq)+log_dev_No-0.5*square(cv_penal(2)));}

 // Recruitments by year starting in t=1
  for (i=2;i<=20;i++){
    N(i,1)=mfexp(log_Rmed1+log_dev_Rt(i-1)-0.5*square(cv_penal(3)));
    }
    for (i=21;i<=nyear;i++){
      N(i,1)=mfexp(log_Rmed2+log_dev_Rt(i-1)-0.5*square(cv_penal(3)));
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
  //alfa=4*h*mfexp(log_Ro)/(5*h-1);
  //beta=So*(1-h)/(5*h-1);
  alfa=(1-h)*(So/mfexp(log_Rmed1))/(4*h);
  beta1=(5*h-1)/(4*h*mfexp(log_Rmed1));
  beta2=(5*h-1)/(4*h*mfexp(log_Rmed2));
  Rpred(1)=mfexp(log_Rmed1);
  for(int i=2;i<=20;i++){
	  Rpred(i)=SSB(i-1)/(alfa+beta1*SSB(i-1));
  }
  for(int i=21;i<=nyear;i++){
	  Rpred(i)=SSB(i-1)/(alfa+beta2*SSB(i-1));
  }
  
  //Rpred=elem_div(alfa*SSB,beta+SSB);

FUNCTION Eval_indexes
 
 int i;
   qacus=mfexp(log_qacus);
   qps=mfexp(log_qps);
   qst=mfexp(log_qst);
   qsa=mfexp(log_qsa);
   
   pred_Acous=qacus*Bio_hb;
   pred_CPUEps=qps*Bio_ps;
   pred_CPUEst=qst(1)*Bio_trw; 
   Bio_sa=qsa*sum(elem_prod(N_sa,Win(1)));// Swept-area

   for (i=1;i<=nyear;i++)
  {
     if (years(i)>=indq(1)){
     pred_CPUEst(i)=qst(2)*Bio_trw(i);}

     if (years(i)>=indq(2)){
     pred_CPUEst(i)=qst(3)*Bio_trw(i);}
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
   likeval(7)=-nmus(1)*sum(elem_prod(pobs_ps,log(ppred_ps)));
   likeval(8)=-nmus(2)*sum(elem_prod(pobs_st,log(ppred_st)));
   likeval(9)=-nmus(3)*sum(elem_prod(pobs_ct,log(ppred_ct)));
   likeval(10)=-nmus(4)*sum(elem_prod(pobs_hb,log(ppred_hb)));
   likeval(11)=-nmus(5)*sum(elem_prod(pobs_psl,log(ppred_psl)));
   likeval(12)=-nmus(6)*sum(elem_prod(pobs_sa,log(ppred_sa)));

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
  rep(alfa)
  rep(beta1)
  rep(beta2)
  rep(exp(log_qacus))
  rep(sum(likeval))
  //rep(mu40)

GLOBALS_SECTION
  //const double pi = 3.141592654;
  #include <admodel.h> 
  #include <time.h>
  #include <string.h>
  adstring simname;
  
 //ESCRIBE ARCHIVOS DE SALIDA DEL MOP
 ofstream byr_out("byr_mcmc_op",ios::app);
 ofstream syr_out("syr_mcmc_op",ios::app);
 ofstream ryr_out("ryr_mcmc_op",ios::app);
 ofstream mps_out("mps_mcmc_op",ios::app);
 ofstream mct_out("mct_mcmc_op",ios::app);
 ofstream mst_out("mst_mcmc_op",ios::app);
 ofstream mt_out("mtot_mcmc_op",ios::app);
 ofstream yps_out("yps_mcmc_op",ios::app);
 ofstream yct_out("yct_mcmc_op",ios::app);
 ofstream yst_out("yst_mcmc_op",ios::app);
 ofstream ytot_out("ytot_mcmc_op",ios::app);
 
 //ESCRIBE ARCHIVOS DE SALIDA DEL ESTIMADOR
 ofstream ebyr_out("byr_mcmc_es",ios::app);
 ofstream esyr_out("syr_mcmc_es",ios::app);
 ofstream eryr_out("ryr_mcmc_es",ios::app);
 ofstream emps_out("mps_mcmc_es",ios::app);
 ofstream emct_out("mct_mcmc_es",ios::app);
 ofstream emst_out("mst_mcmc_es",ios::app);
 ofstream emt_out("mtot_mcmc_es",ios::app);
 //ofstream eyps_out("myps_mcmc_es");
 //ofstream eyct_out("myct_mcmc_es");
 //ofstream eyst_out("myst_mcmc_es");
 ofstream eytot_out("mytot_mcmc_es",ios::app);
  //Razon de potencial reproductivo dinamico
 ofstream rpr_out("rpr_mcmc",ios::app);
 //Razon de potencial reproductivo estatico
 ofstream rpr2_out("rpr2_mcmc",ios::app);
 ofstream yrfav_out("yrfav_mcmc",ios::app);
  

  #undef rep
  #define rep(object) report << #object "\n" << object << endl;

  #undef depur
  #define depur(object) cout << #object "\n" << object << endl; exit(1);
