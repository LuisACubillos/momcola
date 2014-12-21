# Sexual	maturity
#1	2	3	4	5	6	7	8	9	10	11	12	13	14																	
0.0449	0.15605	0.362	0.58693	0.7508	0.8479	0.90221	0.93305	0.9513	0.96263	0.96997	0.97492	0.97838	0.98085
#	Biomass	in swept-area	survey	1972 (thounsands) (Bsa)
720																
#	Length	comp.	in	survey	(swept-area)	of	1972					
#20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80
0	0	6	39	396	833	1026	538	299	219	524	541	421	427	335	267	205	173	150	143	177	179	169	102	87	87	52	41	23	24	21
# Sample size
# (nmus) Nm para:
# Cps & Ctrw & Ctrw_cs  &Csurvey & Clength_ps & CNlength_1972_sa
30    50       50        100        30         0
# dt between operation months
# central   southern  purse-seine
0.25        0.42      0.17
# Growth parameters (Loo, k)
101.3   0.176
# Year blocks in selectivity
#  (indTrw)  Trawl southen area (ST)
1991 2008
#  (indTrw_cs) Trawl central area (CT)
2006 2020
#  (indSrv) acoustic survey 
2008 2020
#  (indq) trawling (ST) catchability
1997 2002
# Numbers of initial years on equilibrium condition (Nyeq)
1
# Phases and estimation options
#--------------------------------------------------
#  (phase1) Selectivity
1
#   (phase2) Bo (si está en equilibrio!!)
1
#  (phase3) Rt deviates 
1
#  (phase4) No deviates (<0 means in steady-state) 
1
#  (phase5) F 
-1
#   (phase6) q CPUE 
2
#   (phase7) q Acoustic 
2
#   (phase8) q swept-area 
-2
#   (phase9) L1 growth 
-2
#   (phase10) cv growth 
2
#   (phase11) M estimation 
-1
#   (phase12) h estimation 
-2
#   Priors on M and h (lognormal)
0.35  0.75
#   penalties (cv_penal)
# Initial N  & S/R deviates  & Ro devs_yrs  & cv_M  &  cv_h  & S_st logistic (0.05) 
0.6          0.6          0.01            0.1      0.1       0.05
# SrType Tipo de relacion S-R (1:S-R; 2: Trigonometrico)
1
# ErrObsType (1:Sin error de observacion; 2: con error de observacion)
1
# SrErrProc (1: Con error de proceso del reclutamiento; 2: Sin error de proceso del reclutamiento)
2
#Numero de años con simulacion
40


