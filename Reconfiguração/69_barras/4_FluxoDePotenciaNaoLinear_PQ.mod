# conjuntos

set Ob; # conjunto de barras
set Ol within Ob cross Ob;
set Ot; #conjunto de periodos estimados

#-------------------------------------------------------------------------------------------
# parametros
param Vini{Ob};
param L{Ol};

#barras 
param Tb{Ob};
param PD{Ob};
param QD{Ob};
param QC{Ob};

param Pd{Ob, Ot};
param Qd{Ob, Ot};
 
#lineas
param R{Ol};
param X{Ol};
param Imax{Ol};
param Z2{Ol};

#para reconfiguracion
param linea{Ol};  #los numeros de las lineas
param N; # el numero de barras de la red

#datos de tension
param Vnom;
param Vmin;
param Vmax;

#para los periodos
param Ki{Ot};
param Si{Ot};
param Ti{Ot};

#parametros para la linealizacion
param Y;
param ms{Ol, Ot, y in 1..Y};
param ds{Ol, Ot};

#parametros para la linealizacion Vsqr y Isqr
param S;
param DeltaV = (Vmax^2-Vmin^2)/(S+1);

#-------------------------------------------------------------------------------------------
#variables

#variables para la linealizacion Vsqr y Isqr
var xv{j in Ob, t in Ot, s in 1..S} binary;
var Pc{(i,j) in Ol, t in Ot, s in 1..S} >=0;


#barras
var PS{Ob,Ot};	#potencia activa generada en la SE
var QS{Ob,Ot};	#potencia reactiva generada en la SE
var Vsqr{Ob,Ot};	#tensi n en las barras de la red

#lineas
var Isqr{Ol,Ot};	#corriente en las lineas
var P{Ol,Ot};	#potencia activa en la linea
var Q{Ol,Ot};	#potencia reactiva en la linea

#variables para la linealizacion

var DP{Ol,Ot, y in 1..Y} >= 0;
var DQ{Ol,Ot, y in 1..Y} >= 0;
var Pp{Ol,Ot} >= 0;
var Pn{Ol,Ot} >= 0;
var Qp{Ol,Ot} >= 0;
var Qn{Ol,Ot} >= 0;

#variables para la reconfiguracion
var b{Ol,Ot};
var y_positivo{Ol,Ot}, binary;
var y_negativo{Ol,Ot}, binary;

#-------------------------------------------------------------------------------------------
#Funcion objetivo
#1
minimize FuncionObjetivo:
	sum {t in Ot}(Ki[t]*Ti[t]*sum{(i,j) in Ol}(R[i,j]*Isqr[i,j,t]));


#-------------------------------------------------------------------------------------------

#Restricciones
#2
subject to PotenciaActiva{i in Ob, t in Ot}:
	sum{(k,i) in Ol}(P[k,i,t])- sum{(i,j) in Ol}(P[i,j,t]+R[i,j]*Isqr[i,j,t])+PS[i,t]=Pd[i,t];
	
#3
subject to PotenciaReactiva{i in Ob, t in Ot}:
	sum{(k,i) in Ol}(Q[k,i,t])- sum{(i,j) in Ol}(Q[i,j,t]+X[i,j]*Isqr[i,j,t])+QC[i]+ QS[i,t]=Qd[i,t];

#4 en 12

#6
subject to LimiteTension{i in Ob, t in Ot}:
	Vmin^2 <= Vsqr[i,t] <= Vmax^2;

#7	
subject to LimiteCorriente{(i,j) in Ol, t in Ot}:	
	Isqr[i,j,t] <= Imax[i,j]^2 * (y_positivo[i,j,t] + y_negativo[i,j,t]);

#8
subject to RestriccionLineal_p1{(i,j) in Ol, t in Ot}:
	Pp[i,j,t] - Pn[i,j,t]=P[i,j,t];
	
#9 var definition 0 <= Pp
subject to Restriccion_Recon_1{(i,j)in Ol, t in Ot}:
	Pp[i,j,t]<=Vmax*Imax[i,j]*y_positivo[i,j,t];
	
#10 var definition 0 <= Pn
subject to Restriccion_Recon_2{(i,j)in Ol, t in Ot}:
	Pn[i,j,t]<=Vmax*Imax[i,j]*y_negativo[i,j,t];
		
#11
subject to Restriccion_Recon_3{(i,j) in Ol, t in Ot}:
	-Vmax*Imax[i,j]*(y_positivo[i,j,t]+y_negativo[i,j,t]) <= Q[i,j,t];

subject to Restriccion_Recon_4{(i,j) in Ol, t in Ot}:
	Q[i,j,t] <= Vmax*Imax[i,j]*(y_positivo[i,j,t]+y_negativo[i,j,t]);
	
#12
subject to queda_de_tensao_1 {(i,j) in Ol, t in Ot}:
	Vsqr[i,t] - 2 * (R[i,j] * P[i,j,t] + X[i,j] * Q[i,j,t]) - Z2[i,j] * Isqr[i,j,t] - Vsqr[j,t] 
	>= - (Vmax^2 - Vmin^2) * (1 - (y_positivo[i,j,t]+y_negativo[i,j,t]));
	
subject to queda_de_tensao_2 {(i,j) in Ol, t in Ot}:
	Vsqr[i,t] - 2 * (R[i,j] * P[i,j,t] + X[i,j] * Q[i,j,t]) - Z2[i,j] * Isqr[i,j,t] - Vsqr[j,t] 
	<= (Vmax^2 - Vmin^2) * (1 - (y_positivo[i,j,t]+y_negativo[i,j,t]));
	
#14
subject to Restriccion_Recon_8{(i,j) in Ol, t in Ot}:
	(y_positivo[i,j,t] + y_negativo[i,j,t]) <=1;
	
#15
subject to jabra1{i in Ob, t in Ot: Tb[i] == 1}:		###################################
	sum{(i,j) in Ol}(y_positivo[i,j,t]) + sum{(j,i) in Ol}(y_negativo[j,i,t]) >= 1;
	
#16
subject to jabra2{i in Ob, t in Ot: Tb[i] == 0}:        ####################################
	sum{(i,j) in Ol}(y_negativo[i,j,t]) + sum{(j,i) in Ol}(y_positivo[j,i,t]) = 1;
	
#18
subject to PotenciaAparente{(i,j) in Ol, t in Ot}:
	(Vmin^2+(1/2)*DeltaV)*Isqr[i,j,t] + sum{s in 1..S}(Pc[i,j,t,s])= sum{y in 1..Y}(ms[i,j,t,y]*DP[i,j,t,y]) + sum{y in 1..Y}(ms[i,j,t,y]*DQ[i,j,t,y]);
	
#19a
subject to RestriccionLineal_q1{(i,j) in Ol, t in Ot}:
	Qp[i,j,t] - Qn[i,j,t]=Q[i,j,t];

#19b
subject to RestriccionLineal_p2{(i,j) in Ol, t in Ot}:
	Pp[i,j,t] + Pn[i,j,t]=sum{y in 1..Y}(DP[i,j,t,y]);
	
#19c
subject to RestriccionLineal_q2{(i,j) in Ol, t in Ot}:
	Qp[i,j,t] + Qn[i,j,t]=sum{y in 1..Y}(DQ[i,j,t,y]);
	
#19d 0 <= DP
subject to RestriccionLineal_p3{(i,j) in Ol, t in Ot, y in 1..Y}:
	DP[i,j,t,y] <= ds[i,j,t];
	
#19e 0 <= DQ
subject to RestriccionLineal_q3{(i,j) in Ol, t in Ot, y in 1..Y}:
	DQ[i,j,t,y] <= ds[i,j,t];
	
#19f e 19g var definition
# 20 e 21 system data	

#-------------------------------------------------------------------------------------------
#Restricciones linealizacion Vsqr[j] * Isqr[i,j]
	
subject to RestriccionVI_1{i in Ob, t in Ot}:
	Vmin^2+DeltaV*sum{s in 1..S}(xv[i,t,s])<=Vsqr[i,t];

subject to RestriccionVI_2{i in Ob, t in Ot}:
	Vsqr[i,t]<=Vmin^2+ DeltaV + DeltaV*sum{s in 1..S}(xv[i,t,s]);
		
subject to RestriccionVI_3{i in Ob, t in Ot,s in 2..S}:
	xv[i,t,s] <= xv[i,t,s-1];
	
subject to RestriccionVI_4{(i,j) in Ol, t in Ot,s in 1..S}:
	-DeltaV * (Vmax^2)*(Imax[i,j]^2)*(1-xv[j,t,s])<= Pc[i,j,t,s] - DeltaV*Isqr[i,j,t];

subject to RestriccionVI_5{(i,j) in Ol, t in Ot,s in 1..S}:
	Pc[i,j,t,s] - DeltaV * Isqr[i,j,t] <= DeltaV *(Vmax^2)* Imax[i,j]^2 * (1-xv[j,t,s]);

subject to RestriccionVI_6{(i,j) in Ol, t in Ot, s in 1..S}:
	Pc[i,j,t,s] <= DeltaV*(Vmax^2)*(Imax[i,j]^2)*xv[j,t,s];

		

	